#!/usr/bin/env python

"""Resolve tip and sequence names to species names.

Gene trees and alignments often name one species in several ways, for
example `A_pumilus_sample1` and `AmaranthusPumilus-Sample2`. A `NameMapper`
resolves every raw name to a species by trying, in order, explicit map rows,
optional name parsing by delimiter or regex, and then exact, normalized,
token and genus-abbreviation matches against a vocabulary (the tips of the
species tree). Normalization is used only to compare names. Raw names are
never changed.

This module does not import toytree.
"""

from __future__ import annotations

import difflib
import re
import unicodedata
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable, Mapping, Sequence, TextIO

from twig.utils import TwigError
from twig.labeling.core.batch import atomic_write
from twig.labeling.core.constants import DEFAULT_NORMALIZE, NORMALIZE_STEPS

TABLE_COLUMNS = (
    "name",
    "species",
    "specimen",
    "method",
    "status",
    "leftover",
    "candidates",
    "n_files",
)

# a genus abbreviation can have at most this many characters
ABBREV_MAX = 3

# number of close vocabulary entries listed for names that cannot be resolved
N_SUGGESTIONS = 3


###########################################################################
# normalization
###########################################################################

def normalize_steps(values: Sequence[str] | None) -> tuple[str, ...]:
    """Return a validated tuple of normalization steps.

    None gives the default steps and ["none"] gives no normalization.
    """
    if values is None:
        return DEFAULT_NORMALIZE
    values = tuple(values)
    if "none" in values:
        if len(values) > 1:
            raise TwigError("normalize step 'none' cannot be combined with other steps")
        return ()
    bad = [i for i in values if i not in NORMALIZE_STEPS]
    if bad:
        raise TwigError(f"unknown normalize step(s) {bad}; choose from {list(NORMALIZE_STEPS)} or 'none'")
    return tuple(i for i in NORMALIZE_STEPS if i in values)


def _strip_accents(text: str) -> str:
    """Return text with accents removed (NFKD, combining marks dropped)."""
    decomposed = unicodedata.normalize("NFKD", text)
    return "".join(ch for ch in decomposed if not unicodedata.combining(ch))


def _split_camel(text: str) -> str:
    """Return text with '_' inserted at lower-to-upper and acronym boundaries."""
    out = []
    last = len(text) - 1
    for idx, ch in enumerate(text):
        if idx:
            prev = text[idx - 1]
            following = text[idx + 1] if idx < last else ""
            lower_to_upper = prev.islower() and ch.isupper()
            acronym_end = prev.isupper() and ch.isupper() and following.islower()
            if lower_to_upper or acronym_end:
                out.append("_")
        out.append(ch)
    return "".join(out)


def _split_digits(text: str) -> str:
    """Return text with '_' inserted at letter-digit boundaries."""
    return re.sub(r"(?<=[^\W\d_])(?=\d)|(?<=\d)(?=[^\W\d_])", "_", text)


def _unify_separators(text: str) -> str:
    """Return text with each run of non-alphanumerics replaced by one '_'."""
    return re.sub(r"[\W_]+", "_", text).strip("_")


def canonical(text: str, steps: Sequence[str] = DEFAULT_NORMALIZE) -> str:
    """Return the canonical form of a name used only for comparisons.

    Steps are always applied in the order accents, camel, digits, sep, case,
    whatever the order given.

    Examples
    --------
    >>> canonical("AmaranthusPumilus-Sample1")
    'amaranthus_pumilus_sample1'
    >>> canonical("A.pumilus-s1", steps=("sep",))
    'A_pumilus_s1'
    """
    out = text
    if "accents" in steps:
        out = _strip_accents(out)
    if "camel" in steps:
        out = _split_camel(out)
    if "digits" in steps:
        out = _split_digits(out)
    if "sep" in steps:
        out = _unify_separators(out)
    if "case" in steps:
        out = out.casefold()
    return out


def tokens(canon: str) -> tuple[str, ...]:
    """Return the alphanumeric tokens of a canonical name."""
    return tuple(re.findall(r"[^\W_]+", canon))


def _remove_run(items: Sequence[str], run: Sequence[str]) -> tuple[str, ...]:
    """Return items without the first contiguous run, or without those values."""
    items = tuple(items)
    run = tuple(run)
    size = len(run)
    if size:
        for start in range(len(items) - size + 1):
            if items[start:start + size] == run:
                return items[:start] + items[start + size:]
    remaining = list(items)
    for item in run:
        if item in remaining:
            remaining.remove(item)
    return tuple(remaining)


###########################################################################
# map files
###########################################################################

@dataclass(frozen=True)
class MapRow:
    """One row of a map file."""
    sample: str
    species: str
    specimen: str | None
    is_regex: bool
    source: str
    pattern: re.Pattern | None = None


def _strip_comment(line: str) -> str:
    """Return a map file line without its comment (a # that starts a token)."""
    match = re.search(r"(^|\s)#", line)
    return line[:match.start()] if match else line


def parse_map_file(path: str | Path) -> list[MapRow]:
    """Return the rows of a map file.

    Rows have one to three whitespace-separated columns, `<sample> [<species>
    [<specimen>]]`. A one-column row maps a name to itself, and a sample that
    starts with `~` is a regular expression applied with `re.search`. A row
    whose sample holds a comma is refused when the row has one column or the
    sample ends with the comma, because the file is then most likely
    comma-separated and every row would map a name to itself.
    """
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8-sig")       # utf-8-sig drops a byte order mark
    except OSError as exc:
        raise TwigError(f"cannot read map file {path}: {exc}") from exc
    except UnicodeDecodeError as exc:
        raise TwigError(f"map file {path} is not valid UTF-8 text") from exc

    rows: list[MapRow] = []
    first_seen: dict[str, int] = {}
    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = _strip_comment(raw).strip()
        if not line:
            continue
        where = f"{path}:{lineno}"
        parts = line.split()
        if len(parts) > 3:
            raise TwigError(f"{where}: expected 1 to 3 columns (<sample> [<species> [<specimen>]]), found {len(parts)}")
        sample = parts[0]
        species = parts[1] if len(parts) > 1 else None
        specimen = parts[2] if len(parts) > 2 else None
        if not sample.startswith("~") and "," in sample and (species is None or sample.endswith(",")):
            raise TwigError(f"{where}: the row looks comma-separated, but the columns of a map file are separated by spaces or tabs")

        if sample.startswith("~"):
            expression = sample[1:]
            if not expression:
                raise TwigError(f"{where}: empty regular expression")
            if species is None:
                raise TwigError(f"{where}: regular expression rows need a species column")
            try:
                compiled = re.compile(expression)
            except re.error as exc:
                raise TwigError(f"{where}: invalid regular expression {expression!r} ({exc})") from exc
            rows.append(MapRow(expression, species, specimen, True, where, compiled))
            continue

        if sample in first_seen:
            raise TwigError(f"{where}: duplicate map row for {sample!r} (first seen at line {first_seen[sample]})")
        first_seen[sample] = lineno
        rows.append(MapRow(sample, species or sample, specimen, False, where))

    if not rows:
        raise TwigError(f"{path}: no map rows found")
    return rows


###########################################################################
# vocabulary
###########################################################################

def read_vocab_list(path: str | Path) -> list[str]:
    """Return the species names listed in a file, one per line.

    Blank lines and lines that start with # are skipped, and only the first
    word of a line is used.
    """
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8-sig")       # utf-8-sig drops a byte order mark
    except OSError as exc:
        raise TwigError(f"cannot read species list {path}: {exc}") from exc
    except UnicodeDecodeError as exc:
        raise TwigError(f"species list {path} is not valid UTF-8 text") from exc
    names = [line.split()[0] for line in text.splitlines() if line.strip() and not line.strip().startswith("#")]
    if not names:
        raise TwigError(f"{path}: no species names found")
    return names


class Vocabulary:
    """The species names that raw names are matched to.

    `match` tries exact, canonical, token and abbreviation matches in that
    order and returns the first step that finds any species, keeping only the
    species that cover the most tokens.
    """

    def __init__(self, species: Iterable[str], steps: Sequence[str] = DEFAULT_NORMALIZE):
        names: list[str] = []
        seen: set[str] = set()
        for item in species:
            if item not in seen:
                seen.add(item)
                names.append(item)
        if not names:
            raise TwigError("the species vocabulary is empty")
        self.names = tuple(names)
        self.steps = tuple(steps)
        self.canon_of = {i: canonical(i, self.steps) for i in self.names}
        self.tokens_of = {i: tokens(self.canon_of[i]) for i in self.names}

        self.by_canon: dict[str, list[str]] = {}
        self._by_first: dict[str, list[str]] = {}
        self._by_prefix: dict[str, list[str]] = {}
        for name in self.names:
            self.by_canon.setdefault(self.canon_of[name], []).append(name)
            parts = self.tokens_of[name]
            if parts:
                self._by_first.setdefault(parts[0], []).append(name)
            if len(parts) > 1:
                for size in range(1, min(ABBREV_MAX, len(parts[0])) + 1):
                    self._by_prefix.setdefault(parts[0][:size], []).append(name)

    def __contains__(self, name: str) -> bool:
        return name in self.canon_of

    def match(self, text: str) -> tuple[str, list[str], tuple[str, ...]]:
        """Return (method, species, used tokens) for the best match of text.

        The method is exact, canon, token, abbrev or none. More than one
        species is returned only when several tie, and the used tokens are
        the tokens of `text` that matched the species name.
        """
        if text in self.canon_of:
            return "exact", [text], tokens(canonical(text, self.steps))

        canon = canonical(text, self.steps)
        found = self.by_canon.get(canon)
        if found:
            return "canon", list(found), tokens(canon)

        parts = tokens(canon)
        if not parts:
            return "none", [], ()

        hits = self._token_hits(parts)
        if hits:
            return "token", *self._best(hits, parts)

        hits = self._abbrev_hits(parts)
        if hits:
            return "abbrev", *self._best(hits, parts)
        return "none", [], ()

    def _token_hits(self, parts: tuple[str, ...]) -> list[tuple[int, str, int, tuple[str, ...]]]:
        """Return (n tokens, species, start, used) for token-run matches."""
        hits = []
        for start, token in enumerate(parts):
            for name in self._by_first.get(token, ()):
                run = self.tokens_of[name]
                if parts[start:start + len(run)] == run:
                    hits.append((len(run), name, start, run))
        return hits

    def _abbrev_hits(self, parts: tuple[str, ...]) -> list[tuple[int, str, int, tuple[str, ...]]]:
        """Return (n tokens, species, start, used) for genus-abbreviation matches."""
        hits = []
        for start, token in enumerate(parts):
            if len(token) > ABBREV_MAX:
                continue
            for name in self._by_prefix.get(token, ()):
                run = self.tokens_of[name]
                if parts[start + 1:start + len(run)] == run[1:]:
                    hits.append((len(run), name, start, parts[start:start + len(run)]))
        return hits

    @staticmethod
    def _best(hits: list[tuple[int, str, int, tuple[str, ...]]], parts: tuple[str, ...]) -> tuple[list[str], tuple[str, ...]]:
        """Return the species covering the most tokens, and the tokens used."""
        longest = max(i[0] for i in hits)
        best = [i for i in hits if i[0] == longest]
        species: list[str] = []
        for _, name, _, _ in best:
            if name not in species:
                species.append(name)
        first = next(i for i in best if i[1] == species[0])
        return species, first[3]

    def suggest(self, text: str, n: int = N_SUGGESTIONS) -> tuple[tuple[str, float], ...]:
        """Return the n vocabulary entries most similar to text, with scores.

        Suggestions are only shown to users and are never used to assign a
        species.
        """
        canon = canonical(text, self.steps)
        scored = [
            (name, difflib.SequenceMatcher(None, canon, self.canon_of[name]).ratio())
            for name in self.names
        ]
        scored.sort(key=lambda x: (-x[1], x[0]))
        return tuple((name, round(score, 3)) for name, score in scored[:n])


###########################################################################
# mapping
###########################################################################

@dataclass(frozen=True)
class NameHit:
    """The result of resolving one raw name.

    The status is ok, inferred (an abbreviation match that has not been
    accepted), ambiguous or unmapped. The specimen is set only when a map row
    or a regex group gave one explicitly, and `leftover` lists the tokens of
    the name that are not part of the species.
    """
    name: str
    species: str | None
    specimen: str | None
    method: str
    status: str
    leftover: tuple[str, ...] = ()
    candidates: tuple[tuple[str, float], ...] = ()

    @property
    def resolved(self) -> bool:
        """True if the name was assigned a species that can be used."""
        return self.status == "ok"

    @property
    def specimen_key(self) -> str:
        """The explicit specimen, else the leftover tokens joined by '_'."""
        return self.specimen or "_".join(self.leftover)


class NameMapper:
    """Resolve raw names to species names.

    Precedence is exact map row, regex map row, name parsing followed by
    vocabulary matching, vocabulary matching on the whole name, and finally
    unmapped. Without a vocabulary the parsed key (or the raw name) is the
    species.

    Parameters
    ----------
    map_rows: Sequence[MapRow]
        Rows from `parse_map_file`.
    regex: str
        Regex with a named group `species` and an optional group `specimen`.
        Names that do not match fall through to the next step.
    delim, delim_idxs, delim_join:
        Parse a species key by splitting the name on a delimiter, keeping the
        listed indices, and joining them (by the delimiter if no join is set).
        Names that cannot be split as asked, including names without the
        delimiter, fall through.
    vocabulary: Iterable[str]
        Species names that names are matched to.
    normalize: Sequence[str]
        Normalization steps used when comparing names.
    accept_inferred: bool
        If True, genus-abbreviation matches are accepted with status ok.
    """

    def __init__(
        self,
        map_rows: Sequence[MapRow] = (),
        regex: str | None = None,
        delim: str | None = None,
        delim_idxs: Sequence[int] = (0,),
        delim_join: str | None = None,
        vocabulary: Iterable[str] | None = None,
        normalize: Sequence[str] | None = DEFAULT_NORMALIZE,
        accept_inferred: bool = False,
    ):
        self.steps = normalize_steps(normalize)
        self.accept_inferred = accept_inferred
        self.vocab = Vocabulary(vocabulary, self.steps) if vocabulary is not None else None

        # parsing rules
        if regex is not None and delim is not None:
            raise TwigError("use either a delimiter (--delim) or a name regex (--name-regex), not both")
        self.regex = None
        if regex is not None:
            try:
                self.regex = re.compile(regex)
            except re.error as exc:
                raise TwigError(f"invalid name regex {regex!r} ({exc})") from exc
            if "species" not in self.regex.groupindex:
                raise TwigError("the name regex must contain a named group (?P<species>...)")
        self.delim = delim
        self.delim_idxs = tuple(delim_idxs)
        self.delim_join = delim if delim_join is None else delim_join
        self._extract_method = "regex" if self.regex is not None else "delim"

        # map rows, checked against the vocabulary if there is one
        self._exact_rows: dict[str, MapRow] = {}
        self._regex_rows: list[MapRow] = []
        for row in map_rows:
            if self.vocab is not None:
                row = self._check_row(row)
            if row.is_regex:
                self._regex_rows.append(row)
            else:
                self._exact_rows[row.sample] = row

        self._cache: dict[str, NameHit] = {}
        self._row_vocabs: dict[str, Vocabulary] = {}

    @classmethod
    def from_args(cls, args, vocabulary: Iterable[str] | None = None) -> "NameMapper":
        """Return a NameMapper configured from parsed command line arguments."""
        imap = getattr(args, "imap", None)
        return cls(
            map_rows=parse_map_file(imap) if imap else (),
            regex=getattr(args, "name_regex", None),
            delim=getattr(args, "delim", None),
            delim_idxs=tuple(getattr(args, "delim_idxs", None) or (0,)),
            delim_join=getattr(args, "delim_join", None),
            vocabulary=vocabulary,
            normalize=getattr(args, "normalize", None),
            accept_inferred=bool(getattr(args, "accept_inferred", False)),
        )

    def _check_row(self, row: MapRow) -> MapRow:
        """Return the row with its species spelled as in the vocabulary."""
        method, found, _ = self.vocab.match(row.species)
        if method in ("exact", "canon") and len(found) == 1:
            return replace(row, species=found[0])
        near = ", ".join(i[0] for i in self.vocab.suggest(row.species))
        raise TwigError(f"{row.source}: species {row.species!r} is not in the species tree (closest are {near})")

    def resolve(self, name: str) -> NameHit:
        """Return the NameHit for a raw name (results are memoized)."""
        hit = self._cache.get(name)
        if hit is None:
            hit = self._resolve(name)
            self._cache[name] = hit
        return hit

    def resolve_all(self, names: Iterable[str]) -> dict[str, NameHit]:
        """Return a dict mapping each raw name to its NameHit."""
        return {name: self.resolve(name) for name in names}

    def hits(self) -> list[NameHit]:
        """Return all hits resolved so far, in the order first resolved."""
        return list(self._cache.values())

    def _extract(self, name: str) -> tuple[str | None, str | None]:
        """Return the species key and specimen parsed from a name, if any."""
        if self.regex is not None:
            found = self.regex.search(name)
            if not found:
                return None, None
            return (found.group("species") or None), (found.groupdict().get("specimen") or None)
        if self.delim is not None:
            if self.delim not in name:
                return None, None
            parts = name.split(self.delim)
            try:
                kept = [parts[i] for i in self.delim_idxs]
            except IndexError:
                return None, None
            return (self.delim_join.join(kept) or None), None
        return None, None

    def _leftover(self, name: str, used: Sequence[str]) -> tuple[str, ...]:
        """Return the canonical tokens of name that are not in the species part."""
        return _remove_run(tokens(canonical(name, self.steps)), used)

    def _from_row(self, name: str, row: MapRow, method: str, specimen: str | None) -> NameHit:
        """Return the hit for a name that a map row resolves.

        The leftover words are found by matching the name to the species of
        the row alone, so a sample keeps the same specimen key when its species
        comes from a frozen map instead of from the vocabulary.
        """
        single = self._row_vocabs.get(row.species)
        if single is None:
            single = self._row_vocabs[row.species] = Vocabulary([row.species], self.steps)
        _, found, used = single.match(name)
        if not found:
            used = tokens(canonical(row.species, self.steps))
        return NameHit(name, row.species, row.specimen or specimen, method, "ok", self._leftover(name, used))

    def _from_vocab(self, name: str, text: str, specimen: str | None) -> NameHit:
        """Match text (the name or a key parsed from it) against the vocabulary."""
        method, found, used = self.vocab.match(text)
        if not found:
            return NameHit(name, None, specimen, "none", "unmapped", (), self.vocab.suggest(name))
        if len(found) > 1:
            scored = self.vocab.suggest(name, n=len(self.vocab.names))
            ties = tuple(i for i in scored if i[0] in found)
            return NameHit(name, None, specimen, method, "ambiguous", (), ties)
        status = "ok"
        if method == "abbrev" and not self.accept_inferred:
            status = "inferred"
        return NameHit(name, found[0], specimen, method, status, self._leftover(name, used))

    def _resolve(self, name: str) -> NameHit:
        key, specimen = self._extract(name)

        # 1. exact map rows, by raw name then by parsed key
        for probe in (name, key):
            row = self._exact_rows.get(probe) if probe is not None else None
            if row is not None:
                return self._from_row(name, row, "map", specimen)

        # 2. regex map rows
        for row in self._regex_rows:
            if row.pattern.search(name):
                return self._from_row(name, row, "map-regex", specimen)

        # 3. parsed key, matched to the vocabulary if there is one
        if key is not None:
            if self.vocab is None:
                leftover = self._leftover(name, tokens(canonical(key, self.steps)))
                return NameHit(name, key, specimen, self._extract_method, "ok", leftover)
            hit = self._from_vocab(name, key, specimen)
            if hit.status != "unmapped":
                return hit

        # 4. the whole name, matched to the vocabulary
        if self.vocab is not None:
            return self._from_vocab(name, name, specimen)

        # 5. nothing to match against, so the name is the species
        return NameHit(name, name, specimen, "self", "ok")


###########################################################################
# writers
###########################################################################

def _format_candidates(candidates: Sequence[tuple[str, float]]) -> str:
    return ";".join(f"{name}({score:g})" for name, score in candidates)


def write_hit_table(
    hits: Iterable[NameHit],
    handle: TextIO,
    counts: Mapping[str, int] | None = None,
    header: bool = True,
) -> None:
    """Write hits as a tab-separated table, using '-' for empty cells.

    The n_files column gives the number of input files in which each name
    occurs, if counts are provided.
    """
    if header:
        handle.write("\t".join(TABLE_COLUMNS) + "\n")
    for hit in hits:
        row = (
            hit.name,
            hit.species or "-",
            hit.specimen or "-",
            hit.method,
            hit.status,
            ",".join(hit.leftover) or "-",
            _format_candidates(hit.candidates) or "-",
            str(counts[hit.name]) if counts and hit.name in counts else "-",
        )
        handle.write("\t".join(row) + "\n")


def write_frozen_map(hits: Iterable[NameHit], path: str | Path, force: bool = False) -> None:
    """Write a map file of resolved names that can be reused with --imap.

    Names with status ok or inferred are written. The specimen column is
    written only when a specimen was given explicitly, because otherwise a
    reloaded map would make every sample its own specimen and silently turn off
    the detection of duplicate samples. Unresolved names are listed in comments.
    """
    hits = list(hits)
    lines = ["# map file written by twig name-map", "# <sample> <species> [<specimen>]"]
    unresolved = []
    rows = []
    for hit in hits:
        if hit.status in ("ok", "inferred") and hit.species:
            if re.search(r"\s", hit.name) or hit.name[0] in "~#":
                raise TwigError(f"name {hit.name!r} cannot be written to a map file")
            cols = [hit.name, hit.species] + ([hit.specimen] if hit.specimen else [])
            rows.append(" ".join(cols) + ("  # inferred by genus abbreviation, review" if hit.status == "inferred" else ""))
        else:
            unresolved.append(hit)
    lines.extend(rows)
    for hit in unresolved:
        suggestions = _format_candidates(hit.candidates)
        lines.append(f"# {hit.status.upper()} {hit.name}" + (f" {suggestions}" if suggestions else ""))
    atomic_write(path, "\n".join(lines) + "\n", force=force)
