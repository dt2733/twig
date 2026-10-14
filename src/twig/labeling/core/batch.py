#!/usr/bin/env python

"""Find and pair input files, name loci, and write outputs safely.

Inputs can be files, directories or quoted globs, which are expanded here so
that thousands of files never reach the shell's argument limit. Outputs are
written to a temporary file and renamed, so that a run that is killed never
leaves a file that looks finished.

This module does not import toytree.
"""

from __future__ import annotations

import glob
import hashlib
import json
import os
import re
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence

from twig.utils import TwigError
from twig.labeling.core.constants import DEFAULT_LOCUS_REGEX
from twig.labeling.core.treeio import iter_trees

RUN_RECORD = "twig-run.json"

# settings that do not change outputs, so they are ignored when resuming
RESUME_IGNORE = ("jobs", "log_level", "force", "resume", "check")

_GLOB_CHARS = re.compile(r"[*?\[]")


###########################################################################
# input files
###########################################################################

def _listed(name: str) -> bool:
    """Return False for a hidden file or a temporary file, which directories and globs leave out."""
    return not name.startswith(".") and not name.endswith(".tmp")


def expand_inputs(items: Iterable[str | Path]) -> list[Path]:
    """Return the files named by paths, directories and globs.

    A directory gives its regular files and a value with glob characters gives
    the regular files it matches, in both cases without hidden files and the
    temporary files that `atomic_write` leaves when a run is killed. Any other
    value must be an existing file. Results are sorted within each item and
    duplicates are dropped, keeping the first occurrence.
    """
    found: list[Path] = []
    seen: set[Path] = set()
    for item in items:
        text = os.fspath(item)
        if os.path.isdir(text):
            paths = [
                Path(text) / name for name in sorted(os.listdir(text))
                if _listed(name) and os.path.isfile(os.path.join(text, name))
            ]
            if not paths:
                raise TwigError(f"directory {text} contains no files")
        elif _GLOB_CHARS.search(text) and not os.path.exists(text):
            paths = [
                Path(i) for i in sorted(glob.glob(text, recursive=True))
                if _listed(os.path.basename(i)) and os.path.isfile(i)
            ]
            if not paths:
                raise TwigError(f"no files match {text!r}")
        elif os.path.isfile(text):
            paths = [Path(text)]
        else:
            raise TwigError(f"input file not found: {text}")
        for path in paths:
            resolved = path.resolve()
            if resolved not in seen:
                seen.add(resolved)
                found.append(path)
    return found


_TREE_START = re.compile(r"^(\[[^\]]*\]\s*)*\(")
_TREE_WITH_LOCUS = re.compile(r"^[^\s(\[;]+\t+(\[[^\]]*\]\s*)*\(")
_COMMENTS_ONLY = re.compile(r"^(\[[^\]]*\]\s*)+$")
_PAML_HEADER = re.compile(r"^\d+\s+\d+$")


def sniff_kind(path: str | Path) -> str:
    """Return `fasta`, `tree` or `unknown` from the first lines of a file.

    A line of two numbers starts a PAML tree file and also a PHYLIP alignment,
    so the line after it decides. A line made only of [...] comments is
    skipped.
    """
    header = False
    try:
        with open(path, "r", encoding="utf-8-sig", errors="replace") as data:      # utf-8-sig drops a byte order mark
            for _, line in zip(range(200), data):
                text = line.strip()
                if not text:
                    continue
                if text.upper().startswith("#NEXUS") and not header:
                    return "tree"
                if text.startswith("#") or _COMMENTS_ONLY.match(text):
                    continue
                if text.startswith(">") and not header:
                    return "fasta"
                if _TREE_START.match(text) or _TREE_WITH_LOCUS.match(line.lstrip()):
                    return "tree"
                if _PAML_HEADER.match(text) and not header:
                    header = True
                    continue
                return "unknown"
    except OSError as exc:
        raise TwigError(f"cannot read {path}: {exc}") from exc
    return "unknown"


def locus_name(path: str | Path, regex: str = DEFAULT_LOCUS_REGEX) -> str:
    """Return the locus name taken from a file name.

    The name is the first group of the regex (or the whole match if it has no
    groups) applied with `re.search` to the file name, so with the default
    regex `OG0001.nt.fa` gives `OG0001`.
    """
    name = Path(path).name
    try:
        found = re.search(regex, name)
    except re.error as exc:
        raise TwigError(f"invalid locus regex {regex!r} ({exc})") from exc
    if not found:
        raise TwigError(f"locus regex {regex!r} does not match file name {name!r}")
    locus = found.group(1) if found.re.groups else found.group(0)
    if not locus:
        raise TwigError(f"locus regex {regex!r} gives an empty name for {name!r}")
    return locus


def tree_records(paths: Iterable[str | Path], regex: str = DEFAULT_LOCUS_REGEX) -> Iterator[tuple[str, str, Path]]:
    """Yield (locus, newick, path) for every tree of every tree file.

    A file with one tree is named by `locus_name`. A file with several trees
    gives each tree the name `<locus>_<number>`, counting from 1, unless the
    tree has its own name in the form `locus<TAB>newick`. Locus names must be
    unique across all files.
    """
    used: dict[str, Path] = {}
    for path in paths:
        path = Path(path)
        records = list(iter_trees(path))
        if not records:
            raise TwigError(f"{path}: no trees found")
        base = None
        for number, (named, newick) in enumerate(records, start=1):
            if named:
                locus = named
            else:
                base = base or locus_name(path, regex)
                locus = base if len(records) == 1 else f"{base}_{number}"
            if locus in used:
                raise TwigError(f"locus name {locus!r} is used by trees in {used[locus]} and {path}")
            used[locus] = path
            yield locus, newick, path


@dataclass(frozen=True)
class Pairing:
    """Trees paired to alignments by locus name."""
    pairs: tuple[tuple[str, Path, Path | None], ...]
    unpaired_trees: tuple[str, ...]
    orphan_alignments: tuple[str, ...]


def pair_by_locus(
    trees: Sequence[Path],
    alignments: Sequence[Path] = (),
    regex: str = DEFAULT_LOCUS_REGEX,
) -> Pairing:
    """Pair tree files to alignment files that share a locus name.

    Duplicate locus names are an error. Trees without an alignment are paired
    with None and listed in `unpaired_trees` (when alignments were given), and
    alignments without a tree are listed in `orphan_alignments`.
    """
    tree_by_locus = _by_locus(trees, regex, "tree")
    aln_by_locus = _by_locus(alignments, regex, "alignment")
    pairs = []
    for locus, tree in tree_by_locus.items():
        pairs.append((locus, tree, aln_by_locus.get(locus)))
    unpaired = tuple(i for i, _, a in pairs if a is None) if alignments else ()
    orphans = tuple(i for i in aln_by_locus if i not in tree_by_locus)
    return Pairing(tuple(pairs), unpaired, orphans)


def _by_locus(paths: Sequence[Path], regex: str, kind: str) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for path in paths:
        locus = locus_name(path, regex)
        if locus in out:
            raise TwigError(f"locus name {locus!r} is shared by {kind} files {out[locus]} and {path}")
        out[locus] = path
    return out


@dataclass(frozen=True)
class ManifestRow:
    """One row of a manifest."""
    locus: str
    tree: Path
    alignment: Path | None


def read_manifest(path: str | Path) -> list[ManifestRow]:
    """Return the rows of a manifest, `<locus> <tree> [<alignment>]`.

    Columns are separated by whitespace and # starts a comment. A relative path
    is taken relative to the working directory if it exists there, otherwise
    relative to the directory of the manifest.
    """
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8-sig")       # utf-8-sig drops a byte order mark
    except OSError as exc:
        raise TwigError(f"cannot read manifest {path}: {exc}") from exc
    except UnicodeDecodeError as exc:
        raise TwigError(f"manifest {path} is not valid UTF-8 text") from exc

    def locate(value: str, where: str) -> Path:
        candidate = Path(value)
        if candidate.is_file():
            return candidate
        beside = path.parent / value
        if beside.is_file():
            return beside
        raise TwigError(f"{where}: file not found: {value}")

    rows: list[ManifestRow] = []
    seen: dict[str, int] = {}
    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = re.sub(r"(^|\s)#.*$", "", raw).strip()
        if not line:
            continue
        where = f"{path}:{lineno}"
        parts = line.split()
        if len(parts) not in (2, 3):
            raise TwigError(f"{where}: expected 2 or 3 columns (<locus> <tree> [<alignment>]), found {len(parts)}")
        if parts[0] in seen:
            raise TwigError(f"{where}: duplicate locus {parts[0]!r} (first seen at line {seen[parts[0]]})")
        seen[parts[0]] = lineno
        rows.append(ManifestRow(
            parts[0],
            locate(parts[1], where),
            locate(parts[2], where) if len(parts) == 3 else None,
        ))
    if not rows:
        raise TwigError(f"{path}: no manifest rows found")
    return rows


###########################################################################
# output files
###########################################################################

def write_stdout(text: str) -> None:
    """Write text to stdout, and stay quiet if the reader of a pipe has gone away.

    This avoids a traceback for commands such as `twig name-map ... | head`.
    """
    try:
        sys.stdout.write(text)
        sys.stdout.flush()
    except BrokenPipeError:
        # without this Python reports the error again when it flushes at exit
        try:
            os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
        except (OSError, ValueError):
            pass


def check_writable(path: str | Path, force: bool = False) -> None:
    """Raise a TwigError if path exists and overwriting was not requested."""
    if Path(path).exists() and not force:
        raise TwigError(f"output path exists, use --force to overwrite: {path}")


def atomic_write(path: str | Path, data: str | bytes, force: bool = False) -> None:
    """Write data to path through a temporary file that is renamed when done."""
    path = Path(path)
    check_writable(path, force)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    mode = "wb" if isinstance(data, bytes) else "w"
    kwargs = {} if isinstance(data, bytes) else {"encoding": "utf-8"}
    try:
        with open(tmp, mode, **kwargs) as out:
            out.write(data)
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


###########################################################################
# run record and resuming
###########################################################################

def _version(package: str) -> str | None:
    try:
        return metadata.version(package)
    except metadata.PackageNotFoundError:
        return None


def settings_hash(settings: Mapping[str, object], ignore: Sequence[str] = RESUME_IGNORE) -> str:
    """Return a short hash of settings that affect outputs."""
    kept = {k: v for k, v in settings.items() if k not in ignore}
    blob = json.dumps(kept, sort_keys=True, default=str)
    return hashlib.sha256(blob.encode()).hexdigest()[:16]


def text_hash(text: str | None) -> str | None:
    """Return a short hash of text, or None if there is no text."""
    if text is None:
        return None
    return hashlib.sha256(text.encode()).hexdigest()[:16]


def write_run_record(
    outdir: str | Path,
    command: str,
    settings: Mapping[str, object],
    sptree_text: str | None = None,
) -> dict:
    """Write twig-run.json in outdir and return its content."""
    record = {
        "command": command,
        "settings_hash": settings_hash(settings),
        "sptree_hash": text_hash(sptree_text),
        "twig_version": _version("twig"),
        "toytree_version": _version("toytree"),
        "python": sys.version.split()[0],
        "started": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    atomic_write(Path(outdir) / RUN_RECORD, json.dumps(record, indent=2) + "\n", force=True)
    return record


def check_resume(
    outdir: str | Path,
    settings: Mapping[str, object],
    sptree_text: str | None = None,
) -> dict:
    """Return the saved run record, or raise if settings or species tree differ."""
    path = Path(outdir) / RUN_RECORD
    if not path.exists():
        raise TwigError(f"cannot resume, {path} does not exist")
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise TwigError(f"cannot resume, {path} is not readable ({exc})") from exc
    differs = []
    if record.get("settings_hash") != settings_hash(settings):
        differs.append("settings")
    if record.get("sptree_hash") != text_hash(sptree_text):
        differs.append("species tree")
    if differs:
        raise TwigError(f"cannot resume, the {' and '.join(differs)} differ from the run recorded in {path}")
    return record
