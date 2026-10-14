#!/usr/bin/env python

"""Strict FASTA reading and writing, and simple alignment and codon checks.

A header is its first whitespace-separated token, and the rest of the header
line is dropped. The reader is strict about empty and duplicate headers, about
records without a sequence and, for alignments, about rows of unequal length.
Several other twig tools have their own FASTA readers, which are deliberately
left unchanged.

This module does not import toytree.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from twig.utils import TwigError

# characters that mean no data in an alignment column
GAP_CHARS = "-?."

NUCLEOTIDES = frozenset("ACGTU")
IUPAC = frozenset("ACGTURYSWKMBDHVN")

# stop codons (as DNA) by NCBI genetic code table number
STOP_CODONS: dict[int, frozenset[str]] = {
    1: frozenset({"TAA", "TAG", "TGA"}),   # standard
    2: frozenset({"TAA", "TAG", "AGA", "AGG"}),   # vertebrate mitochondrial
    3: frozenset({"TAA", "TAG"}),   # yeast mitochondrial
    4: frozenset({"TAA", "TAG"}),   # mold, protozoan and mycoplasma
    5: frozenset({"TAA", "TAG"}),   # invertebrate mitochondrial
    6: frozenset({"TGA"}),   # ciliate nuclear
    9: frozenset({"TAA", "TAG"}),   # echinoderm and flatworm mitochondrial
    10: frozenset({"TAA", "TAG"}),   # euplotid nuclear
    11: frozenset({"TAA", "TAG", "TGA"}),   # bacterial, archaeal and plastid
    12: frozenset({"TAA", "TAG", "TGA"}),   # alternative yeast nuclear
    13: frozenset({"TAA", "TAG"}),   # ascidian mitochondrial
    14: frozenset({"TAG"}),   # alternative flatworm mitochondrial
}


###########################################################################
# reading
###########################################################################

def _read_text(path: str | Path) -> str:
    path = Path(path)
    try:
        return path.read_text(encoding="utf-8-sig")       # utf-8-sig drops a byte order mark
    except OSError as exc:
        raise TwigError(f"cannot read {path}: {exc}") from exc
    except UnicodeDecodeError as exc:
        raise TwigError(f"{path} is not valid UTF-8 text") from exc


def parse_fasta(text: str, source: str = "<text>", aligned: bool = False) -> dict[str, str]:
    """Return FASTA records as an ordered dict of {header: sequence}.

    A header line starts with `>`, after any leading whitespace. Raises a
    TwigError for sequence before the first header, empty or duplicate
    headers, records without a sequence, no records, and (if aligned) rows of
    unequal length.
    """
    records: dict[str, list[str]] = {}
    starts: dict[str, int] = {}
    name = None
    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            fields = line[1:].split()
            if not fields:
                raise TwigError(f"{source}:{lineno}: empty FASTA header")
            name = fields[0]
            if name in records:
                raise TwigError(f"{source}:{lineno}: duplicate FASTA header {name!r}")
            records[name] = []
            starts[name] = lineno
        else:
            if name is None:
                raise TwigError(f"{source}:{lineno}: sequence found before any FASTA header")
            records[name].append("".join(line.split()))
    if not records:
        raise TwigError(f"{source}: no FASTA records")
    for name, lines in records.items():
        if not lines:
            raise TwigError(f"{source}:{starts[name]}: FASTA record {name!r} has no sequence")
    seqs = {k: "".join(v) for k, v in records.items()}
    if aligned:
        lengths = {len(s) for s in seqs.values()}
        if len(lengths) != 1:
            preview = ", ".join(str(i) for i in sorted(lengths)[:6])
            raise TwigError(f"{source}: aligned sequences have different lengths ({preview})")
    return seqs


def read_fasta(path: str | Path, aligned: bool = False) -> dict[str, str]:
    """Return the records of a FASTA file as an ordered dict {header: sequence}."""
    return parse_fasta(_read_text(path), source=str(path), aligned=aligned)


def read_headers(path: str | Path) -> list[str]:
    """Return the headers of a FASTA file, in order, without keeping sequences.

    Headers are found as in `parse_fasta`, and a record without a sequence is
    an error there too. Duplicates are kept, so that callers can report them.
    """
    path = Path(path)
    names = []
    empty = None        # the line of the last header, until a sequence follows it
    try:
        with open(path, "r", encoding="utf-8-sig") as data:
            for lineno, line in enumerate(data, start=1):
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if empty is not None:
                        break
                    fields = line[1:].split()
                    if not fields:
                        raise TwigError(f"{path}:{lineno}: empty FASTA header")
                    names.append(fields[0])
                    empty = lineno
                else:
                    empty = None
    except OSError as exc:
        raise TwigError(f"cannot read {path}: {exc}") from exc
    except UnicodeDecodeError as exc:
        raise TwigError(f"{path} is not valid UTF-8 text") from exc
    if not names:
        raise TwigError(f"{path}: no FASTA records")
    if empty is not None:
        raise TwigError(f"{path}:{empty}: FASTA record {names[-1]!r} has no sequence")
    return names


###########################################################################
# writing and subsetting
###########################################################################

def format_fasta(records: Mapping[str, str] | Iterable[tuple[str, str]], width: int = 0) -> str:
    """Return FASTA text. A width of 0 writes each sequence on one line."""
    items = records.items() if isinstance(records, Mapping) else records
    out = []
    for name, seq in items:
        out.append(f">{name}")
        if width and width > 0:
            out.extend(seq[i:i + width] for i in range(0, len(seq), width))
        else:
            out.append(seq)
    return "\n".join(out) + "\n" if out else ""


def write_fasta(records: Mapping[str, str], dest: str | Path | TextIO, width: int = 0, force: bool = False) -> None:
    """Write records to a path (atomically, refusing to overwrite) or to a text handle."""
    text = format_fasta(records, width)
    if hasattr(dest, "write"):
        dest.write(text)
        return
    from twig.labeling.core.batch import atomic_write
    atomic_write(dest, text, force=force)


def subset_rows(records: Mapping[str, str], names: Iterable[str]) -> dict[str, str]:
    """Return the records for names, in the order of names."""
    out = {}
    for name in names:
        if name not in records:
            raise TwigError(f"sequence {name!r} is not in the alignment")
        out[name] = records[name]
    return out


def drop_empty_columns(records: Mapping[str, str], mode: str = "site") -> dict[str, str]:
    """Return records without columns in which every row has a gap.

    With mode `site` single columns are removed. With mode `codon` a codon
    (three columns) is removed only if all three columns are empty, which keeps
    the reading frame. Mode `none` returns the records unchanged.
    """
    if mode not in ("none", "site", "codon"):
        raise TwigError(f"unknown mode {mode!r} for dropping empty columns (use none, site or codon)")
    if mode == "none" or not records:
        return dict(records)
    names = list(records)
    seqs = [records[i] for i in names]
    length = len(seqs[0])
    if any(len(s) != length for s in seqs):
        raise TwigError("cannot drop empty columns, sequences have different lengths")
    empty = [all(ch in GAP_CHARS for ch in column) for column in zip(*seqs)]
    if mode == "site":
        keep = [i for i, flag in enumerate(empty) if not flag]
    else:
        if length % 3:
            raise TwigError(f"cannot drop empty codons, alignment length {length} is not a multiple of 3")
        keep = []
        for start in range(0, length, 3):
            if not all(empty[start:start + 3]):
                keep.extend(range(start, start + 3))
    return {name: "".join(seq[i] for i in keep) for name, seq in zip(names, seqs)}


###########################################################################
# sequence statistics
###########################################################################

def occupancy(seq: str) -> int:
    """Return the number of A, C, G, T and U characters in a sequence."""
    return sum(1 for ch in seq.upper() if ch in NUCLEOTIDES)


def nongap_length(seq: str) -> int:
    """Return the number of characters that are not gaps (ambiguity codes count)."""
    return sum(1 for ch in seq if ch not in GAP_CHARS)


###########################################################################
# codon checks
###########################################################################

@dataclass(frozen=True)
class CodonReport:
    """Result of checking a nucleotide sequence for use as codons."""
    length: int
    bad_length: bool
    stops: tuple[int, ...]
    bad_chars: tuple[str, ...]

    @property
    def ok(self) -> bool:
        """True if the sequence has no problems."""
        return not (self.bad_length or self.stops or self.bad_chars)


def stop_codons(code: int = 1) -> frozenset[str]:
    """Return the stop codons (as DNA) for an NCBI genetic code table number."""
    try:
        return STOP_CODONS[code]
    except KeyError:
        raise TwigError(f"genetic code {code} is not supported, choose from {sorted(STOP_CODONS)}") from None


def find_stop_codons(seq: str, stops: Iterable[str]) -> tuple[int, ...]:
    """Return the indices (0-based, in codons) of in-frame stop codons."""
    stops = frozenset(stops)
    seq = seq.upper().replace("U", "T")
    return tuple(i // 3 for i in range(0, len(seq) - 2, 3) if seq[i:i + 3] in stops)


def check_codons(seq: str, stops: Iterable[str]) -> CodonReport:
    """Check length, in-frame stop codons and the nucleotide alphabet."""
    allowed = IUPAC | set(GAP_CHARS)
    bad_chars = tuple(sorted({ch for ch in seq.upper() if ch not in allowed}))
    return CodonReport(
        length=len(seq),
        bad_length=len(seq) % 3 != 0,
        stops=find_stop_codons(seq, stops),
        bad_chars=bad_chars,
    )


def mask_stop_codons(seq: str, stops: Iterable[str], fill: str = "---") -> str:
    """Return the sequence with in-frame stop codons replaced by fill."""
    if len(fill) != 3:
        raise TwigError("the fill for a masked codon must have three characters")
    codons = [seq[i:i + 3] for i in range(0, len(seq), 3)]
    for index in find_stop_codons(seq, stops):
        codons[index] = fill
    return "".join(codons)
