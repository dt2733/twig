#!/usr/bin/env python

"""Resolve tip and sequence names to species and report the mapping.

The names of all tips in the trees and all headers in the FASTA files are
resolved to the species of a species tree (or of a species list). The table of
names is printed so that it can be checked, and the resolved mapping can be
frozen into a map file that other twig commands read with `--imap`.

>>> twig name-map -i 'trees/*.treefile' 'aln/*.nt.fa' -s sp.nwk
>>> twig name-map -i 'trees/*.treefile' -s sp.nwk --show problems
>>> twig name-map -i 'trees/*.treefile' 'aln/*.nt.fa' -s sp.nwk --freeze names.tsv
"""

from __future__ import annotations

import io
import os
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

from loguru import logger

from twig.utils import TwigError, set_log_level
from twig.labeling.core import batch, fasta, treeio
from twig.labeling.core.names import NameHit, NameMapper, write_frozen_map, write_hit_table

# how many items of a long list are shown in a log message
SHOWN = 8


@dataclass
class Unit:
    """The names in one tree or one FASTA file."""
    path: Path
    kind: str
    names: list[str]


def read_vocab_list(path: Path) -> list[str]:
    """Return the species names listed in a file, one per line."""
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise TwigError(f"cannot read species list {path}: {exc}") from exc
    names = [line.split()[0] for line in text.splitlines() if line.strip() and not line.strip().startswith("#")]
    if not names:
        raise TwigError(f"{path}: no species names found")
    return names


def read_units(items: list[Path]) -> list[Unit]:
    """Return the names in every tree and FASTA file named by the inputs."""
    units: list[Unit] = []
    for item in items:
        explicit = os.path.isfile(item)
        for path in batch.expand_inputs([item]):
            kind = batch.sniff_kind(path)
            if kind == "unknown":
                if explicit:
                    raise TwigError(f"{path} is neither a tree file nor a FASTA file")
                logger.warning(f"skipping {path}, which is neither a tree file nor a FASTA file")
            elif kind == "fasta":
                units.append(Unit(path, "fasta", fasta.read_headers(path)))
            else:
                for number, (_, newick) in enumerate(treeio.iter_trees(path), start=1):
                    try:
                        names = treeio.tip_names(newick)
                    except TwigError as exc:
                        raise TwigError(f"{path}, tree {number}: {exc}") from None
                    units.append(Unit(path, "tree", names))
    if not units:
        raise TwigError("no names were found in the input files")
    return units


def _shown(items: list[str], limit: int = SHOWN) -> str:
    """Return a comma separated list that is cut after `limit` items."""
    text = ", ".join(items[:limit])
    return text + f" and {len(items) - limit} more" if len(items) > limit else text


def _count(n: int, noun: str) -> str:
    """Return a count with its noun, such as '1 name' or '2 names'."""
    return f"{n} {noun}" if n == 1 else f"{n} {noun}s"


def needs_attention(hit: NameHit) -> bool:
    """Return True if a name is unresolved, inferred, or has leftover words."""
    return hit.status != "ok" or any(i.isalpha() for i in hit.leftover)


def summarize(hits: list[NameHit], units: list[Unit], vocabulary: list[str] | None) -> None:
    """Write a summary of the mapping to the log."""
    n_files = len({u.path for u in units})
    n_trees = sum(u.kind == "tree" for u in units)
    n_fasta = sum(u.kind == "fasta" for u in units)
    logger.info(f"resolved {len(hits)} unique names from {n_files} files ({n_trees} trees, {n_fasta} FASTA files)")
    logger.info("matched by " + ", ".join(f"{k} {v}" for k, v in sorted(Counter(h.method for h in hits).items())))
    logger.info("status " + ", ".join(f"{k} {v}" for k, v in sorted(Counter(h.status for h in hits).items())))

    for status, text in (
        ("ambiguous", "ambiguous, matching more than one species"),
        ("unmapped", "unmapped, matching no species"),
        ("inferred", "inferred, matched only by genus abbreviation (review and freeze with --freeze)"),
    ):
        rows = [h for h in hits if h.status == status]
        if rows:
            detail = [f"{h.name} ({', '.join(c[0] for c in h.candidates[:3]) or h.species})" for h in rows]
            verb = "is" if len(rows) == 1 else "are"
            logger.warning(f"{_count(len(rows), 'name')} {verb} {text}. {_shown(detail)}")

    if vocabulary:
        used = {h.species for h in hits if h.species}
        missing = [s for s in vocabulary if s not in used]
        if missing:
            verb = "has" if len(missing) == 1 else "have"
            logger.info(f"{_count(len(missing), 'species')} {verb} no samples. {_shown(missing)}")

    leftovers: dict[str, Counter] = defaultdict(Counter)
    for hit in hits:
        if hit.species:
            for token in hit.leftover:
                if token.isalpha():
                    leftovers[hit.species][token] += 1
    for species, counter in sorted(leftovers.items()):
        logger.warning(f"words left after the species name in {species} {_shown(sorted(counter))}")

    by_species: dict[str, list[int]] = defaultdict(list)
    species_of = {h.name: h.species for h in hits}
    for unit in units:
        counts = Counter(species_of[n] for n in unit.names if species_of.get(n))
        for species, k in counts.items():
            by_species[species].append(k)
    several = {s: v for s, v in by_species.items() if max(v) > 1}
    for species, counts in sorted(several.items())[:SHOWN]:
        logger.info(
            f"{species} has up to {max(counts)} samples in one tree or FASTA file, "
            f"which happens in {sum(k > 1 for k in counts)} of the {len(counts)} where it occurs"
        )
    if len(several) > SHOWN:
        logger.info(f"{len(several) - SHOWN} more species have several samples in a tree or FASTA file")


def run_name_map(args) -> None:
    """Run name-map subcommand."""
    set_log_level(args.log_level)
    logger.info("CMD: " + " ".join(sys.argv))

    for path in (args.out, args.freeze):
        if path:
            batch.check_writable(path, args.force)
    if args.out and args.freeze and Path(args.out).resolve() == Path(args.freeze).resolve():
        raise TwigError("--out and --freeze must be different files")

    vocabulary = None
    if args.sptree:
        vocabulary = treeio.species_names(args.sptree)
    if args.vocab:
        listed = read_vocab_list(args.vocab)
        vocabulary = list(dict.fromkeys((vocabulary or []) + listed))
    if vocabulary is None:
        logger.warning("no species tree or species list was given, so names are not matched to species")
    else:
        logger.info(f"the species vocabulary has {len(vocabulary)} species")

    units = read_units(args.input)
    mapper = NameMapper.from_args(args, vocabulary=vocabulary)
    order = list(dict.fromkeys(name for unit in units for name in unit.names))
    mapper.resolve_all(order)
    hits = [mapper.resolve(name) for name in order]

    files_of: dict[str, set[Path]] = defaultdict(set)
    for unit in units:
        for name in unit.names:
            files_of[name].add(unit.path)
    n_files = {name: len(paths) for name, paths in files_of.items()}

    rows = [h for h in hits if needs_attention(h)] if args.show == "problems" else hits
    table = io.StringIO()
    write_hit_table(rows, table, counts=n_files)
    if args.out:
        batch.atomic_write(args.out, table.getvalue(), force=args.force)
        logger.info(f"wrote the table of {len(rows)} names to {args.out}")
    else:
        batch.write_stdout(table.getvalue())

    summarize(hits, units, vocabulary)

    problems = [h for h in hits if needs_attention(h)]
    if args.strict and problems:
        raise TwigError(f"{len(problems)} names need attention, see the table and the messages above")
    if args.freeze:
        write_frozen_map(hits, args.freeze, force=args.force)
        logger.info(f"wrote the resolved map to {args.freeze}, use it with -I")
