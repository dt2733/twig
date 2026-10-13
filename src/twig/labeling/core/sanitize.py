#!/usr/bin/env python

"""Make tip names and label names safe for the programs that read them.

HyPhy renames any tip or sequence name that is not a valid identifier, and its
tip names are matched to sequence names without regard to case. To keep tree
and alignment names identical in every file, and so that HyPhy never renames
anything, names are sanitized once for the whole run with a `NameTable`, which
is reversible and is written to a table.

The rules below were checked against the HyPhy source (`_String::ConvertToAnIdent`
and the built-in function table in `parser.cpp`) and the PAML source (`LSPNAME`).

This module does not import toytree.
"""

from __future__ import annotations

import re
from collections.abc import Iterable
from pathlib import Path

from twig.utils import TwigError
from twig.labeling.core.batch import atomic_write
from twig.labeling.core.constants import STATUS_NAMES

# characters HyPhy accepts in an identifier (a leading digit is allowed for tips)
_ID_CHARS = frozenset("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_")

# names of HyPhy built-in functions, which cannot be used as tip names
HYPHY_RESERVED = frozenset({
    "Abs", "Arctan", "Beta", "BranchCount", "BranchLength", "BranchName",
    "CChi2", "CGammaDist", "Call", "Columns", "Cos", "Differentiate",
    "Eigensystem", "Erf", "Eval", "Exp", "Format", "Gamma", "GammaDist",
    "IBeta", "IGamma", "InvChi2", "Inverse", "Join", "LUDecompose", "LUSolve",
    "LnGamma", "Log", "MAccess", "MCoord", "Max", "Min", "PSTreeString",
    "Polynomial", "Random", "RerootTree", "Rows", "Simplex", "Simplify", "Sin",
    "Sqrt", "TEXTreeString", "Tan", "Time", "TipCount", "TipName",
    "Transpose", "Type", "ZCDF",
})

# codeml and baseml read species names of at most this many characters (LSPNAME)
PAML_MAX_NAME = 96

# branch-set options of HyPhy analyses, which cannot be used as label names
HYPHY_BRANCH_SETS = frozenset({"all", "internal", "leaves", "unlabeled", "fg"})

# words with a meaning in twig masks and reports, which cannot be label names
TWIG_WORDS = frozenset(
    {"keep", "unlabel", "unlabeled", "prune", "relabel"}
    | {i for i in STATUS_NAMES if re.fullmatch(r"[a-z]+", i)}
)

LABEL_PATTERN = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")


###########################################################################
# tip names
###########################################################################

def hyphy_identifier(name: str) -> str:
    """Return name as HyPhy would convert it to an identifier.

    Every character that is not a letter, digit or underscore becomes an
    underscore, and a run of such characters (or one that follows an
    underscore) gives a single underscore.
    """
    out: list[str] = []
    for ch in name:
        if ch in _ID_CHARS:
            out.append(ch)
        elif not out or out[-1] != "_":
            out.append("_")
    return "".join(out)


def sanitize_hyphy(name: str) -> str:
    """Return a tip name that HyPhy leaves unchanged and that is not reserved."""
    out = hyphy_identifier(name)
    if not out:
        raise TwigError("cannot make a HyPhy name from an empty name")
    if out in HYPHY_RESERVED:
        out += "_"
    return out


def sanitize_paml(name: str) -> str:
    """Return a tip name that is safe for PAML (letters, digits and underscores)."""
    out = hyphy_identifier(name)
    if not out:
        raise TwigError("cannot make a PAML name from an empty name")
    if len(out) > PAML_MAX_NAME:
        raise TwigError(f"name {name!r} is longer than the {PAML_MAX_NAME} characters PAML reads, rename it first")
    return out


RULES = {"hyphy": sanitize_hyphy, "paml": sanitize_paml}


def sanitize_name(name: str, rule: str = "hyphy") -> str:
    """Return name sanitized by a rule (hyphy or paml)."""
    try:
        return RULES[rule](name)
    except KeyError:
        raise TwigError(f"unknown name rule {rule!r}, choose from {sorted(RULES)}") from None


class NameTable:
    """A reversible table from raw names to sanitized names.

    Add every raw name from all trees and alignments, then call `build` once.
    Names that already satisfy the rule keep their spelling. Names that
    collide (compared without regard to case) either stop the run with a list
    or, with suffix=True, get `_2`, `_3` and so on, in sorted order of the raw
    names so that the result does not depend on the order files were read.
    """

    def __init__(self, rule: str = "hyphy", suffix: bool = False):
        if rule not in RULES:
            raise TwigError(f"unknown name rule {rule!r}, choose from {sorted(RULES)}")
        self.rule = rule
        self.suffix = suffix
        self._raw: set[str] = set()
        self._written: dict[str, str] = {}
        self._reverse: dict[str, str] = {}

    def add(self, name: str) -> None:
        """Register a raw name."""
        self._raw.add(name)

    def update(self, names: Iterable[str]) -> None:
        """Register many raw names."""
        self._raw.update(names)

    def build(self) -> dict[str, str]:
        """Compute the sanitized names and return the table {raw: written}."""
        base = {name: sanitize_name(name, self.rule) for name in sorted(self._raw)}
        taken: dict[str, str] = {}   # lowercase written name -> raw name
        written: dict[str, str] = {}
        collisions: dict[str, list[str]] = {}

        # names that need no change are placed first, so they keep their spelling
        unchanged = [n for n in base if base[n] == n]
        changed = [n for n in base if base[n] != n]
        for raw in unchanged + changed:
            candidate = base[raw]
            if candidate.lower() in taken:
                if not self.suffix:
                    group = collisions.setdefault(candidate.lower(), [taken[candidate.lower()]])
                    group.append(raw)
                    continue
                size = 2
                while f"{base[raw]}_{size}".lower() in taken:
                    size += 1
                candidate = f"{base[raw]}_{size}"
            taken[candidate.lower()] = raw
            written[raw] = candidate

        if collisions:
            lines = [f"  {' and '.join(group)}" for group in collisions.values()]
            raise TwigError(
                f"names collide after sanitizing for {self.rule} (compared without regard to case)\n"
                + "\n".join(lines)
                + "\nuse --sanitize-suffix to add _2, _3 to colliding names, or rename the samples"
            )
        self._written = written
        self._reverse = {v: k for k, v in written.items()}
        return dict(written)

    def written(self, raw: str) -> str:
        """Return the sanitized name for a raw name."""
        try:
            return self._written[raw]
        except KeyError:
            raise TwigError(f"name {raw!r} was not added to the name table before it was built") from None

    def raw(self, written: str) -> str:
        """Return the raw name for a sanitized name."""
        try:
            return self._reverse[written]
        except KeyError:
            raise TwigError(f"name {written!r} is not in the name table") from None

    def changed(self) -> dict[str, str]:
        """Return the names whose spelling was changed."""
        return {k: v for k, v in self._written.items() if k != v}

    def write_tsv(self, path: str | Path, force: bool = False) -> None:
        """Write the table (raw, written, changed) as a tab-separated file."""
        lines = ["raw\twritten\tchanged"]
        for raw in sorted(self._written):
            written = self._written[raw]
            lines.append(f"{raw}\t{written}\t{'yes' if raw != written else 'no'}")
        atomic_write(path, "\n".join(lines) + "\n", force=force)


###########################################################################
# label names
###########################################################################

def label_name_problem(label: str) -> str | None:
    """Return why a label name is not allowed, or None if it is fine."""
    if not LABEL_PATTERN.match(label):
        return "must start with a letter and contain only letters, digits and underscores"
    lowered = label.lower()
    if lowered in HYPHY_BRANCH_SETS:
        return "collides with a HyPhy branch-set option (All, Internal, Leaves, Unlabeled, fg)"
    if lowered in TWIG_WORDS:
        return "is a word that twig uses for mask actions and branch statuses"
    return None


def check_label_names(labels: Iterable[str]) -> None:
    """Raise a TwigError listing every label name that is not allowed."""
    problems = []
    for label in labels:
        reason = label_name_problem(label)
        if reason:
            problems.append(f"  {label!r} {reason}")
    if problems:
        raise TwigError("invalid label name(s)\n" + "\n".join(problems))
