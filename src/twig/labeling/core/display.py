#!/usr/bin/env python

"""A text listing of the labeled branches of a tree, for checking by eye.

toytree's text drawing shows tip names only, so `--show` uses this listing. Each
branch is described by the tips below it, joined by `+`, so `A+B+C` is the
branch above the common ancestor of A, B and C.

This module does not import toytree.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Callable, Mapping

# how many tips of a clade are named, and how many unlabeled branches are listed
CLADE_SHOWN = 8
UNLABELED_SHOWN = 30


def clade_text(tree, species_of: Callable[[str], str] | None = None, limit: int = CLADE_SHOWN) -> dict[int, str]:
    """Return the tips below every node as text, in tree order, cut after `limit` names."""
    keep = limit + 1
    names: dict[int, list[str]] = {}
    for node in tree:
        if node.is_leaf():
            names[node.idx] = [species_of(node.name) if species_of else node.name]
        else:
            joined: list[str] = []
            for child in node.children:
                if len(joined) >= keep:
                    break
                joined.extend(names[child.idx])
            names[node.idx] = list(dict.fromkeys(joined))[:keep]
    out = {}
    for idx, items in names.items():
        out[idx] = "+".join(items[:limit]) + ("+..." if len(items) > limit else "")
    return out


def role_of(node, kinds: Mapping[int, str], sources: Mapping[int, str]) -> str:
    """Return the role of a branch, which is tip, stem, loss, ambiguous, crown or internal."""
    source = sources.get(node.idx, "")
    if source in ("loss", "ambiguous"):
        return source
    return kinds.get(node.idx) or ("tip" if node.is_leaf() else "internal")


def _table(rows: list[tuple[str, ...]]) -> list[str]:
    """Return rows as lines with the columns padded to equal width."""
    if not rows:
        return []
    widths = [max(len(r[i]) for r in rows) for i in range(len(rows[0]) - 1)]
    return ["  " + "  ".join(cell.ljust(w) for cell, w in zip(row, widths)) + "  " + row[-1] for row in rows]


def branch_listing(
    tree,
    labels: Mapping[int, str | None],
    kinds: Mapping[int, str] | None = None,
    sources: Mapping[int, str] | None = None,
    species_of: Callable[[str], str] | None = None,
    limit: int = CLADE_SHOWN,
) -> list[str]:
    """Return the lines that describe the labeled branches of a tree.

    The labeled branches come first, clade by clade from the top of the tree, and
    then the unlabeled ones, so that stems left without a label are not missed.
    """
    kinds, sources = kinds or {}, sources or {}
    clades = clade_text(tree, species_of, limit)
    nonroot = [n for n in tree.traverse("preorder") if n.up is not None]
    labeled = [n for n in nonroot if labels.get(n.idx)]
    unlabeled = [n for n in nonroot if not labels.get(n.idx)]

    lines = ["labeled branches, clade by clade (X+Y is the branch above the common ancestor of X and Y)"]
    rows = [(labels[n.idx], role_of(n, kinds, sources), clades[n.idx]) for n in labeled]
    lines += _table(rows) if rows else ["  none"]

    lines.append("unlabeled branches")
    rows = [("-", role_of(n, kinds, sources), clades[n.idx]) for n in unlabeled[:UNLABELED_SHOWN]]
    lines += _table(rows) if rows else ["  none"]
    if len(unlabeled) > UNLABELED_SHOWN:
        lines.append(f"  and {len(unlabeled) - UNLABELED_SHOWN} more unlabeled branches")

    counts = Counter(labels[n.idx] for n in labeled)
    parts = [f"{k} {label}" for label, k in counts.most_common()]
    parts.append(f"{len(unlabeled)} unlabeled")
    lines.append("summary  " + ", ".join(parts))
    return lines
