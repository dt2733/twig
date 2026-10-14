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

from twig.labeling.core.model import RootMerge

# how many tips of a clade are named, and how many unlabeled branches are listed
CLADE_SHOWN = 8
UNLABELED_SHOWN = 30


def names_below(tree, name_of: Callable[[str], str], limit: int) -> dict[int, list[str]]:
    """Return the names of the tips below every node, each once, in tree order.

    `name_of` turns a tip name into the name listed, such as its species. A
    list is cut after `limit` + 1 names, so a list longer than `limit` means
    that there are more, and repeated names do not count toward the limit.
    """
    keep = limit + 1
    names: dict[int, list[str]] = {}
    for node in tree:
        if node.is_leaf():
            names[node.idx] = [name_of(node.name)]
            continue
        # a child with more than limit names was cut, but then this node is cut as well
        joined: dict[str, None] = {}
        for child in node.children:
            for name in names[child.idx]:
                joined.setdefault(name)
                if len(joined) >= keep:
                    break
            if len(joined) >= keep:
                break
        names[node.idx] = list(joined)
    return names


def clade_text(tree, species_of: Callable[[str], str] | None = None, limit: int = CLADE_SHOWN) -> dict[int, str]:
    """Return the tips below every node as text, in tree order, cut after `limit` names.

    With `species_of` each species is named once, and repeated species do not
    count toward the limit.
    """
    names = names_below(tree, species_of or (lambda name: name), limit)
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
    merge: RootMerge | None = None,
) -> list[str]:
    """Return the lines that describe the labeled branches of a tree.

    The labeled branches come first, clade by clade from the top of the tree, and
    then the unlabeled ones, so that stems left without a label are not missed.
    With a root merge the lines describe the tree as it is written unrooted, so
    the dissolved branch at the root is left out and the kept one has the
    merged label.
    """
    kinds, sources = kinds or {}, sources or {}
    clades = clade_text(tree, species_of, limit)
    if merge is not None:
        labels = {**labels, merge.kept: merge.label}
    skip = {merge.dissolved} if merge is not None else set()
    nonroot = [n for n in tree.traverse("preorder") if n.up is not None and n.idx not in skip]
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
    if merge is not None:
        kept, dissolved = clades[merge.kept], clades[merge.dissolved]
        lines.append(f"root     the branches {kept} and {dissolved} at the root are written as one branch, listed as {kept}")

    counts = Counter(labels[n.idx] for n in labeled)
    parts = [f"{k} {label}" for label, k in counts.most_common()]
    parts.append(f"{len(unlabeled)} unlabeled")
    lines.append("summary  " + ", ".join(parts))
    return lines
