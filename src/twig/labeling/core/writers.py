#!/usr/bin/env python

"""Writers that turn a labeled tree into text.

A writer has a name, a file suffix and two methods. `check` lists problems
(errors make the output unusable for the program it is meant for, warnings are
worth a look) and `render` returns the text, raising a TwigError for a tree that
cannot be written. All tree writers use the serializer in `newick.py`, so the
tokens always come in the order `name{Label}:len`.

The writers are `braces` (twig's own format, rooted as given), `hyphy` (names
made safe for HyPhy and the tree written unrooted by an explicit root-edge
policy), `newick` (labels removed) and `table` (one row per branch).

This module does not import toytree.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Callable, Iterable, Iterator, Mapping
from dataclasses import dataclass, field
from typing import ClassVar

from twig.utils import TwigError
from twig.labeling.core import rooting, treeio
from twig.labeling.core.constants import INTERNAL_NAME_CHOICES
from twig.labeling.core.display import names_below
from twig.labeling.core.model import RootMerge
from twig.labeling.core.newick import braces_label, to_newick
from twig.labeling.core.sanitize import NameTable, label_name_problem, sanitize_name
from twig.labeling.core.validate import ERROR, WARNING, Issue

# the table lists at most this many tip names in the tips column, and this many species in the species column
TIPS_SHOWN = 200


@dataclass(frozen=True)
class Plan:
    """How one tree is written.

    `merge` says how the two branches at the root become one (None if the root
    is written as it is), `tip_names` and `internal_names` are the names that
    are written, by idx of the node.
    """
    merge: RootMerge | None
    tip_names: Mapping[int, str]
    internal_names: Mapping[int, str]


@dataclass
class WriteContext:
    """Settings and facts about one tree that the writers use.

    `names` is a table of sanitized tip names built from all trees of a run, so
    that a sample is written with the same name everywhere. `kinds` and
    `sources` describe the branches for the table, and `plan` is the plan of
    the tree writer that is used alongside the table.
    """
    locus: str = ""
    dists: bool = True
    dist_fmt: str = "%.12g"
    root_edges: str = "unlabel"
    keep_root: bool = False
    internal_names: str = "none"
    names: NameTable | None = None
    outgroup_side: int | None = None
    kinds: Mapping[int, str] = field(default_factory=dict)
    sources: Mapping[int, str] = field(default_factory=dict)
    species_of: Callable[[str], str] | None = None
    short_len: float = 1e-6
    header: bool = True
    plan: Plan | None = None


def write_dists(tree, ctx: WriteContext) -> bool:
    """Return True if branch lengths are written, which needs lengths in the tree.

    toytree gives 1.0 to a missing branch length, so a tree read without lengths
    is written without them rather than with invented ones, and so is each
    branch that had no length in a tree with lengths elsewhere.
    """
    return ctx.dists and treeio.has_dists(tree)


def branch_length(tree, node, merge: RootMerge | None = None) -> float | None:
    """Return the length of the branch above a node as it is written, or None if it has none.

    The kept child of a root merge has the summed length of the two branches at
    the root, and no length if either of them had none in the input.
    """
    if not treeio.has_dists(tree) or not treeio.has_dist(node):
        return None
    if merge is not None and node.idx == merge.kept:
        dissolved = tree[merge.dissolved]
        return node.dist + dissolved.dist if treeio.has_dist(dissolved) else None
    return node.dist


def written_branches(tree, merge: RootMerge | None) -> Iterator[tuple[object, str | None, float | None]]:
    """Yield (node, label, length) for every branch that is written.

    The dissolved child of the root is not written, and the kept child has the
    merged label and the summed length. The length is None for a branch that
    had no length in the input.
    """
    for node in tree:
        if node.up is None or (merge is not None and node.idx == merge.dissolved):
            continue
        label = merge.label if merge is not None and node.idx == merge.kept else treeio.node_label(node)
        yield node, label, branch_length(tree, node, merge)


def internal_texts(
    tree,
    mode: str,
    tip_texts: Iterable[str],
    skip: Iterable[int] = (),
    sanitize: Callable[[str], str] | None = None,
) -> dict[int, str]:
    """Return {idx: text} for the internal nodes that are written with a name.

    The modes are none (no names), auto (`twn1`, `twn2` and so on in postorder,
    leaving out numbers that would equal a tip name), support (the support
    value of each node that has one) and keep (the existing names, sanitized
    if a function is given). Names must differ from tip names and from each
    other without regard to case, because HyPhy matches names that way.
    """
    if mode not in INTERNAL_NAME_CHOICES:
        raise TwigError(f"unknown choice {mode!r} for internal names, use one of {list(INTERNAL_NAME_CHOICES)}")
    if mode == "none":
        return {}
    skip = set(skip)
    nodes = [n for n in tree.traverse("postorder") if not n.is_leaf() and n.up is not None and n.idx not in skip]
    taken = {t.lower() for t in tip_texts}
    out: dict[int, str] = {}
    if mode == "auto":
        count = 0
        for node in nodes:
            while True:
                count += 1
                name = f"twn{count}"
                if name.lower() not in taken:
                    break
            out[node.idx] = name
    elif mode == "support":
        for node in nodes:
            value = getattr(node, "support", None)
            if value is not None and value == value:
                out[node.idx] = "%.12g" % value
    else:
        seen: dict[str, str] = {}
        for node in nodes:
            if not node.name:
                continue
            text = sanitize(node.name) if sanitize else node.name
            clash = taken | set(seen)
            if text.lower() in clash:
                raise TwigError(
                    f"the internal name {node.name!r} is used twice or equals a tip name, "
                    "use --internal-names auto or support instead"
                )
            seen[text.lower()] = node.name
            out[node.idx] = text
    return out


###########################################################################
# the writer classes
###########################################################################

class Writer:
    """Base class of the writers."""
    name: ClassVar[str] = ""
    suffix: ClassVar[str] = ""
    unrooted: ClassVar[bool] = False
    sanitizer: ClassVar[str | None] = None
    help: ClassVar[str] = ""

    def check(self, tree, ctx: WriteContext) -> list[Issue]:
        """Return the problems found in a tree for this output."""
        return []

    def render(self, tree, ctx: WriteContext) -> str:
        """Return the text for a tree, without a final newline."""
        raise NotImplementedError


WRITERS: dict[str, type[Writer]] = {}


def register_writer(cls: type[Writer]) -> type[Writer]:
    """Register a writer class under its name."""
    WRITERS[cls.name] = cls
    return cls


def get_writer(name: str) -> Writer:
    """Return a writer object for a name."""
    try:
        return WRITERS[name]()
    except KeyError:
        raise TwigError(f"unknown output format {name!r}, choose from {list(WRITERS)}") from None


class TreeWriter(Writer):
    """A writer that produces Newick text."""
    labels: ClassVar[bool] = True

    def label_problem(self, label: str) -> str | None:
        """Return why a label cannot be written, or None."""
        try:
            braces_label(label)
        except TwigError as exc:
            return str(exc)
        return None

    def merge_for(self, tree, ctx: WriteContext) -> RootMerge | None:
        """Return the root merge, which is planned for unrooted writers only."""
        if not self.unrooted or ctx.keep_root:
            return None
        return rooting.plan_root_merge(tree, ctx.root_edges, ctx.outgroup_side)

    def tip_name(self, node, ctx: WriteContext) -> str:
        if self.sanitizer is None:
            return node.name
        if ctx.names is not None:
            return ctx.names.written(node.name)
        return sanitize_name(node.name, self.sanitizer)

    def plan(self, tree, ctx: WriteContext) -> Plan:
        """Return how the tree is written, or raise a TwigError if it cannot be."""
        merge = self.merge_for(tree, ctx)
        if self.sanitizer is not None:
            repeated = sorted(n for n, k in Counter(leaf.name for leaf in tree if leaf.is_leaf()).items() if k > 1)
            if repeated:
                raise TwigError(
                    f"tip names are not unique ({', '.join(repeated[:5])}), and {self.sanitizer} needs a different name for every tip"
                )
        tips = {n.idx: self.tip_name(n, ctx) for n in tree if n.is_leaf()}
        if self.sanitizer is not None and ctx.names is None:
            groups: dict[str, list[str]] = {}
            for node in tree:
                if node.is_leaf():
                    groups.setdefault(tips[node.idx].lower(), []).append(node.name)
            clashes = [" and ".join(g) for g in groups.values() if len(g) > 1]
            if clashes:
                raise TwigError(
                    f"tip names collide after sanitizing for {self.sanitizer} (compared without regard to case), "
                    + "; ".join(clashes[:5])
                )
        sanitize = (lambda name: sanitize_name(name, self.sanitizer)) if self.sanitizer else None
        internal = internal_texts(
            tree, ctx.internal_names, tips.values(), {merge.dissolved} if merge else (), sanitize
        )
        return Plan(merge, tips, internal)

    def check(self, tree, ctx: WriteContext) -> list[Issue]:
        locus = ctx.locus or None
        issues: list[Issue] = []
        try:
            plan = self.plan(tree, ctx)
        except TwigError as exc:
            return [Issue(ERROR, "write", str(exc), locus)]
        merge = plan.merge

        if merge is not None and merge.conflict:
            first, second = (i or "no label" for i in merge.labels)
            issues.append(Issue(
                WARNING, "root-merge",
                f"the two branches at the root have the labels {first} and {second}, so the merged branch gets "
                f"{merge.label or 'no label'} (--root-edges {ctx.root_edges})",
                locus,
            ))
        if self.unrooted and ctx.keep_root:
            kids = tree.treenode.children
            if len(kids) == 2 and treeio.node_label(kids[0]) != treeio.node_label(kids[1]):
                issues.append(Issue(
                    WARNING, "root-merge",
                    "the tree is written rooted, and the program that reads it merges the two branches at the root "
                    "itself and keeps only one of their labels",
                    locus,
                ))

        if self.labels:
            labels = {label for _, label, _ in written_branches(tree, merge) if label}
            for label in sorted(labels):
                problem = self.label_problem(label)
                if problem:
                    issues.append(Issue(ERROR, "label-name", f"label {label!r} {problem}", locus))

        if write_dists(tree, ctx):
            nan = [node.idx for node, _, dist in written_branches(tree, merge) if dist is not None and dist != dist]
            if nan:
                issues.append(Issue(ERROR, "branch-length", f"{len(nan)} branch lengths are not numbers", locus))
            elif self.sanitizer == "hyphy":
                short = [
                    node for node, label, dist in written_branches(tree, merge)
                    if label and not node.is_leaf() and dist is not None and dist <= ctx.short_len
                ]
                if short:
                    issues.append(Issue(
                        WARNING, "short-labeled",
                        f"{len(short)} labeled internal branches have lengths of at most {ctx.short_len:g}, "
                        "which HyPhy may remove with --kill-zero-lengths Yes, so pass Constrain",
                        locus,
                    ))
        return issues

    def render(self, tree, ctx: WriteContext) -> str:
        plan = self.plan(tree, ctx)
        return to_newick(
            tree,
            tip_text=lambda node: plan.tip_names[node.idx],
            internal_text=lambda node: plan.internal_names.get(node.idx, ""),
            label_text=braces_label if self.labels else (lambda label: ""),
            dists=write_dists(tree, ctx),
            dist_fmt=ctx.dist_fmt,
            root_merge=plan.merge,
        )


@register_writer
class BracesWriter(TreeWriter):
    """twig's own format, `name{Label}:len`, rooted as given."""
    name = "braces"
    suffix = ".labeled.nwk"
    help = "labels in braces after the name, as HyPhy reads them (rooted as given)"


@register_writer
class HyphyWriter(TreeWriter):
    """The HyPhy format, with sanitized names and the tree written unrooted."""
    name = "hyphy"
    suffix = ".hyphy.nwk"
    unrooted = True
    sanitizer = "hyphy"
    help = "braces format with names made safe for HyPhy, written unrooted"

    def label_problem(self, label: str) -> str | None:
        return label_name_problem(label)


@register_writer
class NewickWriter(TreeWriter):
    """A plain Newick tree without labels."""
    name = "newick"
    suffix = ".nwk"
    labels = False
    help = "plain Newick with the labels removed"


def _cut(names: list[str]) -> str:
    """Return names joined by semicolons, cut after TIPS_SHOWN with ;... at the end."""
    return ";".join(names[:TIPS_SHOWN]) + (";..." if len(names) > TIPS_SHOWN else "")


def row_ids(tree) -> dict[int, str]:
    """Return an id for every non-root node for the branch table.

    Tips are identified by their names, and internal nodes by their names if
    these are unique, otherwise by `S` followed by the idx.
    """
    counts = Counter(node.name for node in tree if node.name)
    ids: dict[int, str] = {}
    for node in tree:
        if node.up is None:
            continue
        if node.is_leaf() or (node.name and counts[node.name] == 1):
            ids[node.idx] = node.name
        else:
            ids[node.idx] = f"S{node.idx}"
    return ids


@register_writer
class TableWriter(Writer):
    """A table with one row per branch."""
    name = "table"
    suffix = ".branches.tsv"
    help = "a table with one row per branch"

    COLUMNS: ClassVar[tuple[str, ...]] = (
        "locus", "idx", "node", "is_tip", "n_tips", "tips", "label", "kind", "source",
        "dist", "support", "written_name", "written", "merged_into", "note",
    )

    def columns(self, ctx: WriteContext) -> tuple[str, ...]:
        cols = list(self.COLUMNS)
        if ctx.species_of is not None:
            cols.insert(cols.index("tips") + 1, "species")
        return tuple(cols)

    def render(self, tree, ctx: WriteContext) -> str:
        plan = ctx.plan
        merge = plan.merge if plan else None
        ids = row_ids(tree)
        cols = self.columns(ctx)

        # tip names below every node, in tree order, cut at TIPS_SHOWN
        tips: dict[int, list[str]] = {}
        counts: dict[int, int] = {}
        for node in tree:
            if node.is_leaf():
                tips[node.idx], counts[node.idx] = [node.name], 1
            else:
                kids = node.children
                counts[node.idx] = sum(counts[c.idx] for c in kids)
                joined: list[str] = []
                for child in kids:
                    if len(joined) > TIPS_SHOWN:
                        break
                    joined.extend(tips[child.idx])
                tips[node.idx] = joined[: TIPS_SHOWN + 1]
        # species below every node, each once, found from all tips and not only from the tips listed
        species = names_below(tree, ctx.species_of, TIPS_SHOWN) if ctx.species_of else {}

        lines = ["\t".join(cols)] if ctx.header else []
        for node in tree:
            if node.up is None:
                continue
            dissolved = merge is not None and node.idx == merge.dissolved
            kept = merge is not None and node.idx == merge.kept
            label = merge.label if kept else treeio.node_label(node)
            dist = branch_length(tree, node, merge)
            support = getattr(node, "support", None)
            written_name = ""
            if plan is not None:
                written_name = plan.tip_names.get(node.idx) or plan.internal_names.get(node.idx, "")
            note = ""
            if kept:
                note = "root-merge"
                if merge.conflict:
                    note += " (labels differed)"
            elif dissolved:
                note = "root-merge"
            values = {
                "locus": ctx.locus,
                "idx": str(node.idx),
                "node": ids[node.idx],
                "is_tip": "yes" if node.is_leaf() else "no",
                "n_tips": str(counts[node.idx]),
                "tips": _cut(tips[node.idx]),
                "species": _cut(species[node.idx]) if ctx.species_of else "",
                "label": label or "",
                "kind": ctx.kinds.get(node.idx, ""),
                "source": ctx.sources.get(node.idx, ""),
                "dist": ("%.12g" % dist) if dist is not None else "",
                "support": ("%.12g" % support) if support is not None and support == support else "",
                "written_name": written_name,
                "written": "no" if dissolved else "yes",
                "merged_into": ids[merge.kept] if dissolved else "",
                "note": note,
            }
            lines.append("\t".join(values[c] for c in cols))
        return "\n".join(lines)
