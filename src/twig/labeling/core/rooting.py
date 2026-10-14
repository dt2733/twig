#!/usr/bin/env python

"""Roots, outgroups and the two branches at the root.

A rooted tree has two branches at its root. HyPhy and PAML read every tree as
unrooted, so they see those two branches as one. HyPhy keeps only the label of
the first internal root child, which depends on the order in which the children
happen to be written. twig makes the choice itself with an explicit policy
(`plan_root_merge`) and writes the tree as unrooted.

Rerooting an already rooted tree has the same problem. toytree first unroots
the tree. The first internal root child becomes the root, and the branch of the
other child stands for both old root branches with that child's label. The
labels are then moved along the path to the new root, so the label that
survives depends on the order of the children and on where the new root is.
twig merges the two labels by the policy first and gives the result to both
old root children, so it survives whichever branch toytree keeps.

Rooting is decided with sets of tips and not with toytree's monophyly test,
which raises on unrooted trees.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass

from twig.utils import TwigError
from twig.labeling.core import treeio
from twig.labeling.core.constants import ROOT_EDGE_POLICIES
from twig.labeling.core.model import RootMerge
from twig.labeling.core.sanitize import label_name_problem


def check_root_edge_policy(policy: str) -> str:
    """Return the policy if it is a keyword or a valid label name, else raise."""
    if policy in ROOT_EDGE_POLICIES:
        return policy
    problem = label_name_problem(policy)
    if problem:
        raise TwigError(
            f"unknown root edge policy {policy!r}, use one of {list(ROOT_EDGE_POLICIES)} or the name of a label"
        )
    return policy


def merge_root_labels(
    label_a: str | None,
    label_b: str | None,
    policy: str = "unlabel",
    a_is_outgroup: bool | None = None,
) -> tuple[str | None, bool]:
    """Return (label, conflict) for two root branches that become one.

    The policy is `unlabel` (a merged branch of two different labels is
    unlabeled), `labeled` (if only one of the two is labeled its label is used),
    `outgroup` or `ingroup` (the label of the outgroup side or the ingroup
    side), `error`, or the name of a label that wins if either branch carries it.
    The conflict is True if the two labels differ.
    """
    check_root_edge_policy(policy)
    if label_a == label_b:
        return label_a, False
    if policy == "unlabel":
        return None, True
    if policy == "labeled":
        if label_a and not label_b:
            return label_a, True
        if label_b and not label_a:
            return label_b, True
        return None, True
    if policy in ("outgroup", "ingroup"):
        if a_is_outgroup is None:
            raise TwigError(
                "cannot tell which branch at the root leads to the outgroup because the two sides are equally "
                "large, give the outgroups with -r"
            )
        take_a = (policy == "outgroup") == a_is_outgroup
        return (label_a if take_a else label_b), True
    if policy == "error":
        raise TwigError(f"the two branches at the root carry different labels ({label_a!r} and {label_b!r})")
    if policy in (label_a, label_b):
        return policy, True
    return None, True


def plan_root_merge(tree, policy: str = "unlabel", outgroup_side: int | None = None) -> RootMerge | None:
    """Plan how a rooted tree is written as unrooted.

    Returns None if the root does not have exactly two children, or if both are
    tips. Otherwise the root child that is dissolved is a bifurcating internal
    child if there is one, so that PAML never gets more than three children at
    the root. The other child is kept, and gets the merged label and the summed
    branch length when the tree is written.

    `outgroup_side` is the idx of the root child on the side of the outgroup. If
    it is not given, the smaller side is taken as the outgroup side.
    """
    root = tree.treenode
    kids = list(root.children)
    if len(kids) != 2:
        return None
    internal = [k for k in kids if not k.is_leaf()]
    if not internal:
        return None
    dissolved = ([k for k in internal if len(k.children) == 2] or internal)[0]
    kept = kids[1] if kids[0] is dissolved else kids[0]

    a_is_outgroup: bool | None
    if outgroup_side is not None:
        a_is_outgroup = dissolved.idx == outgroup_side
    else:
        n_dissolved, n_kept = len(dissolved), len(kept)
        a_is_outgroup = None if n_dissolved == n_kept else n_dissolved < n_kept
    label_d, label_k = treeio.node_label(dissolved), treeio.node_label(kept)
    label, conflict = merge_root_labels(label_d, label_k, policy, a_is_outgroup)
    return RootMerge(dissolved.idx, kept.idx, tuple(treeio.node_label(k) for k in kids), label, conflict)


###########################################################################
# rerooting
###########################################################################

@dataclass(frozen=True)
class Rooted:
    """A tree rooted on an outgroup, and what happened while doing it.

    `outgroup_side` is the idx of the root child on the side of the outgroup,
    which the `outgroup` and `ingroup` root-edge policies need.
    """
    tree: object
    description: str
    notes: tuple[str, ...]
    outgroup_side: int | None


def tip_sets(tree) -> dict[int, frozenset[int]]:
    """Return the idx of the tips below every node."""
    below: dict[int, frozenset[int]] = {}
    for node in tree:
        if node.is_leaf():
            below[node.idx] = frozenset({node.idx})
        else:
            below[node.idx] = frozenset().union(*(below[c.idx] for c in node.children))
    return below


def tip_name_sets(tree) -> dict[int, frozenset[str]]:
    """Return the names of the tips below every node."""
    below: dict[int, frozenset[str]] = {}
    for node in tree:
        if node.is_leaf():
            below[node.idx] = frozenset({node.name})
        else:
            below[node.idx] = frozenset().union(*(below[c.idx] for c in node.children))
    return below


def parse_outgroups(values: Iterable[str]) -> list[tuple[str, ...]]:
    """Return ordered outgroup clades from values whose tips are joined by commas."""
    clades = []
    for value in values:
        tips = tuple(i.strip() for i in value.split(",") if i.strip())
        if not tips:
            raise TwigError(f"empty outgroup {value!r}")
        clades.append(tips)
    if not clades:
        raise TwigError("no outgroups were given")
    return clades


def _side_of(node, root):
    """Return the child of the root that is node or contains it."""
    while node.up is not root:
        node = node.up
    return node


def root_on(tree, target_idx: int, policy: str = "unlabel") -> tuple["object", list[str]]:
    """Return a copy of the tree rooted on the branch above a node, and notes.

    Labels stay on their branches. If the tree already had a root with two
    children, their branches become one, and the label of that branch is chosen
    by the policy before the tree is rerooted. Otherwise toytree would keep the
    label of one of them, chosen by the order of the children.
    """
    out = treeio.copy_tree(tree)
    target = out[target_idx]
    root = out.treenode
    kids = list(root.children)
    notes: list[str] = []
    if target.up is root and len(kids) == 2:
        return treeio.carry_over_dists(tree, treeio.ensure_labels(out)), notes   # already rooted on this split
    if len(kids) == 2:
        near = _side_of(target, root)
        far = kids[1] if kids[0] is near else kids[0]
        label, conflict = merge_root_labels(treeio.node_label(near), treeio.node_label(far), policy, True)
        if conflict:
            first, second = (i or "no label" for i in (treeio.node_label(near), treeio.node_label(far)))
            notes.append(
                f"the two branches at the old root had the labels {first} and {second}, "
                f"and the merged branch has {label or 'no label'}"
            )
        # toytree keeps the branch of one child and replaces the label of the other while rerooting
        near.label = far.label = label
    with treeio.deep_recursion(out):
        rooted = out.root(target.idx, edge_features=["label"])
    treeio.ensure_labels(rooted)
    return treeio.carry_over_dists(tree, rooted), notes


def reroot(
    tree,
    outgroups: Sequence[str],
    policy: str = "unlabel",
    species_of: Callable[[str], str] | None = None,
) -> Rooted:
    """Return the tree rooted on the first outgroup that fits.

    Each outgroup is a clade of tips joined by commas, and is used if some of its
    tips are in the tree and they are one side of a branch. Outgroups are tried in
    the order given. With `species_of` the outgroup names are species and not tip
    names, so one outgroup works for trees with several samples per species.
    """
    check_root_edge_policy(policy)
    clades = parse_outgroups(outgroups)
    key = species_of or (lambda name: name)
    tips = [leaf for leaf in tree if leaf.is_leaf()]
    below = tip_sets(tree)
    everything = frozenset(leaf.idx for leaf in tips)
    for clade in clades:
        wanted_keys = set(clade)
        present = frozenset(leaf.idx for leaf in tips if key(leaf.name) in wanted_keys)
        if not present or present == everything:
            continue
        for wanted in (present, everything - present):
            target = next((n for n in tree if n.up is not None and below[n.idx] == wanted), None)
            if target is None:
                continue
            rooted, notes = root_on(tree, target.idx, policy)
            names = {tree[i].name for i in present}
            side = None
            sets = tip_name_sets(rooted)
            for child in rooted.treenode.children:
                if names <= sets[child.idx]:
                    side = child.idx
            return Rooted(rooted, "rooted on " + ",".join(clade), tuple(notes), side)
    raise TwigError(
        "none of the outgroups (" + "; ".join(",".join(c) for c in clades) + ") can root the tree, because none "
        "has tips in the tree that form one side of a branch"
    )
