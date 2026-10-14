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
from twig.labeling.core.constants import HAS_DIST, LABEL, ROOT_EDGE_POLICIES
from twig.labeling.core.model import RootMerge
from twig.labeling.core.sanitize import label_name_problem


def _one_edit_apart(a: str, b: str) -> bool:
    """Return True if one inserted, removed or changed letter, or two swapped neighbors, turn a into b."""
    if a == b or abs(len(a) - len(b)) > 1:
        return False
    if len(a) > len(b):
        a, b = b, a
    start = next((i for i, (x, y) in enumerate(zip(a, b)) if x != y), len(a))
    if len(a) < len(b):
        return a[start:] == b[start + 1:]
    if a[start + 1:] == b[start + 1:]:
        return True
    return a[start:start + 2] == b[start:start + 2][::-1] and a[start + 2:] == b[start + 2:]


def check_root_edge_policy(policy: str) -> str:
    """Return the policy if it is a keyword or a valid label name, else raise.

    A label name that is a keyword written in another case, or a keyword with
    one letter wrong, is refused, because as a label that no branch has it
    would silently give the merged branch no label.
    """
    if policy in ROOT_EDGE_POLICIES:
        return policy
    problem = label_name_problem(policy)
    for keyword in ROOT_EDGE_POLICIES:
        if policy.lower() == keyword:
            raise TwigError(
                f"the root edge policy {policy!r} would be read as a label name, but it looks like {keyword!r}. "
                f"Policies are written in lower case, so give {keyword} or choose another label name"
            )
        if not problem and _one_edit_apart(policy.lower(), keyword):
            raise TwigError(
                f"the root edge policy {policy!r} would be read as a label name, but it looks like a misspelling "
                f"of {keyword!r}. Give {keyword} or choose another label name"
            )
    if problem:
        raise TwigError(
            f"unknown root edge policy {policy!r}, use one of {list(ROOT_EDGE_POLICIES)} or the name of a label"
        )
    return policy


def missing_policy_label(policy: str, labels: Iterable[str]) -> str | None:
    """Return a warning if the policy names a label that no branch has, else None.

    Such a policy gives no label to a merged branch whose two halves differ, as
    unlabel does, which is most likely not what was meant.
    """
    labels = set(labels)
    if policy in ROOT_EDGE_POLICIES or policy in labels or not labels:
        return None
    message = (
        f"--root-edges names the label {policy}, which no branch has, so a merged root branch whose two "
        "halves differ gets no label, as with unlabel. Check the spelling"
    )
    alike = sorted(i for i in labels if i.lower() == policy.lower())
    if alike:
        message += f", or give {alike[0]}, which differs only in case"
    return message


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


def _tips_below(tree) -> dict[int, frozenset[int]]:
    """Return the ids of the tip objects below every node.

    Unlike idx these stay the same when a tree is rerooted in place, and unlike
    names they differ for tips that share a name.
    """
    below: dict[int, frozenset[int]] = {}
    for node in tree:
        if node.is_leaf():
            below[node.idx] = frozenset({id(node)})
        else:
            below[node.idx] = frozenset().union(*(below[c.idx] for c in node.children))
    return below


def _named_branches(tree) -> list[tuple[frozenset[int], str]]:
    """Return (tips below, name) for every internal branch that has a name."""
    named = [n for n in tree if not n.is_leaf() and n.up is not None and n.name]
    if not named:
        return []
    below = _tips_below(tree)
    return [(below[n.idx], n.name) for n in named]


def _restore_names(tree, named: list[tuple[frozenset[int], str]]) -> list[str]:
    """Give each internal name back to the branch it named, and return notes on names that are dropped.

    toytree leaves names on their nodes when it reroots, so on the path to the
    new root every name would move to the next branch. A branch is found by the
    tips on one side of it. The two old root branches are now one, so a name of
    one of them is dropped if the other is a tip or has a name of its own.
    """
    if not named:
        return []
    below = _tips_below(tree)
    everything = below[tree.treenode.idx]
    by_tips = {below[n.idx]: n for n in tree if n.up is not None}
    for node in tree:
        if not node.is_leaf() and node.up is not None:
            node.name = ""
    # a branch that still has the same tips below it takes its name first
    late = []
    for tips, name in named:
        node = by_tips.get(tips)
        if node is not None and not node.is_leaf():
            node.name = name
        else:
            late.append((tips, name))
    notes = []
    for tips, name in late:
        node = by_tips.get(everything - tips)
        if node is not None and not node.is_leaf() and not node.name:
            node.name = name
            continue
        where = f"the branch of {node.name}" if node is not None and node.name else "another branch"
        notes.append(f"the internal name {name} is dropped, because its branch at the old root was joined with {where}")
    return notes


def root_on(tree, target_idx: int, policy: str = "unlabel") -> tuple["object", list[str]]:
    """Return a copy of the tree rooted on the branch above a node, and notes.

    Labels and internal names stay on their branches. If the tree already had a
    root with two children, their branches become one, and the label of that
    branch is chosen by the policy before the tree is rerooted. Otherwise
    toytree would keep the label of one of them, chosen by the order of the
    children. The merged branch has the summed length, and no length if either
    of the two had none in the input.
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
        # the merged branch has a length only if both of its halves had one
        if not (treeio.has_dist(near) and treeio.has_dist(far)):
            setattr(near, HAS_DIST, False)
            setattr(far, HAS_DIST, False)
    named = _named_branches(out)
    with treeio.deep_recursion(out):
        rooted = out.root(target.idx, edge_features=[LABEL, HAS_DIST], inplace=True)
    treeio.ensure_labels(rooted)
    notes += _restore_names(rooted, named)
    return treeio.carry_over_dists(tree, rooted), notes


def reroot(
    tree,
    outgroups: Sequence[str],
    policy: str = "unlabel",
    species_of: Callable[[str], str | None] | None = None,
) -> Rooted:
    """Return the tree rooted on the first outgroup that fits.

    Each outgroup is a clade of tips joined by commas, and is used if some of its
    tips are in the tree and they are one side of a branch. Outgroups are tried in
    the order given. With `species_of` the outgroup names are species and not tip
    names, so one outgroup works for trees with several samples per species, and
    a tip whose species is None belongs to no outgroup.
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
