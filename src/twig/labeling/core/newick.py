#!/usr/bin/env python

"""Write labeled Newick strings with one iterative serializer.

All writers share `to_newick`, so the order of tokens is always
`name{Label}:len` for tips and `)name{Label}:len` for internal nodes. toytree's
own writer is not used for labeled output because it cannot place `{...}` text
on some nodes only, and the serializer can also write a rooted tree as
unrooted in one pass, as HyPhy and PAML read it.

This module does not import toytree.
"""

from __future__ import annotations

import re
from collections.abc import Callable

from twig.utils import TwigError
from twig.labeling.core.constants import LABEL
from twig.labeling.core.model import RootMerge
from twig.labeling.core.treeio import has_dist

# characters that cannot appear in an unquoted Newick name
_BAD_NAME = re.compile(r"[()\[\]{},:;'\"\s]")


def braces_label(label: str | None) -> str:
    """Return the HyPhy form of a label, `{Test}`, or '' if there is no label."""
    if not label:
        return ""
    bad = _BAD_NAME.search(label)
    if bad:
        raise TwigError(f"cannot write label {label!r}, it contains {bad.group()!r}")
    return "{%s}" % label


def _default_tip_text(node) -> str:
    return node.name


def _default_internal_text(node) -> str:
    return ""


def _check_text(text: str, node) -> str:
    bad = _BAD_NAME.search(text)
    if bad:
        raise TwigError(f"cannot write name {text!r} of node {node.idx}, it contains {bad.group()!r}")
    return text


def to_newick(
    tree,
    *,
    tip_text: Callable = _default_tip_text,
    internal_text: Callable = _default_internal_text,
    label_text: Callable = braces_label,
    dists: bool = True,
    dist_fmt: str = "%.12g",
    root_merge: RootMerge | None = None,
) -> str:
    """Return a Newick string for a tree whose nodes have a `label` feature.

    Parameters
    ----------
    tree: ToyTree
        Nodes may have a `label` attribute, which is written after the name.
        The root has no branch, so its label and length are never written.
    tip_text: Callable
        Function of a node that returns the text written for a tip.
    internal_text: Callable
        Function of a node that returns the text written after the `)` of an
        internal node (a name or support value), or ''.
    label_text: Callable
        Function of a label (or None) that returns the text written for it,
        such as `{Test}`. The default is the HyPhy form.
    dists: bool
        Write branch lengths. Trees read without lengths should be written
        without them, since toytree fills in 1.0 for missing lengths, and a
        branch whose length was missing in the input is written without one.
    dist_fmt: str
        Format string for branch lengths.
    root_merge: RootMerge
        Write the tree as unrooted. The `dissolved` child of the root is
        removed and its children move up to the root, and the `kept` child is
        written with the merged label and the summed branch length, or without
        a length if either of the two had none.
    """
    root = tree.treenode
    dissolved = kept = None
    if root_merge is not None:
        dissolved, kept = tree[root_merge.dissolved], tree[root_merge.kept]
        if dissolved.up is not root or kept.up is not root or len(root.children) != 2:
            raise TwigError("a root merge needs a root with exactly two children")
        if dissolved.is_leaf():
            raise TwigError("cannot dissolve a root child that is a tip")

    def children_of(node):
        if node is root and dissolved is not None:
            out = []
            for child in node.children:
                out.extend(dissolved.children if child is dissolved else (child,))
            return out
        return node.children

    def suffix(node) -> str:
        if node is root:
            return ""
        label = root_merge.label if node is kept else getattr(node, LABEL, None)
        text = label_text(label)
        if dists and has_dist(node) and (node is not kept or has_dist(dissolved)):
            dist = node.dist + (dissolved.dist if node is kept else 0.0)
            if dist != dist:
                raise TwigError(f"branch length of node {node.idx} is NaN")
            text += ":" + dist_fmt % dist
        return text

    out: list[str] = []
    stack: list[tuple[str, object]] = [("node", root)]
    while stack:
        kind, item = stack.pop()
        if kind == "text":
            out.append(item)
            continue
        node = item
        if kind == "close":
            out.append(")" + _check_text(internal_text(node), node) + suffix(node))
            continue
        kids = children_of(node)
        if not kids:
            out.append(_check_text(tip_text(node), node) + suffix(node))
            continue
        out.append("(")
        stack.append(("close", node))
        for idx in range(len(kids) - 1, -1, -1):
            stack.append(("node", kids[idx]))
            if idx:
                stack.append(("text", ","))
    return "".join(out) + ";"
