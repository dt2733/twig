#!/usr/bin/env python

"""Read trees with branch labels into toytree.

Labels are rewritten at the string level before toytree parses anything. This
avoids several problems in toytree itself. Its distance parser raises on
`A:0.1#1`, `{...}` ends up in node names (or in a `trait` feature in toytree
master), and a single `95{Test}` internal name stops toytree from recognizing
the other internal names as support values.

Labels are stored in the node feature `label` (None means unlabeled) on the
node below the branch, and `label` is registered as an edge feature so that it
follows its branch when a tree is rerooted. The root has no branch and never
carries a label.

Recognized label syntaxes are HyPhy braces `A{Test}`, PAML marks `A #1` and
`$1`, and the `label` key of NHX `[&&NHX:label=Test]` and FigTree
`[&label="Test"]` annotations.

Only functions that build trees import toytree, so scanning names with
`tip_names` is fast and does not need toytree.
"""

from __future__ import annotations

import re
import sys
from collections import Counter
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from loguru import logger

from twig.utils import TwigError
from twig.labeling.core.constants import LABEL, LABEL_SYNTAXES

if TYPE_CHECKING:  # pragma: no cover
    from toytree import ToyTree

# node features used to carry labels through toytree's parser
_TMP_LABEL = "twiglab"
_TMP_PAML = "twigpaml"
_TMP_CLADE = "twigpamlc"

_TOKEN = re.compile(
    r"(?P<comment>\[[^\]]*\])"
    r"|(?P<brace>\{[^}]*\})"
    r"|(?P<mark>[#$]\d+)"
    r"|(?P<punct>[(),;:])"
    r"|(?P<space>\s+)"
    r"|(?P<text>[^()\[\]{},:;#$\s]+)"
)
_PAML_HEADER = re.compile(r"^\s*\d+\s+\d+\s*$")
_LOCUS_PREFIX = re.compile(r"^([^\s(\[;]+)\t+(\(.*)$")
_NUMBER = r"[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?"
_IS_NUMBER = re.compile(rf"^{_NUMBER}$")
_IS_DOUBLE_SUPPORT = re.compile(rf"^({_NUMBER})/({_NUMBER})$")

_toytree = None


def import_toytree():
    """Import toytree, silencing its log messages, and return the module."""
    global _toytree
    if _toytree is None:
        import toytree
        logger.disable("toytree")
        _toytree = toytree
    return _toytree


def toytree_errors() -> tuple[type[Exception], ...]:
    """Return the exception types raised by toytree that callers should catch."""
    import_toytree()
    from toytree.utils import NodeDataError, ToytreeError, ToytreeRegexError
    return (ToytreeError, NodeDataError, ToytreeRegexError)


###########################################################################
# splitting files into trees
###########################################################################

def _read_text(path: str | Path) -> str:
    path = Path(path)
    try:
        return path.read_text(encoding="utf-8")
    except OSError as exc:
        raise TwigError(f"cannot read tree file {path}: {exc}") from exc
    except UnicodeDecodeError as exc:
        raise TwigError(f"tree file {path} is not valid UTF-8 text") from exc


def _compact(text: str, source: str, number: int) -> str:
    """Remove whitespace outside [...] comments, and refuse quoted names."""
    out = []
    depth = 0
    for ch in text:
        if ch == "[":
            depth += 1
        elif ch == "]":
            depth = max(0, depth - 1)
        if depth:
            out.append(" " if ch in "\n\r\t" else ch)
        elif not ch.isspace():
            if ch in "'\"":
                raise TwigError(
                    f"{source}: tree {number} has quoted labels, which are not supported; remove the quotes"
                )
            out.append(ch)
    return "".join(out)


def iter_trees_text(text: str, source: str = "<text>") -> Iterator[tuple[str | None, str]]:
    """Yield (locus, newick) for each tree in text.

    A tree ends at the first `;` outside [...] comments, so trees can span
    lines and several trees can share a line. Blank lines, `#` comment lines
    between trees and a PAML `ns ntree` header are skipped, and a line can have
    the form `locus<TAB>newick`. The locus is None otherwise. NEXUS input is not
    supported.
    """
    if text.lstrip()[:6].upper() == "#NEXUS":
        raise TwigError(f"{source}: NEXUS format is not supported, convert the trees to plain Newick")
    pieces: list[str] = []
    in_tree = False     # the tree in progress has some text
    locus: str | None = None
    depth = 0
    number = 0
    for line in text.replace("\r\n", "\n").replace("\r", "\n").split("\n"):
        stripped = line.strip()
        if not in_tree:
            if not stripped or stripped.startswith("#"):
                continue
            if number == 0 and _PAML_HEADER.match(stripped):
                continue
            found = _LOCUS_PREFIX.match(line.lstrip())
            if found:
                locus, stripped = found.group(1), found.group(2).strip()
        start = 0
        for pos, ch in enumerate(stripped):
            if ch == "[":
                depth += 1
            elif ch == "]":
                depth = max(0, depth - 1)
            elif ch == ";" and depth == 0:
                pieces.append(stripped[start:pos + 1])
                number += 1
                yield locus, _compact("".join(pieces), source, number)
                pieces, in_tree, locus, start = [], False, None, pos + 1
        rest = stripped[start:]
        if rest:
            pieces.append(rest)
            in_tree = True
    # only comments may follow the last tree
    if in_tree and re.sub(r"\[[^\]]*\]", "", "".join(pieces)).strip():
        raise TwigError(f"{source}: the last tree has no final ';'")


def iter_trees(path: str | Path) -> Iterator[tuple[str | None, str]]:
    """Yield (locus, newick) for each tree in a file. See `iter_trees_text`."""
    yield from iter_trees_text(_read_text(path), source=str(path))


###########################################################################
# tokens
###########################################################################

def _tokens(newick: str) -> Iterator[tuple[str, str, int]]:
    """Yield (kind, text, position) for each token, skipping whitespace."""
    pos, size = 0, len(newick)
    while pos < size:
        found = _TOKEN.match(newick, pos)
        if found is None:
            raise TwigError(f"unexpected character {newick[pos]!r} at position {pos} of the tree")
        kind = found.lastgroup
        if kind != "space":
            yield kind, found.group(), pos
        pos = found.end()


def tip_names(newick: str) -> list[str]:
    """Return the tip names of a Newick string, in the order written.

    This scans the string without toytree. Labels, comments, support values and
    branch lengths are skipped. A tip without a name is an error.
    """
    names: list[str] = []
    depth = 0
    pending = True      # a tip name may start here
    after_colon = False
    for kind, token, pos in _tokens(newick):
        if kind == "punct":
            if token == "(":
                depth += 1
                pending = True
            elif token == ",":
                if pending:
                    raise TwigError(f"the tree has a tip without a name before position {pos}")
                pending = True
            elif token == ")":
                if pending:
                    raise TwigError(f"the tree has a tip without a name before position {pos}")
                depth -= 1
                if depth < 0:
                    raise TwigError(f"unbalanced parentheses at position {pos} of the tree")
            elif token == ";":
                break
            after_colon = token == ":"
        elif kind == "text":
            if "'" in token or '"' in token:
                raise TwigError("quoted labels are not supported; remove the quotes")
            if after_colon:
                after_colon = False
            elif pending:
                names.append(token)
                pending = False
    if depth != 0:
        raise TwigError("unbalanced parentheses in the tree")
    if pending and names:
        raise TwigError("the tree has a tip without a name")
    return names


###########################################################################
# rewriting label syntaxes
###########################################################################

@dataclass
class _Rewrite:
    text: str
    labels: list[str]
    syntaxes: set[str]
    has_dists: bool


def _split_pairs(body: str, delim: str) -> list[str]:
    """Split body on delim, but not inside {} or double quotes."""
    pairs, current, depth, quoted = [], [], 0, False
    for ch in body:
        if ch == '"':
            quoted = not quoted
        elif not quoted and ch == "{":
            depth += 1
        elif not quoted and ch == "}":
            depth = max(0, depth - 1)
        if ch == delim and depth == 0 and not quoted:
            pairs.append("".join(current))
            current = []
        else:
            current.append(ch)
    pairs.append("".join(current))
    return [i for i in pairs if i]


def _comment_pairs(comment: str) -> list[tuple[str, str]] | None:
    """Return (key, value) pairs of a [&...] or [&&NHX:...] annotation.

    Plain comments are not annotations and give None. Items without `=` are
    dropped.
    """
    body = comment[1:-1]
    if body.startswith("&&NHX:"):
        items = _split_pairs(body[len("&&NHX:"):], ":")
    elif body.startswith("&"):
        items = _split_pairs(body.lstrip("&"), ",")
    else:
        return None
    return [tuple(i.split("=", 1)) for i in items if "=" in i]


def _rewrite(newick: str, label_syntax: str) -> _Rewrite:
    """Replace label syntaxes by [&twiglab=N] style annotations toytree can parse.

    Labels are stored in a table and referred to by index, so that any label
    text is safe and toytree does not convert it to a number.
    """
    modes = {
        "auto": ("braces", "paml", "nhx"),
        "braces": ("braces",),
        "paml": ("paml",),
        "nhx": ("nhx",),
        "none": (),
    }[label_syntax]
    table: list[str] = []
    index: dict[str, int] = {}
    syntaxes: set[str] = set()
    out: list[str] = []
    in_node_part = 0     # annotations seen on the current node or edge
    has_dists = False

    def slot(label: str) -> int:
        if label not in index:
            index[label] = len(table)
            table.append(label)
        return index[label]

    def annotate(pos: int) -> None:
        nonlocal in_node_part
        in_node_part += 1
        if in_node_part > 1:
            raise TwigError(
                f"the node near position {pos} of the tree has more than one annotation "
                "and twig can read only one per node or edge"
            )

    def refuse(syntax: str, token: str, pos: int) -> None:
        raise TwigError(
            f"found {token!r} at position {pos} of the tree, but the label syntax is {label_syntax!r} "
            f"and this is {syntax} syntax"
        )

    for kind, token, pos in _tokens(newick):
        if kind == "punct":
            in_node_part = 0
            has_dists = has_dists or token == ":"
            out.append(token)
        elif kind == "text":
            if "'" in token or '"' in token:
                raise TwigError("quoted labels are not supported; remove the quotes")
            out.append(token)
        elif kind == "brace":
            if label_syntax == "none":
                continue
            if "braces" not in modes:
                refuse("braces", token, pos)
            label = token[1:-1].strip()
            if not label:
                raise TwigError(f"empty label {token!r} at position {pos} of the tree")
            annotate(pos)
            syntaxes.add("braces")
            out.append(f"[&{_TMP_LABEL}={slot(label)}]")
        elif kind == "mark":
            if label_syntax == "none":
                continue
            if "paml" not in modes:
                refuse("PAML", token, pos)
            annotate(pos)
            syntaxes.add("paml")
            key = _TMP_PAML if token[0] == "#" else _TMP_CLADE
            out.append(f"[&{key}={token[1:]}]")
        else:   # comment
            pairs = _comment_pairs(token)
            if not pairs:
                continue
            kept = []
            for key, value in pairs:
                if key.lower() == "label":
                    if label_syntax == "none":
                        continue
                    if "nhx" not in modes:
                        refuse("NHX or FigTree", token, pos)
                    value = value.strip().strip('"')
                    if not value:
                        raise TwigError(f"empty label in {token!r} at position {pos} of the tree")
                    syntaxes.add("nhx")
                    kept.append(f"{_TMP_LABEL}={slot(value)}")
                else:
                    kept.append(f"{key}={value}")
            if kept:
                annotate(pos)
                out.append("[&" + ",".join(kept) + "]")

    if len(syntaxes) > 1:
        raise TwigError(
            f"the tree mixes label syntaxes ({', '.join(sorted(syntaxes))}); "
            "set the label syntax to read only one of them"
        )
    if not newick.rstrip().endswith(";"):
        out.append(";")
    return _Rewrite("".join(out), table, syntaxes, has_dists)


###########################################################################
# reading trees
###########################################################################

@dataclass
class ParsedTree:
    """A tree read from Newick, and what was found while reading it.

    `syntax` is the label syntax found (braces, paml, nhx) or None, and
    `root_label` is a label that was found on the root, which has no branch
    and so is not kept on the tree.
    """
    tree: "ToyTree"
    syntax: str | None
    has_dists: bool
    root_label: str | None


def _parse_with_toytree(text: str, where: str) -> "ToyTree":
    toytree = import_toytree()
    errors = toytree_errors() + (ValueError, IndexError, KeyError)
    limit = sys.getrecursionlimit()
    try:
        sys.setrecursionlimit(max(limit, 6000))
        return toytree.tree(text, internal_labels="name")
    except RecursionError as exc:
        raise TwigError(f"{where}: the tree is too deep to parse") from exc
    except errors as exc:
        raise TwigError(f"{where}: cannot parse the tree ({exc})") from exc
    finally:
        sys.setrecursionlimit(limit)


def _is_number(text: str) -> bool:
    return bool(_IS_NUMBER.match(text))


def _move_supports(tree: "ToyTree") -> None:
    """Move numeric internal names to `support`.

    This is done only if every internal name is a number (or empty), as toytree
    does, so that names such as Node1 are kept as names. IQ-TREE names such as
    `80.5/95` (SH-aLRT and UFBoot) become `support_alrt` and `support`.
    """
    internal = [n for n in tree if not n.is_leaf()]
    named = [n for n in internal if n.name and not n.is_root()]
    if not named:
        return
    if all(_is_number(n.name) for n in named):
        for node in internal:
            if node.name and _is_number(node.name):
                node.support = float(node.name)
                node.name = ""
    elif all(_IS_DOUBLE_SUPPORT.match(n.name) for n in named):
        for node in tree:
            node.support_alrt = float("nan")
        for node in internal:
            found = _IS_DOUBLE_SUPPORT.match(node.name)
            if found:
                node.support_alrt = float(found.group(1))
                node.support = float(found.group(2))
                node.name = ""
        tree.edge_features.add("support_alrt")


def _paml_labels(tree: "ToyTree", codes: Mapping[int | str, str] | None) -> None:
    """Set `label` from PAML `#` and `$` marks, as codeml does.

    A `$` mark labels its clade, meaning the branch above the node and all
    branches below it. An explicit `#` beats `$`, and an inner `$` beats an
    outer one.
    """
    names = {str(k): v for k, v in (codes or {}).items()}
    clade_of: dict[int, int | None] = {}
    for node in tree.traverse("preorder"):
        own_clade = node.__dict__.pop(_TMP_CLADE, None)
        own_mark = node.__dict__.pop(_TMP_PAML, None)
        inherited = clade_of[node.up.idx] if node.up is not None else None
        clade = int(own_clade) if own_clade is not None else inherited
        clade_of[node.idx] = clade
        mark = int(own_mark) if own_mark is not None else clade
        node.label = None if mark is None else names.get(str(mark), str(mark))


def parse_newick(
    newick: str,
    label_syntax: str = "auto",
    paml_codes: Mapping[int | str, str] | None = None,
    source: str = "tree",
) -> ParsedTree:
    """Return a ParsedTree from one Newick string.

    See `read_tree` for the parameters.
    """
    if label_syntax not in LABEL_SYNTAXES:
        raise TwigError(f"unknown label syntax {label_syntax!r}, choose from {list(LABEL_SYNTAXES)}")
    try:
        rewritten = _rewrite(newick, label_syntax)
    except TwigError as exc:
        raise TwigError(f"{source}: {exc}") from None
    tree = _parse_with_toytree(rewritten.text, source)

    # labels stored by the rewrite -> the `label` feature
    syntax = next(iter(rewritten.syntaxes), None)
    if syntax == "paml":
        _paml_labels(tree, paml_codes)
    else:
        for node in tree:
            slot = node.__dict__.pop(_TMP_LABEL, None)
            node.label = None if slot is None else rewritten.labels[int(slot)]
    root_label = tree.treenode.label
    tree.treenode.label = None

    _move_supports(tree)
    for feature in (_TMP_LABEL, _TMP_PAML, _TMP_CLADE):
        tree.edge_features.discard(feature)
    tree.edge_features.add(LABEL)
    tree._twig_has_dists = rewritten.has_dists
    return ParsedTree(tree, syntax, rewritten.has_dists, root_label)


def _single_newick(source: str | Path) -> tuple[str, str]:
    """Return (newick, description) for one tree given as a string or a file."""
    if isinstance(source, Path) or (
        isinstance(source, str) and ";" not in source and not source.lstrip().startswith("(")
    ):
        path = Path(source)
        if not path.is_file():
            raise TwigError(f"tree file not found: {source}")
        records = list(iter_trees(path))
        where = str(path)
    else:
        records = list(iter_trees_text(str(source)))
        where = "tree"
    if not records:
        raise TwigError(f"{where}: no tree found")
    if len(records) > 1:
        raise TwigError(f"{where}: found {len(records)} trees but expected one")
    return records[0][1], where


def read_tree(
    source: str | Path,
    label_syntax: str = "auto",
    paml_codes: Mapping[int | str, str] | None = None,
) -> "ToyTree":
    """Return a ToyTree read from a Newick string or a tree file.

    Every node has a `label` feature (None if unlabeled) that is registered as
    an edge feature, a label on the root is dropped, numeric internal names are
    moved to `support`, and `has_dists(tree)` says whether the input had branch
    lengths.

    Parameters
    ----------
    source: str or Path
        A Newick string, or the path of a file with a single tree.
    label_syntax: str
        One of auto, braces, paml, nhx or none. With auto, any one of the
        syntaxes is read and a tree that mixes syntaxes is an error. With none,
        labels are removed.
    paml_codes: Mapping
        Names for PAML label numbers, for example {1: "Test", 2: "Reference"}.
        Numbers without a name are used as the label text.
    """
    newick, where = _single_newick(source)
    return parse_newick(newick, label_syntax, paml_codes, where).tree


###########################################################################
# labels on trees
###########################################################################

def node_label(node) -> str | None:
    """Return the label of the branch above a node, or None."""
    return getattr(node, LABEL, None)


def has_dists(tree: "ToyTree") -> bool:
    """Return True unless the tree was read from Newick without branch lengths."""
    return getattr(tree, "_twig_has_dists", True)


def carry_over_dists(old: "ToyTree", new: "ToyTree") -> "ToyTree":
    """Copy the has-lengths flag from old to new (for trees built from old)."""
    new._twig_has_dists = has_dists(old)
    return new


def strip_labels(tree: "ToyTree", inplace: bool = False) -> "ToyTree":
    """Return a tree in which every branch is unlabeled."""
    out = tree if inplace else tree.copy()
    for node in out:
        node.label = None
    out.edge_features.add(LABEL)
    return out


def labeled_branches(tree: "ToyTree") -> dict[int, str]:
    """Return {idx: label} for every labeled node."""
    return {n.idx: n.label for n in tree if node_label(n)}


###########################################################################
# species trees
###########################################################################

def _check_unique(names: list[str], where: str) -> None:
    dups = sorted(n for n, k in Counter(names).items() if k > 1)
    if dups:
        raise TwigError(f"{where}: tip names are not unique ({', '.join(dups[:10])})")


def species_names(source: str | Path) -> list[str]:
    """Return the unique tip names of a species tree, without using toytree."""
    newick, where = _single_newick(source)
    try:
        names = tip_names(newick)
    except TwigError as exc:
        raise TwigError(f"{where}: {exc}") from None
    if not names:
        raise TwigError(f"{where}: the species tree has no tips")
    _check_unique(names, where)
    return names


def read_species_tree(source: str | Path, label_syntax: str = "auto") -> "ToyTree":
    """Return a species tree, possibly labeled, after checking its tip names.

    Tip names must be unique, and a label found on the root is ignored with a
    warning because the root has no branch.
    """
    newick, where = _single_newick(source)
    parsed = parse_newick(newick, label_syntax, source=where)
    _check_unique([n.name for n in parsed.tree if n.is_leaf()], where)
    if parsed.root_label:
        logger.warning(f"{where}: ignoring the label {parsed.root_label!r} on the root, which has no branch")
    return parsed.tree


def branch_ids(tree: "ToyTree") -> dict[int, str]:
    """Return {idx: id} for every non-root node of a species tree.

    The id of an internal branch is the node's name if it has one, otherwise
    `S` followed by its idx. Tips are identified by their names.
    """
    ids: dict[int, str] = {}
    for node in tree:
        if node.is_root():
            continue
        ids[node.idx] = node.name if node.name else f"S{node.idx}"
    _check_unique(list(ids.values()), "species tree branch ids")
    return ids
