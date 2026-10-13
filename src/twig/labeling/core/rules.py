#!/usr/bin/env python

"""Rules that turn a selection of tips into labeled branches.

A rule takes the tips selected by a label set and returns the branches that get
the label of the set. The rules are registered by name, so a new rule is one
function with `register_rule`. HyPhy's names for its own `label-tree` rules are
accepted as aliases, but `Parsimony` means the real parsimony of `fitch` and not
HyPhy's rule of the same name, which behaves exactly like `All descendants`
because the unselected tips never get a state.

Branches are identified by the idx of the node below them. The root has no
branch and is never labeled.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Iterable
from dataclasses import dataclass, field
from fractions import Fraction

from twig.utils import TwigError
from twig.labeling.core.parsimony import sankoff_binary, to_fraction

AMBIGUOUS_MODES = ("deltran", "acctran", "unlabel")


@dataclass(frozen=True)
class RuleOptions:
    """Options of the parsimony rule.

    `ambiguous` says what to do with nodes that have both states in their
    state set, which happens when two histories are equally parsimonious. It is
    acctran, deltran, unlabel, or the name of a label for such nodes.
    """
    gain_cost: Fraction = Fraction(1)
    loss_cost: Fraction = Fraction(1)
    ambiguous: str = "acctran"
    root_state: str = "free"


@dataclass(frozen=True)
class RuleInput:
    """What a rule is given.

    `selected` holds the idx of the selected tips, `unknown` the tips whose
    state is not known (used by parsimony), and `groups` the tips chosen by each
    selector value separately (used by `clade` and `stem`).
    """
    tree: object
    label: str
    selected: frozenset[int]
    unknown: frozenset[int] = frozenset()
    groups: tuple[frozenset[int], ...] = ()
    options: RuleOptions = field(default_factory=RuleOptions)


@dataclass
class RuleOutput:
    """What a rule returns.

    `labels` maps branches to labels. `direct` is the set of branches labeled
    because their own tip was selected, as opposed to branches labeled because
    of the tips below them. `losses` and `ambiguous` are filled by the
    parsimony rule only.
    """
    labels: dict[int, str]
    direct: frozenset[int] = frozenset()
    losses: frozenset[int] = frozenset()
    ambiguous: frozenset[int] = frozenset()
    cost: Fraction | None = None


Rule = Callable[[RuleInput], RuleOutput]
RULES: dict[str, Rule] = {}
HELP: dict[str, str] = {}
ALIASES: dict[str, str] = {}


def normalize_rule_name(name: str) -> str:
    """Return a rule name in lower case, without commas and extra spaces."""
    return re.sub(r"\s+", " ", name.replace(",", " ").strip().lower())


def register_rule(name: str, *aliases: str, help: str = ""):
    """Register a rule function under a name and optional aliases."""
    def decorator(func: Rule) -> Rule:
        RULES[name] = func
        HELP[name] = help
        for alias in aliases:
            ALIASES[normalize_rule_name(alias)] = name
        return func
    return decorator


def resolve_rule(name: str) -> str:
    """Return the canonical rule name for a name or an alias."""
    key = normalize_rule_name(name)
    if key in RULES:
        return key
    if key in ALIASES:
        return ALIASES[key]
    raise TwigError(f"unknown rule {name!r}, choose from {list(RULES)}")


###########################################################################
# tree helpers
###########################################################################

def tip_counts(tree, selected: Iterable[int]) -> tuple[dict[int, int], dict[int, int]]:
    """Return the number of selected tips and of all tips below every node."""
    selected = set(selected)
    n_selected: dict[int, int] = {}
    n_tips: dict[int, int] = {}
    for node in tree:
        if node.is_leaf():
            n_tips[node.idx] = 1
            n_selected[node.idx] = int(node.idx in selected)
        else:
            n_tips[node.idx] = sum(n_tips[c.idx] for c in node.children)
            n_selected[node.idx] = sum(n_selected[c.idx] for c in node.children)
    return n_selected, n_tips


def mrca(tree, idxs: Iterable[int]):
    """Return the most recent common ancestor of nodes given by idx."""
    shared: list[int] | None = None
    for idx in idxs:
        chain = []
        node = tree[idx]
        while node is not None:
            chain.append(node.idx)
            node = node.up
        if shared is None:
            shared = chain
        else:
            keep = set(chain)
            shared = [i for i in shared if i in keep]
    if not shared:
        raise TwigError("cannot find the common ancestor of no tips")
    return tree[shared[0]]


def below(node) -> list:
    """Return a node and all the nodes below it."""
    return list(node.traverse("preorder"))


def is_root(node) -> bool:
    return node.up is None


###########################################################################
# the rules
###########################################################################

@register_rule("tips", "None", help="the selected tips only")
def rule_tips(inp: RuleInput) -> RuleOutput:
    labels = {idx: inp.label for idx in inp.selected}
    return RuleOutput(labels, direct=frozenset(labels))


@register_rule("all", "All descendants", help="every branch all of whose descendant tips are selected")
def rule_all(inp: RuleInput) -> RuleOutput:
    n_selected, n_tips = tip_counts(inp.tree, inp.selected)
    labels = {
        node.idx: inp.label for node in inp.tree
        if not is_root(node) and n_selected[node.idx] and n_selected[node.idx] == n_tips[node.idx]
    }
    return RuleOutput(labels, direct=frozenset(labels) & inp.selected)


@register_rule("crown", "All descendants, no MRCA", help="the all result without the branches that start a labeled clade")
def rule_crown(inp: RuleInput) -> RuleOutput:
    labeled = rule_all(inp).labels
    labels = {}
    for idx, label in labeled.items():
        node = inp.tree[idx]
        # tips stay labeled, an internal branch stays only below a labeled branch
        if node.is_leaf() or (not is_root(node.up) and node.up.idx in labeled):
            labels[idx] = label
    return RuleOutput(labels, direct=frozenset(labels) & inp.selected)


@register_rule("some", "Some descendants", help="every branch with at least one selected descendant tip")
def rule_some(inp: RuleInput) -> RuleOutput:
    n_selected, _ = tip_counts(inp.tree, inp.selected)
    labels = {node.idx: inp.label for node in inp.tree if not is_root(node) and n_selected[node.idx]}
    return RuleOutput(labels, direct=frozenset(labels) & inp.selected)


@register_rule("fitch", "Parsimony", help="the branches reconstructed in the selected state by two-state parsimony")
def rule_fitch(inp: RuleInput) -> RuleOutput:
    opts = inp.options
    tips = [node for node in inp.tree if node.is_leaf()]
    states = {
        node.idx: None if node.idx in inp.unknown else int(node.idx in inp.selected)
        for node in tips
    }
    explicit = opts.ambiguous in ("acctran", "deltran")
    result = sankoff_binary(
        inp.tree, states, opts.gain_cost, opts.loss_cost, opts.root_state,
        tie=opts.ambiguous if explicit else "acctran",
    )

    # a tip can be ambiguous too, if its state is unknown
    ambiguous = {
        node.idx for node in inp.tree
        if not is_root(node) and len(result.state_sets[node.idx]) == 2
    }
    labels: dict[int, str] = {}
    losses = set()
    for node in inp.tree:
        if is_root(node):
            continue
        idx, parent = node.idx, node.up
        if explicit:
            if result.states[idx] == 1:
                labels[idx] = inp.label
            if result.states[idx] == 0 and result.states[parent.idx] == 1:
                losses.add(idx)
        else:
            sets = result.state_sets
            if sets[idx] == {1}:
                labels[idx] = inp.label
            elif idx in ambiguous and opts.ambiguous != "unlabel":
                labels[idx] = opts.ambiguous
            if sets[idx] == {0} and sets[parent.idx] == {1}:
                losses.add(idx)
    return RuleOutput(
        labels,
        direct=frozenset(i for i in labels if i in inp.selected),
        losses=frozenset(losses),
        ambiguous=frozenset(ambiguous),
        cost=result.cost,
    )


def _group_nodes(inp: RuleInput):
    """Return the common ancestor of each selector's tips, which must not be the root."""
    if not inp.groups:
        raise TwigError(f"rule needs selectors for the label {inp.label!r}")
    out = []
    for number, group in enumerate(inp.groups, start=1):
        if not group:
            raise TwigError(f"selector {number} of the label {inp.label!r} matches no tips")
        node = mrca(inp.tree, group)
        if is_root(node):
            raise TwigError(
                f"the tips of selector {number} of the label {inp.label!r} have the root as their common "
                "ancestor, and the root has no branch to label"
            )
        out.append(node)
    return out


@register_rule("clade", help="for each selector, the branch above the common ancestor of its tips and everything below it")
def rule_clade(inp: RuleInput) -> RuleOutput:
    labels = {}
    for node in _group_nodes(inp):
        for member in below(node):
            labels[member.idx] = inp.label
    return RuleOutput(labels, direct=frozenset(labels) & inp.selected)


@register_rule("stem", help="for each selector, only the branch above the common ancestor of its tips")
def rule_stem(inp: RuleInput) -> RuleOutput:
    labels = {node.idx: inp.label for node in _group_nodes(inp)}
    return RuleOutput(labels, direct=frozenset(labels) & inp.selected)


def make_options(gain=1, loss=1, ambiguous: str = "acctran", root_state: str = "free") -> RuleOptions:
    """Return RuleOptions from values as they are given on the command line."""
    return RuleOptions(to_fraction(gain), to_fraction(loss), ambiguous, str(root_state))
