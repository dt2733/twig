#!/usr/bin/env python

"""Two-state parsimony by the Sankoff algorithm.

The trait has two states, 0 for absent and 1 for present, and a change costs
`gain` from 0 to 1 and `loss` from 1 to 0. Costs are exact fractions, so ties
between reconstructions are detected without rounding errors, and polytomies
are handled exactly.

The algorithm has three parts. A downpass in idx order, children first, gives
the cost of each subtree for each state of its root. An uppass gives the cost
of everything outside each subtree, so that the sum of the two says which
states occur in some most parsimonious reconstruction (the state set of a
node). A preorder pass then picks one reconstruction. At a node where both
states are equally good given the parent, `deltran` delays the change toward
the tips by keeping the state of the parent, and `acctran` moves the change
toward the root by changing at once.

toytree's own ancestral state code is not used, because in the versions
supported it only runs a downpass.

This module does not import toytree.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction

from twig.utils import TwigError

TIES = ("deltran", "acctran")


@dataclass
class ParsimonyResult:
    """The result of a parsimony reconstruction.

    `states` is one most parsimonious reconstruction as {idx: state}, chosen
    with the tie rule. `state_sets` gives for each idx the states that occur in
    at least one most parsimonious reconstruction, so a node with both states
    is ambiguous. `cost` is the cost of the reconstruction.
    """
    states: dict[int, int]
    state_sets: dict[int, frozenset[int]]
    cost: Fraction

    def ambiguous(self) -> set[int]:
        """Return the idx of nodes that have both states in their state set."""
        return {idx for idx, states in self.state_sets.items() if len(states) == 2}


def to_fraction(value) -> Fraction:
    """Return a cost as an exact fraction, reading floats as decimal text."""
    if isinstance(value, Fraction):
        return value
    try:
        out = Fraction(str(value))
    except (ValueError, ZeroDivisionError):
        raise TwigError(f"invalid change cost {value!r}") from None
    if out < 0:
        raise TwigError(f"change costs cannot be negative, got {value!r}")
    return out


def sankoff_binary(
    tree,
    tip_states: Mapping[int, int | None],
    gain=1,
    loss=1,
    root_state: str | int = "free",
    tie: str = "acctran",
) -> ParsimonyResult:
    """Reconstruct a binary trait on a tree with parsimony.

    Parameters
    ----------
    tree: ToyTree
        Its nodes are visited in idx order, which puts children before parents.
    tip_states: Mapping[int, int or None]
        The state (0 or 1) of every tip by idx. None means unknown, so the tip
        costs nothing in either state.
    gain, loss:
        The cost of a change from 0 to 1 and from 1 to 0.
    root_state: "free", 0 or 1
        Fix the state of the root, or let it be chosen. If both states of a
        free root are equally good, the root is given state 0.
    tie: "deltran" or "acctran"
        How to pick between equally good states of a node, given its parent.
        `deltran` keeps the state of the parent, which delays changes toward
        the tips and favors separate gains over a gain followed by a loss.
        `acctran` changes at once, which moves changes toward the root and
        favors a single gain followed by a loss.
    """
    if tie not in TIES:
        raise TwigError(f"unknown tie rule {tie!r}, choose from {list(TIES)}")
    if str(root_state) not in ("free", "0", "1"):
        raise TwigError(f"root state must be free, 0 or 1, got {root_state!r}")
    cost = ((Fraction(0), to_fraction(gain)), (to_fraction(loss), Fraction(0)))
    root = tree.treenode

    # the cost of an impossible state, more than any history can cost, kept as
    # a fraction because a float would overflow when it meets very large costs
    impossible = (cost[0][1] + cost[1][0] + 1) * (sum(1 for _ in tree) + 1)

    # downpass, d[idx][s] is the cost of the subtree below idx if idx has state s
    down: dict[int, tuple] = {}
    for node in tree:
        if node.is_leaf():
            if node.idx not in tip_states:
                raise TwigError(f"no state was given for tip {node.name!r}")
            state = tip_states[node.idx]
            if state not in (0, 1, None):
                raise TwigError(f"the state of tip {node.name!r} must be 0, 1 or None, got {state!r}")
            down[node.idx] = tuple(Fraction(0) if state in (s, None) else impossible for s in (0, 1))
        else:
            down[node.idx] = tuple(
                sum(min(cost[s][t] + down[child.idx][t] for t in (0, 1)) for child in node.children)
                for s in (0, 1)
            )

    root_cost = tuple(Fraction(0) if str(root_state) in ("free", str(s)) else impossible for s in (0, 1))
    best = min(down[root.idx][s] + root_cost[s] for s in (0, 1))

    # uppass, up[idx][s] is the cost of everything outside idx if idx has state s
    up: dict[int, tuple] = {root.idx: root_cost}
    for node in tree.traverse("preorder"):
        if node.is_leaf():
            continue
        best_child = {
            child.idx: tuple(min(cost[t][r] + down[child.idx][r] for r in (0, 1)) for t in (0, 1))
            for child in node.children
        }
        total = tuple(sum(best_child[child.idx][t] for child in node.children) for t in (0, 1))
        for child in node.children:
            up[child.idx] = tuple(
                min(
                    cost[t][s] + up[node.idx][t] + total[t] - best_child[child.idx][t]
                    for t in (0, 1)
                )
                for s in (0, 1)
            )

    state_sets = {
        node.idx: frozenset(s for s in (0, 1) if down[node.idx][s] + up[node.idx][s] == best)
        for node in tree
    }

    # preorder, choose one reconstruction given the state of the parent
    states: dict[int, int] = {}
    rooted = tuple(s for s in (0, 1) if down[root.idx][s] + root_cost[s] == best)
    states[root.idx] = rooted[0]
    for node in tree.traverse("preorder"):
        if node is root:
            continue
        parent = states[node.up.idx]
        options = tuple(cost[parent][t] + down[node.idx][t] for t in (0, 1))
        good = [t for t in (0, 1) if options[t] == min(options)]
        if len(good) == 1:
            states[node.idx] = good[0]
        elif node.is_leaf() or tie == "deltran":
            states[node.idx] = parent
        else:
            states[node.idx] = 1 - parent
    return ParsimonyResult(states, state_sets, Fraction(best))
