"""Tests for two-state parsimony (twig.labeling.core.parsimony)."""

import itertools
import random
from fractions import Fraction

import pytest

pytest.importorskip("toytree")

from twig.utils import TwigError
from twig.labeling.core import treeio
from twig.labeling.core.parsimony import ParsimonyResult, sankoff_binary, to_fraction

# the species tree of the design, with names on the internal nodes
NAMED = "(O,((M1,(A,(B,C)BC)ABC)X,(M2,(D,E)DE)Y)I)R;"


def species():
    return treeio.read_tree(NAMED)


def tip_states(tree, ones=(), unknown=()):
    return {n.idx: (None if n.name in unknown else int(n.name in ones)) for n in tree if n.is_leaf()}


def names(tree, idxs):
    return sorted(tree[i].name for i in idxs)


def ones(tree, result):
    return names(tree, [i for i, s in result.states.items() if s == 1])


class TestSpeciesTreeOfTheDesign:
    """Selecting A, B, D and E on (O,((M1,(A,(B,C))),(M2,(D,E)))); has cost 3."""

    def test_acctran_puts_one_gain_on_the_stem_and_a_loss_on_c(self):
        tree = species()
        result = sankoff_binary(tree, tip_states(tree, "ABDE"), tie="acctran")
        assert result.cost == 3
        assert ones(tree, result) == ["A", "ABC", "B", "BC", "D", "DE", "E"]

    def test_deltran_puts_separate_gains_on_a_and_b(self):
        tree = species()
        result = sankoff_binary(tree, tip_states(tree, "ABDE"), tie="deltran")
        assert result.cost == 3
        assert ones(tree, result) == ["A", "B", "D", "DE", "E"]

    def test_abc_and_bc_are_ambiguous_but_de_is_not(self):
        tree = species()
        result = sankoff_binary(tree, tip_states(tree, "ABDE"))
        assert names(tree, result.ambiguous()) == ["ABC", "BC"]
        de = next(n.idx for n in tree if n.name == "DE")
        assert result.state_sets[de] == {1}

    def test_a_tip_with_unknown_state_removes_the_tie(self):
        # C is a nuisance species whose own state is not used, so one gain on the stem is cheapest
        tree = species()
        for tie in ("acctran", "deltran"):
            result = sankoff_binary(tree, tip_states(tree, "ABDE", unknown="C"), tie=tie)
            assert result.cost == 2
            assert ones(tree, result) == ["A", "ABC", "B", "BC", "C", "D", "DE", "E"]
            assert result.ambiguous() == set()

    def test_a_loss_costs_more_than_a_gain_then_the_separate_gains_win(self):
        tree = species()
        result = sankoff_binary(tree, tip_states(tree, "ABDE"), gain=1, loss=2)
        assert ones(tree, result) == ["A", "B", "D", "DE", "E"] and result.cost == 3
        assert result.ambiguous() == set()

    def test_a_gain_costs_more_than_a_loss_then_the_stem_wins(self):
        tree = species()
        result = sankoff_binary(tree, tip_states(tree, "ABDE"), gain=2, loss=1, root_state=0, tie="deltran")
        assert ones(tree, result) == ["A", "ABC", "B", "BC", "D", "DE", "E"] and result.cost == 5

    def test_with_cheap_losses_a_free_root_prefers_the_trait_at_the_root(self):
        # losing the trait in O, M1, M2 and C costs 4, which is less than gaining it twice and losing it once
        tree = species()
        result = sankoff_binary(tree, tip_states(tree, "ABDE"), gain=2, loss=1)
        assert result.cost == 4 and result.states[tree.treenode.idx] == 1

    def test_fixing_the_root_state(self):
        tree = species()
        low = sankoff_binary(tree, tip_states(tree, "ABDE"), root_state=0)
        high = sankoff_binary(tree, tip_states(tree, "ABDE"), root_state=1)
        assert low.cost == 3 and high.cost == 4
        assert low.states[tree.treenode.idx] == 0 and high.states[tree.treenode.idx] == 1

    def test_all_tips_in_one_state_costs_nothing(self):
        tree = species()
        assert sankoff_binary(tree, tip_states(tree, "")).cost == 0
        result = sankoff_binary(tree, {n.idx: 1 for n in tree if n.is_leaf()})
        assert result.cost == 0 and set(result.states.values()) == {1}


class TestExactness:
    def test_costs_are_read_as_decimals(self):
        assert to_fraction(0.1) == Fraction(1, 10)
        assert to_fraction("1/3") == Fraction(1, 3)
        assert to_fraction(2) == 2
        with pytest.raises(TwigError, match="negative"):
            to_fraction(-1)
        with pytest.raises(TwigError, match="invalid"):
            to_fraction("cheap")

    def test_ties_are_found_with_fractional_costs(self):
        # three gains of 0.1 equal one gain of 0.1 plus two losses of 0.1
        tree = treeio.read_tree("((A,B)AB,C)R;")
        result = sankoff_binary(tree, {0: 1, 1: 1, 2: 0}, gain="0.1", loss="0.1")
        assert result.cost == Fraction(1, 10)
        assert isinstance(result.cost, Fraction)

    def test_very_large_costs_are_exact(self):
        # an impossible state used to cost a float infinity, which overflowed when it met a cost of 1e400
        tree = species()
        result = sankoff_binary(tree, tip_states(tree, "ABDE"), gain="1e400", root_state=0)
        assert result.cost == Fraction(10) ** 400 + 3
        assert ones(tree, result) == ["A", "ABC", "B", "BC", "D", "DE", "E", "I", "X", "Y"]
        result = sankoff_binary(tree, tip_states(tree, "ABDE"), gain="1e400")
        assert result.cost == 4 and result.states[tree.treenode.idx] == 1

    def test_polytomies_are_handled(self):
        tree = treeio.read_tree("(A,B,C,D)R;")
        result = sankoff_binary(tree, {0: 1, 1: 1, 2: 1, 3: 0})
        assert result.cost == 1 and result.states[tree.treenode.idx] == 1


class TestErrors:
    def test_bad_arguments(self):
        tree = species()
        states = tip_states(tree, "AB")
        with pytest.raises(TwigError, match="unknown tie"):
            sankoff_binary(tree, states, tie="random")
        with pytest.raises(TwigError, match="root state"):
            sankoff_binary(tree, states, root_state=2)
        with pytest.raises(TwigError, match="no state was given"):
            sankoff_binary(tree, {0: 1})
        bad = dict(states)
        bad[0] = 5
        with pytest.raises(TwigError, match="must be 0, 1 or None"):
            sankoff_binary(tree, bad)


###########################################################################
# comparison with a search over every reconstruction
###########################################################################

def random_newick(rng, n_tips):
    """Return a random tree of n_tips tips with some polytomies."""
    nodes = [f"t{i}" for i in range(n_tips)]
    while len(nodes) > 1:
        k = min(len(nodes), rng.choice([2, 2, 2, 3]))
        picked = rng.sample(range(len(nodes)), k)
        group = [nodes[i] for i in picked]
        nodes = [x for i, x in enumerate(nodes) if i not in picked]
        nodes.append("(" + ",".join(group) + ")")
    return nodes[0] + ";"


def brute_force(tree, tip_state, cost, root_state):
    """Return (best cost, optimal assignments) over all states of the internal nodes."""
    internal = [n.idx for n in tree if not n.is_leaf()]
    best, optimal = None, []
    for combo in itertools.product((0, 1), repeat=len(internal)):
        state = dict(zip(internal, combo))
        for n in tree:
            if n.is_leaf():
                s = tip_state[n.idx]
                state[n.idx] = s
        total = Fraction(0)
        free_tips = [n.idx for n in tree if n.is_leaf() and tip_state[n.idx] is None]
        for fill in itertools.product((0, 1), repeat=len(free_tips)):
            for idx, value in zip(free_tips, fill):
                state[idx] = value
            if root_state != "free" and state[tree.treenode.idx] != int(root_state):
                continue
            total = sum(
                (cost[state[n.up.idx]][state[n.idx]] for n in tree if n.up is not None), Fraction(0)
            )
            if best is None or total < best:
                best, optimal = total, [dict(state)]
            elif total == best:
                optimal.append(dict(state))
    return best, optimal


@pytest.mark.parametrize("seed", range(40))
def test_agrees_with_a_search_over_all_reconstructions(seed):
    rng = random.Random(seed)
    tree = treeio.read_tree(random_newick(rng, rng.randint(3, 7)))
    gain, loss = rng.choice([(1, 1), (1, 2), (2, 1), ("0.5", 1), (0, 1), (1, 0)])
    root_state = rng.choice(["free", "free", 0, 1])
    states = {n.idx: rng.choice([0, 1, 1, 0, None]) for n in tree if n.is_leaf()}
    cost = ((Fraction(0), to_fraction(gain)), (to_fraction(loss), Fraction(0)))

    best, optimal = brute_force(tree, states, cost, root_state)
    for tie in ("acctran", "deltran"):
        result = sankoff_binary(tree, states, gain, loss, root_state, tie)
        assert isinstance(result, ParsimonyResult)
        assert result.cost == best
        # the chosen reconstruction is one of the most parsimonious ones
        chosen_cost = sum((cost[result.states[n.up.idx]][result.states[n.idx]] for n in tree if n.up is not None), Fraction(0))
        assert chosen_cost == best
        for idx, value in states.items():
            if value is not None:
                assert result.states[idx] == value
        if root_state != "free":
            assert result.states[tree.treenode.idx] == int(root_state)
    # the state sets are exactly the states found in some most parsimonious reconstruction
    expected = {n.idx: frozenset(o[n.idx] for o in optimal) for n in tree}
    result = sankoff_binary(tree, states, gain, loss, root_state)
    for n in tree:
        if n.is_leaf() and states[n.idx] is not None:
            assert result.state_sets[n.idx] == expected[n.idx] == {states[n.idx]}
        else:
            assert result.state_sets[n.idx] == expected[n.idx], (seed, n.idx)
