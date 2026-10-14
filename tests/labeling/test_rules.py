"""Tests for the label rules (twig.labeling.core.rules).

The expected results are the table of the design for the species tree
(O,((M1,(A,(B,C))),(M2,(D,E)))); with the tips A, B, D and E selected, where
the internal nodes are named BC, ABC, X (M1 and ABC), DE, Y (M2 and DE) and I
(the whole ingroup).
"""

from fractions import Fraction

import pytest

pytest.importorskip("toytree")

from twig.utils import TwigError
from twig.labeling.core import constants, treeio
from twig.labeling.core.rules import (
    ALIASES,
    HELP,
    RULES,
    RuleInput,
    make_options,
    mrca,
    resolve_rule,
    tip_counts,
)

NAMED = "(O,((M1,(A,(B,C)BC)ABC)X,(M2,(D,E)DE)Y)I)R;"


@pytest.fixture
def tree():
    return treeio.read_tree(NAMED)


def idxs(tree, names):
    return frozenset(n.idx for n in tree if n.name in set(names))


def run(tree, rule, select="ABDE", unknown="", groups=(), **options):
    inp = RuleInput(tree, "Test", idxs(tree, select), idxs(tree, unknown), tuple(idxs(tree, g) for g in groups),
                    make_options(**options))
    return RULES[resolve_rule(rule)](inp)


def where(tree, output):
    return sorted(tree[i].name for i in output.labels)


class TestDesignTable:
    def test_tips(self, tree):
        assert where(tree, run(tree, "tips")) == ["A", "B", "D", "E"]

    def test_all(self, tree):
        assert where(tree, run(tree, "all")) == ["A", "B", "D", "DE", "E"]

    def test_crown(self, tree):
        assert where(tree, run(tree, "crown")) == ["A", "B", "D", "E"]

    def test_some(self, tree):
        # every branch with a selected descendant, including the ingroup branch but not the root
        assert where(tree, run(tree, "some")) == ["A", "ABC", "B", "BC", "D", "DE", "E", "I", "X", "Y"]

    def test_fitch(self, tree):
        out = run(tree, "fitch")
        assert where(tree, out) == ["A", "ABC", "B", "BC", "D", "DE", "E"]
        assert out.cost == 3
        assert where(tree, type("R", (), {"labels": dict.fromkeys(out.losses)})) == ["C"]
        assert sorted(tree[i].name for i in out.ambiguous) == ["ABC", "BC"]

    def test_crown_and_fitch_say_where_the_label_begins(self, tree):
        assert names_of(tree, run(tree, "crown").stems) == ["DE"]
        assert names_of(tree, run(tree, "fitch").stems) == ["ABC", "DE"]
        assert names_of(tree, run(tree, "fitch", ambiguous="deltran").stems) == ["DE"]
        # an ambiguous branch is not a gain for certain, and DE is
        assert names_of(tree, run(tree, "fitch", ambiguous="Maybe").stems) == ["DE"]
        for rule in ("tips", "all", "some"):
            assert run(tree, rule).stems is None, rule

    def test_hyphys_parsimony_gives_the_same_as_all_but_fitch_does_not(self, tree):
        # the regression that shows why real parsimony matters
        assert where(tree, run(tree, "all")) == ["A", "B", "D", "DE", "E"]
        assert where(tree, run(tree, "fitch")) != where(tree, run(tree, "all"))

    def test_fitch_with_deltran_gives_separate_gains(self, tree):
        assert where(tree, run(tree, "fitch", ambiguous="deltran")) == ["A", "B", "D", "DE", "E"]

    def test_fitch_with_an_unknown_tip_has_a_single_origin_and_no_loss(self, tree):
        out = run(tree, "fitch", unknown="C")
        assert where(tree, out) == ["A", "ABC", "B", "BC", "C", "D", "DE", "E"] or where(tree, out) == ["A", "ABC", "B", "BC", "D", "DE", "E"]
        assert out.losses == frozenset() and out.ambiguous == frozenset() and out.cost == 2

    def test_ambiguous_nodes_can_be_left_unlabeled_or_get_their_own_label(self, tree):
        out = run(tree, "fitch", ambiguous="unlabel")
        assert where(tree, out) == ["A", "B", "D", "DE", "E"]
        out = run(tree, "fitch", ambiguous="Maybe")
        assert where(tree, out) == ["A", "ABC", "B", "BC", "D", "DE", "E"]
        assert {tree[i].name: v for i, v in out.labels.items() if tree[i].name in ("ABC", "BC", "A", "DE")} == {
            "ABC": "Maybe", "BC": "Maybe", "A": "Test", "DE": "Test"}

    def test_an_unknown_tip_can_be_ambiguous(self, tree):
        # M1 has either state in the equally good histories, which used to be ignored for tips
        out = run(tree, "fitch", unknown=["M1"], ambiguous="unlabel")
        assert "M1" in names_of(tree, out.ambiguous) and "M1" not in where(tree, out)
        out = run(tree, "fitch", unknown=["M1"], ambiguous="Maybe")
        assert {tree[i].name: v for i, v in out.labels.items()}["M1"] == "Maybe"
        assert "M1" in where(tree, run(tree, "fitch", unknown=["M1"], ambiguous="acctran"))
        assert "M1" not in where(tree, run(tree, "fitch", unknown=["M1"], ambiguous="deltran"))

    def test_losses_are_only_certain_when_ambiguous_nodes_are_not_resolved(self, tree):
        # with the stem ambiguous the loss on C is not certain
        assert run(tree, "fitch", ambiguous="unlabel").losses == frozenset()
        assert names_of(tree, run(tree, "fitch", ambiguous="acctran").losses) == ["C"]
        assert run(tree, "fitch", ambiguous="deltran").losses == frozenset()


def names_of(tree, idxs):
    return sorted(tree[i].name for i in idxs)


class TestDirectLabels:
    def test_direct_labels_are_the_selected_tips(self, tree):
        for rule in ("tips", "all", "crown", "some", "fitch"):
            assert names_of(tree, run(tree, rule).direct) == ["A", "B", "D", "E"], rule


class TestCladeAndStem:
    def test_clade_labels_the_branch_and_everything_below(self, tree):
        out = run(tree, "clade", select="BC", groups=["BC"])
        assert where(tree, out) == ["B", "BC", "C"]
        out = run(tree, "clade", select="ABD", groups=["AB", "D"])      # two selectors, two clades
        assert where(tree, out) == ["A", "ABC", "B", "BC", "C", "D"]

    def test_stem_labels_only_the_branch_above(self, tree):
        out = run(tree, "stem", select="BC", groups=["BC"])
        assert where(tree, out) == ["C"] or where(tree, out) == ["BC"]
        # the common ancestor of B and C is BC
        assert where(tree, out) == ["BC"]
        out = run(tree, "stem", select="ABDE", groups=["AB", "DE"])
        assert where(tree, out) == ["ABC", "DE"]

    def test_a_single_tip_is_its_own_common_ancestor(self, tree):
        assert where(tree, run(tree, "stem", select="A", groups=["A"])) == ["A"]
        assert where(tree, run(tree, "clade", select="A", groups=["A"])) == ["A"]

    def test_the_root_has_no_branch(self, tree):
        with pytest.raises(TwigError, match="root as their common ancestor"):
            run(tree, "stem", select="OA", groups=["OA"])
        with pytest.raises(TwigError, match="root as their common ancestor"):
            run(tree, "clade", select="OE", groups=["OE"])

    def test_empty_groups_are_an_error(self, tree):
        with pytest.raises(TwigError, match="matches no tips"):
            run(tree, "clade", select="A", groups=["A", ""])
        with pytest.raises(TwigError, match="needs selectors"):
            run(tree, "stem", select="A", groups=[])


class TestAliases:
    @pytest.mark.parametrize("alias, rule", [
        ("None", "tips"), ("none", "tips"), ("All descendants", "all"), ("ALL DESCENDANTS", "all"),
        ("All descendants, no MRCA", "crown"), ("all descendants no mrca", "crown"),
        ("Some descendants", "some"), ("Parsimony", "fitch"), ("fitch", "fitch"), ("Clade", "clade"),
    ])
    def test_hyphy_names_and_case(self, alias, rule):
        assert resolve_rule(alias) == rule

    def test_unknown_rule(self):
        with pytest.raises(TwigError, match="unknown rule 'magic'"):
            resolve_rule("magic")

    def test_constants_list_the_registered_rules(self):
        assert constants.RULE_NAMES == tuple(RULES)
        assert set(HELP) == set(RULES) and all(HELP.values())
        assert {k.title() for k in ALIASES} >= {"None", "All Descendants"}
        assert set(ALIASES.values()) <= set(RULES)
        assert len(constants.RULE_ALIASES) == len(ALIASES)


class TestHelpers:
    def test_tip_counts(self, tree):
        selected, total = tip_counts(tree, idxs(tree, "ABDE"))
        by_name = {n.name: (selected[n.idx], total[n.idx]) for n in tree if n.name}
        assert by_name["ABC"] == (2, 3) and by_name["DE"] == (2, 2) and by_name["R"] == (4, 8)

    def test_mrca(self, tree):
        assert mrca(tree, idxs(tree, "BC")).name == "BC"
        assert mrca(tree, idxs(tree, "AC")).name == "ABC"
        assert mrca(tree, idxs(tree, "AD")).name == "I"
        assert mrca(tree, idxs(tree, "OA")).name == "R"
        assert mrca(tree, idxs(tree, "A")).name == "A"
        with pytest.raises(TwigError):
            mrca(tree, [])

    def test_make_options_reads_costs_exactly(self):
        opts = make_options(gain="0.1", loss=2, ambiguous="deltran", root_state=0)
        assert opts.gain_cost == Fraction(1, 10) and opts.loss_cost == 2 and opts.root_state == "0"
