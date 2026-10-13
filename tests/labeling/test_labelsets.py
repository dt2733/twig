"""Tests for applying label sets to a tree (twig.labeling.core.labelsets)."""

from fractions import Fraction

import pytest

pytest.importorskip("toytree")

from twig.utils import TwigError
from twig.labeling.core import treeio
from twig.labeling.core.labelsets import (
    LabelingConfig,
    LabelSet,
    apply_label_sets,
    nuisance_set,
    parse_label_args,
    read_spec,
    select_tips,
)
from twig.labeling.core.newick import to_newick
from conftest import LABELED, SPECIES

NAMED = "(O,((M1,(A,(B,C)BC)ABC)X,(M2,(D,E)DE)Y)I)R;"
DIOECIOUS = ("A", "B", "D", "E")


def fresh(newick=SPECIES):
    return treeio.read_tree(newick)


def written(tree, result):
    """Return the Newick of the tree with the labels of a result."""
    labeled = treeio.set_labels(tree, result.labels)
    return to_newick(labeled, dists=False)


def exact(*names, rule="fitch", label="Test"):
    return LabelSet(label, rule, tuple(names))


def cfg(**kwargs):
    kwargs.setdefault("exact", True)
    return LabelingConfig(**kwargs)


def by_name(tree, result, field="labels"):
    """Return {node name: value} for the named nodes of a tree read with NAMED."""
    data = getattr(result, field)
    return {n.name: data.get(n.idx) for n in tree if n.name and n.up is not None}


###########################################################################
# the labeled species tree of the design
###########################################################################

class TestDesignExample:
    def test_fitch_with_stems_unlabeled_losses_named_and_background(self):
        tree = fresh()
        config = cfg(stems="unlabel", losses="Rev", background="Reference")
        result = apply_label_sets(tree, [exact(*DIOECIOUS)], config)
        assert written(tree, result) == LABELED

    def test_the_same_tree_when_c_is_a_nuisance_species(self):
        # C is marked directly instead of being found as a reversion
        tree = fresh()
        config = cfg(stems="unlabel", background="Reference")
        result = apply_label_sets(tree, [exact(*DIOECIOUS), nuisance_set(["C"], "Rev")], config)
        assert written(tree, result) == LABELED

    def test_without_stems_switch_the_stems_are_labeled_test(self):
        tree = fresh(NAMED)
        result = apply_label_sets(tree, [exact(*DIOECIOUS)], cfg())
        assert {k: v for k, v in by_name(tree, result).items() if v} == {
            "A": "Test", "B": "Test", "BC": "Test", "ABC": "Test", "D": "Test", "E": "Test", "DE": "Test"}

    def test_stems_unlabeled(self):
        tree = fresh(NAMED)
        result = apply_label_sets(tree, [exact(*DIOECIOUS)], cfg(stems="unlabel"))
        assert sorted(k for k, v in by_name(tree, result).items() if v) == ["A", "B", "BC", "D", "E"]

    def test_stems_get_a_label_of_their_own(self):
        tree = fresh(NAMED)
        result = apply_label_sets(tree, [exact(*DIOECIOUS)], cfg(stems="Stem"))
        labels = by_name(tree, result)
        assert labels["ABC"] == labels["DE"] == "Stem" and labels["BC"] == "Test" and labels["A"] == "Test"
        assert by_name(tree, result, "source")["ABC"] == "stem"

    def test_kinds_and_sources(self):
        tree = fresh(NAMED)
        result = apply_label_sets(tree, [exact(*DIOECIOUS)], cfg(stems="unlabel", losses="Rev", background="Reference"))
        kinds, sources = by_name(tree, result, "kind"), by_name(tree, result, "source")
        assert kinds["A"] == "tip" and kinds["ABC"] == "stem" and kinds["DE"] == "stem" and kinds["BC"] == "crown"
        assert kinds["X"] == "internal"
        assert sources["BC"] == "Test:fitch" and sources["C"] == "loss" and sources["O"] == "background"
        assert sources["ABC"] == ""                  # an unlabeled stem
        assert by_name(tree, result)["C"] == "Rev"

    def test_the_counts_describe_what_was_done(self):
        tree = fresh()
        result = apply_label_sets(tree, [exact(*DIOECIOUS)], cfg(stems="unlabel", losses="Rev", background="Reference"))
        assert result.counts == {"stems": 2, "losses": 1, "ambiguous": 0, "background": 6}
        report = result.reports[0]
        assert (report.selected, report.claimed, report.won, report.cost) == (4, 7, 7, Fraction(3))

    def test_deltran_gives_separate_gains(self):
        tree = fresh(NAMED)
        result = apply_label_sets(tree, [exact(*DIOECIOUS)], cfg(ambiguous="deltran"))
        assert sorted(k for k, v in by_name(tree, result).items() if v) == ["A", "B", "D", "DE", "E"]

    def test_ambiguous_nodes_can_be_left_out_of_the_background(self):
        tree = fresh(NAMED)
        result = apply_label_sets(tree, [exact(*DIOECIOUS)], cfg(ambiguous="unlabel", background="Reference"))
        labels = by_name(tree, result)
        assert labels["ABC"] is None and labels["BC"] is None        # uncertain, so neither test nor reference
        assert labels["X"] == "Reference" and labels["A"] == "Test" and labels["DE"] == "Test"
        assert result.counts["ambiguous"] == 2

    def test_other_rules(self):
        tree = fresh(NAMED)
        result = apply_label_sets(tree, [exact(*DIOECIOUS, rule="all")], cfg())
        assert sorted(k for k, v in by_name(tree, result).items() if v) == ["A", "B", "D", "DE", "E"]
        result = apply_label_sets(tree, [exact(*DIOECIOUS, rule="crown")], cfg())
        assert sorted(k for k, v in by_name(tree, result).items() if v) == ["A", "B", "D", "E"]
        result = apply_label_sets(tree, [exact(*DIOECIOUS, rule="some")], cfg(stems="unlabel"))
        # the stems are the labeled branches whose parent is not labeled, here the ingroup branch
        assert by_name(tree, result)["I"] is None and by_name(tree, result)["X"] == "Test"


###########################################################################
# nuisance species
###########################################################################

class TestNuisance:
    def test_a_nuisance_tip_is_unknown_to_parsimony_so_the_stem_is_a_single_gain(self):
        tree = fresh(NAMED)
        with_nuisance = apply_label_sets(tree, [exact(*DIOECIOUS), nuisance_set(["C"])], cfg(ambiguous="deltran"))
        # without the nuisance designation deltran would have given separate gains in A and B
        labels = by_name(tree, with_nuisance)
        assert labels["ABC"] == labels["BC"] == "Test" and labels["C"] == "Nuisance"
        assert with_nuisance.reports[0].cost == 2

    def test_the_nuisance_label_beats_every_other_label_and_causes_no_warning(self):
        tree = fresh(NAMED)
        sets = [exact(*DIOECIOUS), nuisance_set(["D"], "Skip")]
        for precedence in ("specific", "first", "last"):
            result = apply_label_sets(tree, sets, cfg(precedence=precedence))
            labels = by_name(tree, result)
            assert labels["D"] == "Skip" and labels["E"] == "Test", precedence
            assert not [w for w in result.warnings if "went to another label" in w]

    def test_a_clade_of_nuisance_species_is_labeled_as_a_whole(self):
        tree = fresh(NAMED)
        result = apply_label_sets(tree, [nuisance_set(["D", "E"])], cfg())
        assert sorted(k for k, v in by_name(tree, result).items() if v) == ["D", "DE", "E"]
        assert set(by_name(tree, result).values()) <= {"Nuisance", None}

    def test_the_stems_switch_does_not_touch_nuisance_clades(self):
        tree = fresh(NAMED)
        result = apply_label_sets(tree, [nuisance_set(["D", "E"])], cfg(stems="unlabel", background="Reference"))
        assert by_name(tree, result)["DE"] == "Nuisance"

    def test_nuisance_tips_keep_their_label_when_tip_labels_are_skipped(self):
        tree = fresh(NAMED)
        sets = [exact(*DIOECIOUS, rule="all"), nuisance_set(["C"], "Rev")]
        result = apply_label_sets(tree, sets, cfg(leaf_nodes="skip"))
        labels = by_name(tree, result)
        assert labels["C"] == "Rev" and labels["A"] is None and labels["DE"] == "Test"

    def test_nuisance_species_are_left_out_of_the_background(self):
        tree = fresh(NAMED)
        result = apply_label_sets(tree, [nuisance_set(["C"])], cfg(background="Reference"))
        assert by_name(tree, result)["C"] == "Nuisance" and by_name(tree, result)["BC"] == "Reference"

    def test_a_custom_name_and_the_default_name(self):
        assert nuisance_set(["C"]).name == "Nuisance" and nuisance_set(["C"], "Rev").name == "Rev"
        assert nuisance_set(["C"]).nuisance and nuisance_set(["C"]).rule == "all"


###########################################################################
# precedence, existing labels and other switches
###########################################################################

class TestPrecedence:
    SETS = [
        LabelSet("Test", "clade", ("^[AB]$",)),     # one selector, so the clade of A and B, which includes C
        LabelSet("Rev", "tips", ("^C$",)),
    ]

    def test_a_label_selected_for_a_tip_beats_one_inherited_from_a_clade(self):
        tree = fresh(NAMED)
        result = apply_label_sets(tree, self.SETS, cfg(exact=False))
        assert by_name(tree, result)["C"] == "Rev" and by_name(tree, result)["B"] == "Test"
        assert any("Test: 1 of" in w for w in result.warnings)       # the overridden set is not silent
        assert result.reports[0].claimed > result.reports[0].won

    def test_first_and_last(self):
        tree = fresh(NAMED)
        assert by_name(tree, apply_label_sets(tree, self.SETS, cfg(exact=False, precedence="first")))["C"] == "Test"
        assert by_name(tree, apply_label_sets(tree, self.SETS, cfg(exact=False, precedence="last")))["C"] == "Rev"

    def test_error(self):
        tree = fresh(NAMED)
        with pytest.raises(TwigError, match="--precedence error forbids"):
            apply_label_sets(tree, self.SETS, cfg(exact=False, precedence="error"))

    def test_existing_labels_have_priority_unless_cleared(self):
        tree = treeio.read_tree("(O,((M1,(A{Old},(B,C))),(M2,(D,E))));")
        sets = [exact("A", "B", rule="tips")]
        kept = apply_label_sets(tree, sets, cfg())
        assert kept.labels[next(n.idx for n in tree if n.name == "A")] == "Old"
        assert kept.source[next(n.idx for n in tree if n.name == "A")] == "existing"
        cleared = apply_label_sets(tree, sets, cfg(existing="clear"))
        assert cleared.labels[next(n.idx for n in tree if n.name == "A")] == "Test"

    def test_skipping_tip_labels(self):
        tree = fresh(NAMED)
        result = apply_label_sets(tree, [exact(*DIOECIOUS, rule="all")], cfg(leaf_nodes="skip"))
        assert [k for k, v in by_name(tree, result).items() if v] == ["DE"]

    def test_the_tips_rule_labels_nothing_when_tip_labels_are_skipped(self):
        with pytest.raises(TwigError, match="labels nothing"):
            apply_label_sets(fresh(), [exact("A", rule="tips")], cfg(leaf_nodes="skip"))

    def test_invert_labels_the_other_tips(self):
        tree = fresh(NAMED)
        result = apply_label_sets(tree, [exact("O", "M1", "M2", rule="tips")], cfg(invert=True))
        assert sorted(k for k, v in by_name(tree, result).items() if v) == ["A", "B", "C", "D", "E"]

    def test_losses_can_be_unlabeled_and_kept_out_of_the_background(self):
        tree = fresh(NAMED)
        result = apply_label_sets(tree, [exact(*DIOECIOUS)], cfg(losses="unlabel", background="Reference"))
        assert by_name(tree, result)["C"] is None
        result = apply_label_sets(tree, [exact(*DIOECIOUS)], cfg(losses="keep", background="Reference"))
        assert by_name(tree, result)["C"] == "Reference"

    def test_selecting_on_species(self):
        tree = treeio.read_tree("((a_1,a_2),(b_1,c_1));")
        species = lambda name: name.split("_")[0]
        chosen, unused = select_tips(tree, ["a"], exact=True, select_on="species", species_of=species)
        assert sorted(tree[i].name for i in chosen) == ["a_1", "a_2"] and unused == []
        with pytest.raises(TwigError, match="way to map names to species"):
            select_tips(tree, ["a"], select_on="species")
        result = apply_label_sets(tree, [LabelSet("Test", "all", ("a",))], cfg(select_on="species"), species_of=species)
        assert sorted(tree[i].name for i, v in result.labels.items() if v and tree[i].name) == ["a_1", "a_2"]


class TestWarningsAndErrors:
    def test_selectors_that_match_nothing_are_reported(self):
        tree = fresh()
        result = apply_label_sets(tree, [exact("A", "ZZZ", rule="tips")], cfg())
        assert result.reports[0].unused == ["ZZZ"] and any("ZZZ" in w for w in result.warnings)
        result = apply_label_sets(tree, [exact("ZZZ", rule="tips")], cfg())
        assert any("choose no tips" in w for w in result.warnings)

    def test_unequal_costs_on_an_unrooted_tree(self):
        unrooted = fresh("((M1,(A,(B,C))),(M2,(D,E)),O);")
        result = apply_label_sets(unrooted, [exact(*DIOECIOUS)], cfg(gain_cost=Fraction(2)))
        assert any("depends on the position of the root" in w for w in result.warnings)
        result = apply_label_sets(unrooted, [exact(*DIOECIOUS)], cfg())
        assert not result.warnings

    def test_label_names_are_checked(self):
        with pytest.raises(TwigError, match="invalid label name"):
            apply_label_sets(fresh(), [LabelSet("All", "tips", ("A",))], cfg())
        with pytest.raises(TwigError, match="invalid label name"):
            apply_label_sets(fresh(), [exact("A")], cfg(background="1ref"))
        with pytest.raises(TwigError, match="invalid label name"):
            apply_label_sets(fresh(), [exact("A")], cfg(stems="Stem-1"))

    def test_costs_cannot_both_be_zero(self):
        # every history would cost nothing, and acctran used to give an arbitrary alternating labeling
        with pytest.raises(TwigError, match="costs nothing to gain or to lose"):
            apply_label_sets(fresh(), [exact(*DIOECIOUS)], cfg(gain_cost=Fraction(0), loss_cost=Fraction(0)))
        free = LabelSet("Test", "fitch", DIOECIOUS, gain_cost=Fraction(0), loss_cost=Fraction(0))
        with pytest.raises(TwigError, match="costs nothing to gain or to lose"):
            apply_label_sets(fresh(), [free], cfg())
        # one free change is allowed, and so are zero costs for a rule that does not use them
        apply_label_sets(fresh(), [exact(*DIOECIOUS)], cfg(gain_cost=Fraction(0)))
        apply_label_sets(fresh(), [exact(*DIOECIOUS, rule="all")], cfg(gain_cost=Fraction(0), loss_cost=Fraction(0)))

    def test_per_set_costs_override_the_defaults(self):
        tree = fresh(NAMED)
        loss_costly = LabelSet("Test", "fitch", DIOECIOUS, loss_cost=Fraction(2))
        result = apply_label_sets(tree, [loss_costly], cfg())
        assert sorted(k for k, v in by_name(tree, result).items() if v) == ["A", "B", "D", "DE", "E"]


###########################################################################
# reading label sets
###########################################################################

class TestReadingSets:
    def test_label_arguments(self):
        sets = parse_label_args([["Test:fitch", "dioecious.txt"], ["Out", "^O$", "^Z$"]], default_rule="tips")
        assert sets[0] == LabelSet("Test", "fitch", ("dioecious.txt",))
        assert sets[1] == LabelSet("Out", "tips", ("^O$", "^Z$"))

    def test_hyphy_rule_names_work_in_label_arguments(self):
        assert parse_label_args([["Test:Parsimony", "x"]])[0].rule == "fitch"

    def test_label_argument_errors(self):
        with pytest.raises(TwigError, match="needs a name and at least one selector"):
            parse_label_args([["Test"]])
        with pytest.raises(TwigError, match="unknown rule"):
            parse_label_args([["Test:magic", "A"]])

    def test_spec_file(self, tmp_path):
        path = tmp_path / "labels.tsv"
        path.write_text(
            "# label  rule  selectors  options\n"
            "Test\tfitch\tdioecious.txt\tambiguous=deltran gain-cost=2 loss-cost=1\n"
            "Out\ttips\t^O$ ^Z$\n"
            "Rest\t-\tM1 M2\n"
        )
        sets = read_spec(path, default_rule="all")
        assert sets[0] == LabelSet("Test", "fitch", ("dioecious.txt",), "deltran", Fraction(2), Fraction(1))
        assert sets[1] == LabelSet("Out", "tips", ("^O$", "^Z$"))
        assert sets[2] == LabelSet("Rest", "all", ("M1", "M2"))

    @pytest.mark.parametrize("text, message", [
        ("Test\tfitch\n", r"labels.tsv:1: expected 3 or 4"),
        ("Test\tmagic\tA\n", r"labels.tsv:1: unknown rule"),
        ("Test\tfitch\tA\tcolor=red\n", r"labels.tsv:1: unknown option"),
        ("\tfitch\tA\n", r"cannot be empty"),
        ("# nothing\n", r"no label sets found"),
    ])
    def test_spec_errors(self, tmp_path, text, message):
        path = tmp_path / "labels.tsv"
        path.write_text(text)
        with pytest.raises(TwigError, match=message):
            read_spec(path)
