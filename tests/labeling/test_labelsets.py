"""Tests for applying label sets to a tree (twig.labeling.core.labelsets)."""

from fractions import Fraction

import pytest

pytest.importorskip("toytree")

from twig.utils import TwigError
from twig.labeling.core import rooting, treeio
from twig.labeling.core.labelsets import (
    LabelingConfig,
    LabelSet,
    apply_label_sets,
    nuisance_set,
    parse_label_args,
    read_spec,
    select_tips,
    unused_settings,
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

    def test_an_unknown_tip_with_an_ambiguous_state_is_treated_like_an_ambiguous_branch(self):
        # M1 has either state in the equally good histories, and it used to get the background label
        tree = fresh(NAMED)
        result = apply_label_sets(tree, [exact(*DIOECIOUS)], cfg(unknown=("M1",), ambiguous="unlabel", background="Reference"))
        labels = by_name(tree, result)
        assert labels["M1"] is None and labels["X"] is None and labels["M2"] == "Reference"
        result = apply_label_sets(tree, [exact(*DIOECIOUS)], cfg(unknown=("M1",), ambiguous="Unsure", background="Reference"))
        assert by_name(tree, result)["M1"] == "Unsure"

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
# stems, the branches on which a trait arose
###########################################################################

class TestStems:
    CROWN = "(O,(M,((D1,D2)D,E)DE)Y)R;"

    def test_with_crown_the_branch_left_out_is_the_stem(self):
        # the first crown branch D used to be the stem, and the real stem DE got the background label
        tree = fresh(self.CROWN)
        sets = [exact("D1", "D2", "E", rule="crown")]
        result = apply_label_sets(tree, sets, cfg(stems="Stem"))
        labels, kinds = by_name(tree, result), by_name(tree, result, "kind")
        assert labels["DE"] == "Stem" and labels["D"] == "Test"
        assert kinds["DE"] == "stem" and kinds["D"] == "crown"
        result = apply_label_sets(tree, sets, cfg(stems="unlabel", background="Reference"))
        labels = by_name(tree, result)
        assert labels["D"] == "Test" and labels["DE"] is None and labels["Y"] == "Reference"
        result = apply_label_sets(tree, sets, cfg(background="Reference"))
        assert by_name(tree, result)["DE"] is None and result.counts["stems"] == 1

    def test_with_crown_on_the_species_tree(self):
        tree = fresh(NAMED)
        result = apply_label_sets(tree, [exact(*DIOECIOUS, rule="crown")], cfg(stems="Stem"))
        assert {k: v for k, v in by_name(tree, result).items() if v} == {
            "A": "Test", "B": "Test", "D": "Test", "E": "Test", "DE": "Stem"}

    def test_the_label_of_ambiguous_branches_is_not_a_trait_label(self):
        # ABC, which is only ambiguous, used to become a stem
        tree = fresh(NAMED)
        result = apply_label_sets(tree, [exact(*DIOECIOUS)], cfg(ambiguous="Unsure", stems="Stem"))
        labels = by_name(tree, result)
        assert labels["ABC"] == labels["BC"] == "Unsure" and labels["DE"] == "Stem"
        assert by_name(tree, result, "source")["ABC"] == "ambiguous"
        result = apply_label_sets(tree, [exact(*DIOECIOUS)], cfg(ambiguous="Unsure", stems="unlabel"))
        labels = by_name(tree, result)
        assert labels["ABC"] == "Unsure" and labels["DE"] is None and result.counts["stems"] == 1

    def test_with_fitch_a_branch_below_a_root_that_has_the_trait_is_not_a_stem(self):
        # the root has the trait and C and E lost it, so no gain happens anywhere
        tree = fresh("((A,B)AB,(C,(D,E)DE)CDE)R;")
        result = apply_label_sets(tree, [exact("A", "B", "D")], cfg(gain_cost=Fraction(2), stems="unlabel"))
        labels = by_name(tree, result)
        assert labels["AB"] == labels["CDE"] == labels["DE"] == "Test" and result.counts["stems"] == 0
        result = apply_label_sets(tree, [exact("A", "B", "C", "D")], cfg(stems="Stem"))
        assert by_name(tree, result)["AB"] == by_name(tree, result)["CDE"] == "Test"
        # with the root fixed to lack the trait, the gains are stems again
        result = apply_label_sets(tree, [exact("A", "B", "C", "D")], cfg(stems="Stem", root_state="0"))
        assert by_name(tree, result)["AB"] == "Stem"

    def test_with_other_rules_a_label_on_every_child_of_the_root_has_no_stem_there(self):
        tree = fresh("((A,B)AB,(C,D)CD)R;")
        result = apply_label_sets(tree, [exact("A", "B", "C", "D", rule="all")], cfg(stems="Stem"))
        assert by_name(tree, result)["AB"] == by_name(tree, result)["CD"] == "Test"
        result = apply_label_sets(tree, [exact("A", "B", rule="all")], cfg(stems="Stem"))
        assert by_name(tree, result)["AB"] == "Stem"


class TestStemsOfLabelsInTheTree:
    TEST_ONLY = "(O,((M1,(A{Test},(B{Test},C)BC{Test})ABC{Test})X,(M2,(D{Test},E{Test})DE{Test})Y)I)R;"

    def test_stems_act_on_labels_kept_from_the_input_tree(self):
        # --stems used to do nothing to labels that the tree already had
        tree = fresh(self.TEST_ONLY)
        result = apply_label_sets(tree, [], cfg(stems="Stem"))
        labels = by_name(tree, result)
        assert labels["ABC"] == labels["DE"] == "Stem" and labels["BC"] == "Test"
        assert any("--stems treats 2 branches labeled in the input tree as stems (A+B+C, D+E)" in w for w in result.warnings)
        result = apply_label_sets(tree, [], cfg(stems="unlabel", background="Reference"))
        labels = by_name(tree, result)
        assert labels["ABC"] is None and labels["DE"] is None and labels["X"] == "Reference" and labels["BC"] == "Test"

    def test_labels_in_the_tree_are_left_alone_with_keep_or_clear(self):
        tree = fresh(self.TEST_ONLY)
        result = apply_label_sets(tree, [], cfg())
        assert by_name(tree, result)["ABC"] == "Test" and not result.warnings and result.counts["stems"] == 0
        result = apply_label_sets(tree, [exact("A", "B", rule="all")], cfg(stems="Stem", existing="clear"))
        assert by_name(tree, result)["DE"] is None and not result.warnings

    def test_labels_that_the_run_gives_to_the_background_are_not_traits(self):
        # the loss clade below ABC carries the background label of an earlier run
        tree = fresh("(O{Ref},((A{Test},(B{Ref},C{Ref})BC{Ref})ABC{Test},D{Ref})X{Ref})R;")
        result = apply_label_sets(tree, [], cfg(stems="Stem", background="Ref"))
        labels = by_name(tree, result)
        assert labels["BC"] == "Ref" and labels["ABC"] == "Stem"
        result = apply_label_sets(tree, [], cfg(stems="Stem"))
        assert by_name(tree, result)["BC"] == "Stem"      # without --background, Ref looks like a trait, and the warning says so
        assert any("A+B+C, B+C" in w for w in result.warnings)


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

    def test_a_branch_that_goes_to_the_same_label_is_not_lost(self):
        # B is claimed by two sets with the label Test, and a warning used to say that it went to another label
        tree = fresh(NAMED)
        result = apply_label_sets(tree, [exact("A", "B", rule="tips"), exact("B", "C", rule="tips")], cfg())
        assert by_name(tree, result)["B"] == "Test"
        assert not [w for w in result.warnings if "went to another label" in w]
        assert [r.kept for r in result.reports] == [2, 1]

    def test_a_branch_lost_to_a_label_in_the_tree_names_the_existing_option(self):
        tree = fresh("(O,((M1,(A,(B,C)BC{Old})ABC)X,(M2,(D,E)DE)Y)I)R;")
        result = apply_label_sets(tree, [exact(*DIOECIOUS)], cfg())
        assert by_name(tree, result)["BC"] == "Old"
        assert any("label Test: 1 of 7 branches kept the label they had in the input tree because of --existing keep" in w for w in result.warnings)
        assert not [w for w in result.warnings if "--precedence" in w]

    def test_the_error_names_the_branch_by_its_tips(self):
        # the branch used to be named by its internal number, such as branch 3
        tree = fresh(NAMED)
        sets = [LabelSet("Outer", "stem", ("^[BC]$",)), LabelSet("Inner", "stem", ("^[BC]$",))]
        with pytest.raises(TwigError, match=r"the labels 'Inner' and 'Outer' claim the same branch above B\+C, which --precedence error forbids"):
            apply_label_sets(tree, sets, cfg(exact=False, precedence="error"))

    def test_what_a_set_keeps_is_counted_after_skipped_tips_and_stems(self):
        # the log used to say that Test keeps 7 branches when only BC was left
        tree = fresh(NAMED)
        result = apply_label_sets(tree, [exact(*DIOECIOUS)], cfg(stems="unlabel", leaf_nodes="skip"))
        report = result.reports[0]
        assert (report.claimed, report.won, report.kept) == (7, 7, 1)
        result = apply_label_sets(tree, [exact(*DIOECIOUS)], cfg(ambiguous="Unsure", stems="Stem"))
        assert result.reports[0].kept == 6          # A, B, BC, ABC, D and E, but not DE, which became Stem

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

    def test_skipped_tips_are_left_out_of_the_background(self):
        # the skipped tips of the selected species used to get the background label
        tree = fresh(NAMED)
        result = apply_label_sets(tree, [exact(*DIOECIOUS)], cfg(stems="unlabel", leaf_nodes="skip", background="Reference"))
        labels = by_name(tree, result)
        assert [labels[n] for n in "ABDE"] == [None] * 4
        assert labels["BC"] == "Test" and labels["C"] == "Reference" and labels["M1"] == "Reference"
        result = apply_label_sets(tree, [exact(*DIOECIOUS, rule="all")], cfg(leaf_nodes="skip", background="Reference"))
        labels = by_name(tree, result)
        assert [labels[n] for n in "ABDE"] == [None] * 4 and labels["DE"] == "Test" and labels["C"] == "Reference"

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

    def test_a_tip_without_a_species_is_never_chosen(self):
        # tree-label -U ignore gives tips that cannot be resolved the species None
        tree = treeio.read_tree("((a_1,a_2),(b_1,c_1));")
        species = {"a_1": "a", "a_2": None, "b_1": "b", "c_1": "c"}.get
        chosen, unused = select_tips(tree, ["a_2", "a"], exact=True, select_on="species", species_of=species)
        assert sorted(tree[i].name for i in chosen) == ["a_1"] and unused == ["a_2"]
        chosen, _ = select_tips(tree, ["."], select_on="species", species_of=species)
        assert sorted(tree[i].name for i in chosen) == ["a_1", "b_1", "c_1"]
        config = cfg(select_on="species", invert=True, unknown=("a_2",))
        result = apply_label_sets(tree, [LabelSet("Test", "tips", ("a",))], config, species_of=species)
        labeled = sorted(tree[i].name for i, v in result.labels.items() if v)
        assert labeled == ["b_1", "c_1"]


class TestWarningsAndErrors:
    def test_selectors_that_match_nothing_are_reported(self):
        tree = fresh()
        result = apply_label_sets(tree, [exact("A", "ZZZ", rule="tips")], cfg())
        assert result.reports[0].unused == ["ZZZ"] and any("ZZZ" in w for w in result.warnings)
        result = apply_label_sets(tree, [exact("ZZZ", rule="tips")], cfg())
        assert any("choose no tips" in w for w in result.warnings)

    def test_every_rule_but_tips_warns_on_an_unrooted_tree(self):
        unrooted = fresh("((M1,(A,(B,C))),(M2,(D,E)),O);")
        result = apply_label_sets(unrooted, [exact(*DIOECIOUS)], cfg(gain_cost=Fraction(2)))
        assert any("depends on the position of the root" in w for w in result.warnings)
        # equal costs used to give no warning, although the labels depend on the root as much
        result = apply_label_sets(unrooted, [exact(*DIOECIOUS)], cfg())
        assert [w for w in result.warnings if "depends on the position of the root" in w and "-r with the outgroup" in w]
        result = apply_label_sets(unrooted, [exact(*DIOECIOUS, rule="all")], cfg())
        assert any("the all rule gives a result that depends on the position of the root" in w for w in result.warnings)
        assert not apply_label_sets(unrooted, [exact(*DIOECIOUS, rule="tips")], cfg()).warnings
        assert not apply_label_sets(fresh(), [exact(*DIOECIOUS)], cfg()).warnings

    def test_the_labels_of_an_unrooted_tree_depend_on_where_its_basal_node_lies(self):
        # the same unrooted tree as above, written with the basal node next to A
        tree = fresh("(A,(B,C),(M1,(O,(M2,(D,E)))));")
        config = cfg(stems="unlabel", losses="Rev", background="Reference")
        result = apply_label_sets(tree, [exact(*DIOECIOUS)], config)
        assert any("Give -r with the outgroup" in w for w in result.warnings)
        rooted = rooting.reroot(tree, ["O"]).tree
        result = apply_label_sets(rooted, [exact(*DIOECIOUS)], config)
        assert not result.warnings
        assert written(rooted, result) == (
            "(O{Reference},((M2{Reference},(D{Test},E{Test})){Reference},"
            "(M1{Reference},(A{Test},(B{Test},C{Rev}){Test})){Reference}){Reference});"
        )

    def test_stems_of_labels_in_an_unrooted_tree_warn_too(self):
        tree = fresh("((M1,(A{Test},(B{Test},C){Test}){Test}),(M2,(D{Test},E{Test}){Test}),O);")
        result = apply_label_sets(tree, [], cfg(stems="Stem"))
        assert any("the stems of the labels in the input tree depend on the position of the root" in w for w in result.warnings)

    def test_label_names_are_checked(self):
        with pytest.raises(TwigError, match="invalid label name"):
            apply_label_sets(fresh(), [LabelSet("All", "tips", ("A",))], cfg())
        with pytest.raises(TwigError, match="invalid label name"):
            apply_label_sets(fresh(), [exact("A")], cfg(background="1ref"))
        with pytest.raises(TwigError, match="invalid label name"):
            apply_label_sets(fresh(), [exact("A")], cfg(stems="Stem-1"))

    @pytest.mark.parametrize("settings", [
        dict(ambiguous="Deltran"), dict(ambiguous="ACCTRAN"), dict(ambiguous="Unlabel"),
        dict(stems="Keep"), dict(losses="UNLABEL"),
    ])
    def test_modes_in_another_case_are_refused(self, settings):
        # Deltran used to become a label name that was put on the ambiguous branches
        with pytest.raises(TwigError, match="looks like the mode"):
            apply_label_sets(fresh(), [exact(*DIOECIOUS)], cfg(**settings))

    def test_a_mode_in_another_case_is_refused_in_a_spec_option(self):
        item = LabelSet("Test", "fitch", DIOECIOUS, ambiguous="Deltran")
        with pytest.raises(TwigError, match="ambiguous=Deltran of the label 'Test' would be read as a label name"):
            apply_label_sets(fresh(), [item], cfg())

    def test_parsimony_options_without_a_fitch_label_are_reported(self):
        # they used to be accepted and ignored without a word
        config = cfg(losses="Rev", ambiguous="Amb", gain_cost=Fraction(5), unknown=("C",))
        assert unused_settings([exact(*DIOECIOUS, rule="all")], config) == [
            "no label uses the fitch rule, so --losses, --ambiguous, --gain-cost and --unknown have no effect"
        ]
        assert unused_settings([exact(*DIOECIOUS, rule="all"), exact("A")], config) == []
        assert unused_settings([exact(*DIOECIOUS, rule="all")], cfg()) == []
        assert unused_settings([], cfg(root_state="1")) == ["no label uses the fitch rule, so --root-state has no effect"]
        with_options = LabelSet("Test", "all", DIOECIOUS, ambiguous="deltran", gain_cost=Fraction(2))
        assert unused_settings([with_options], cfg()) == [
            "label Test: the all rule is not fitch, so ambiguous= and gain-cost= have no effect"
        ]

    def test_unknown_tips_that_match_nothing_are_reported(self):
        result = apply_label_sets(fresh(), [exact(*DIOECIOUS)], cfg(unknown=("Cx",)))
        assert any("--unknown has 1 listed name or pattern that matches no tip (Cx)" in w for w in result.warnings)

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

    def test_spec_file_with_a_byte_order_mark(self, tmp_path):
        path = tmp_path / "labels.tsv"
        path.write_bytes("\ufeffTest\ttips\tA\n".encode("utf-8"))
        assert read_spec(path) == [LabelSet("Test", "tips", ("A",))]
        path.write_bytes(b"Test\ttips\t\xe9\n")
        with pytest.raises(TwigError, match="not valid UTF-8"):
            read_spec(path)

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
