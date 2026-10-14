"""Tests for roots and outgroups (twig.labeling.core.rooting)."""

import pytest

pytest.importorskip("toytree")

from twig.utils import TwigError
from twig.labeling.core import rooting, treeio
from conftest import LABELED, SPECIES

# every branch has its own label and a length of 1, and the two branches at the root differ
DISTINCT = "(O{X}:1,((M1{R}:1,(A{T}:1,(B{T}:1,C{V}:1){T}:1){S}:1){R}:1,(M2{R}:1,(D{T}:1,E{T}:1){W}:1){R}:1){Z}:1);"


def splits(tree):
    """Return {split: [(label, length), ...]} where a split is the smaller side of a branch."""
    below = rooting.tip_name_sets(tree)
    everything = frozenset(n.name for n in tree if n.is_leaf())
    out = {}
    for node in tree:
        if node.up is None:
            continue
        side = below[node.idx]
        key = tuple(sorted(min(side, everything - side, key=lambda s: (len(s), sorted(s)))))
        out.setdefault(key, []).append((node.label, node.dist))
    return out


def labels_of(tree):
    """Return {split: label} for splits whose edges all have the same label."""
    out = {}
    for key, edges in splits(tree).items():
        if len({label for label, _ in edges}) == 1:
            out[key] = edges[0][0]
    return out


class TestMergeRootLabels:
    @pytest.mark.parametrize("policy", ["unlabel", "labeled", "outgroup", "ingroup", "error", "Test"])
    def test_equal_labels_are_kept_by_every_policy(self, policy):
        assert rooting.merge_root_labels("Test", "Test", policy, True) == ("Test", False)
        assert rooting.merge_root_labels(None, None, policy, False) == (None, False)

    def test_unlabel_removes_the_label_of_a_branch_that_merges_different_labels(self):
        assert rooting.merge_root_labels("A", "B", "unlabel") == (None, True)
        assert rooting.merge_root_labels("A", None, "unlabel") == (None, True)

    def test_labeled_uses_the_only_label_but_not_a_choice_between_two(self):
        assert rooting.merge_root_labels("A", None, "labeled") == ("A", True)
        assert rooting.merge_root_labels(None, "B", "labeled") == ("B", True)
        assert rooting.merge_root_labels("A", "B", "labeled") == (None, True)

    def test_outgroup_and_ingroup_take_the_label_of_a_side(self):
        assert rooting.merge_root_labels("A", "B", "outgroup", True) == ("A", True)
        assert rooting.merge_root_labels("A", "B", "outgroup", False) == ("B", True)
        assert rooting.merge_root_labels("A", "B", "ingroup", True) == ("B", True)
        assert rooting.merge_root_labels("A", "B", "ingroup", False) == ("A", True)

    def test_outgroup_needs_to_know_the_side(self):
        with pytest.raises(TwigError, match="give the outgroups with -r"):
            rooting.merge_root_labels("A", "B", "outgroup", None)

    def test_error_policy_stops_on_different_labels(self):
        with pytest.raises(TwigError, match="different labels"):
            rooting.merge_root_labels("A", "B", "error")

    def test_a_label_name_wins_if_either_branch_has_it(self):
        assert rooting.merge_root_labels("Test", "Ref", "Test") == ("Test", True)
        assert rooting.merge_root_labels("Ref", "Test", "Test") == ("Test", True)
        assert rooting.merge_root_labels("Ref", "Rev", "Test") == (None, True)

    @pytest.mark.parametrize("policy", ["bad name", "All", "9x", "prune"])
    def test_other_policies_are_rejected(self, policy):
        with pytest.raises(TwigError, match="unknown root edge policy"):
            rooting.check_root_edge_policy(policy)


class TestPlanRootMerge:
    def test_equal_labels_at_the_root_of_the_designed_tree(self):
        tree = treeio.read_tree(LABELED)
        plan = rooting.plan_root_merge(tree)
        assert tree[plan.kept].name == "O" and not tree[plan.dissolved].is_leaf()
        assert plan.label == "Reference" and not plan.conflict
        assert plan.labels == ("Reference", "Reference")

    def test_an_unrooted_tree_needs_no_merge(self):
        assert rooting.plan_root_merge(treeio.read_tree("((M1,(A,(B,C))),(M2,(D,E)),O);")) is None

    def test_two_tips_at_the_root_need_no_merge(self):
        assert rooting.plan_root_merge(treeio.read_tree("(A,B);")) is None

    def test_a_bifurcating_child_is_dissolved_first(self):
        tree = treeio.read_tree("((A,B,C),(D,E));")
        plan = rooting.plan_root_merge(tree)
        assert sorted(n.name for n in tree[plan.dissolved].children) == ["D", "E"]
        assert sorted(l.name for l in tree[plan.kept].children) == ["A", "B", "C"]

    def test_the_only_internal_child_is_dissolved(self):
        tree = treeio.read_tree("(O,(A,B,C,D));")
        plan = rooting.plan_root_merge(tree)
        assert tree[plan.kept].name == "O" and len(tree[plan.dissolved].children) == 4

    def test_conflicting_labels_follow_the_policy(self):
        tree = treeio.read_tree("(A{X},(B{Y},C{Y}){Z});")
        assert rooting.plan_root_merge(tree, "unlabel").label is None
        assert rooting.plan_root_merge(tree, "unlabel").conflict
        assert rooting.plan_root_merge(tree, "labeled").label is None
        assert rooting.plan_root_merge(tree, "Z").label == "Z"

    def test_without_outgroups_the_smaller_side_is_the_outgroup(self):
        tree = treeio.read_tree("(A{X},(B{Y},C{Y}){Z});")
        assert rooting.plan_root_merge(tree, "outgroup").label == "X"
        assert rooting.plan_root_merge(tree, "ingroup").label == "Z"

    def test_given_the_outgroup_side_it_decides(self):
        tree = treeio.read_tree("(A{X},(B{Y},C{Y}){Z});")
        bc = next(n for n in tree if not n.is_leaf() and n.up is not None)
        assert rooting.plan_root_merge(tree, "outgroup", outgroup_side=bc.idx).label == "Z"
        assert rooting.plan_root_merge(tree, "ingroup", outgroup_side=bc.idx).label == "X"

    def test_equal_sides_need_the_outgroup(self):
        tree = treeio.read_tree("((A{X},B{X}){X},(C{Y},D{Y}){Y});")
        with pytest.raises(TwigError, match="equally large"):
            rooting.plan_root_merge(tree, "outgroup")

    def test_error_policy(self):
        with pytest.raises(TwigError, match="different labels"):
            rooting.plan_root_merge(treeio.read_tree("(A{X},(B,C){Z});"), "error")


class TestRootOn:
    def test_labels_stay_with_their_branches(self):
        tree = treeio.read_tree(DISTINCT)
        before = labels_of(tree)
        bc = next(n for n in tree if n.is_leaf() is False and {l.name for l in n.traverse() if l.is_leaf()} == {"B", "C"})
        rooted, notes = rooting.root_on(tree, bc.idx)
        after = labels_of(rooted)
        for key, label in before.items():
            if key != ("O",):
                assert after[key] == label, key

    def test_the_branch_that_is_cut_has_its_label_on_both_halves(self):
        tree = treeio.read_tree(DISTINCT)
        bc = next(n for n in tree if not n.is_leaf() and {l.name for l in n.traverse() if l.is_leaf()} == {"B", "C"})
        rooted, _ = rooting.root_on(tree, bc.idx)
        halves = splits(rooted)[("B", "C")]
        assert halves == [("T", 0.5), ("T", 0.5)]
        assert sorted(c.label for c in rooted.treenode.children) == ["T", "T"]

    def test_the_two_old_root_branches_become_one_with_the_label_of_the_policy(self):
        tree = treeio.read_tree(DISTINCT)
        c = next(n for n in tree if n.name == "C")
        rooted, notes = rooting.root_on(tree, c.idx, "unlabel")
        assert splits(rooted)[("O",)] == [(None, 2.0)]
        assert len(notes) == 1 and "X" in notes[0] and "Z" in notes[0] and "no label" in notes[0]
        rooted, notes = rooting.root_on(tree, c.idx, "X")
        assert splits(rooted)[("O",)] == [("X", 2.0)]
        rooted, notes = rooting.root_on(tree, c.idx, "ingroup")
        assert splits(rooted)[("O",)] == [("X", 2.0)]       # the outgroup is on the side of C, so the ingroup side is O
        rooted, notes = rooting.root_on(tree, c.idx, "outgroup")
        assert splits(rooted)[("O",)] == [("Z", 2.0)]

    # two internal children at the root, written in both orders, because toytree keeps the branch of one of them by order
    TWO_SIDES = ("((A{Test}:1,B:1){T}:2,(C{Ref}:1,D:1){R}:3);", "((C{Ref}:1,D:1){R}:3,(A{Test}:1,B:1){T}:2);")

    @pytest.mark.parametrize("newick", TWO_SIDES)
    @pytest.mark.parametrize("target, near, far", [("A", "T", "R"), ("B", "T", "R"), ("C", "R", "T"), ("D", "R", "T")])
    @pytest.mark.parametrize("policy", ["unlabel", "labeled", "outgroup", "ingroup", "T", "R"])
    def test_the_merged_branch_follows_the_policy_whatever_the_order_of_the_root_children(self, newick, target, near, far, policy):
        tree = treeio.read_tree(newick)
        node = next(n for n in tree if n.name == target)
        rooted, notes = rooting.root_on(tree, node.idx, policy)
        expected = {"unlabel": None, "labeled": None, "outgroup": near, "ingroup": far}.get(policy, policy)
        assert splits(rooted)[("A", "B")] == [(expected, 5.0)]
        assert len(notes) == 1 and notes[0].endswith(f"the merged branch has {expected or 'no label'}")
        before, after = labels_of(tree), labels_of(rooted)
        assert all(after[key] == label for key, label in before.items()), (before, after)

    @pytest.mark.parametrize("newick", [
        "((A{Test}:1,B:1):2,(C{Ref}:1,D:1){Ref}:3);",
        "((C{Ref}:1,D:1){Ref}:3,(A{Test}:1,B:1):2);",
    ])
    def test_a_label_that_the_policy_removes_does_not_come_back(self, newick):
        # rooting on C kept Ref on the merged branch when the side of C came second, although the note said no label
        tree = treeio.read_tree(newick)
        c = next(n for n in tree if n.name == "C")
        rooted, notes = rooting.root_on(tree, c.idx, "unlabel")
        assert splits(rooted)[("A", "B")] == [(None, 5.0)]
        assert notes == ["the two branches at the old root had the labels Ref and no label, and the merged branch has no label"]
        rooted, _ = rooting.root_on(tree, c.idx, "labeled")
        assert splits(rooted)[("A", "B")] == [("Ref", 5.0)]
        assert labels_of(rooted)[("A",)] == "Test" and labels_of(rooted)[("C",)] == "Ref"

    def test_equal_labels_at_the_old_root_give_no_note(self):
        tree = treeio.read_tree(LABELED)
        c = next(n for n in tree if n.name == "C")
        rooted, notes = rooting.root_on(tree, c.idx)
        assert notes == [] and labels_of(rooted)[("O",)] == "Reference"

    def test_total_branch_length_does_not_change(self):
        tree = treeio.read_tree("((A:1,B:2):3,(C:4,D:5):6);")
        a = next(n for n in tree if n.name == "A")
        rooted, _ = rooting.root_on(tree, a.idx)
        assert sum(n.dist for n in rooted) == pytest.approx(21)
        assert treeio.has_dists(rooted)

    def test_rooting_where_the_tree_is_already_rooted_changes_nothing(self):
        tree = treeio.read_tree(DISTINCT)
        o = next(n for n in tree if n.name == "O")
        rooted, notes = rooting.root_on(tree, o.idx)
        assert notes == [] and splits(rooted) == splits(tree)
        assert [n.name for n in rooted] == [n.name for n in tree]

    def test_the_input_tree_is_not_changed(self):
        tree = treeio.read_tree(DISTINCT)
        c = next(n for n in tree if n.name == "C")
        rooting.root_on(tree, c.idx)
        assert splits(tree)[("O",)] == [("X", 1.0), ("Z", 1.0)]

    def test_nodes_made_by_toytree_have_the_label_feature(self):
        tree = treeio.read_tree(DISTINCT)
        c = next(n for n in tree if n.name == "C")
        rooted, _ = rooting.root_on(tree, c.idx)
        assert all(hasattr(n, "label") for n in rooted) and "label" in rooted.edge_features


class TestReroot:
    def test_roots_on_a_clade_and_says_which_side_it_is(self):
        tree = treeio.read_tree(DISTINCT)
        result = rooting.reroot(tree, ["B,C"])
        assert result.description == "rooted on B,C"
        side = result.tree[result.outgroup_side]
        assert {l.name for l in side.traverse() if l.is_leaf()} == {"B", "C"}
        assert {c.idx for c in result.tree.treenode.children} == {result.outgroup_side, next(
            c.idx for c in result.tree.treenode.children if c.idx != result.outgroup_side)}

    def test_the_outgroup_can_be_the_larger_side(self):
        tree = treeio.read_tree("((A,B),((C,D),E));")
        result = rooting.reroot(tree, ["C,D,E"])
        side = result.tree[result.outgroup_side]
        assert {l.name for l in side.traverse() if l.is_leaf()} == {"C", "D", "E"}

    def test_an_outgroup_that_is_not_a_clade_is_skipped_for_the_next_one(self):
        tree = treeio.read_tree(DISTINCT)
        result = rooting.reroot(tree, ["A,C", "D,E"])
        assert result.description == "rooted on D,E"

    def test_tips_missing_from_the_tree_are_ignored(self):
        tree = treeio.read_tree(DISTINCT)
        result = rooting.reroot(tree, ["Q,B,C"])
        side = result.tree[result.outgroup_side]
        assert {l.name for l in side.traverse() if l.is_leaf()} == {"B", "C"}

    def test_no_usable_outgroup_is_an_error(self):
        tree = treeio.read_tree(DISTINCT)
        with pytest.raises(TwigError, match="none of the outgroups"):
            rooting.reroot(tree, ["Q", "A,C"])
        with pytest.raises(TwigError, match="none of the outgroups"):
            rooting.reroot(tree, ["O,M1,A,B,C,M2,D,E"])

    def test_empty_outgroups_are_an_error(self):
        with pytest.raises(TwigError, match="empty outgroup"):
            rooting.parse_outgroups([","])
        with pytest.raises(TwigError, match="no outgroups"):
            rooting.parse_outgroups([])

    def test_species_can_be_the_outgroup_when_a_mapping_is_given(self):
        tree = treeio.read_tree("((A_1,A_2),((B_1,B_2),C_1));")
        species = lambda name: name.split("_")[0]
        result = rooting.reroot(tree, ["A"], species_of=species)
        side = result.tree[result.outgroup_side]
        assert {l.name for l in side.traverse() if l.is_leaf()} == {"A_1", "A_2"}
        with pytest.raises(TwigError, match="none of the outgroups"):
            rooting.reroot(tree, ["A"])

    def test_unrooted_trees_can_be_rooted(self):
        tree = treeio.read_tree("((M1,(A,(B,C))),(M2,(D,E)),O);")
        result = rooting.reroot(tree, ["O"])
        assert [c.name for c in result.tree.treenode.children if c.is_leaf()] == ["O"]
        assert result.tree[result.outgroup_side].name == "O"

    def test_the_policy_is_checked_even_if_nothing_needs_merging(self):
        with pytest.raises(TwigError, match="unknown root edge policy"):
            rooting.reroot(treeio.read_tree(SPECIES), ["O"], policy="bad name")

    def test_the_outgroup_side_makes_the_outgroup_policy_work(self):
        tree = treeio.read_tree(DISTINCT)
        result = rooting.reroot(tree, ["O"])
        assert result.tree[result.outgroup_side].name == "O"
        assert rooting.plan_root_merge(result.tree, "outgroup", result.outgroup_side).label == "X"
        assert rooting.plan_root_merge(result.tree, "ingroup", result.outgroup_side).label == "Z"
        result = rooting.reroot(tree, ["M1,A,B,C,M2,D,E"])
        assert rooting.plan_root_merge(result.tree, "outgroup", result.outgroup_side).label == "Z"


class TestDeepTrees:
    def test_a_deep_tree_can_be_rerooted(self):
        text = "T0:1"
        for i in range(1, 400):
            text = f"({text},T{i}:1):1"
        tree = treeio.read_tree(text.rsplit(":", 1)[0] + ";")
        result = rooting.reroot(tree, ["T200"])
        assert result.description == "rooted on T200"
        assert result.tree[result.outgroup_side].name == "T200"
        assert sum(n.dist for n in result.tree) == pytest.approx(sum(n.dist for n in tree))
