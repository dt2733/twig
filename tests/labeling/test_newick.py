"""Tests for the Newick serializer (twig.labeling.core.newick)."""

import pytest

pytest.importorskip("toytree")

from twig.utils import TwigError
from twig.labeling.core import treeio
from twig.labeling.core.model import RootMerge
from twig.labeling.core.newick import braces_label, to_newick
from conftest import LABELED, SPECIES

# the concordant gene of the design, written as unrooted with HyPhy and no lengths
UNROOTED = (
    "(O{Reference},(M1{Reference},(A{Test},(B{Test},C{Rev}){Test})){Reference},"
    "(M2{Reference},(D{Test},E{Test})){Reference});"
)


def write(newick, **kwargs):
    tree = treeio.read_tree(newick)
    kwargs.setdefault("dists", treeio.has_dists(tree))
    return to_newick(tree, **kwargs)


def root_merge(tree, label):
    """Dissolve the first internal child of the root, keeping the other one."""
    children = tree.treenode.children
    dissolved = next(c for c in children if not c.is_leaf())
    kept = next(c for c in children if c is not dissolved)
    return RootMerge(dissolved.idx, kept.idx, (dissolved.label, kept.label), label, dissolved.label != kept.label)


class TestRoundTrip:
    @pytest.mark.parametrize("newick", [
        SPECIES,
        LABELED,
        "((A:0.1,B:0.2):0.3,C:0.4);",
        "((A{Test}:0.1,B:0.2){Rev}:0.3,C{Test}:0.4);",
        "(A,B,C,D);",
        "(A,B,(C,D,E){Test});",
        "A;",
    ])
    def test_reading_then_writing_gives_the_same_text(self, newick):
        assert write(newick) == newick

    def test_missing_lengths_are_not_invented(self):
        tree = treeio.read_tree(SPECIES)
        assert all(n.dist == 1.0 for n in tree if not n.is_root())      # what toytree fills in
        assert ":" not in to_newick(tree, dists=treeio.has_dists(tree))

    def test_lengths_are_formatted_compactly(self):
        assert write("(A:1.0,B:0.123456789012345);") == "(A:1,B:0.123456789012);"
        assert to_newick(treeio.read_tree("(A:0.5,B:0.25);"), dist_fmt="%.3f") == "(A:0.500,B:0.250);"

    def test_the_root_never_gets_a_label_or_length(self):
        assert write("(A:1,B:1):5;") == "(A:1,B:1);"

    def test_internal_names_are_not_written_unless_asked(self):
        assert write("(A,B,(C,D,E)X);") == "(A,B,(C,D,E));"
        assert write("(A,B,(C,D,E)X);", internal_text=lambda n: n.name) == "(A,B,(C,D,E)X);"

    def test_deep_trees_do_not_need_recursion(self):
        n = 1200      # deeper than Python's default recursion limit of 1000
        newick = "(" * n + "t0" + "".join(f",t{i + 1})" for i in range(n)) + ";"
        tree = treeio.read_tree(newick)
        assert to_newick(tree, dists=False) == newick


class TestHooks:
    def test_tip_text_for_sanitized_names(self):
        tree = treeio.read_tree("(A-b{Test},C.d);")
        text = to_newick(tree, dists=False, tip_text=lambda n: n.name.replace("-", "_").replace(".", "_"))
        assert text == "(A_b{Test},C_d);"

    def test_internal_text_for_support_or_names(self):
        tree = treeio.read_tree("((A,B)95,(C,D)80);")
        text = to_newick(tree, dists=False, internal_text=lambda n: "" if n.support != n.support else f"{n.support:g}")
        assert text == "((A,B)95,(C,D)80);"
        text = to_newick(tree, dists=False, internal_text=lambda n: f"twn{n.idx}")
        assert text == "((A,B)twn4,(C,D)twn5)twn6;"

    def test_label_text_for_other_programs(self):
        tree = treeio.read_tree("(A{Test},B{Reference},C);")
        codes = {"Test": " #1", "Reference": " #0"}
        text = to_newick(tree, dists=False, label_text=lambda lab: codes.get(lab, ""))
        assert text == "(A #1,B #0,C);"

    def test_braces_label(self):
        assert braces_label("Test") == "{Test}" and braces_label(None) == "" and braces_label("") == ""


class TestRootMerge:
    def test_concordant_gene_is_written_as_the_design_expects(self):
        tree = treeio.read_tree(LABELED)
        assert to_newick(tree, dists=False, root_merge=root_merge(tree, "Reference")) == UNROOTED

    def test_the_merged_label_and_summed_length_go_to_the_kept_branch(self):
        tree = treeio.read_tree("(O{Reference}:1,((A{Test}:1,B{Test}:1){Test}:2,C{Reference}:3){Rev}:4);")
        merge = root_merge(tree, None)
        assert merge.conflict
        assert to_newick(tree, root_merge=merge) == "(O:5,(A{Test}:1,B{Test}:1){Test}:2,C{Reference}:3);"
        merge = RootMerge(merge.dissolved, merge.kept, merge.labels, "Rev", True)
        assert to_newick(tree, root_merge=merge) == "(O{Rev}:5,(A{Test}:1,B{Test}:1){Test}:2,C{Reference}:3);"

    def test_the_dissolved_child_is_replaced_in_place(self):
        tree = treeio.read_tree("(((A,B),(C,D)),E);")
        children = tree.treenode.children
        merge = RootMerge(children[0].idx, children[1].idx, (None, None), None, False)
        assert to_newick(tree, dists=False, root_merge=merge) == "((A,B),(C,D),E);"
        tree = treeio.read_tree("(E,((A,B),(C,D)));")
        children = tree.treenode.children
        merge = RootMerge(children[1].idx, children[0].idx, (None, None), None, False)
        assert to_newick(tree, dists=False, root_merge=merge) == "(E,(A,B),(C,D));"

    def test_the_kept_branch_can_be_internal(self):
        tree = treeio.read_tree("(((A:1,B:1){Test}:1,C:1){Rev}:2,((D:1,E:1):1,F:1){Test}:3);")
        merge = root_merge(tree, "Test")
        text = to_newick(tree, root_merge=merge)
        # the first internal child is dissolved, and the second keeps its children and gets the summed length
        assert text == "((A:1,B:1){Test}:1,C:1,((D:1,E:1):1,F:1){Test}:5);"

    def test_errors(self):
        tree = treeio.read_tree("(A,B,C);")
        with pytest.raises(TwigError, match="exactly two children"):
            to_newick(tree, root_merge=RootMerge(0, 1, (None, None), None, False))
        tree = treeio.read_tree("(A,(B,C));")
        a = tree.get_nodes("A")[0]
        other = next(c for c in tree.treenode.children if c is not a)
        with pytest.raises(TwigError, match="cannot dissolve a root child that is a tip"):
            to_newick(tree, root_merge=RootMerge(a.idx, other.idx, (None, None), None, False))


class TestSafety:
    @pytest.mark.parametrize("name", ["A b", "A,b", "A(b", "A:b", "A;b", "A[b", "A{b", "A'b"])
    def test_names_that_cannot_be_written_are_an_error(self, name):
        tree = treeio.read_tree("(A,B);")
        tree[0].name = name
        with pytest.raises(TwigError, match="cannot write name"):
            to_newick(tree, dists=False)

    @pytest.mark.parametrize("label", ["a b", "a,b", "a}b", "a:b", "a(b"])
    def test_labels_that_cannot_be_written_are_an_error(self, label):
        tree = treeio.read_tree("(A,B);")
        tree[0].label = label
        with pytest.raises(TwigError, match="cannot write label"):
            to_newick(tree, dists=False)

    def test_nan_lengths(self):
        tree = treeio.read_tree("(A:1,B:1);")
        tree[0]._dist = float("nan")
        with pytest.raises(TwigError, match="NaN"):
            to_newick(tree)
