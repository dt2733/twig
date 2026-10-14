"""Tests for reading trees and labels (twig.labeling.core.treeio)."""

import random

import pytest

pytest.importorskip("toytree")

from twig.utils import TwigError
from twig.labeling.core import treeio
from conftest import LABELED, SPECIES


def labels(tree):
    """Return {name or idx: label} for the labeled nodes of a tree."""
    return {(n.name or f"#{n.idx}"): n.label for n in tree if n.label}


###########################################################################
# splitting files into trees
###########################################################################

class TestIterTrees:
    def test_single_line_multi_line_and_several_per_line(self):
        text = "(A,B);\n(C,\n  D\n);\n(E,F);(G,H);\n"
        trees = [t for _, t in treeio.iter_trees_text(text)]
        assert trees == ["(A,B);", "(C,D);", "(E,F);", "(G,H);"]

    def test_blank_lines_comment_lines_and_crlf(self):
        text = "\r\n# a comment\r\n(A,B);\r\n\r\n#another\r\n(C,D);\r\n"
        assert [t for _, t in treeio.iter_trees_text(text)] == ["(A,B);", "(C,D);"]

    def test_locus_prefix_and_paml_header(self):
        text = "2 2\nOG0001\t(A,B);\nOG0002\t(C,D);\n"
        assert list(treeio.iter_trees_text(text)) == [("OG0001", "(A,B);"), ("OG0002", "(C,D);")]
        assert list(treeio.iter_trees_text("(A,B);")) == [(None, "(A,B);")]

    def test_semicolons_inside_comments_do_not_end_a_tree(self):
        text = "(A[&note=a;b],B);\n"
        assert [t for _, t in treeio.iter_trees_text(text)] == ["(A[&note=a;b],B);"]

    def test_paml_marks_are_not_comment_lines_inside_a_tree(self):
        text = "(A,\n#1\nB);\n"
        assert [t for _, t in treeio.iter_trees_text(text)] == ["(A,#1B);"]

    def test_a_comment_after_the_last_tree_is_fine(self):
        assert [t for _, t in treeio.iter_trees_text("(A,B);\n[end of file]\n")] == ["(A,B);"]

    def test_a_tree_over_many_lines_is_read_in_linear_time(self):
        import time
        n = 20000
        text = "(\n" + ",\n".join(f"t{i}:1" for i in range(n)) + "\n);\n"
        start = time.perf_counter()
        trees = list(treeio.iter_trees_text(text))
        assert len(treeio.tip_names(trees[0][1])) == n
        assert time.perf_counter() - start < 5

    def test_errors(self):
        with pytest.raises(TwigError, match="NEXUS"):
            list(treeio.iter_trees_text("#NEXUS\nbegin trees;\n"))
        with pytest.raises(TwigError, match="no final ';'"):
            list(treeio.iter_trees_text("(A,B)"))
        with pytest.raises(TwigError, match="tree 2 has quoted labels"):
            list(treeio.iter_trees_text("(A,B);\n('C d',E);\n"))

    def test_file(self, tmp_path):
        path = tmp_path / "trees.nwk"
        path.write_text("(A,B);\n(C,D);\n")
        assert [t for _, t in treeio.iter_trees(path)] == ["(A,B);", "(C,D);"]
        with pytest.raises(TwigError, match="cannot read tree file"):
            list(treeio.iter_trees(tmp_path / "missing.nwk"))
        bad = tmp_path / "bad.nwk"
        bad.write_bytes(b"\xff\xfe(A,B);")
        with pytest.raises(TwigError, match="not valid UTF-8"):
            list(treeio.iter_trees(bad))

    def test_a_byte_order_mark_is_ignored(self, tmp_path):
        # the mark stayed in the text, so toytree read the whole tree as one tip
        path = tmp_path / "bom.nwk"
        path.write_bytes(b"\xef\xbb\xbf((A,B),C);\n")
        assert list(treeio.iter_trees(path)) == [(None, "((A,B),C);")]
        assert list(treeio.iter_trees_text("\ufeff((A,B),C);\n")) == [(None, "((A,B),C);")]
        assert [n.name for n in treeio.read_tree(path) if n.is_leaf()] == ["A", "B", "C"]
        assert treeio.tip_names("\ufeff((A,B),C);") == ["A", "B", "C"]

    @pytest.mark.parametrize("text", [
        "#NEXUS\nbegin trees;\n",
        "\ufeff#NEXUS\nbegin trees;\ntree t1 = ((A,B),C);\nend;\n",
        "begin trees;\ntree t1 = ((A,B),C);\nend;\n",
        "tree t1 = ((A,B),C);\n",
        "TREE * t1 = [&U] ((A,B),C);\n",
    ])
    def test_nexus_is_refused_with_a_byte_order_mark_or_without_its_header(self, text):
        # these were split into trees such as 'treet1=((A,B),C);', which toytree read as a single tip
        with pytest.raises(TwigError, match="NEXUS format is not supported"):
            list(treeio.iter_trees_text(text))

    def test_whitespace_between_tokens_is_removed(self):
        text = "( A , B ) ;\n(A #1 , B);\n((A:1, B:2) 95 : 0.3, C { Test } :1);\n(A[a\nnote],\nB);\n"
        trees = [t for _, t in treeio.iter_trees_text(text)]
        assert trees == ["(A,B);", "(A#1,B);", "((A:1,B:2)95:0.3,C{Test}:1);", "(A[a note],B);"]

    @pytest.mark.parametrize("text, words", [
        ("((Homo sapiens:1,B:1):1,C:1);", "Homo sapiens"),
        ("((Homo\tsapiens,B),C);", "Homo sapiens"),
        ("((Homo\nsapiens,B),C);", "Homo sapiens"),
        ("((A{Te st},B),C);", "Te st"),
        ("((A,B)9 5,C);", "9 5"),
    ])
    def test_a_space_inside_a_name_or_label_is_an_error(self, text, words):
        # the whitespace was removed, so Homo sapiens became Homosapiens without a word
        with pytest.raises(TwigError, match=f"tree 1 has a space in '{words}'. Names and labels cannot contain spaces"):
            list(treeio.iter_trees_text(text))

    def test_a_locus_name_can_be_followed_by_a_comment_before_the_tree(self):
        assert list(treeio.iter_trees_text("L1\t[&R] ((A,B),C);\n")) == [("L1", "[&R]((A,B),C);")]
        # a locus named tree is not taken for NEXUS
        assert list(treeio.iter_trees_text("tree\t((A[&&NHX:S=x],B),C);\n")) == [("tree", "((A[&&NHX:S=x],B),C);")]


###########################################################################
# tip names without toytree
###########################################################################

class TestTipNames:
    def test_names_in_order_skipping_labels_support_and_lengths(self):
        assert treeio.tip_names(LABELED) == ["O", "M1", "A", "B", "C", "M2", "D", "E"]
        assert treeio.tip_names("((A:0.1,B:0.2)95:0.3[&x=1],(C,D)Node7);") == ["A", "B", "C", "D"]
        assert treeio.tip_names("((A #1,B $2) #3,C);") == ["A", "B", "C"]
        assert treeio.tip_names("A;") == ["A"]

    @pytest.mark.parametrize("newick, message", [
        ("(,B);", "tip without a name"),
        ("(A,);", "tip without a name"),
        ("(A:0.1,:0.2);", "tip without a name"),
        ("((A,B);", "unbalanced"),
        ("(A,B));", "unbalanced"),
        ("('A b',C);", "quoted"),
        ("(A,B)[unclosed;", "unexpected character"),
        ("(A,B){unclosed;", "unexpected character"),
        ("(A#x,B);", "unexpected character"),
        ("A,", "tip without a name"),
        ("treet1=((A,B),C);", "unexpected text 'treet1=' at position 0 of the tree, before its first"),
        ("locus1((A,B),C);", "unexpected text 'locus1' at position 0 of the tree, before its first"),
        ("((A{Test}x,B),C);", "unexpected text 'x' at position 9 of the tree, after a label or comment"),
        ("((Homo sapiens,B),C);", "space in 'Homo sapiens' at position 7"),
        ("((A{Te st},B),C);", "the label '{Te st}' at position 3 of the tree has a space"),
    ])
    def test_errors(self, newick, message):
        with pytest.raises(TwigError, match=message):
            treeio.tip_names(newick)

    def test_agrees_with_toytree_on_random_trees(self):
        import toytree
        rng = random.Random(7)
        for _ in range(25):
            tree = toytree.rtree.unittree(rng.randint(2, 40), seed=rng.randint(0, 10**6))
            newick = tree.write(internal_labels=None)
            assert sorted(treeio.tip_names(newick)) == sorted(tree.get_tip_labels())

    def test_does_not_import_toytree(self):
        import subprocess, sys
        code = "import sys; from twig.labeling.core import treeio; treeio.tip_names('(A,B);'); assert 'toytree' not in sys.modules"
        assert subprocess.run([sys.executable, "-c", code]).returncode == 0


###########################################################################
# reading labels
###########################################################################

class TestReadTree:
    def test_hyphy_braces_on_tips_and_internal_nodes(self):
        tree = treeio.read_tree(LABELED)
        assert [n.label for n in tree if n.is_leaf()] == ["Reference", "Reference", "Test", "Test", "Rev", "Reference", "Test", "Test"]
        internal = {n.idx: n.label for n in tree if not n.is_leaf()}
        assert sorted(v for v in internal.values() if v) == ["Reference", "Reference", "Reference", "Test"]
        assert tree.treenode.label is None            # the root has no branch

    def test_every_node_has_a_label_feature_registered_on_edges(self):
        tree = treeio.read_tree(SPECIES)
        assert all(hasattr(n, "label") and n.label is None for n in tree)
        assert "label" in tree.edge_features

    def test_copy_keeps_labels_and_the_edge_feature(self):
        tree = treeio.read_tree(LABELED)
        copy = tree.copy()
        assert labels(copy) == labels(tree) and "label" in copy.edge_features

    def test_labels_follow_their_branches_when_the_tree_is_rerooted(self):
        # the tree was rooted on E before, which is where it is rooted already, so nothing moved
        from twig.labeling.core import rooting
        tree = treeio.read_tree("(((A{Test}:1,B:1):1,(C:1,D:1){Rev}:1):1,E:2);")
        a = next(n for n in tree if n.name == "A")
        rooted, _ = rooting.root_on(tree, a.idx)
        below = rooting.tip_name_sets(rooted)
        assert sorted(",".join(sorted(below[c.idx])) for c in rooted.treenode.children) == ["A", "B,C,D,E"]
        # A's branch is cut by the new root and keeps its label on both halves, and Rev stays above C and D
        found = sorted((",".join(sorted(below[n.idx])), n.label) for n in rooted if n.label)
        assert found == [("A", "Test"), ("B,C,D,E", "Test"), ("C,D", "Rev")]

    def test_label_text_is_kept_exactly_and_never_converted_to_numbers(self):
        tree = treeio.read_tree("(A{nan}:1,B{1e5}:1,C{Inf}:1);")
        assert [n.label for n in tree if n.is_leaf()] == ["nan", "1e5", "Inf"]

    def test_internal_name_with_a_label_still_gives_support(self):
        # toytree alone treats the node as a name because of the braces, and then loses all supports
        tree = treeio.read_tree("((A:1,B:1)95{Test}:1,(C:1,D:1)80:1);")
        supports = sorted(n.support for n in tree if not n.is_leaf() and n.support == n.support)
        assert supports == [80.0, 95.0]
        assert [n.label for n in tree if n.label] == ["Test"]

    def test_named_internal_nodes_are_kept_as_names(self):
        tree = treeio.read_tree("((A,B)Node3{Test},C);")
        assert [(n.name, n.label) for n in tree if not n.is_leaf() and n.name] == [("Node3", "Test")]

    def test_iqtree_double_support(self):
        tree = treeio.read_tree("((A:1,B:1)80.5/95:1,(C:1,D:1)99.1/100:1);")
        values = sorted((n.support_alrt, n.support) for n in tree if not n.is_leaf() and n.name == "")
        assert [v for v in values if v[0] == v[0]] == [(80.5, 95.0), (99.1, 100.0)]
        assert {"support", "support_alrt"} <= tree.edge_features

    def test_numeric_support_is_not_confused_with_names_when_mixed(self):
        tree = treeio.read_tree("((A,B)95,(C,D)Node2);")
        assert sorted(n.name for n in tree if not n.is_leaf() and n.name) == ["95", "Node2"]

    def test_lengths_flag(self):
        assert treeio.has_dists(treeio.read_tree("((A:1,B:2):1,C:3);")) is True
        assert treeio.has_dists(treeio.read_tree(SPECIES)) is False
        assert treeio.has_dists(treeio.read_tree(LABELED)) is False
        assert treeio.has_dists(treeio.read_tree(LABELED).copy()) is False
        # the root has no branch, so a length on it alone does not count
        assert treeio.has_dists(treeio.read_tree("((A,B),C):0;")) is False

    def test_branches_without_a_length_are_marked(self):
        # as in an ASTRAL species tree, which has internal lengths and none for the tips
        tree = treeio.read_tree("(O,((M1,(A,(B,C)1:0.5)0.9:0.3)1:1.2,(M2,(D,E)1:0.8)1:0.4));")
        assert treeio.has_dists(tree) is True
        below = {n.idx: "".join(sorted(l.name for l in n.iter_leaves())) for n in tree}
        given = sorted(below[n.idx] for n in tree if n.up is not None and treeio.has_dist(n))
        assert given == ["ABC", "ABCM1", "BC", "DE", "DEM2"]
        assert "has_dist" in tree.edge_features
        assert all(treeio.has_dist(n) for n in treeio.read_tree("((A:1,B:2):1,C:3);"))

    def test_root_label_is_dropped(self):
        parsed = treeio.parse_newick("(A{Test},B{Test}){Test};")
        assert parsed.tree.treenode.label is None and parsed.root_label == "Test"

    def test_from_file_and_errors(self, tmp_path):
        path = tmp_path / "sp.nwk"
        path.write_text(LABELED + "\n")
        assert labels(treeio.read_tree(path)) == labels(treeio.read_tree(LABELED))
        assert labels(treeio.read_tree(str(path))) == labels(treeio.read_tree(LABELED))
        path.write_text("(A,B);\n(C,D);\n")
        with pytest.raises(TwigError, match="found 2 trees but expected one"):
            treeio.read_tree(path)
        with pytest.raises(TwigError, match="tree file not found"):
            treeio.read_tree(tmp_path / "missing.nwk")
        with pytest.raises(TwigError, match="parse"):
            treeio.read_tree("(A,B;")

    def test_very_deep_tree(self):
        n = 1500
        newick = "(" * n + "t0" + "".join(f",t{i + 1})" for i in range(n)) + ";"
        tree = treeio.read_tree(newick)
        assert tree.ntips == n + 1


class TestLabelSyntaxes:
    def test_paml_branch_marks_with_and_without_lengths(self):
        tree = treeio.read_tree("((A #1,B:0.2 #2):0.3,C:0.4#1);", paml_codes={1: "Test", 2: "Reference"})
        by_name = {n.name: n.label for n in tree if n.is_leaf()}
        assert by_name == {"A": "Test", "B": "Reference", "C": "Test"}
        assert treeio.has_dists(tree)

    def test_paml_marks_without_names_use_the_number(self):
        tree = treeio.read_tree("((A #1,B),C);")
        assert [n.label for n in tree if n.is_leaf()] == ["1", None, None]

    def test_paml_clade_marks_expand_like_codeml(self):
        # $1 labels the branch above the node and every branch below it; # beats $; an inner $ beats an outer one
        tree = treeio.read_tree("(((A,B #2) $1,(C,D $3) $4) $5,E);", paml_codes={1: "One", 2: "Two", 3: "Three", 4: "Four", 5: "Five"})
        by_name = {n.name: n.label for n in tree if n.is_leaf()}
        assert by_name == {"A": "One", "B": "Two", "C": "Four", "D": "Three", "E": None}
        clade_ab = tree.get_mrca_node("A", "B")
        assert clade_ab.label == "One"
        assert tree.get_mrca_node("A", "C").label == "Five"      # the branch above the $5 clade

    def test_nhx_and_figtree_label_keys(self):
        nhx = treeio.read_tree("(A[&&NHX:S=human:label=Test]:0.1,B[&&NHX:S=mouse]:0.2);")
        assert [(n.name, n.label) for n in nhx if n.is_leaf()] == [("A", "Test"), ("B", None)]
        assert nhx[0].S == "human"                 # other NHX attributes are kept as features
        fig = treeio.read_tree('(A[&label="Test"]:0.1,B[&rate=2.5]:0.2);')
        assert [(n.name, n.label) for n in fig if n.is_leaf()] == [("A", "Test"), ("B", None)]
        assert fig[1].rate == 2.5

    def test_beast_style_annotations_with_braces_and_commas_are_kept(self):
        tree = treeio.read_tree("(A[&height_95%_HPD={1.0,2.0},label=Test,rate=0.5]:1,B:1);")
        assert (tree[0].label, tree[0].rate) == ("Test", 0.5)
        assert tree[0].__dict__["height_95%_HPD"] == "{1.0,2.0}"

    def test_nhx_labels_are_refused_by_other_explicit_syntaxes(self):
        with pytest.raises(TwigError, match="NHX or FigTree syntax"):
            treeio.read_tree("(A[&label=Test],B);", label_syntax="braces")
        with pytest.raises(TwigError, match="empty label"):
            treeio.read_tree('(A[&label=""],B);')

    def test_a_missing_final_semicolon_and_quotes_in_direct_use(self):
        assert treeio.parse_newick("(A{Test},B)").tree[0].label == "Test"
        with pytest.raises(TwigError, match="quoted"):
            treeio.parse_newick("('A',B);")

    def test_trees_that_are_too_deep_are_a_clear_error(self, monkeypatch):
        import toytree

        def too_deep(*args, **kwargs):
            raise RecursionError()

        monkeypatch.setattr(toytree, "tree", too_deep)
        with pytest.raises(TwigError, match="too deep to parse"):
            treeio.read_tree("(A,B);")

    def test_plain_comments_are_ignored(self):
        tree = treeio.read_tree("(A[a comment]{Test},B[another]);")
        assert [(n.name, n.label) for n in tree if n.is_leaf()] == [("A", "Test"), ("B", None)]

    def test_edge_annotation_after_a_length_can_carry_a_label(self):
        tree = treeio.read_tree("(A:0.1[&label=Test],B:0.2);")
        assert [(n.name, n.label) for n in tree if n.is_leaf()] == [("A", "Test"), ("B", None)]

    def test_mixing_syntaxes_is_an_error_in_auto_mode(self):
        with pytest.raises(TwigError, match="mixes label syntaxes"):
            treeio.read_tree("(A{Test},B #1);")

    def test_explicit_syntax_refuses_the_others(self):
        assert treeio.read_tree("(A{Test},B);", label_syntax="braces")[0].label == "Test"
        with pytest.raises(TwigError, match="label syntax is 'paml'"):
            treeio.read_tree("(A{Test},B);", label_syntax="paml")
        with pytest.raises(TwigError, match="label syntax is 'braces'"):
            treeio.read_tree("(A #1,B);", label_syntax="braces")
        with pytest.raises(TwigError, match="unknown label syntax"):
            treeio.read_tree("(A,B);", label_syntax="beast")

    def test_none_removes_all_labels(self):
        tree = treeio.read_tree("(A{Test}:1,(B #2,C[&label=x]):1);", label_syntax="none")
        assert not any(n.label for n in tree)

    def test_a_node_with_two_annotations_is_a_clear_error(self):
        with pytest.raises(TwigError, match="more than one annotation"):
            treeio.read_tree("(A[&x=1]{Test},B);")
        with pytest.raises(TwigError, match="more than one annotation"):
            treeio.read_tree("(A{Test}[&x=1],B);")
        with pytest.raises(TwigError, match="more than one annotation"):
            treeio.read_tree("(A #1 $2,B);")
        # one annotation on the node and one on its edge is fine
        assert treeio.read_tree("(A{Test}:0.1[&x=1],B);")[0].label == "Test"

    @pytest.mark.parametrize("newick", [
        "((A{Test}:0.1,B:0.1){Test}:0.2{Ref},C:0.1);",
        "((A#1:0.1#2,B:0.1),C:0.1);",
        "((A[&label=Test]:0.1[&label=Ref],B:0.1),C:0.1);",
        "((A,B)$1:0.1#2,C);",
        "(A[&label=Test,label=Ref],B);",
    ])
    def test_two_labels_on_one_branch_are_a_clear_error(self, newick):
        # the label after the length won without a word
        with pytest.raises(TwigError, match="has more than one label"):
            treeio.read_tree(newick)

    def test_a_label_and_another_annotation_on_one_branch_are_fine(self):
        assert treeio.read_tree("(A[&x=1]:0.1{Test},B);")[0].label == "Test"
        assert treeio.read_tree("((A,B)95{Test}:0.2[&x=1],C);").get_mrca_node(0, 1).label == "Test"

    @pytest.mark.parametrize("newick, text", [
        ("((A{Test},B){Test}95,C);", "95"),
        ("((A{Test}x,B),C);", "x"),
        ("((A#1x,B),C);", "x"),
        ("((Pop#12a,B),C);", "a"),
        ("((A[&label=Test]x,B),C);", "x"),
        ("((A,B)[&x=1]95,C);", "95"),
    ])
    def test_text_right_after_a_label_or_comment_is_a_clear_error(self, newick, text):
        # toytree glued the text to the annotation, which ended in a ValueError or in a lost support value
        with pytest.raises(TwigError, match=f"unexpected text '{text}' at position \\d+ of the tree, after a label or comment"):
            treeio.parse_newick(newick)

    @pytest.mark.parametrize("newick", ["tree t1=((A,B),C);", "locus1 ((A,B),C);"])
    def test_text_before_the_tree_is_a_clear_error(self, newick):
        with pytest.raises(TwigError, match="before its first '\\('"):
            treeio.parse_newick(newick)

    def test_a_label_with_a_space_is_refused_from_text_as_from_files(self):
        # a file gave the label Test and a string gave 'Te st'
        with pytest.raises(TwigError, match="has a space"):
            treeio.parse_newick("((A{Te st},B),C);")
        with pytest.raises(TwigError, match="has a space"):
            treeio.read_tree("((A{Te st},B),C);")
        assert treeio.parse_newick("((A{ Test },B),C);").tree[0].label == "Test"

    def test_a_single_tip_and_a_leading_comment_are_fine(self):
        assert [n.name for n in treeio.parse_newick("A;").tree] == ["A"]
        assert [n.name for n in treeio.parse_newick("[&R]((A,B),C);").tree if n.is_leaf()] == ["A", "B", "C"]

    def test_a_label_number_that_toytree_changed_is_a_clear_error(self):
        tree = treeio.read_tree("((A,B),C);")
        tree[0].twiglab = "0]9"
        with pytest.raises(TwigError, match="cannot read the label of the node 'A'"):
            treeio._table_labels(tree, ["Test"])
        tree[0].twigpaml = "1]"
        with pytest.raises(TwigError, match="cannot read the label of the node 'A'"):
            treeio._paml_labels(tree, None)

    def test_empty_labels_and_quotes(self):
        with pytest.raises(TwigError, match="empty label"):
            treeio.read_tree("(A{},B);")
        with pytest.raises(TwigError, match="quoted"):
            treeio.read_tree("('A',B);")


class TestStripAndHelpers:
    def test_strip_labels_returns_a_copy_by_default(self):
        tree = treeio.read_tree(LABELED)
        stripped = treeio.strip_labels(tree)
        assert not labels(stripped) and labels(tree)
        assert "label" in stripped.edge_features
        treeio.strip_labels(tree, inplace=True)
        assert not labels(tree)

    def test_labeled_branches_and_node_label(self):
        tree = treeio.read_tree("(A{Test},(B,C){Rev});")
        found = treeio.labeled_branches(tree)
        assert sorted(found.values()) == ["Rev", "Test"]
        assert treeio.node_label(tree[0]) == "Test" and treeio.node_label(object()) is None

    def test_carry_over_dists(self):
        old = treeio.read_tree(SPECIES)
        new = treeio.carry_over_dists(old, old.mod.prune(0, 1, 2))
        assert treeio.has_dists(new) is False


class TestExampleData:
    """The example species trees of the design are shared by the tests and the documentation."""

    def test_species_tree_file(self, data_dir):
        assert (data_dir / "species.nwk").read_text().strip() == SPECIES
        assert treeio.species_names(data_dir / "species.nwk") == ["O", "M1", "A", "B", "C", "M2", "D", "E"]

    def test_labeled_species_tree_file_round_trips(self, data_dir):
        from twig.labeling.core.newick import to_newick
        path = data_dir / "species.labeled.nwk"
        assert path.read_text().strip() == LABELED
        tree = treeio.read_tree(path)
        assert to_newick(tree, dists=treeio.has_dists(tree)) == LABELED

    def test_names_data_species_tree(self, names_dir):
        assert len(treeio.species_names(names_dir / "species.nwk")) == 5


class TestSpeciesTrees:
    def test_species_names_without_toytree(self):
        assert treeio.species_names(LABELED) == ["O", "M1", "A", "B", "C", "M2", "D", "E"]

    def test_duplicate_tips_are_an_error(self):
        with pytest.raises(TwigError, match="tip names are not unique"):
            treeio.species_names("((A,B),(A,C));")
        with pytest.raises(TwigError, match="tip names are not unique"):
            treeio.read_species_tree("((A,B),(A,C));")

    def test_label_on_the_root_is_ignored_with_a_warning(self, capsys):
        from twig.utils import set_log_level
        set_log_level("WARNING")
        tree = treeio.read_species_tree("(A{Test},B{Test}){Test};")
        assert tree.treenode.label is None
        assert "ignoring the label 'Test' on the root" in capsys.readouterr().err

    def test_branch_ids_use_names_or_idx(self):
        tree = treeio.read_species_tree("((A,B)AB,(C,D));")
        ids = treeio.branch_ids(tree)
        assert ids[tree.get_mrca_node("A", "B").idx] == "AB"
        assert ids[tree.get_mrca_node("C", "D").idx] == f"S{tree.get_mrca_node('C', 'D').idx}"
        assert ids[0] == "A"
        assert tree.treenode.idx not in ids      # the root has no branch

    def test_duplicate_branch_ids_are_an_error(self):
        tree = treeio.read_tree("((A,B)X,(C,D)X);")
        with pytest.raises(TwigError, match="not unique"):
            treeio.branch_ids(tree)


def ladder(n):
    """Return a Newick string for a tree that is n levels deep."""
    text = "T0:1"
    for i in range(1, n):
        text = f"({text},T{i}:1):1"
    return text.rsplit(":", 1)[0] + ";"


class TestDeepTrees:
    """toytree copies and reroots trees with deepcopy, which fails on trees more than about 140 levels deep."""

    def test_depth_is_the_number_of_branches_on_the_longest_path(self):
        assert treeio.tree_depth(treeio.read_tree("A;")) == 0
        assert treeio.tree_depth(treeio.read_tree("(A,B);")) == 1
        assert treeio.tree_depth(treeio.read_tree(SPECIES)) == 5
        assert treeio.tree_depth(treeio.read_tree(ladder(50))) == 49

    def test_copying_a_deep_tree_works_and_the_limit_is_restored(self):
        import sys
        tree = treeio.read_tree(ladder(400))
        before = sys.getrecursionlimit()
        with pytest.raises(RecursionError):
            tree.copy()
        assert sys.getrecursionlimit() == before
        copy = treeio.copy_tree(tree)
        assert [n.name for n in copy] == [n.name for n in tree]
        assert sys.getrecursionlimit() == before

    def test_the_limit_is_restored_after_an_error_inside_the_block(self):
        import sys
        tree = treeio.read_tree(ladder(400))
        before = sys.getrecursionlimit()
        with pytest.raises(ValueError):
            with treeio.deep_recursion(tree):
                assert sys.getrecursionlimit() > before
                raise ValueError("inside")
        assert sys.getrecursionlimit() == before

    def test_labels_can_be_set_on_a_deep_tree(self):
        tree = treeio.read_tree(ladder(400))
        labeled = treeio.set_labels(tree, {0: "Test"})
        assert [n.label for n in labeled if n.label] == ["Test"] and not any(n.label for n in tree)
        assert not any(n.label for n in treeio.strip_labels(labeled))

    def test_a_tree_that_is_too_deep_is_an_error_and_not_a_crash(self, monkeypatch):
        monkeypatch.setattr(treeio, "MAX_RECURSION", 2000)
        tree = treeio.read_tree(ladder(400))
        with pytest.raises(TwigError, match="400 levels deep|399 levels deep"):
            treeio.copy_tree(tree)

    def test_shallow_trees_do_not_change_the_limit(self):
        import sys
        tree = treeio.read_tree(SPECIES)
        before = sys.getrecursionlimit()
        with treeio.deep_recursion(tree):
            assert sys.getrecursionlimit() == before
