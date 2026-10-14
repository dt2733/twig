"""Tests for the output writers (twig.labeling.core.writers)."""

import pytest

pytest.importorskip("toytree")

from twig.utils import TwigError
from twig.labeling.core import treeio
from twig.labeling.core.constants import FORMAT_NAMES
from twig.labeling.core.sanitize import NameTable
from twig.labeling.core.validate import WARNING, errors, warnings
from twig.labeling.core.writers import (
    WRITERS,
    WriteContext,
    get_writer,
    internal_texts,
    row_ids,
    written_branches,
)
from conftest import LABELED, SPECIES

# the concordant gene of the design, written as unrooted by HyPhy and without lengths
UNROOTED = (
    "(O{Reference},(M1{Reference},(A{Test},(B{Test},C{Rev}){Test})){Reference},"
    "(M2{Reference},(D{Test},E{Test})){Reference});"
)
DISTINCT = "(O{X}:1,((M1{R}:1,(A{T}:1,(B{T}:1,C{V}:1){T}:1){S}:1){R}:1,(M2{R}:1,(D{T}:1,E{T}:1){W}:1){R}:1){Z}:1);"


def render(name, newick, **kwargs):
    return get_writer(name).render(treeio.read_tree(newick), WriteContext(**kwargs))


def check(name, newick, **kwargs):
    return get_writer(name).check(treeio.read_tree(newick), WriteContext(**kwargs))


class TestRegistry:
    def test_the_writers_are_the_formats_named_in_the_constants(self):
        assert tuple(WRITERS) == FORMAT_NAMES

    def test_every_writer_has_a_suffix_and_help(self):
        for cls in WRITERS.values():
            assert cls.suffix.startswith(".") and cls.help

    def test_unknown_names_are_an_error(self):
        with pytest.raises(TwigError, match="unknown output format 'fasta'"):
            get_writer("fasta")


class TestBraces:
    def test_the_labeled_tree_is_written_as_it_was_read(self):
        assert render("braces", LABELED) == LABELED

    def test_lengths_are_written_when_the_input_had_them(self):
        newick = "((A{Test}:0.1,B:0.25){Rev}:0.3,C:0.4);"
        assert render("braces", newick) == newick
        assert render("braces", newick, dists=False) == "((A{Test},B){Rev},C);"

    def test_the_length_format_can_be_changed(self):
        assert render("braces", "(A:0.123456,B:2);", dist_fmt="%.2f") == "(A:0.12,B:2.00);"

    def test_the_root_is_written_as_it_is(self):
        assert render("braces", DISTINCT) == DISTINCT

    def test_a_label_with_a_character_that_cannot_be_written_is_an_error(self):
        tree = treeio.read_tree(SPECIES)
        tree[0].label = "a b"
        found = errors(get_writer("braces").check(tree, WriteContext()))
        assert len(found) == 1 and found[0].code == "label-name"
        with pytest.raises(TwigError, match="cannot write label"):
            get_writer("braces").render(tree, WriteContext())

    def test_a_clean_tree_has_no_issues(self):
        assert check("braces", LABELED) == []


class TestNewick:
    def test_labels_are_removed(self):
        assert render("newick", LABELED) == SPECIES

    def test_lengths_stay(self):
        assert render("newick", "((A{T}:0.1,B:0.2):0.3,C:0.4);") == "((A:0.1,B:0.2):0.3,C:0.4);"

    def test_labels_are_not_checked(self):
        tree = treeio.read_tree(SPECIES)
        tree[0].label = "a b"
        assert get_writer("newick").check(tree, WriteContext()) == []


class TestHyphy:
    def test_the_designed_tree_is_written_unrooted(self):
        assert render("hyphy", LABELED) == UNROOTED

    def test_an_unrooted_tree_is_written_as_it_is(self):
        newick = "(O{Reference},(M1{Reference},(A{Test},B{Test})),(M2{Reference},(D{Test},E{Test})));"
        assert render("hyphy", newick) == newick

    def test_the_summed_length_goes_to_the_kept_branch(self):
        out = render("hyphy", "(A{X}:1,(B{Y}:2,C{Y}:3){X}:4);")
        assert out == "(A{X}:5,B{Y}:2,C{Y}:3);"

    def test_keep_root_writes_the_rooted_tree(self):
        assert render("hyphy", LABELED, keep_root=True) == LABELED

    def test_lengths_that_the_input_did_not_have_are_not_invented(self):
        # the tips and the merged root branch were written with toytree's 1.0, and O with 2
        out = render("hyphy", "(O,((M1,(A,(B,C)1:0.5)0.9:0.3)1:1.2,(M2,(D,E)1:0.8)1:0.4));")
        assert out == "(O,(M1,(A,(B,C):0.5):0.3):1.2,(M2,(D,E):0.8):0.4);"
        assert render("hyphy", "((A,B):1,(C:1,D:1):2);") == "(A,B,(C:1,D:1):3);"
        assert render("hyphy", "((A:1,B:1),(C:1,D:1):2);") == "(A:1,B:1,(C:1,D:1));"

    def test_short_branches_without_a_length_are_not_reported(self):
        assert check("hyphy", "((A:1,B:1){T},(C:1,D:1):0.5,E:1);") == []

    @pytest.mark.parametrize("policy, expected", [
        ("unlabel", "(O:2,"),
        ("labeled", "(O:2,"),
        ("X", "(O{X}:2,"),
        ("Z", "(O{Z}:2,"),
        ("outgroup", "(O{X}:2,"),
        ("ingroup", "(O{Z}:2,"),
    ])
    def test_root_edge_policies(self, policy, expected):
        assert render("hyphy", DISTINCT, root_edges=policy).startswith(expected)

    def test_the_error_policy_stops_the_write(self):
        with pytest.raises(TwigError, match="different labels"):
            render("hyphy", DISTINCT, root_edges="error")
        found = errors(check("hyphy", DISTINCT, root_edges="error"))
        assert [i.code for i in found] == ["write"]

    def test_equal_sizes_need_an_outgroup_side_for_the_outgroup_policy(self):
        newick = "((A{X},B{X}){X},(C{Y},D{Y}){Y});"
        with pytest.raises(TwigError, match="equally large"):
            render("hyphy", newick, root_edges="outgroup")
        tree = treeio.read_tree(newick)
        first, second = (c.idx for c in tree.treenode.children)
        write = lambda side, policy: get_writer("hyphy").render(tree, WriteContext(root_edges=policy, outgroup_side=side))
        assert write(second, "outgroup") == "(A{X},B{X},(C{Y},D{Y}){Y});"
        assert write(first, "outgroup") == "(A{X},B{X},(C{Y},D{Y}){X});"
        assert write(second, "ingroup") == "(A{X},B{X},(C{Y},D{Y}){X});"

    def test_a_label_that_differs_at_the_root_gives_a_warning(self):
        found = warnings(check("hyphy", DISTINCT))
        assert [i.code for i in found] == ["root-merge"]
        assert "X and Z" in str(found[0]) and "no label" in str(found[0])

    def test_no_warning_when_the_labels_agree(self):
        assert check("hyphy", LABELED) == []

    def test_writing_rooted_warns_when_the_root_labels_differ(self):
        found = warnings(check("hyphy", DISTINCT, keep_root=True))
        assert len(found) == 1 and "merges the two branches at the root itself" in found[0].message
        assert check("hyphy", LABELED, keep_root=True) == []

    @pytest.mark.parametrize("label", ["Unlabeled", "All", "1abc", "a-b"])
    def test_label_names_are_checked_for_hyphy(self, label):
        found = errors(check("hyphy", f"(A{{{label}}},B,C);"))
        assert [i.code for i in found] == ["label-name"]

    def test_names_are_made_safe(self):
        out = render("hyphy", "(AmaranthusPumilus-Sample1{T},Beta{T},(C.1,D:1.5));".replace("D:1.5", "D"))
        assert out == "(AmaranthusPumilus_Sample1{T},Beta_{T},(C_1,D));"

    def test_names_that_collide_after_sanitizing_are_an_error(self):
        tree = treeio.read_tree("(A-1,A.1,(B,C));")
        found = errors(get_writer("hyphy").check(tree, WriteContext()))
        assert len(found) == 1 and "collide" in found[0].message and "A-1 and A.1" in found[0].message
        with pytest.raises(TwigError, match="collide"):
            get_writer("hyphy").render(tree, WriteContext())

    def test_names_that_differ_only_in_case_collide(self):
        with pytest.raises(TwigError, match="collide"):
            render("hyphy", "(A_b,a_B,(C,D));")

    def test_a_name_table_gives_the_same_written_name_in_every_tree(self):
        table = NameTable("hyphy", suffix=True)
        table.update(["A-1", "A.1", "B", "C"])
        table.build()
        out = render("hyphy", "(A-1,A.1,(B,C));", names=table)
        assert out == "(A_1,A_1_2,(B,C));"
        assert render("hyphy", "(A.1,(B,C),A-1);", names=table) == "(A_1_2,(B,C),A_1);"

    def test_labels_are_not_changed_by_sanitizing(self):
        assert render("hyphy", "(A-1{Test},(B,C){Rev},D);") == "(A_1{Test},(B,C){Rev},D);"

    def test_short_labeled_internal_branches_are_reported(self):
        found = warnings(check("hyphy", "((A:1,B:1){T}:0.0000001,(C:1,D:1){T}:0.5,E:1);"))
        assert [i.code for i in found] == ["short-labeled"]
        assert "1 labeled internal branch has a length of at most" in found[0].message and "Constrain" in found[0].message

    def test_short_tips_and_unlabeled_branches_are_not_reported(self):
        assert check("hyphy", "((A:0.0000001{T},B:1):0.0000001,(C:1,D:1){T}:0.5,E:1);") == []

    def test_short_lengths_are_not_checked_without_lengths(self):
        assert check("hyphy", "((A,B){T},(C,D),E);") == []
        assert check("hyphy", "((A:1,B:1){T}:0.0000001,C:1,D:1);", dists=False) == []

    def test_branch_lengths_that_are_not_numbers_are_an_error(self):
        tree = treeio.read_tree("((A:1,B:1):1,C:1,D:1);")
        tree[0]._dist = float("nan")
        found = errors(get_writer("hyphy").check(tree, WriteContext()))
        assert [i.code for i in found] == ["branch-length"]
        assert found[0].message == "1 branch length is not a number"

    def test_issues_name_the_locus(self):
        found = get_writer("hyphy").check(treeio.read_tree(DISTINCT), WriteContext(locus="OG0001"))
        assert str(found[0]).startswith("OG0001: ") and found[0].level == WARNING


class TestInternalNames:
    NAMED = "((A:1,B:1)n1:1,(C:1,D:1)n2:1)root;"
    SUPPORT = "((A:1,B:1)95:1,(C:1,D:1)0.8:1,E:1);"

    def test_none_writes_no_names(self):
        assert render("braces", self.NAMED) == "((A:1,B:1):1,(C:1,D:1):1);"

    def test_keep_writes_existing_names(self):
        assert render("braces", self.NAMED, internal_names="keep") == "((A:1,B:1)n1:1,(C:1,D:1)n2:1);"

    def test_the_root_name_is_never_written(self):
        assert "root" not in render("braces", self.NAMED, internal_names="keep")

    def test_auto_names_nodes_in_postorder(self):
        assert render("braces", "(((A,B),C),(D,E));", internal_names="auto") == "(((A,B)twn1,C)twn2,(D,E)twn3);"

    def test_auto_leaves_out_names_that_equal_a_tip(self):
        assert render("braces", "((twn1,B),(C,D));", internal_names="auto") == "((twn1,B)twn2,(C,D)twn3);"

    def test_auto_skips_the_dissolved_root_child(self):
        out = render("hyphy", "(O,((A,B),(C,D)));", internal_names="auto")
        assert out == "(O,(A,B)twn1,(C,D)twn2);"

    def test_support_writes_the_support_values(self):
        assert render("braces", self.SUPPORT, internal_names="support") == "((A:1,B:1)95:1,(C:1,D:1)0.8:1,E:1);"

    def test_support_skips_nodes_without_a_value(self):
        assert render("braces", "((A,B)95,(C,D),E);", internal_names="support") == "((A,B)95,(C,D),E);"

    def test_keep_sanitizes_names_for_hyphy(self):
        out = render("hyphy", "((A,B)n-1,(C,D)n.2,E);", internal_names="keep")
        assert out == "((A,B)n_1,(C,D)n_2,E);"

    def test_keep_rejects_a_name_used_twice(self):
        with pytest.raises(TwigError, match="used twice or equals a tip name"):
            render("braces", "((A,B)n,(C,D)n,E);", internal_names="keep")

    def test_keep_rejects_a_name_equal_to_a_tip_even_in_another_case(self):
        with pytest.raises(TwigError, match="used twice or equals a tip name"):
            render("braces", "((A,B)c,(C,D)x,E);", internal_names="keep")

    def test_unknown_choices_are_an_error(self):
        with pytest.raises(TwigError, match="unknown choice 'all'"):
            render("braces", self.NAMED, internal_names="all")

    def test_the_function_returns_texts_by_idx(self):
        tree = treeio.read_tree("((A,B),(C,D));")
        texts = internal_texts(tree, "auto", ["A", "B", "C", "D"])
        assert sorted(texts.values()) == ["twn1", "twn2"]
        assert all(not tree[i].is_leaf() and tree[i].up is not None for i in texts)


class TestTable:
    def rows(self, newick, writer="braces", **kwargs):
        tree = treeio.read_tree(newick)
        ctx = WriteContext(locus="sp", **kwargs)
        if writer != "table":
            ctx.plan = get_writer(writer).plan(tree, ctx)
        lines = get_writer("table").render(tree, ctx).split("\n")
        header = lines[0].split("\t")
        return header, [dict(zip(header, line.split("\t"))) for line in lines[1:]]

    def test_one_row_per_branch_and_no_row_for_the_root(self):
        header, rows = self.rows(LABELED)
        assert header[:4] == ["locus", "idx", "node", "is_tip"]
        assert len(rows) == 14
        assert {r["locus"] for r in rows} == {"sp"}

    def test_rows_describe_the_branches(self):
        _, rows = self.rows(LABELED)
        by_node = {r["node"]: r for r in rows}
        assert by_node["C"]["label"] == "Rev" and by_node["C"]["is_tip"] == "yes" and by_node["C"]["n_tips"] == "1"
        bc = next(r for r in rows if r["tips"] == "B;C")
        assert bc["label"] == "Test" and bc["is_tip"] == "no" and bc["n_tips"] == "2" and bc["node"].startswith("S")
        abc = next(r for r in rows if r["tips"] == "A;B;C")
        assert abc["label"] == ""

    def test_kinds_and_sources_come_from_the_context(self):
        tree = treeio.read_tree(LABELED)
        ctx = WriteContext(locus="sp", kinds={2: "tip"}, sources={2: "Test:fitch"})
        lines = get_writer("table").render(tree, ctx).split("\n")
        row = dict(zip(lines[0].split("\t"), lines[3].split("\t")))
        assert row["node"] == "A" and row["kind"] == "tip" and row["source"] == "Test:fitch"

    def test_lengths_and_support(self):
        _, rows = self.rows("((A:1,B:2)95:3,C:4);")
        by_node = {r["node"]: r for r in rows}
        assert by_node["B"]["dist"] == "2"
        abc = next(r for r in rows if r["tips"] == "A;B")
        assert abc["dist"] == "3" and abc["support"] == "95"

    def test_the_species_column_names_every_species_of_a_large_clade(self):
        # the column was built from the first 200 tip names only, so B and C were missing without a mark
        names = ",".join(f"A_{i}" for i in range(1, 211))
        _, rows = self.rows(f"(O_1,((({names}),B_1),C_1));", writer="table", species_of=lambda n: n.split("_")[0])
        by_size = {int(r["n_tips"]): r["species"] for r in rows}
        assert by_size[210] == "A" and by_size[211] == "A;B" and by_size[212] == "A;B;C"
        assert by_size[1] in {"A", "B", "C", "O"}

    def test_the_species_column_is_cut_after_200_species(self):
        names = ",".join(f"S{i}_1" for i in range(1, 206))
        _, rows = self.rows(f"(O_1,({names}));", writer="table", species_of=lambda n: n.split("_")[0])
        big = next(r for r in rows if r["n_tips"] == "205")
        assert big["species"].endswith(";...") and len(big["species"].split(";")) == 201
        assert big["tips"].endswith(";...") and len(big["tips"].split(";")) == 201

    def test_no_lengths_for_a_tree_read_without_them(self):
        _, rows = self.rows(LABELED)
        assert {r["dist"] for r in rows} == {""}

    def test_no_length_for_a_branch_that_had_none(self):
        _, rows = self.rows("((A,B:2):1,C);", writer="hyphy")
        assert {r["node"]: r["dist"] for r in rows if r["is_tip"] == "yes"} == {"A": "", "B": "2", "C": ""}
        assert [r["dist"] for r in rows if r["is_tip"] == "no"] == ["1"]

    def test_the_root_merge_is_marked(self):
        _, rows = self.rows(DISTINCT, writer="hyphy")
        by_node = {r["node"]: r for r in rows}
        dissolved = next(r for r in rows if r["written"] == "no")
        assert dissolved["merged_into"] == "O" and dissolved["note"] == "root-merge"
        assert dissolved["label"] == "Z"
        kept = by_node["O"]
        assert kept["written"] == "yes" and kept["label"] == "" and kept["dist"] == "2"
        assert kept["note"] == "root-merge (labels differed)"

    def test_written_names_are_listed(self):
        _, rows = self.rows("(A-1{T},(B,C)),", writer="hyphy") if False else self.rows("(A-1{T},(B,C));", writer="hyphy")
        by_node = {r["node"]: r for r in rows}
        assert by_node["A-1"]["written_name"] == "A_1"

    def test_auto_names_join_to_the_table(self):
        _, rows = self.rows("(((A,B),C),(D,E));", writer="braces", internal_names="auto")
        assert sorted(r["written_name"] for r in rows if r["is_tip"] == "no") == ["twn1", "twn2", "twn3"]

    def test_species_column_is_added_when_a_mapping_is_given(self):
        header, rows = self.rows("((A_1,A_2),(B_1,C_1));", species_of=lambda n: n.split("_")[0])
        assert header[header.index("tips") + 1] == "species"
        top = next(r for r in rows if r["n_tips"] == "2" and r["tips"] == "A_1;A_2")
        assert top["species"] == "A"

    def test_header_can_be_left_out(self):
        tree = treeio.read_tree(SPECIES)
        text = get_writer("table").render(tree, WriteContext(header=False))
        assert not text.startswith("locus")

    def test_long_tip_lists_are_cut(self):
        names = [f"t{i}" for i in range(300)]
        newick = "(" + ",".join(names) + ");"
        tree = treeio.read_tree(newick)
        # the root is not a branch, so use a tree with one clade
        tree = treeio.read_tree("((" + ",".join(names) + "),x);")
        lines = get_writer("table").render(tree, WriteContext()).split("\n")
        big = next(dict(zip(lines[0].split("\t"), line.split("\t"))) for line in lines[1:] if line.split("\t")[4] == "300")
        assert big["tips"].endswith(";...") and big["tips"].count(";") == 200

    def test_row_ids_use_unique_internal_names(self):
        tree = treeio.read_tree("((A,B)n1,(C,D)n1,E);")
        ids = row_ids(tree)
        assert sorted(v for v in ids.values() if v not in "ABCDE") == [f"S{i}" for i in sorted(
            n.idx for n in tree if not n.is_leaf() and n.up is not None)]
        tree = treeio.read_tree("((A,B)n1,(C,D)n2,E);")
        assert {"n1", "n2"} <= set(row_ids(tree).values())

    def test_written_branches_skip_the_dissolved_child(self):
        tree = treeio.read_tree(DISTINCT)
        plan = get_writer("hyphy").plan(tree, WriteContext())
        rows = list(written_branches(tree, plan.merge))
        assert len(rows) == 13 and tree[plan.merge.dissolved] not in [n for n, _, _ in rows]
        kept = next(r for r in rows if r[0].idx == plan.merge.kept)
        assert kept[1] is None and kept[2] == 2.0


class TestRepeatedTips:
    def test_a_tree_with_repeated_tip_names_can_be_written_in_braces_format(self):
        assert render("braces", "((A,A),B);") == "((A,A),B);"

    def test_hyphy_needs_a_different_name_for_every_tip(self):
        with pytest.raises(TwigError, match="tip names are not unique"):
            render("hyphy", "((A,A),B,C);")
        assert errors(check("hyphy", "((A,A),B,C);"))[0].code == "write"
