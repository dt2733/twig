"""Tests for the tree-label command (twig.labeling.tree_label)."""

import os

import pytest

pytest.importorskip("toytree")

from twig.labeling.core import treeio
from twig.labeling.tree_label import read_paml_codes
from twig.utils import TwigError
from conftest import DATA, LABELED, SPECIES

NAMED = DATA / "species.named.nwk"          # (O,((M1,(A,(B,C)BC)ABC)X,(M2,(D,E)DE)Y)Z);
DIOECIOUS = DATA / "dioecious.txt"          # A, B, D and E

# the labeled species tree written for HyPhy: the two branches at the root become one
UNROOTED = (
    "(O{Reference},(M1{Reference},(A{Test},(B{Test},C{Rev}){Test})){Reference},"
    "(M2{Reference},(D{Test},E{Test})){Reference});"
)
# step 1 of the analysis in the design
STEP_1 = f"-i {DATA / 'species.nwk'} -L Test {DIOECIOUS} --rule fitch --stems unlabel --losses Rev --background Reference"


def labeled_names(out):
    """Return {name: label} for the labeled branches of a tree, named by tip or internal name."""
    tree = treeio.read_tree(out.strip())
    return {node.name: node.label for node in tree if node.label}


def trees_of(out):
    return [newick for _, newick in treeio.iter_trees_text(out)]


def table_of(text):
    lines = text.rstrip("\n").split("\n")
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


@pytest.fixture
def run(run_twig):
    """Run tree-label with arguments and return (code, stdout, stderr)."""
    return lambda *args: run_twig("tree-label", *args)


class TestRulesOnTheDesignedSpeciesTree:
    @pytest.mark.parametrize("options, expected", [
        ("--rule tips", {"A", "B", "D", "E"}),
        ("--rule all", {"A", "B", "D", "E", "DE"}),
        ("--rule crown", {"A", "B", "D", "E"}),
        ("--rule some", {"A", "B", "D", "E", "BC", "ABC", "X", "DE", "Y", "Z"}),
        ("--rule fitch", {"A", "B", "BC", "ABC", "D", "E", "DE"}),
        ("--rule Parsimony", {"A", "B", "BC", "ABC", "D", "E", "DE"}),
        ("--rule fitch --stems unlabel", {"A", "B", "BC", "D", "E"}),
        ("--rule fitch --ambiguous unlabel", {"A", "B", "D", "E", "DE"}),
        ("--rule fitch --ambiguous deltran", {"A", "B", "D", "E", "DE"}),
        ("--internal-nodes all", {"A", "B", "D", "E", "DE"}),
    ])
    def test_labeled_branches(self, run, options, expected):
        code, out, err = run(f"-i {NAMED} -L Test {DIOECIOUS} {options} --internal-names keep")
        assert code == 0, err
        assert set(labeled_names(out)) == expected
        assert set(labeled_names(out).values()) == {"Test"}

    def test_fitch_reports_its_cost(self, run):
        _, _, err = run(f"-i {NAMED} -L Test {DIOECIOUS} --rule fitch")
        assert "label Test (fitch) selected 4 tips, covers 7 branches and keeps 7, parsimony cost 3" in err

    def test_the_log_counts_what_a_label_keeps_in_the_end(self, run):
        _, _, err = run(f"-i {NAMED} -L Test {DIOECIOUS} --rule fitch --stems unlabel --leaf-nodes skip")
        assert "label Test (fitch) selected 4 tips, covers 7 branches and keeps 1, parsimony cost 3" in err

    def test_the_rule_can_be_given_to_one_label(self, run):
        _, out, _ = run(f"-i {NAMED} -L Test:fitch {DIOECIOUS} --internal-names keep")
        assert set(labeled_names(out)) == {"A", "B", "BC", "ABC", "D", "E", "DE"}

    def test_regular_expressions_select_tips(self, run):
        _, out, _ = run(f"-i {NAMED} -L Test ^[AB]$ --rule tips")
        assert set(labeled_names(out)) == {"A", "B"}

    def test_exact_names_do_not_match_parts_of_other_names(self, run, tmp_path):
        tree = tmp_path / "t.nwk"
        tree.write_text("(MA,(A,B));")
        assert set(labeled_names(run(f"-i {tree} -L Test A --rule tips")[1])) == {"MA", "A"}
        assert set(labeled_names(run(f"-i {tree} -L Test A --rule tips --exact")[1])) == {"A"}

    def test_invert_labels_the_other_tips(self, run):
        _, out, _ = run(f"-i {NAMED} -L Test {DIOECIOUS} --rule tips --invert")
        assert set(labeled_names(out)) == {"O", "M1", "C", "M2"}

    def test_leaf_nodes_can_be_skipped(self, run):
        _, out, _ = run(f"-i {NAMED} -L Test {DIOECIOUS} --rule all --leaf-nodes skip --internal-names keep")
        assert set(labeled_names(out)) == {"DE"}

    def test_skipped_tips_do_not_get_the_background_label(self, run):
        _, out, _ = run(f"-i {NAMED} -L Test {DIOECIOUS} --rule fitch --stems unlabel --leaf-nodes skip --background Reference --internal-names keep")
        names = labeled_names(out)
        assert not {"A", "B", "D", "E"} & set(names) and names["BC"] == "Test" and names["C"] == "Reference"

    def test_the_tips_rule_cannot_skip_tips(self, run):
        code, _, err = run(f"-i {NAMED} -L Test {DIOECIOUS} --rule tips --leaf-nodes skip")
        assert code == 1 and "labels nothing" in err and "Traceback" not in err


class TestTheDesignedAnalysis:
    def test_step_one_writes_the_labeled_species_tree(self, run):
        code, out, err = run(STEP_1)
        assert code == 0, err
        assert out == LABELED + "\n"

    def test_a_nuisance_species_gives_the_same_tree(self, run):
        code, out, err = run(
            f"-i {DATA / 'species.nwk'} -L Test {DIOECIOUS} --rule fitch --nuisance C --nuisance-label Rev "
            "--stems unlabel --background Reference")
        assert code == 0, err
        assert out == LABELED + "\n"

    def test_the_nuisance_label_is_called_nuisance_by_default(self, run):
        _, out, _ = run(f"-i {DATA / 'species.nwk'} -L Test {DIOECIOUS} --rule fitch --nuisance C --stems unlabel --background Reference")
        assert "C{Nuisance}" in out

    def test_nuisance_species_do_not_count_as_evidence_against_the_trait(self, run):
        # without the states of M1 and M2 one gain above all of the ingroup and one loss on C is the best history
        _, out, _ = run(f"-i {NAMED} -L Test {DIOECIOUS} --rule fitch --nuisance ^M[12]$ --internal-names keep")
        names = labeled_names(out)
        assert names["M1"] == names["M2"] == "Nuisance"
        assert {n for n, label in names.items() if label == "Test"} == {"A", "B", "D", "E", "BC", "ABC", "DE", "X", "Y", "Z"}

    def test_the_listing_checks_the_result_by_eye(self, run):
        code, out, err = run(STEP_1 + " --show")
        assert code == 0 and out == LABELED + "\n"
        assert "labeled branches, clade by clade" in err
        assert "  Rev        loss      C" in err
        assert "  -  stem  A+B+C" in err and "  -  stem  D+E" in err
        assert "summary  6 Reference, 5 Test, 1 Rev, 2 unlabeled" in err

    def test_the_log_says_what_each_label_did(self, run):
        _, _, err = run(STEP_1)
        assert "read 1 tree from" in err
        assert "label Test (fitch) selected 4 tips" in err
        assert "14 branches, labeled Reference 6, Test 5, Rev 1, unlabeled 2, 2 stems" in err

    def test_the_labeled_tree_is_the_input_for_hyphy(self, run, tmp_path):
        labeled = tmp_path / "sp.labeled.nwk"
        labeled.write_text(LABELED + "\n")
        code, out, err = run(f"-i {labeled} -F hyphy")
        assert code == 0, err
        assert out == UNROOTED + "\n"
        assert "no label sets were given, so the labels already in the tree are kept" in err

    def test_labeling_and_writing_for_hyphy_in_one_step(self, run):
        _, out, _ = run(STEP_1 + " -F hyphy")
        assert out == UNROOTED + "\n"

    def test_stems_can_get_their_own_label(self, run):
        _, out, _ = run(f"-i {NAMED} -L Test {DIOECIOUS} --rule fitch --stems Stem --internal-names keep")
        names = labeled_names(out)
        assert names["ABC"] == names["DE"] == "Stem" and names["BC"] == "Test"

    def test_the_stem_of_crown_is_the_branch_it_leaves_out(self, run):
        _, out, err = run(f"-i {NAMED} -L Test {DIOECIOUS} --rule crown --stems Stem --background Reference --internal-names keep --show")
        names = labeled_names(out)
        assert names["DE"] == "Stem" and names["Y"] == "Reference" and names["D"] == "Test"
        assert "  Stem       stem      D+E" in err

    def test_an_ambiguous_branch_is_not_a_stem(self, run):
        _, out, err = run(f"-i {NAMED} -L Test {DIOECIOUS} --rule fitch --ambiguous Unsure --stems Stem --internal-names keep --show")
        names = labeled_names(out)
        assert names["ABC"] == names["BC"] == "Unsure" and names["DE"] == "Stem"
        assert "  Unsure  ambiguous  A+B+C" in err

    def test_no_stem_below_a_root_that_has_the_trait(self, run, tmp_path):
        tree = tmp_path / "t.nwk"
        tree.write_text("((A,B)AB,(C,(D,E)DE)CDE);")
        _, out, _ = run(f"-i {tree} -L Test ^[ABD]$ --rule fitch --gain-cost 2 --stems unlabel --internal-names keep")
        assert set(labeled_names(out)) == {"A", "B", "AB", "CDE", "D", "DE"}

    def test_ambiguous_branches_can_get_their_own_label(self, run):
        _, out, _ = run(f"-i {NAMED} -L Test {DIOECIOUS} --rule fitch --ambiguous Amb --internal-names keep")
        assert {n for n, label in labeled_names(out).items() if label == "Amb"} == {"ABC", "BC"}


class TestOutputs:
    def test_output_file(self, run, tmp_path):
        out = tmp_path / "out.nwk"
        code, stdout, err = run(f"{STEP_1} -o {out}")
        assert code == 0 and stdout == ""
        assert out.read_text() == LABELED + "\n"
        assert f"wrote 1 tree to {out}" in err

    def test_an_existing_output_needs_force(self, run, tmp_path):
        out = tmp_path / "out.nwk"
        out.write_text("keep me")
        code, _, err = run(f"{STEP_1} -o {out}")
        assert code == 1 and "use --force" in err and out.read_text() == "keep me"
        assert run(f"{STEP_1} -o {out} -f")[0] == 0
        assert out.read_text() == LABELED + "\n"

    def test_outputs_must_differ_from_each_other_and_from_the_input(self, run, tmp_path):
        both = tmp_path / "both.txt"
        code, _, err = run(f"{STEP_1} -o {both} -B {both}")
        assert code == 1 and "must be different files" in err
        tree = tmp_path / "t.nwk"
        tree.write_text(SPECIES)
        code, _, err = run(f"-i {tree} -L Test A -o {tree} -f")
        assert code == 1 and "same as the input file" in err and tree.read_text() == SPECIES

    def test_an_output_that_is_not_a_regular_file_is_refused_even_with_force(self, run, tmp_path):
        folder = tmp_path / "folder"
        folder.mkdir()
        link = tmp_path / "null"
        link.symlink_to(os.devnull)
        for option in ("-o", "-B"):
            for path in (folder, link):
                code, out, err = run(f"{STEP_1} {option} {path} -f")
                assert code == 1 and "not a regular file, so it is never replaced" in err, (option, path)
                assert out == "" and "Traceback" not in err
        assert folder.is_dir() and link.is_symlink() and sorted(p.name for p in tmp_path.iterdir()) == ["folder", "null"]
        # without -o the tree still goes to stdout
        code, out, _ = run(STEP_1)
        assert code == 0 and out == LABELED + "\n"

    def test_plain_newick_has_no_labels(self, run):
        _, out, _ = run(STEP_1 + " -F newick")
        assert out == SPECIES + "\n"

    def test_branch_table_next_to_the_tree(self, run, tmp_path):
        table = tmp_path / "branches.tsv"
        code, out, err = run(f"{STEP_1} -B {table}")
        assert code == 0 and out == LABELED + "\n"
        rows = table_of(table.read_text())
        assert len(rows) == 14 and {r["locus"] for r in rows} == {"species"}
        by_node = {r["node"]: r for r in rows}
        assert by_node["C"]["label"] == "Rev" and by_node["C"]["kind"] == "tip" and by_node["C"]["source"] == "loss"
        assert by_node["A"]["source"] == "Test:fitch"
        assert by_node["O"]["source"] == "background"
        stems = [r for r in rows if r["kind"] == "stem"]
        assert sorted(r["tips"] for r in stems) == ["A;B;C", "D;E"] and {r["label"] for r in stems} == {""}

    def test_table_as_the_output(self, run):
        code, out, _ = run(STEP_1 + " -F table")
        assert code == 0
        assert out.startswith("locus\tidx\tnode\tis_tip")
        assert len(table_of(out)) == 14

    def test_the_table_shows_what_hyphy_will_see(self, run, tmp_path):
        table = tmp_path / "branches.tsv"
        run(f"{STEP_1} -F hyphy -B {table}")
        rows = table_of(table.read_text())
        dissolved = [r for r in rows if r["written"] == "no"]
        assert len(dissolved) == 1 and dissolved[0]["merged_into"] == "O"

    def test_internal_names_for_hyphy(self, run):
        _, out, _ = run(STEP_1 + " -F hyphy --internal-names auto")
        assert out == (
            "(O{Reference},(M1{Reference},(A{Test},(B{Test},C{Rev})twn1{Test})twn2)twn3{Reference},"
            "(M2{Reference},(D{Test},E{Test})twn4)twn5{Reference});\n"
        )

    def test_branch_lengths(self, run, tmp_path):
        tree = tmp_path / "t.nwk"
        tree.write_text("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6,E:0.7);")
        _, out, _ = run(f"-i {tree} -L Test ^[AB]$")
        assert out == "((A{Test}:0.1,B{Test}:0.2){Test}:0.3,(C:0.4,D:0.5):0.6,E:0.7);\n"
        _, out, _ = run(f"-i {tree} -L Test ^[AB]$ --no-dists")
        assert out == "((A{Test},B{Test}){Test},(C,D),E);\n"

    def test_lengths_that_the_input_did_not_have_are_not_invented(self, run, tmp_path):
        # an ASTRAL species tree has internal lengths only, and every tip was written with toytree's 1.0
        tree = tmp_path / "astral.nwk"
        tree.write_text("(O,((M1,(A,(B,C)1:0.5)0.9:0.3)1:1.2,(M2,(D,E)1:0.8)1:0.4));")
        code, out, err = run(f"-i {tree} -L Test ^[ABDE]$ --rule fitch -F hyphy")
        assert code == 0, err
        assert out == "(O,(M1,(A{Test},(B{Test},C){Test}:0.5){Test}:0.3):1.2,(M2,(D{Test},E{Test}){Test}:0.8):0.4);\n"
        table = tmp_path / "t.tsv"
        run(f"-i {tree} -B {table}")
        assert {row["node"]: row["dist"] for row in table_of(table.read_text())}["A"] == ""
        tree.write_text("((A,B),C):0;")
        _, out, _ = run(f"-i {tree}")
        assert out == "((A,B),C);\n"

    def test_keep_root_writes_the_tree_rooted(self, run):
        _, out, _ = run(STEP_1 + " -F hyphy --keep-root")
        assert out == LABELED + "\n"

    def test_several_trees_in_one_file(self, run, tmp_path):
        trees = tmp_path / "trees.nwk"
        trees.write_text("((A,B),(C,D),E);\n((A,C),(B,D),E);\n")
        table = tmp_path / "b.tsv"
        code, out, err = run(f"-i {trees} -L Test ^[AB]$ -B {table}")
        assert code == 0, err
        assert trees_of(out) == ["((A{Test},B{Test}){Test},(C,D),E);", "((A{Test},C),(B{Test},D),E);"]
        assert {r["locus"] for r in table_of(table.read_text())} == {"trees_1", "trees_2"}
        assert "read 2 trees from" in err

    def test_the_listing_names_each_tree_when_there_are_several(self, run, tmp_path):
        trees = tmp_path / "trees.nwk"
        trees.write_text("((A,B),(C,D),E);\n((A,C),(B,D),E);\n")
        _, _, err = run(f"-i {trees} -L Test ^[AB]$ --show")
        assert "tree trees_1" in err and "tree trees_2" in err

    def test_the_listing_describes_the_tree_written_for_hyphy(self, run):
        # the listing counted 6 Reference branches while the HyPhy tree has 5
        _, out, err = run(STEP_1 + " -F hyphy --show")
        assert out == UNROOTED + "\n"
        assert "root     the branches O and M1+A+B+C+M2+D+E at the root are written as one branch, listed as O" in err
        assert "summary  5 Reference, 5 Test, 1 Rev, 2 unlabeled" in err
        _, _, err = run(STEP_1 + " --show")
        assert "summary  6 Reference, 5 Test, 1 Rev, 2 unlabeled" in err and "root     " not in err

    def test_the_listing_does_not_name_a_single_tree(self, run):
        _, _, err = run(STEP_1 + " --show")
        assert "tree species" not in err


class TestRootEdgesAndRooting:
    DISTINCT = "(O{X}:1,((M1{R}:1,(A{T}:1,(B{T}:1,C{V}:1){T}:1){S}:1){R}:1,(M2{R}:1,(D{T}:1,E{T}:1){W}:1){R}:1){Z}:1);"

    @pytest.fixture
    def distinct(self, tmp_path):
        path = tmp_path / "distinct.nwk"
        path.write_text(self.DISTINCT)
        return path

    @pytest.mark.parametrize("policy, label", [("unlabel", ""), ("X", "{X}"), ("Z", "{Z}"), ("outgroup", "{X}"), ("ingroup", "{Z}")])
    def test_the_label_of_the_joined_root_branch(self, run, distinct, policy, label):
        code, out, err = run(f"-i {distinct} -F hyphy --root-edges {policy}")
        assert code == 0, err
        assert out.startswith(f"(O{label}:2,")

    def test_a_warning_says_what_happened_at_the_root(self, run, distinct):
        _, _, err = run(f"-i {distinct} -F hyphy")
        assert "the two branches at the root have the labels X and Z, so the merged branch gets no label" in err

    def test_the_error_policy_stops(self, run, distinct):
        code, out, err = run(f"-i {distinct} -F hyphy --root-edges error")
        assert code == 1 and out == "" and "different labels" in err and "Traceback" not in err

    def test_an_unknown_policy_is_an_error(self, run, distinct):
        code, _, err = run(f"-i {distinct} --root-edges bad-name")
        assert code == 1 and "unknown root edge policy" in err

    def test_a_misspelled_policy_is_refused_before_any_tree_is_read(self, run, tmp_path):
        # outgrup was taken as a label name that no branch has, and the merged branch silently got no label
        for policy in ("outgrup", "Outgroup"):
            code, out, err = run(f"-i {tmp_path / 'nothing.nwk'} -F hyphy --root-edges {policy}")
            assert code == 1 and out == "" and "looks like" in err and "'outgroup'" in err
            assert "nothing.nwk not found" not in err and "Traceback" not in err

    def test_a_policy_label_that_no_branch_has_gives_a_warning(self, run, distinct):
        code, out, err = run(f"-i {distinct} -F hyphy --root-edges Q")
        assert code == 0 and out.startswith("(O:2,")
        assert "--root-edges names the label Q, which no branch has" in err
        _, _, err = run(f"-i {distinct} -F hyphy --root-edges Z")
        assert "which no branch has" not in err

    def test_reroot_on_an_outgroup_before_labeling(self, run, tmp_path):
        tree = tmp_path / "unrooted.nwk"
        tree.write_text("((M1,(A,(B,C))),(M2,(D,E)),O);")
        code, out, err = run(f"-i {tree} -r O -L Test {DIOECIOUS} --rule fitch --stems unlabel --losses Rev --background Reference -F hyphy")
        assert code == 0, err
        assert out == UNROOTED + "\n"
        assert "rooted on O" in err

    def test_an_unrooted_tree_without_r_gives_a_warning_whatever_the_costs(self, run, tmp_path):
        # with the basal node next to A the stem and the loss came out wrong, and nothing was said
        tree = tmp_path / "unrooted.nwk"
        tree.write_text("(A,(B,C),(M1,(O,(M2,(D,E)))));")
        options = f"-L Test {DIOECIOUS} --rule fitch --stems unlabel --losses Rev --background Reference -F hyphy"
        code, _, err = run(f"-i {tree} {options}")
        assert code == 0 and "the root has 3 children, so the tree looks unrooted" in err and "Give -r with the outgroup" in err
        code, out, err = run(f"-i {tree} -r O {options}")
        assert code == 0 and "looks unrooted" not in err
        assert labeled_names(out)["C"] == "Rev"

    def test_the_outgroup_policy_uses_the_side_of_the_outgroup_given_with_r(self, run, distinct):
        # here the outgroup is the whole ingroup, so the outgroup side is the branch above it
        _, out, _ = run(f"-i {distinct} -r M1,A,B,C,M2,D,E -F hyphy --root-edges outgroup")
        assert out.startswith("(O{Z}:2,")
        _, out, _ = run(f"-i {distinct} -r M1,A,B,C,M2,D,E -F hyphy --root-edges ingroup")
        assert out.startswith("(O{X}:2,")

    def test_a_rerooting_that_is_impossible_is_an_error(self, run, distinct):
        code, _, err = run(f"-i {distinct} -r Q")
        assert code == 1 and "none of the outgroups" in err and "Traceback" not in err

    def test_the_labels_of_the_old_root_are_merged_by_the_policy_when_rerooting(self, run, distinct):
        _, _, err = run(f"-i {distinct} -r B,C")
        assert "the two branches at the old root had the labels" in err

    def test_the_merged_old_root_branch_does_not_depend_on_the_order_of_the_children(self, run, tmp_path):
        # with the side of C written second, Ref stayed on the merged branch although the warning said no label
        for newick in ("((A{Test}:1,B:1):2,(C{Ref}:1,D:1){Ref}:3);", "((C{Ref}:1,D:1){Ref}:3,(A{Test}:1,B:1):2);"):
            tree = tmp_path / "t.nwk"
            tree.write_text(newick)
            code, out, err = run(f"-i {tree} -r C")
            assert code == 0 and out == "(C{Ref}:0.5,(D:1,(A{Test}:1,B:1):5){Ref}:0.5);\n"
            assert "the merged branch has no label" in err
            _, out, _ = run(f"-i {tree} -r C --root-edges labeled")
            assert out == "(C{Ref}:0.5,(D:1,(A{Test}:1,B:1){Ref}:5){Ref}:0.5);\n"
            _, out, _ = run(f"-i {tree} -r C -L Test ^A$ -F hyphy")
            assert out == "(C{Ref}:1,D:1,(A{Test}:1,B:1):5);\n"

    def test_internal_names_stay_on_their_branches_when_rerooting(self, run, tmp_path):
        # Z named the branch of Y after rerooting, and Y the tip branch of M2
        code, out, err = run(f"-i {NAMED} -r M2 --internal-names keep -L Test:stem ^[DE]$")
        assert code == 0 and out == "(M2,((D,E)DE{Test},((M1,(A,(B,C)BC)ABC)X,O)Y));\n"
        assert "the internal name Z is dropped, because its branch at the old root was joined with the branch of O" in err
        table = tmp_path / "t.tsv"
        run(f"-i {NAMED} -r M2 -B {table}")
        rows = {row["node"]: row["tips"] for row in table_of(table.read_text())}
        assert rows["Y"] == "M1;A;B;C;O" and rows["X"] == "M1;A;B;C" and "Z" not in rows


class TestExistingLabels:
    def test_existing_labels_are_kept_by_default(self, run, tmp_path):
        tree = tmp_path / "t.nwk"
        tree.write_text("((A{Old},B),(C,D),E);")
        _, out, _ = run(f"-i {tree} -L Test ^[AB]$ --rule tips")
        assert out == "((A{Old},B{Test}),(C,D),E);\n"

    def test_existing_labels_can_be_cleared(self, run, tmp_path):
        tree = tmp_path / "t.nwk"
        tree.write_text("((A{Old},B),(C,D),E);")
        _, out, _ = run(f"-i {tree} -L Test ^[AB]$ --rule tips --existing clear")
        assert out == "((A{Test},B{Test}),(C,D),E);\n"

    def test_clearing_without_new_labels_removes_all_labels(self, run):
        _, out, err = run(f"-i {DATA / 'species.labeled.nwk'} --existing clear")
        assert out == SPECIES + "\n"
        assert "the labels already in the tree are cleared" in err and "no branch is labeled" in err

    def test_a_label_on_the_root_is_ignored(self, run, tmp_path):
        tree = tmp_path / "t.nwk"
        tree.write_text("((A{T},B{T}),C,D){X};")
        code, out, err = run(f"-i {tree}")
        assert code == 0 and out == "((A{T},B{T}),C,D);\n"
        assert "ignoring the label 'X' on the root" in err

    def test_stems_act_on_labels_already_in_the_tree(self, run, tmp_path):
        tree = tmp_path / "t.nwk"
        tree.write_text("(O,((M1,(A{Test},(B{Test},C){Test}){Test}),(M2,(D{Test},E{Test}){Test})));")
        code, out, err = run(f"-i {tree} --stems Stem")
        assert code == 0, err
        assert out == "(O,((M1,(A{Test},(B{Test},C){Test}){Stem}),(M2,(D{Test},E{Test}){Stem})));\n"
        assert "--stems treats 2 branches labeled in the input tree as stems (A+B+C, D+E)" in err
        _, out, _ = run(f"-i {tree} --stems unlabel --background Reference")
        assert out == (
            "(O{Reference},((M1{Reference},(A{Test},(B{Test},C{Reference}){Test})){Reference},"
            "(M2{Reference},(D{Test},E{Test})){Reference}){Reference});\n"
        )

    def test_paml_labels_are_read_and_named(self, run, tmp_path):
        tree = tmp_path / "t.nwk"
        tree.write_text("(A #1,B,(C,D) #2);")
        code, out, err = run(f"-i {tree} --paml-codes 1=Test 2=Reference")
        assert code == 0, err
        assert out == "(A{Test},B,(C,D){Reference});\n"

    def test_paml_codes_must_look_like_number_equals_name(self):
        assert read_paml_codes(["1=Test", "2=Ref"]) == {"1": "Test", "2": "Ref"}
        assert read_paml_codes(["01=Test"]) == {"1": "Test"}
        assert read_paml_codes(None) is None
        for bad in ("1", "=Test", "1=", "x=Test"):
            with pytest.raises(TwigError, match="--paml-codes expects"):
                read_paml_codes([bad])

    @pytest.mark.parametrize("code, message", [
        ("0=unlabel", "--paml-codes gives 0 the name 'unlabel', which is a word that twig uses"),
        ("1=9x", "--paml-codes gives 1 the name '9x', which must start with a letter"),
        ("1=All", "--paml-codes gives 1 the name 'All', which collides with a HyPhy branch-set option"),
    ])
    def test_paml_code_names_are_checked_like_other_label_names(self, run, tmp_path, code, message):
        # 0=unlabel wrote a label named unlabel
        code_, out, err = run(f"-i {tmp_path / 'nothing.nwk'} --paml-codes {code}")
        assert code_ == 1 and out == "" and message in err and "nothing.nwk not found" not in err

    def test_paml_zero_is_the_background(self, run, tmp_path):
        # #0 was written as the label 0, and -F hyphy stopped because 0 is not a label name
        tree = tmp_path / "t.nwk"
        tree.write_text("((A #1,B #0)$1,(C,D) #0,E);")
        code, out, err = run(f"-i {tree} --paml-codes 1=Test")
        assert code == 0 and out == "((A{Test},B){Test},(C,D),E);\n"
        code, out, err = run(f"-i {tree} --paml-codes 1=Test -F hyphy")
        assert code == 0 and out == "((A{Test},B){Test},(C,D),E);\n"
        _, out, _ = run(f"-i {tree} --paml-codes 1=Test 0=Reference")
        assert out == "((A{Test},B{Reference}){Test},(C,D){Reference},E);\n"

    def test_paml_numbers_without_names_give_a_warning(self, run, tmp_path):
        # {1} and {2} were written without a word, and only -F hyphy said that they are not label names
        tree = tmp_path / "t.nwk"
        tree.write_text("(O #2,((M1 #2,(A #1,(B #1,C) #1)),(M2,(D #1,E #1))));\n(O #2,(A #1,B));\n")
        code, out, err = run(f"-i {tree}")
        assert code == 0 and trees_of(out)[0] == "(O{2},((M1{2},(A{1},(B{1},C){1})),(M2,(D{1},E{1}))));"
        message = "the PAML label numbers 1 and 2 have no names, and HyPhy cannot use a number as a label. Name them with --paml-codes"
        assert err.count(message) == 1                  # once for the run, not once for each tree
        _, _, err = run(f"-i {tree} --paml-codes 1=Test 2=Reference")
        assert "have no names" not in err
        _, _, err = run(f"-i {tree} --existing clear")
        assert "have no names" not in err

    def test_labels_in_the_tree_that_hyphy_cannot_use_give_a_warning(self, run, tmp_path):
        tree = tmp_path / "t.nwk"
        tree.write_text("((A{Internal},B{Test}),C);")
        code, out, err = run(f"-i {tree}")
        assert code == 0 and out == "((A{Internal},B{Test}),C);\n"
        assert "the label 'Internal' read from the tree collides with a HyPhy branch-set option" in err
        assert "'Test' read from the tree" not in err


class TestSeveralLabels:
    def test_a_label_for_each_selection(self, run):
        _, out, _ = run(f"-i {DATA / 'species.nwk'} -L Test:fitch {DIOECIOUS} -L Out:tips ^O$ --background Reference --stems unlabel")
        assert out.startswith("(O{Out},")
        assert "M1{Reference}" in out and "A{Test}" in out

    def test_a_spec_file(self, run, tmp_path):
        spec = tmp_path / "labels.tsv"
        spec.write_text(f"# label\trule\tselectors\toptions\nTest\tfitch\t{DIOECIOUS}\tambiguous=deltran\nOut\ttips\t^O$\n")
        code, out, err = run(f"-i {NAMED} --spec {spec} --internal-names keep")
        assert code == 0, err
        names = labeled_names(out)
        assert names["O"] == "Out" and set(n for n, label in names.items() if label == "Test") == {"A", "B", "D", "E", "DE"}

    def test_the_first_label_wins_with_first_precedence(self, run):
        _, out, err = run(f"-i {NAMED} -L One ^A$ --rule tips -L Two ^[AB]$ --rule tips --precedence first")
        assert labeled_names(out) == {"A": "One", "B": "Two"}
        assert "label Two: 1 of 2 branches went to another label" in err

    def test_error_precedence_stops_when_two_labels_claim_a_branch(self, run):
        code, _, err = run(f"-i {NAMED} -L One ^A$ --rule tips -L Two ^[AB]$ --rule tips --precedence error")
        assert code == 1 and "forbids" in err

    def test_selectors_that_match_nothing_are_reported(self, run):
        code, out, err = run(f"-i {NAMED} -L Test ^Q --rule tips")
        assert code == 0 and "the selectors choose no tips" in err and "no branch is labeled" in err
        assert "label Test: 1 listed name or pattern matches no tip (^Q)" in err
        assert labeled_names(out) == {}
        _, _, err = run(f"-i {NAMED} -L Test ^Q ^R ^A --rule tips")
        assert "label Test: 2 listed names or patterns match no tip (^Q, ^R)" in err

    def test_species_selections_work_on_trees_with_several_samples(self, run, tmp_path):
        tree = tmp_path / "gene.nwk"
        tree.write_text("((A_1,A_2),(B_1,C_1),D_1);")
        code, out, err = run(f"-i {tree} -L Test A B --select-on species -d _ --rule all")
        assert code == 0, err
        assert out == "((A_1{Test},A_2{Test}){Test},(B_1{Test},C_1),D_1);\n"

    def test_a_map_file_with_a_byte_order_mark_maps_its_first_row(self, run, tmp_path):
        tree = tmp_path / "gene.nwk"
        tree.write_text("((s1,s3),(s2,s4));")
        mapping = tmp_path / "map.txt"
        mapping.write_bytes("\ufeffs1 Amaranthus_pumilus\ns3 Amaranthus_pumilus\ns2 Beta_vulgaris\ns4 Beta_vulgaris\n".encode("utf-8"))
        code, out, err = run(f"-i {tree} -L Test Amaranthus_pumilus --select-on species -I {mapping} --rule all --background Reference --exact")
        assert code == 0, err
        assert out == "((s1{Test},s3{Test}){Test},(s2{Reference},s4{Reference}){Reference});\n"

    def test_species_selection_without_a_map_gives_a_warning(self, run):
        _, _, err = run(f"-i {NAMED} -L Test A --select-on species --rule tips")
        assert "every tip is its own species" in err


class TestSpeciesMapping:
    @pytest.fixture
    def gene(self, tmp_path):
        tree = tmp_path / "g.nwk"
        tree.write_text("((s1,s3),(s2,s4));")
        return tree

    def test_a_tip_missing_from_the_map_stops_the_run(self, run, gene, tmp_path):
        # s3 used to become a species of its own, so it was silently left out of Test
        mapping = tmp_path / "map.txt"
        mapping.write_text("s1 Amaranthus_pumilus\ns2 Beta_vulgaris\ns4 Beta_vulgaris\n")
        command = f"-i {gene} -L Test Amaranthus_pumilus --select-on species -I {mapping} --rule all --background Reference --exact"
        code, out, err = run(command)
        assert code == 1 and out == ""
        assert "1 tip name cannot be resolved to a species, so -U error stops the run\n  s3 (it is not resolved by the map file)\n" in err
        assert "give -U ignore" in err

    def test_ignore_leaves_the_tip_without_a_species(self, run, gene, tmp_path):
        mapping = tmp_path / "map.txt"
        mapping.write_text("s1 Amaranthus_pumilus\ns2 Beta_vulgaris\ns4 Beta_vulgaris\n")
        code, out, err = run(f"-i {gene} -L Test Amaranthus_pumilus --select-on species -I {mapping} --rule all --background Reference --exact -U ignore")
        assert code == 0, err
        assert out == "((s1{Test},s3{Reference}){Reference},(s2{Reference},s4{Reference}){Reference});\n"
        assert "1 tip name cannot be resolved to a species, so -U ignore leaves it without one and no selector chooses it. s3 (it is not resolved by the map file)" in err
        # its own name does not make it a species that a selector or --invert can choose
        _, out, _ = run(f"-i {gene} -L Test s3 --select-on species -I {mapping} --rule tips -U ignore")
        assert out == "((s1,s3),(s2,s4));\n"
        _, out, _ = run(f"-i {gene} -L Test Beta_vulgaris --invert --select-on species -I {mapping} --rule tips -U ignore")
        assert out == "((s1{Test},s3),(s2,s4));\n"

    def test_a_name_without_the_delimiter_is_not_resolved(self, run, tmp_path):
        tree = tmp_path / "g.nwk"
        tree.write_text("(O_1,((M1_1,((A_1,A_2),(B_1,C_1))),(M2_1,(D_1,Ee))));")
        code, _, err = run(f"-i {tree} -L Test ^[ABDE]$ --select-on species -d _ --rule all")
        assert code == 1 and "  Ee (it is not resolved by the delimiter)" in err
        code, out, err = run(f"-i {tree} -L Test ^[ABDE]$ --select-on species -d _ --rule all -U ignore")
        assert code == 0 and "Ee (it is not resolved by the delimiter)" in err
        assert "(D_1{Test},Ee)" in out

    def test_every_name_missing_from_a_map_is_listed(self, run, data_dir, tmp_path):
        mapping = tmp_path / "map.txt"
        mapping.write_text("A_1 A\nB_1 B\n")
        code, _, err = run(f"-i {data_dir / 'gene.nwk'} -L Test ^A$ --select-on species -I {mapping} --rule tips")
        assert code == 1 and "7 tip names cannot be resolved to a species" in err
        listed = [line.split()[0] for line in err.splitlines() if line.startswith("  ")]
        assert listed == ["O_1", "M1_1", "A_2", "C_1", "M2_1", "D_1", "E_1"]

    def test_a_species_list_resolves_names_written_in_several_ways(self, run, tmp_path):
        tree = tmp_path / "mixed.nwk"
        tree.write_text("((Amaranthus_pumilus_s1,AmaranthusPumilus-s2),(Beta_vulgaris_s1,A_hybridus_s1));")
        species = tmp_path / "species.txt"
        species.write_text("Amaranthus_pumilus\nAmaranthus_hybridus\nBeta_vulgaris\n")
        command = f"-i {tree} -L Test ^Amaranthus --select-on species -V {species} --rule tips"
        code, _, err = run(command)
        assert code == 1 and "A_hybridus_s1 (only an abbreviated genus matches it to Amaranthus_hybridus, which --accept-inferred accepts)" in err
        code, out, err = run(command + " --accept-inferred")
        assert code == 0, err
        assert out == "((Amaranthus_pumilus_s1{Test},AmaranthusPumilus-s2{Test}),(Beta_vulgaris_s1,A_hybridus_s1{Test}));\n"
        # without case folding AmaranthusPumilus-s2 is not resolved
        code, out, err = run(command + " --accept-inferred -N none -U ignore")
        assert out == "((Amaranthus_pumilus_s1{Test},AmaranthusPumilus-s2),(Beta_vulgaris_s1,A_hybridus_s1{Test}));\n"
        assert "AmaranthusPumilus-s2 (it matches no species of the species list, the closest are Amaranthus_pumilus" in err

    def test_a_species_list_that_is_not_text_is_a_clear_error(self, run, gene, tmp_path):
        species = tmp_path / "species.txt"
        species.write_bytes(b"Amaranthus_pumilus\n\xff\xfe\nBeta_vulgaris\n")
        code, out, err = run(f"-i {gene} -L Test A --select-on species -V {species} --rule tips")
        assert code == 1 and out == ""
        assert "species.txt is not valid UTF-8 text" in err and "Traceback" not in err

    def test_ambiguous_names_and_map_rows_are_checked_against_the_species_list(self, run, tmp_path):
        tree = tmp_path / "g.nwk"
        tree.write_text("((A_pumilus_s1,Beta_vulgaris_s1),x1);")
        species = tmp_path / "species.txt"
        species.write_text("Amaranthus_pumilus\nAtriplex_pumilus\nBeta_vulgaris\n")
        code, _, err = run(f"-i {tree} -L Test Beta_vulgaris --select-on species -V {species} --accept-inferred --rule tips")
        assert code == 1 and "A_pumilus_s1 (it matches Atriplex_pumilus and Amaranthus_pumilus equally well)" in err
        mapping = tmp_path / "map.txt"
        mapping.write_text("A_pumilus_s1 Amaranthus_pumilus\nx1 Beta_vulgarus\n")
        code, _, err = run(f"-i {tree} -L Test Beta_vulgaris --select-on species -V {species} -I {mapping} --rule tips")
        assert code == 1 and "species 'Beta_vulgarus' is not in the species tree" in err
        mapping.write_text("A_pumilus_s1 Amaranthus_pumilus\nx1 Beta_vulgaris\n")
        code, out, err = run(f"-i {tree} -L Test Beta_vulgaris --select-on species -V {species} -I {mapping} --rule tips")
        assert code == 0 and out == "((A_pumilus_s1,Beta_vulgaris_s1{Test}),x1{Test});\n"

    def test_mapping_options_that_have_no_effect_are_reported(self, run, data_dir, tmp_path):
        gene = data_dir / "gene.nwk"
        _, _, err = run(f"-i {gene} -L Test ^A$ -d _ -N none --rule tips")
        assert "-d and -N have no effect without --select-on species" in err
        for extra in ("--accept-inferred", "-N none", "-N digits"):
            _, _, err = run(f"-i {gene} -L Test A --select-on species -d _ --rule tips {extra}")
            assert "no effect without a species list given with -V" in err, extra
        _, _, err = run(f"-i {gene} -L Test A --select-on species -x ^(?P<species>[^_]+)_ -di 1 --rule tips")
        assert "-di has no effect without -d" in err
        species = tmp_path / "species.txt"
        species.write_text("A\nB\nC\nD\nE\nO\nM1\nM2\n")
        _, _, err = run(f"-i {gene} -L Test A --select-on species -d _ -V {species} -N none --rule tips")
        assert "no effect" not in err

    def test_parsimony_options_that_change_nothing_are_reported(self, run):
        code, _, err = run(f"-i {NAMED} -L Test:all {DIOECIOUS} --losses Rev --ambiguous Amb --gain-cost 5")
        assert code == 0 and "no label uses the fitch rule, so --losses, --ambiguous and --gain-cost have no effect" in err
        _, _, err = run(f"-i {NAMED} -L Test:fitch {DIOECIOUS} --unknown Cx")
        assert "--unknown has 1 listed name or pattern that matches no tip" in err and "no effect" not in err

    def test_an_unknown_tip_takes_the_state_of_the_history_inferred_from_the_other_tips(self, run):
        _, out, _ = run(f"-i {NAMED} -L Test {DIOECIOUS} --rule fitch --unknown C --internal-names keep")
        assert labeled_names(out)["BC"] == "Test" and labeled_names(out)["C"] == "Test"
        _, out, _ = run(f"-i {NAMED} -L Test {DIOECIOUS} --rule fitch --internal-names keep")
        assert "C" not in labeled_names(out)


class TestHyphyNames:
    def test_names_are_made_safe_and_the_table_can_be_saved(self, run, tmp_path):
        tree = tmp_path / "t.nwk"
        tree.write_text("(AmaranthusPumilus-Sample1,Beta,(C.1,D));")
        table = tmp_path / "names.tsv"
        code, out, err = run(f"-i {tree} -F hyphy --names-table {table}")
        assert code == 0, err
        assert out == "(AmaranthusPumilus_Sample1,Beta_,(C_1,D));\n"
        rows = table_of(table.read_text())
        assert {r["raw"]: r["written"] for r in rows}["AmaranthusPumilus-Sample1"] == "AmaranthusPumilus_Sample1"
        assert "3 tip names changed to be safe for hyphy" in err

    def test_colliding_names_stop_the_run_unless_a_suffix_is_added(self, run, tmp_path):
        tree = tmp_path / "t.nwk"
        tree.write_text("(A-1,A.1,(B,C));")
        code, _, err = run(f"-i {tree} -F hyphy")
        assert code == 1 and "collide" in err and "--sanitize-suffix" in err
        code, out, _ = run(f"-i {tree} -F hyphy --sanitize-suffix")
        assert code == 0 and out == "(A_1,A_1_2,(B,C));\n"

    def test_the_names_table_is_only_for_hyphy(self, run, tmp_path):
        table = tmp_path / "names.tsv"
        _, _, err = run(f"-i {NAMED} --names-table {table}")
        assert "--names-table is only used with -F hyphy" in err and not table.exists()

    def test_label_names_must_suit_hyphy(self, run):
        code, _, err = run(f"-i {NAMED} -L Unlabeled A --rule tips")
        assert code == 1 and "invalid label name" in err and "Unlabeled" in err


class TestErrors:
    def test_unknown_rule(self, run):
        code, out, err = run(f"-i {NAMED} -L Test A --rule bogus")
        assert code == 1 and out == "" and "unknown rule 'bogus'" in err and "Traceback" not in err

    def test_a_label_needs_selectors(self, run):
        code, _, err = run(f"-i {NAMED} -L Test")
        assert code == 1 and "needs a name and at least one selector" in err

    def test_label_names_are_checked_before_any_tree_is_read(self, run, tmp_path):
        code, _, err = run(f"-i {tmp_path / 'nothing.nwk'} -L Test {DIOECIOUS} --stems 9bad")
        assert code == 1 and "invalid label name" in err and "'9bad'" in err and "nothing.nwk not found" not in err

    def test_missing_input(self, run, tmp_path):
        code, _, err = run(f"-i {tmp_path / 'nothing.nwk'}")
        assert code == 1 and "nothing.nwk not found" in err and "Traceback" not in err
        assert "no label sets were given" not in err

    def test_a_byte_order_mark_is_ignored(self, run, tmp_path):
        # the mark made the whole tree one tip, so no tip was selected and the tree could not be written
        path = tmp_path / "bom.nwk"
        path.write_bytes(b"\xef\xbb\xbf" + SPECIES.encode() + b"\n")
        code, out, err = run(f"-i {path} -L Test ^[ABDE]$")
        assert code == 0, err
        assert out == "(O,((M1,(A{Test},(B{Test},C))),(M2,(D{Test},E{Test}){Test})));\n"
        code, out, _ = run(f"-i {path} -L Test ^[ABDE]$ -F hyphy")
        assert code == 0 and "{Test}" in out

    @pytest.mark.parametrize("text, message", [
        ("((A{Test},B){Test}95,C);\n", "unexpected text '95' at position 18 of the tree, after a label or comment"),
        ("((Pop#12a,B),C);\n", "unexpected text 'a' at position 8 of the tree, after a label or comment"),
        ("locus1 ((A,B),C);\n", "unexpected text 'locus1' at position 0 of the tree, before its first '('"),
        ("begin trees;\ntree t1 = ((A,B),C);\nend;\n", "NEXUS format is not supported"),
        ("((Homo sapiens:1,B:1):1,C:1);\n", "tree 1 has a space in 'Homo sapiens'"),
    ])
    def test_text_in_the_wrong_place_is_a_clear_error(self, run, tmp_path, text, message):
        path = tmp_path / "x.nwk"
        path.write_text(text)
        for options in ("", "-F hyphy", "-F table"):
            code, out, err = run(f"-i {path} {options}")
            assert code == 1 and out == "" and message in err and "Traceback" not in err

    def test_input_that_is_not_a_tree(self, run, tmp_path):
        path = tmp_path / "x.nwk"
        path.write_text("this is not a tree\n")
        code, _, err = run(f"-i {path}")
        assert code == 1 and "has no final ';'" in err and "Traceback" not in err

    def test_a_broken_tree_names_the_file(self, run, tmp_path):
        path = tmp_path / "x.nwk"
        path.write_text("((A,B),C;\n")
        code, _, err = run(f"-i {path}")
        assert code == 1 and "x.nwk" in err and "Traceback" not in err

    def test_bad_costs(self, run):
        code, _, err = run(f"-i {NAMED} -L Test {DIOECIOUS} --rule fitch --gain-cost abc")
        assert code == 1 and "invalid change cost" in err
        code, _, err = run(f"-i {NAMED} -L Test {DIOECIOUS} --rule fitch --loss-cost -1")
        assert code == 1 and "cannot be negative" in err

    def test_a_mode_written_in_capitals_is_not_taken_as_a_label_name(self, run):
        code, out, err = run(f"-i {NAMED} -L Test {DIOECIOUS} --rule fitch --ambiguous Deltran")
        assert code == 1 and out == "" and "looks like the mode 'deltran'" in err and "Traceback" not in err

    def test_very_large_costs_work_and_zero_costs_are_refused(self, run):
        code, out, err = run(f"-i {NAMED} -L Test {DIOECIOUS} --rule fitch --gain-cost 1e400 --internal-names keep")
        assert code == 0 and "Traceback" not in err and "parsimony cost 4" in err
        assert set(labeled_names(out)) == {"A", "B", "BC", "ABC", "D", "E", "DE", "X", "Y", "Z"}
        code, out, err = run(f"-i {NAMED} -L Test {DIOECIOUS} --rule fitch --gain-cost 0 --loss-cost 0")
        assert code == 1 and out == "" and "costs nothing to gain or to lose" in err

    def test_expensive_gains_are_avoided_by_gaining_once_above_the_ingroup(self, run):
        _, out, _ = run(f"-i {NAMED} -L Test {DIOECIOUS} --rule fitch --gain-cost 3 --internal-names keep")
        assert set(labeled_names(out)) == {"A", "B", "BC", "ABC", "D", "E", "DE", "X", "Y", "Z"}

    def test_the_state_of_the_root_can_be_fixed(self, run, tmp_path):
        tree = tmp_path / "t.nwk"
        tree.write_text("((A,B)AB,C);")
        _, out, _ = run(f"-i {tree} -L Test ^[AC]$ --rule fitch --internal-names keep")
        assert set(labeled_names(out)) == {"A", "AB", "C"}          # the trait is at the root, and B lost it
        _, out, _ = run(f"-i {tree} -L Test ^[AC]$ --rule fitch --root-state 0 --internal-names keep")
        assert set(labeled_names(out)) == {"A", "C"}


class TestChoices:
    """The parser enforces every option with a fixed set of choices, so run_tree_label does not check them again."""

    @pytest.mark.parametrize("bad", [
        "-F fasta", "--select-on tip", "--leaf-nodes drop", "--precedence some", "--existing erase",
        "--root-state 2", "--internal-names all", "-U skip",
    ])
    def test_a_choice_that_does_not_exist_is_refused_by_the_parser(self, run, bad):
        with pytest.raises(SystemExit) as info:
            run(f"-i {NAMED} {bad}")
        assert info.value.code == 2

    def test_every_choice_that_exists_is_accepted(self, run):
        for good in ("-F newick", "--select-on name", "--leaf-nodes label", "--precedence last", "--existing clear", "--root-state 1", "--internal-names auto", "-U ignore"):
            code, _, err = run(f"-i {NAMED} {good}")
            assert code == 0, (good, err)


class TestDeepTrees:
    def test_a_deep_tree_is_labeled_and_written(self, run, tmp_path):
        text = "T0:1"
        for i in range(1, 400):
            text = f"({text},T{i}:1):1"
        tree = tmp_path / "deep.nwk"
        tree.write_text(text.rsplit(":", 1)[0] + ";")
        table = tmp_path / "deep.tsv"
        code, out, err = run(f"-i {tree} -L Test ^T[0-9]*[05]$ --rule fitch -r T200 -F hyphy -B {table}")
        assert code == 0, err
        # 80 tips are selected, and the label of T200 is dropped when the two branches at the root become one
        assert out.count("{Test}") == 79
        assert len(table_of(table.read_text())) == 2 * 400 - 2
