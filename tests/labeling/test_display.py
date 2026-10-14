"""Tests for the listing of labeled branches and for issues (twig.labeling.core.display, validate)."""

import pytest

pytest.importorskip("toytree")

from twig.labeling.core import display, rooting, treeio
from twig.labeling.core.labelsets import LabelingConfig, LabelSet, apply_label_sets
from twig.labeling.core.validate import ERROR, WARNING, Issue, errors, warnings
from conftest import LABELED, SPECIES


def listing(newick=SPECIES, sets=None, config=None, **kwargs):
    tree = treeio.read_tree(newick)
    result = apply_label_sets(tree, sets or [], config or LabelingConfig())
    return display.branch_listing(tree, result.labels, result.kind, result.source, **kwargs)


def designed():
    sets = [LabelSet("Test", "fitch", ("^[ABDE]$",))]
    config = LabelingConfig(stems="unlabel", losses="Rev", background="Reference")
    return listing(SPECIES, sets, config)


class TestCount:
    @pytest.mark.parametrize("n, noun, expected", [
        (1, "tip", "1 tip"),
        (2, "tip", "2 tips"),
        (0, "file", "0 files"),
        (2, "branch", "2 branches"),
        (1, "species", "1 species"),
        (2, "species", "2 species"),
        (3, "FASTA file", "3 FASTA files"),
    ])
    def test_nouns(self, n, noun, expected):
        assert display.count(n, noun) == expected

    def test_a_plural_can_be_given(self):
        assert display.count(2, "listed name or pattern", "listed names or patterns") == "2 listed names or patterns"
        assert display.count(1, "listed name or pattern", "listed names or patterns") == "1 listed name or pattern"


class TestBranchListing:
    def test_labeled_branches_come_first_clade_by_clade(self):
        lines = designed()
        assert lines[0].startswith("labeled branches, clade by clade")
        rows = [line.split() for line in lines[1:lines.index("unlabeled branches")]]
        assert rows[0] == ["Reference", "tip", "O"]
        assert ["Test", "crown", "B+C"] in rows
        assert ["Rev", "loss", "C"] in rows
        assert len(rows) == 12

    def test_unlabeled_stems_are_listed(self):
        lines = designed()
        unlabeled = [line.split() for line in lines[lines.index("unlabeled branches") + 1:-1]]
        assert unlabeled == [["-", "stem", "A+B+C"], ["-", "stem", "D+E"]]

    def test_summary_counts_each_label(self):
        assert designed()[-1] == "summary  6 Reference, 5 Test, 1 Rev, 2 unlabeled"

    def test_nothing_labeled(self):
        lines = listing()
        assert lines[1] == "  none"
        assert lines[-1] == "summary  14 unlabeled"

    def test_columns_line_up(self):
        lines = designed()
        body = lines[1:lines.index("unlabeled branches")]
        assert len({len(line) - len(line.split()[-1]) for line in body}) == 1

    def test_long_clades_are_cut(self):
        names = ",".join(f"t{i}" for i in range(20))
        lines = listing(f"(({names}),x);", limit=5)
        assert any(line.endswith("t0+t1+t2+t3+t4+...") for line in lines)

    def test_species_are_listed_once_each(self):
        tree = treeio.read_tree("((A_1,A_2),(B_1,C_1));")
        text = display.clade_text(tree, lambda n: n.split("_")[0])
        assert sorted(text.values()) == ["A", "A", "A", "A+B+C", "B", "B+C", "C"]

    def test_repeated_species_do_not_count_toward_the_cut(self):
        # the cut counted the samples of A to E twice, so F, G and H were left out without +...
        tree = treeio.read_tree("(O_1,((A_1,B_1,C_1,D_1,E_1),(A_2,B_2,C_2,D_2,E_2),(F_1,G_1,H_1)));")
        species = lambda n: n.split("_")[0]
        stem = tree.treenode.children[1].idx
        assert display.clade_text(tree, species)[stem] == "A+B+C+D+E+F+G+H"
        assert display.clade_text(tree, species, limit=7)[stem] == "A+B+C+D+E+F+G+..."
        assert display.clade_text(tree, species, limit=3)[stem] == "A+B+C+..."
        assert display.clade_text(tree, species, limit=5)[tree.treenode.children[1].children[0].idx] == "A+B+C+D+E"

    def test_with_a_root_merge_the_tree_is_described_as_it_is_written(self):
        # the dissolved ingroup branch was listed as a Reference branch of its own, which HyPhy never sees
        tree = treeio.read_tree(SPECIES)
        sets = [LabelSet("Test", "fitch", ("^[ABDE]$",))]
        result = apply_label_sets(tree, sets, LabelingConfig(stems="unlabel", losses="Rev", background="Reference"))
        labeled = treeio.set_labels(tree, result.labels)
        merge = rooting.plan_root_merge(labeled)
        lines = display.branch_listing(labeled, result.labels, result.kind, result.source, merge=merge)
        assert "M1+A+B+C+M2+D+E" not in " ".join(line for line in lines if line.startswith("  "))
        assert lines[-2] == "root     the branches O and M1+A+B+C+M2+D+E at the root are written as one branch, listed as O"
        assert lines[-1] == "summary  5 Reference, 5 Test, 1 Rev, 2 unlabeled"

    def test_the_kept_root_branch_has_the_merged_label(self):
        tree = treeio.read_tree("(A{X},(B{Y},C{Y}){Z});")
        lines = display.branch_listing(tree, treeio.labeled_branches(tree), merge=rooting.plan_root_merge(tree))
        assert ["-", "tip", "A"] in [line.split() for line in lines]
        assert lines[-1] == "summary  2 Y, 1 unlabeled"

    def test_a_long_list_of_unlabeled_branches_is_cut(self):
        names = ",".join(f"t{i}" for i in range(50))
        lines = listing(f"({names});")
        assert lines[-2] == "  and 20 more unlabeled branches"
        lines = listing("(" + ",".join(f"t{i}" for i in range(31)) + ");")
        assert lines[-2] == "  and 1 more unlabeled branch"

    def test_roles_of_labeled_branches_read_from_a_tree_with_labels(self):
        tree = treeio.read_tree(LABELED)
        lines = display.branch_listing(tree, treeio.labeled_branches(tree))
        assert lines[-1] == "summary  6 Reference, 5 Test, 1 Rev, 2 unlabeled"
        assert lines[1].split() == ["Reference", "tip", "O"]


class TestIssue:
    def test_str_names_the_locus(self):
        assert str(Issue(ERROR, "write", "bad tree", "OG1")) == "OG1: bad tree"
        assert str(Issue(WARNING, "root-merge", "careful")) == "careful"

    def test_filters(self):
        issues = [Issue(ERROR, "a", "x"), Issue(WARNING, "b", "y"), Issue(ERROR, "c", "z")]
        assert [i.code for i in errors(issues)] == ["a", "c"]
        assert [i.code for i in warnings(issues)] == ["b"]

    def test_issues_cannot_be_changed(self):
        with pytest.raises(Exception):
            Issue(ERROR, "a", "x").message = "y"
