"""Tests for the listing of labeled branches and for issues (twig.labeling.core.display, validate)."""

import pytest

pytest.importorskip("toytree")

from twig.labeling.core import display, treeio
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

    def test_a_long_list_of_unlabeled_branches_is_cut(self):
        names = ",".join(f"t{i}" for i in range(50))
        lines = listing(f"({names});")
        assert lines[-2] == "  and 20 more unlabeled branches"

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
