"""Tests for tip name and label name rules (twig.labeling.core.sanitize)."""

import pytest

from twig.utils import TwigError
from twig.labeling.core import sanitize
from twig.labeling.core.sanitize import NameTable


class TestHyphyIdentifier:
    @pytest.mark.parametrize("raw, expected", [
        ("AmaranthusPumilus-Sample1", "AmaranthusPumilus_Sample1"),
        ("A.pumilus s1", "A_pumilus_s1"),
        ("a--b", "a_b"),                    # a run of invalid characters gives one underscore
        ("a_-b", "a_b"),                    # an invalid character after an underscore is dropped
        ("a-_b", "a__b"),
        ("-A", "_A"),
        ("1abc", "1abc"),                   # HyPhy allows a leading digit in tip names
        ("OG0001|Beta_vulgaris", "OG0001_Beta_vulgaris"),
        ("Ñandú", "_and_"),
    ])
    def test_conversion_matches_hyphy(self, raw, expected):
        assert sanitize.hyphy_identifier(raw) == expected

    def test_valid_names_are_unchanged_and_conversion_is_idempotent(self):
        for name in ("Beta_vulgaris", "A1", "_x", "OG0001_Amaranthus_pumilus_s1"):
            assert sanitize.hyphy_identifier(name) == name
        once = sanitize.sanitize_name("A.b-c d", "hyphy")
        assert sanitize.sanitize_name(once, "hyphy") == once

    def test_reserved_function_names_get_an_underscore(self):
        assert sanitize.sanitize_name("Beta", "hyphy") == "Beta_"
        assert sanitize.sanitize_name("beta", "hyphy") == "beta"      # HyPhy compares names with case
        assert sanitize.sanitize_name("Beta_vulgaris", "hyphy") == "Beta_vulgaris"
        assert "Beta" in sanitize.HYPHY_RESERVED and "Time" in sanitize.HYPHY_RESERVED

    def test_empty_name_is_an_error(self):
        with pytest.raises(TwigError, match="empty name"):
            sanitize.sanitize_name("", "hyphy")

    def test_paml_rule_has_a_length_limit(self):
        assert sanitize.sanitize_name("A-b", "paml") == "A_b"
        with pytest.raises(TwigError, match="longer than the 96 characters"):
            sanitize.sanitize_name("x" * 97, "paml")

    def test_unknown_rule(self):
        with pytest.raises(TwigError, match="unknown name rule"):
            sanitize.sanitize_name("A", "mrbayes")
        with pytest.raises(TwigError, match="unknown name rule"):
            NameTable("mrbayes")


class TestNameTable:
    def build(self, names, **kwargs):
        table = NameTable("hyphy", **kwargs)
        table.update(names)
        table.build()
        return table

    def test_tree_and_alignment_names_are_identical(self):
        table = self.build(["AmaranthusPumilus-Sample1", "Beta_vulgaris", "Beta"])
        assert table.written("AmaranthusPumilus-Sample1") == "AmaranthusPumilus_Sample1"
        assert table.written("Beta_vulgaris") == "Beta_vulgaris"
        assert table.written("Beta") == "Beta_"
        assert table.changed() == {"AmaranthusPumilus-Sample1": "AmaranthusPumilus_Sample1", "Beta": "Beta_"}

    def test_the_table_is_reversible(self):
        table = self.build(["A-b", "C"])
        assert table.raw(table.written("A-b")) == "A-b"
        with pytest.raises(TwigError, match="not in the name table"):
            table.raw("zzz")

    def test_collisions_stop_the_run_with_a_list(self):
        names = ["A-b", "A_b", "X.y", "X_y", "ok"]
        table = NameTable("hyphy")
        table.update(names)
        with pytest.raises(TwigError) as info:
            table.build()
        message = str(info.value)
        assert "A_b and A-b" in message and "X_y and X.y" in message
        assert "ok" not in message
        assert "--sanitize-suffix" in message

    def test_collisions_are_compared_without_regard_to_case(self):
        table = NameTable("hyphy")
        table.update(["Abc", "abc"])
        with pytest.raises(TwigError, match="without regard to case"):
            table.build()

    def test_suffix_keeps_valid_names_and_numbers_the_rest(self):
        table = self.build(["A_b", "A-b", "A.b"], suffix=True)
        assert table.written("A_b") == "A_b"          # already valid, so it keeps its spelling
        assert table.written("A-b") == "A_b_2"
        assert table.written("A.b") == "A_b_3"
        # the result does not depend on the order in which names were added
        again = self.build(["A.b", "A-b", "A_b"], suffix=True)
        assert again.build() == table.build()

    def test_suffix_skips_names_that_are_already_taken(self):
        table = self.build(["A_b", "A-b", "A_b_2"], suffix=True)
        assert len({v.lower() for v in table.build().values()}) == 3
        assert table.written("A_b_2") == "A_b_2"

    def test_names_must_be_added_before_building(self):
        table = self.build(["A"])
        with pytest.raises(TwigError, match="was not added"):
            table.written("B")

    def test_write_tsv(self, tmp_path):
        table = self.build(["A-b", "C"])
        path = tmp_path / "sanitized_names.tsv"
        table.write_tsv(path)
        assert path.read_text().splitlines() == ["raw\twritten\tchanged", "A-b\tA_b\tyes", "C\tC\tno"]
        with pytest.raises(TwigError, match="use --force"):
            table.write_tsv(path)


class TestLabelNames:
    @pytest.mark.parametrize("label", ["Test", "Reference", "Rev", "Nuisance", "Test2", "Stem_1"])
    def test_valid(self, label):
        assert sanitize.label_name_problem(label) is None
        sanitize.check_label_names([label])

    @pytest.mark.parametrize("label", ["1Test", "_Test", "Te-st", "Te st", "", "Tëst"])
    def test_pattern(self, label):
        assert "must start with a letter" in sanitize.label_name_problem(label)

    @pytest.mark.parametrize("label", ["All", "all", "INTERNAL", "Leaves", "unlabeled", "Unlabeled", "fg"])
    def test_hyphy_branch_set_options_in_any_case(self, label):
        assert "HyPhy branch-set option" in sanitize.label_name_problem(label)

    @pytest.mark.parametrize("label", ["keep", "Unlabel", "prune", "Relabel", "mixed", "Conflict", "paralog", "Masked", "unmapped"])
    def test_twig_words_in_any_case(self, label):
        assert "twig uses" in sanitize.label_name_problem(label)

    def test_check_lists_every_problem(self):
        with pytest.raises(TwigError) as info:
            sanitize.check_label_names(["Test", "All", "1x", "keep"])
        message = str(info.value)
        assert "'All'" in message and "'1x'" in message and "'keep'" in message and "'Test'" not in message
