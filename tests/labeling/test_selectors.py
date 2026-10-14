"""Tests for tip selectors (twig.labeling.core.selectors)."""

import pytest

from twig.utils import TwigError
from twig.labeling.core.selectors import Selector, read_name_list

NAMES = ["A", "MA", "B", "C", "Beta_vulgaris", "Amaranthus_pumilus_s1"]


class TestSelector:
    def test_regex_values_are_searched(self):
        selector = Selector(["^A"])
        assert selector.select(NAMES) == ["A", "Amaranthus_pumilus_s1"]

    def test_unanchored_regex_also_matches_inside_names(self):
        assert Selector(["A"]).select(NAMES) == ["A", "MA", "Amaranthus_pumilus_s1"]

    def test_exact_values_do_not_match_substrings(self):
        assert Selector(["A"], exact=True).select(NAMES) == ["A"]

    def test_several_values_are_combined(self):
        assert Selector(["^B", "^C"]).select(NAMES) == ["B", "C", "Beta_vulgaris"]

    def test_list_file_with_names_regexes_and_comments(self, tmp_path):
        path = tmp_path / "dioecious.txt"
        path.write_text("# dioecious species\nA\n\n  B  \n~^Beta\n")
        names, patterns = read_name_list(path)
        assert names == ["A", "B"] and patterns == ["^Beta"]
        selector = Selector([path])
        assert selector.select(NAMES) == ["A", "B", "Beta_vulgaris"]
        assert not selector.matches("MA")        # list entries are exact names

    def test_a_byte_order_mark_is_not_part_of_the_first_name(self, tmp_path):
        path = tmp_path / "list.txt"
        path.write_bytes("\ufeffA\nB\n".encode("utf-8"))
        assert read_name_list(path) == (["A", "B"], [])
        path.write_bytes("\ufeff~^Beta\n".encode("utf-8"))
        assert read_name_list(path) == ([], ["^Beta"])

    def test_exact_flag_does_not_change_how_lists_are_read(self, tmp_path):
        path = tmp_path / "list.txt"
        path.write_text("~^Beta\n")
        assert Selector([path], exact=True).select(NAMES) == ["Beta_vulgaris"]

    def test_unused_entries_help_find_typos(self, tmp_path):
        path = tmp_path / "list.txt"
        path.write_text("A\nZZZ\n~^Q\n")
        selector = Selector([path])
        assert selector.unused(NAMES) == ["ZZZ", "~^Q"]
        # a pattern from the command line is shown as it was typed, without a ~
        assert Selector(["^Z"]).unused(NAMES) == ["^Z"]

    def test_empty_selector_matches_nothing(self):
        selector = Selector([])
        assert not selector and selector.select(NAMES) == []

    def test_invalid_regex_is_a_clear_error(self, tmp_path):
        with pytest.raises(TwigError, match="command line: invalid regular expression"):
            Selector(["(unclosed"])
        path = tmp_path / "list.txt"
        path.write_text("A\n~(\n")
        with pytest.raises(TwigError, match=r"list.txt:2: invalid regular expression"):
            Selector([path])

    def test_missing_list_is_treated_as_a_regex(self, tmp_path):
        # a value that is not an existing file is a pattern, so a typo in a path is not an error
        assert Selector([str(tmp_path / "missing.txt")]).select(NAMES) == []

    def test_describe(self):
        assert Selector(["^A"]).describe() == "0 names and 1 pattern from ^A"
