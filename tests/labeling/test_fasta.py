"""Tests for FASTA reading and writing (twig.labeling.core.fasta)."""

import io

import pytest

from twig.utils import TwigError
from twig.labeling.core import fasta


class TestParse:
    def test_header_is_the_first_token_and_sequences_can_wrap(self):
        text = ">seq1 some description\nACGT\nAC GT\n\n>seq2\r\nTTTT\r\n"
        records = fasta.parse_fasta(text)
        assert list(records) == ["seq1", "seq2"]
        assert records == {"seq1": "ACGTACGT", "seq2": "TTTT"}

    @pytest.mark.parametrize("text, message", [
        (">a\nAC\n>a\nGT\n", r"<text>:3: duplicate FASTA header 'a'"),
        (">\nACGT\n", r"<text>:1: empty FASTA header"),
        ("ACGT\n>a\nAC\n", r"<text>:1: sequence found before any FASTA header"),
        ("\n\n", r"no FASTA records"),
        (">a\n>b\nACGT\n", r"<text>:1: FASTA record 'a' has no sequence"),
        (">a\nACGT\n>b\n\n", r"<text>:3: FASTA record 'b' has no sequence"),
    ])
    def test_strict_errors(self, text, message):
        with pytest.raises(TwigError, match=message):
            fasta.parse_fasta(text)

    def test_unequal_lengths_are_an_error_only_when_aligned(self):
        text = ">a\nACGT\n>b\nAC\n"
        assert fasta.parse_fasta(text) == {"a": "ACGT", "b": "AC"}
        with pytest.raises(TwigError, match="different lengths"):
            fasta.parse_fasta(text, aligned=True)

    def test_read_fasta_and_headers(self, tmp_path):
        path = tmp_path / "x.fa"
        path.write_text(">a first\nACGT\n>b\nAC--\n>a\nTT\n")
        assert fasta.read_headers(path) == ["a", "b", "a"]      # duplicates are kept for the caller
        with pytest.raises(TwigError, match="duplicate FASTA header"):
            fasta.read_fasta(path)

    def test_a_byte_order_mark_is_ignored(self, tmp_path):
        path = tmp_path / "x.fa"
        path.write_bytes("\ufeff>a\nACGT\n>b\nAC--\n".encode("utf-8"))
        assert fasta.read_headers(path) == ["a", "b"]
        assert fasta.read_fasta(path) == {"a": "ACGT", "b": "AC--"}

    def test_both_readers_find_the_same_headers(self, tmp_path):
        # an indented header was a record for read_fasta but invisible to read_headers
        path = tmp_path / "x.fa"
        path.write_text(" >a\nAC\n  >b desc\nGT\n")
        assert fasta.read_fasta(path) == {"a": "AC", "b": "GT"}
        assert fasta.read_headers(path) == ["a", "b"]

    @pytest.mark.parametrize("text, message", [
        (">a\n>b\nACGT\n", r"x.fa:1: FASTA record 'a' has no sequence"),
        (">a\nACGT\n>b\n", r"x.fa:3: FASTA record 'b' has no sequence"),
        (">a\nACGT\n>b\n\n  \n>c\nAC\n", r"x.fa:3: FASTA record 'b' has no sequence"),
    ])
    def test_both_readers_refuse_a_record_without_a_sequence(self, tmp_path, text, message):
        path = tmp_path / "x.fa"
        path.write_text(text)
        with pytest.raises(TwigError, match=message):
            fasta.read_headers(path)
        with pytest.raises(TwigError, match=message):
            fasta.read_fasta(path)

    def test_header_scan_errors(self, tmp_path):
        empty = tmp_path / "empty.fa"
        empty.write_text("ACGT\n")
        with pytest.raises(TwigError, match="no FASTA records"):
            fasta.read_headers(empty)
        with pytest.raises(TwigError, match="cannot read"):
            fasta.read_headers(tmp_path / "missing.fa")


class TestWrite:
    def test_format(self):
        assert fasta.format_fasta({"a": "ACGTAC", "b": "TT"}) == ">a\nACGTAC\n>b\nTT\n"
        assert fasta.format_fasta({"a": "ACGTAC"}, width=4) == ">a\nACGT\nAC\n"
        assert fasta.format_fasta({}) == ""

    def test_write_to_handle_and_atomically_to_a_path(self, tmp_path):
        out = io.StringIO()
        fasta.write_fasta({"a": "AC"}, out)
        assert out.getvalue() == ">a\nAC\n"
        path = tmp_path / "out" / "x.fa"
        fasta.write_fasta({"a": "AC"}, path)
        assert path.read_text() == ">a\nAC\n"
        assert not (tmp_path / "out" / "x.fa.tmp").exists()
        with pytest.raises(TwigError, match="use --force"):
            fasta.write_fasta({"a": "GG"}, path)
        fasta.write_fasta({"a": "GG"}, path, force=True)
        assert path.read_text() == ">a\nGG\n"

    def test_subset_rows_follows_the_requested_order(self):
        records = {"a": "AC", "b": "GT", "c": "TT"}
        assert list(fasta.subset_rows(records, ["c", "a"])) == ["c", "a"]
        with pytest.raises(TwigError, match="'z' is not in the alignment"):
            fasta.subset_rows(records, ["z"])

    def test_round_trip(self):
        records = {"a": "ATG---CCC", "b": "ATGAAACCC"}
        assert fasta.parse_fasta(fasta.format_fasta(records, width=4), aligned=True) == records


class TestColumns:
    RECORDS = {"a": "AT--GC--", "b": "AT--G---", "c": "A---GC--"}

    def test_drop_empty_sites(self):
        out = fasta.drop_empty_columns(self.RECORDS, "site")
        assert out == {"a": "ATGC", "b": "ATG-", "c": "A-GC"}

    def test_none_returns_a_copy(self):
        out = fasta.drop_empty_columns(self.RECORDS, "none")
        assert out == self.RECORDS and out is not self.RECORDS

    def test_codon_mode_keeps_the_reading_frame(self):
        records = {"a": "ATG---CCC", "b": "ATG---C-C"}
        # the middle codon is empty in every row, but the last codon is not
        assert fasta.drop_empty_columns(records, "codon") == {"a": "ATGCCC", "b": "ATGC-C"}
        # a codon with data in some row keeps all three columns
        records = {"a": "ATG-T-CCC", "b": "ATG---CCC"}
        assert fasta.drop_empty_columns(records, "codon") == records
        # site mode would break the frame
        assert fasta.drop_empty_columns(records, "site")["a"] == "ATGTCCC"

    def test_missing_data_characters_count_as_empty(self):
        assert fasta.drop_empty_columns({"a": "A?.", "b": "A-."}, "site") == {"a": "A", "b": "A"}

    def test_errors(self):
        with pytest.raises(TwigError, match="not a multiple of 3"):
            fasta.drop_empty_columns({"a": "ATGC"}, "codon")
        with pytest.raises(TwigError, match="different lengths"):
            fasta.drop_empty_columns({"a": "AT", "b": "A"}, "site")
        with pytest.raises(TwigError, match="unknown mode"):
            fasta.drop_empty_columns({"a": "AT"}, "all")


class TestSequenceStatistics:
    def test_occupancy_counts_only_unambiguous_bases(self):
        assert fasta.occupancy("ACGTUacgtu-NNRY?") == 10

    def test_length_counts_ambiguity_codes_but_not_gaps(self):
        assert fasta.nongap_length("AC-NNR?.") == 5


class TestCodons:
    STOPS = fasta.stop_codons(1)

    def test_stop_codon_tables(self):
        assert fasta.stop_codons(1) == {"TAA", "TAG", "TGA"}
        assert fasta.stop_codons(2) == {"TAA", "TAG", "AGA", "AGG"}
        with pytest.raises(TwigError, match="not supported"):
            fasta.stop_codons(99)

    def test_find_stops_is_in_frame_and_ignores_gaps(self):
        assert fasta.find_stop_codons("ATGTAAGGGTGA", self.STOPS) == (1, 3)
        assert fasta.find_stop_codons("ATGGTAAGGG", self.STOPS) == ()      # TAA is out of frame
        assert fasta.find_stop_codons("---TAA---", self.STOPS) == (1,)
        assert fasta.find_stop_codons("atgtaa", self.STOPS) == (1,)
        assert fasta.find_stop_codons("AUGUAA", self.STOPS) == (1,)        # RNA

    def test_check_codons(self):
        report = fasta.check_codons("ATGAAATAG", self.STOPS)
        assert report.stops == (2,) and not report.bad_length and not report.ok
        report = fasta.check_codons("ATGAAAC", self.STOPS)
        assert report.bad_length and report.stops == ()
        report = fasta.check_codons("ATGXAA", self.STOPS)
        assert report.bad_chars == ("X",)
        assert fasta.check_codons("ATG---CCC", self.STOPS).ok

    def test_mask_stop_codons(self):
        assert fasta.mask_stop_codons("ATGTAAGGGTGA", self.STOPS) == "ATG---GGG---"
        assert fasta.mask_stop_codons("ATGTAA", self.STOPS, fill="NNN") == "ATGNNN"
        with pytest.raises(TwigError, match="three characters"):
            fasta.mask_stop_codons("ATGTAA", self.STOPS, fill="N")
