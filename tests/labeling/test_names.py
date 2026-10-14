"""Tests for name mapping (twig.labeling.core.names)."""

import io
from types import SimpleNamespace

import pytest

from twig.utils import TwigError
from twig.labeling.core.names import (
    NameMapper,
    Vocabulary,
    canonical,
    normalize_steps,
    parse_map_file,
    read_vocab_list,
    tokens,
    write_frozen_map,
    write_hit_table,
)

SPECIES = ["Amaranthus_pumilus", "Amaranthus_hybridus", "Atriplex_pumilus", "Beta_vulgaris"]


###########################################################################
# species lists
###########################################################################

class TestReadVocabList:
    def test_one_name_per_line_with_comments_and_blank_lines(self, tmp_path):
        path = tmp_path / "species.txt"
        path.write_text("# the species\nAmaranthus_pumilus  extra words\n\n  Beta_vulgaris\n")
        assert read_vocab_list(path) == ["Amaranthus_pumilus", "Beta_vulgaris"]

    def test_a_byte_order_mark_is_not_part_of_the_first_name(self, tmp_path):
        path = tmp_path / "species.txt"
        path.write_bytes("\ufeffAmaranthus_pumilus\nBeta_vulgaris\n".encode("utf-8"))
        assert read_vocab_list(path) == ["Amaranthus_pumilus", "Beta_vulgaris"]

    def test_an_empty_list_is_an_error(self, tmp_path):
        path = tmp_path / "species.txt"
        path.write_text("# nothing here\n\n")
        with pytest.raises(TwigError, match="no species names found"):
            read_vocab_list(path)

    def test_a_missing_file_is_an_error(self, tmp_path):
        with pytest.raises(TwigError, match="cannot read species list"):
            read_vocab_list(tmp_path / "missing.txt")

    def test_a_file_that_is_not_text_is_an_error_and_not_a_traceback(self, tmp_path):
        path = tmp_path / "species.txt"
        path.write_bytes(b"Amaranthus_pumilus\n\xff\xfe\n")
        with pytest.raises(TwigError, match="species.txt is not valid UTF-8 text"):
            read_vocab_list(path)


###########################################################################
# normalization
###########################################################################

class TestCanonical:
    def test_default_steps_unify_both_name_styles(self):
        assert canonical("A_pumilus_sample1") == "a_pumilus_sample1"
        assert canonical("AmaranthusPumilus-Sample1") == "amaranthus_pumilus_sample1"

    @pytest.mark.parametrize("raw, steps, expected", [
        ("Ñandú", ("accents",), "Nandu"),
        ("AmaranthusPumilus", ("camel",), "Amaranthus_Pumilus"),
        ("GCAAmaranthus", ("camel",), "GCA_Amaranthus"),
        ("DNAPolymerase", ("camel",), "DNA_Polymerase"),
        ("ABC", ("camel",), "ABC"),
        ("Sample1", ("digits",), "Sample_1"),
        ("1Sample", ("digits",), "1_Sample"),
        ("A.pumilus-s1", ("sep",), "A_pumilus_s1"),
        ("--A  b__c--", ("sep",), "A_b_c"),
        ("Amaranthus_Pumilus", ("case",), "amaranthus_pumilus"),
        ("Amaranthus-Pumilus", (), "Amaranthus-Pumilus"),
    ])
    def test_each_step(self, raw, steps, expected):
        assert canonical(raw, steps) == expected

    def test_digits_step_is_off_by_default(self):
        assert canonical("M1") == "m1"
        assert canonical("M1", ("digits", "sep", "case")) == "m_1"

    def test_steps_are_applied_in_a_fixed_order(self):
        assert canonical("AmaranthusPumilus", ("case", "camel", "sep")) == canonical("AmaranthusPumilus", ("camel", "sep", "case"))

    def test_normalize_steps(self):
        assert normalize_steps(None) == ("accents", "camel", "sep", "case")
        assert normalize_steps(["none"]) == ()
        assert normalize_steps(["case", "sep"]) == ("sep", "case")
        with pytest.raises(TwigError, match="unknown normalize"):
            normalize_steps(["upper"])
        with pytest.raises(TwigError, match="cannot be combined"):
            normalize_steps(["none", "case"])

    def test_tokens(self):
        assert tokens("amaranthus_pumilus_sample1") == ("amaranthus", "pumilus", "sample1")
        assert tokens("") == ()


###########################################################################
# map files
###########################################################################

class TestMapFile:
    def write(self, tmp_path, text):
        path = tmp_path / "map.tsv"
        path.write_text(text)
        return path

    def test_rows_comments_and_one_column_rows(self, tmp_path):
        path = self.write(tmp_path, (
            "# sample   species   specimen\n"
            "A_pumilus_sample1   Amaranthus_pumilus   s1\n"
            "\n"
            "AmaranthusPumilus-Sample2   Amaranthus_pumilus   # a trailing comment\n"
            "~^Ahyb[_-]   Amaranthus_hybridus\n"
            "Beta_vulgaris\n"
        ))
        rows = parse_map_file(path)
        assert [r.sample for r in rows] == ["A_pumilus_sample1", "AmaranthusPumilus-Sample2", "^Ahyb[_-]", "Beta_vulgaris"]
        assert rows[0].specimen == "s1" and rows[1].specimen is None
        assert rows[2].is_regex and rows[2].pattern.search("Ahyb_x1")
        assert rows[3].species == "Beta_vulgaris"      # a one-column row maps a name to itself
        assert rows[0].source.endswith("map.tsv:2")

    @pytest.mark.parametrize("text, message", [
        ("a b c d\n", r"map.tsv:1: expected 1 to 3 columns"),
        ("a b\na c\n", r"map.tsv:2: duplicate map row for 'a' \(first seen at line 1\)"),
        ("~[unclosed b\n", r"map.tsv:1: invalid regular expression"),
        ("~onlyregex\n", r"regular expression rows need a species column"),
        ("~ Beta_vulgaris\n", r"map.tsv:1: empty regular expression"),
        ("# nothing here\n\n", r"no map rows found"),
        ("s1,Amaranthus_pumilus\n", r"map.tsv:1: the row looks comma-separated"),
        ("s1, Amaranthus_pumilus\n", r"map.tsv:1: the row looks comma-separated"),
    ])
    def test_errors_name_the_line(self, tmp_path, text, message):
        with pytest.raises(TwigError, match=message):
            parse_map_file(self.write(tmp_path, text))

    def test_commas_inside_a_column_are_fine_otherwise(self, tmp_path):
        rows = parse_map_file(self.write(tmp_path, "s1,a Amaranthus_pumilus\n~^A{1,2}_ Beta_vulgaris\n"))
        assert [r.sample for r in rows] == ["s1,a", "^A{1,2}_"]

    def test_missing_file(self, tmp_path):
        with pytest.raises(TwigError, match="cannot read map file"):
            parse_map_file(tmp_path / "missing.tsv")

    def test_file_that_is_not_text(self, tmp_path):
        path = tmp_path / "map.tsv"
        path.write_bytes(b"\xff\xfe\x00")
        with pytest.raises(TwigError, match="not valid UTF-8"):
            parse_map_file(path)

    def test_a_byte_order_mark_is_not_part_of_the_first_sample(self, tmp_path):
        # Windows editors add the mark, which used to make the first row match no name
        path = tmp_path / "map.tsv"
        path.write_bytes("\ufeffAmaranthus_pumilus_s9 Atriplex_pumilus\nx1 Beta_vulgaris\n".encode("utf-8"))
        rows = parse_map_file(path)
        assert rows[0].sample == "Amaranthus_pumilus_s9"
        hit = NameMapper(map_rows=rows, vocabulary=SPECIES).resolve("Amaranthus_pumilus_s9")
        assert (hit.species, hit.method) == ("Atriplex_pumilus", "map")


###########################################################################
# vocabulary matching
###########################################################################

class TestMatching:
    @pytest.fixture
    def mapper(self):
        return NameMapper(vocabulary=SPECIES)

    def test_exact(self, mapper):
        hit = mapper.resolve("Beta_vulgaris")
        assert (hit.species, hit.method, hit.status, hit.leftover) == ("Beta_vulgaris", "exact", "ok", ())

    def test_canonical_form(self, mapper):
        for name in ("beta-vulgaris", "Beta.Vulgaris", "BetaVulgaris"):
            hit = mapper.resolve(name)
            assert (hit.species, hit.method, hit.status) == ("Beta_vulgaris", "canon", "ok"), name

    @pytest.mark.parametrize("name, leftover", [
        ("amaranthus pumilus sample1", ("sample1",)),
        ("s012 amaranthus pumilus", ("s012",)),
        ("og1 amaranthus pumilus g7", ("og1", "g7")),
        ("AmaranthusPumilus-Sample1", ("sample1",)),
        ("Amaranthus_pumilus_x1", ("x1",)),
    ])
    def test_tokens_anywhere_in_the_name(self, mapper, name, leftover):
        hit = mapper.resolve(name)
        assert (hit.species, hit.method, hit.status, hit.leftover) == ("Amaranthus_pumilus", "token", "ok", leftover)

    def test_both_name_styles_of_the_user_resolve_to_one_species(self):
        mapper = NameMapper(vocabulary=["Amaranthus_pumilus", "Beta_vulgaris"])
        a = mapper.resolve("A_pumilus_sample1")
        b = mapper.resolve("AmaranthusPumilus-Sample1")
        assert a.species == b.species == "Amaranthus_pumilus"

    def test_longest_match_wins(self):
        mapper = NameMapper(vocabulary=["A_hybridus", "A_hybridus_quitensis"])
        assert mapper.resolve("A_hybridus_quitensis_s1").species == "A_hybridus_quitensis"
        assert mapper.resolve("A_hybridus_s1").species == "A_hybridus"

    def test_abbreviation_is_inferred_and_not_accepted_by_default(self, mapper):
        hit = mapper.resolve("A_hybridus_s1")
        assert (hit.species, hit.method, hit.status, hit.leftover) == ("Amaranthus_hybridus", "abbrev", "inferred", ("s1",))
        assert not hit.resolved

    def test_abbreviation_can_be_accepted(self):
        mapper = NameMapper(vocabulary=SPECIES, accept_inferred=True)
        hit = mapper.resolve("A_hybridus_s1")
        assert (hit.species, hit.method, hit.status) == ("Amaranthus_hybridus", "abbrev", "ok")

    def test_abbreviation_is_at_most_three_characters(self, mapper):
        assert mapper.resolve("Ama_hybridus").status == "inferred"
        assert mapper.resolve("Amar_hybridus").status == "unmapped"

    def test_homonym_makes_an_abbreviation_ambiguous(self, mapper):
        hit = mapper.resolve("A_pumilus_sample1")
        assert hit.status == "ambiguous" and hit.species is None
        assert {c[0] for c in hit.candidates} == {"Amaranthus_pumilus", "Atriplex_pumilus"}
        # a full name still resolves, even though the abbreviation does not
        assert mapper.resolve("AmaranthusPumilus").status == "ok"

    def test_unmapped_lists_three_close_entries_but_never_assigns(self, mapper):
        hit = mapper.resolve("Beta_vulgarus")
        assert hit.status == "unmapped" and hit.species is None
        assert len(hit.candidates) == 3
        assert hit.candidates[0][0] == "Beta_vulgaris"

    def test_tokens_must_be_contiguous(self, mapper):
        assert mapper.resolve("Amaranthus_x_pumilus").status == "unmapped"

    def test_without_normalization_case_matters_but_words_still_match(self):
        # the page used to say that only exact names match with -N none
        mapper = NameMapper(vocabulary=SPECIES, normalize=["none"])
        assert mapper.resolve("Beta_vulgaris").status == "ok"
        assert mapper.resolve("beta_vulgaris").status == "unmapped"
        assert mapper.resolve("BetaVulgaris").status == "unmapped"
        hits = [mapper.resolve(n) for n in ("Beta_vulgaris_s1", "S012-Beta-vulgaris", "B_vulgaris")]
        assert [(h.species, h.method, h.status) for h in hits] == [
            ("Beta_vulgaris", "token", "ok"), ("Beta_vulgaris", "token", "ok"), ("Beta_vulgaris", "abbrev", "inferred"),
        ]

    def test_steps_given_replace_the_default_steps(self):
        # -N digits alone switches off case folding and the other default steps, as the help says
        mapper = NameMapper(vocabulary=SPECIES, normalize=["digits"])
        assert mapper.steps == ("digits",)
        assert mapper.resolve("beta_vulgaris").status == "unmapped"
        mapper = NameMapper(vocabulary=SPECIES, normalize=["accents", "camel", "digits", "sep", "case"])
        assert mapper.resolve("beta_vulgaris").status == "ok" and mapper.resolve("Beta_vulgaris2").leftover == ("2",)

    def test_species_with_the_same_canonical_form_are_ambiguous_when_not_exact(self):
        vocab = Vocabulary(["A_b", "A-b"])
        assert vocab.match("A_b")[:2] == ("exact", ["A_b"])
        method, found, _ = vocab.match("a b")
        assert method == "canon" and sorted(found) == ["A-b", "A_b"]

    def test_results_are_memoized_in_order(self, mapper):
        first = mapper.resolve("Beta_vulgaris")
        assert mapper.resolve("Beta_vulgaris") is first
        mapper.resolve("Amaranthus_pumilus")
        assert [h.name for h in mapper.hits()] == ["Beta_vulgaris", "Amaranthus_pumilus"]
        assert set(mapper.resolve_all(["Beta_vulgaris", "Zea_mays"])) == {"Beta_vulgaris", "Zea_mays"}

    def test_empty_vocabulary_is_an_error(self):
        with pytest.raises(TwigError, match="vocabulary is empty"):
            NameMapper(vocabulary=[])

    def test_names_without_any_letters_or_digits_are_unmapped(self, mapper):
        hit = mapper.resolve("---")
        assert hit.status == "unmapped" and hit.species is None

    def test_vocabulary_membership(self):
        vocab = Vocabulary(SPECIES)
        assert "Beta_vulgaris" in vocab and "beta_vulgaris" not in vocab


###########################################################################
# precedence and parsing rules
###########################################################################

class TestPrecedence:
    def test_exact_map_row_beats_everything(self, tmp_path):
        path = tmp_path / "map.tsv"
        path.write_text("Beta_vulgaris_x1 Amaranthus_pumilus s9\n")
        mapper = NameMapper(map_rows=parse_map_file(path), vocabulary=SPECIES)
        hit = mapper.resolve("Beta_vulgaris_x1")
        assert (hit.species, hit.specimen, hit.method) == ("Amaranthus_pumilus", "s9", "map")

    def test_regex_row_beats_the_vocabulary_but_not_exact_rows(self, tmp_path):
        path = tmp_path / "map.tsv"
        path.write_text("Ahyb_special Beta_vulgaris\n~^Ahyb[_-] Amaranthus_hybridus\n")
        mapper = NameMapper(map_rows=parse_map_file(path), vocabulary=SPECIES)
        assert mapper.resolve("Ahyb_special").species == "Beta_vulgaris"
        hit = mapper.resolve("Ahyb-1")
        assert (hit.species, hit.method) == ("Amaranthus_hybridus", "map-regex")

    def test_exact_rows_are_matched_by_parsed_key_too(self, tmp_path):
        path = tmp_path / "map.tsv"
        path.write_text("sppA pop1\n")
        mapper = NameMapper(map_rows=parse_map_file(path), delim="-", delim_idxs=[1])
        assert mapper.resolve("OG1-sppA-g1").species == "pop1"

    def test_map_species_must_be_in_the_vocabulary(self, tmp_path):
        path = tmp_path / "map.tsv"
        path.write_text("x1 Zea_mays\n")
        with pytest.raises(TwigError, match=r"map.tsv:1: species 'Zea_mays' is not in the species tree"):
            NameMapper(map_rows=parse_map_file(path), vocabulary=SPECIES)

    def test_map_species_are_spelled_as_in_the_vocabulary(self, tmp_path):
        path = tmp_path / "map.tsv"
        path.write_text("x1 beta-vulgaris\n")
        mapper = NameMapper(map_rows=parse_map_file(path), vocabulary=SPECIES)
        assert mapper.resolve("x1").species == "Beta_vulgaris"

    def test_delimiter_rule_then_vocabulary(self):
        mapper = NameMapper(vocabulary=SPECIES, delim="-", delim_idxs=[1])
        hit = mapper.resolve("OG0001-AmaranthusPumilus-s1")
        assert (hit.species, hit.method, hit.leftover) == ("Amaranthus_pumilus", "canon", ("og0001", "s1"))

    def test_delimiter_rule_that_cannot_split_falls_through(self):
        mapper = NameMapper(vocabulary=SPECIES, delim="|", delim_idxs=[1])
        assert mapper.resolve("Beta_vulgaris").species == "Beta_vulgaris"

    def test_a_name_without_the_delimiter_falls_through(self):
        # with -di 0 the whole name used to become the key, so it looked resolved by the delimiter
        hit = NameMapper(delim="_").resolve("Apum2")
        assert (hit.species, hit.method) == ("Apum2", "self")
        assert NameMapper(delim="_", delim_idxs=[0, -1]).resolve("Apum2").method == "self"
        assert NameMapper(delim="_").resolve("Apum_2").method == "delim"

    def test_delimiter_join_defaults_to_the_delimiter(self):
        mapper = NameMapper(delim="_", delim_idxs=[0, 1])
        assert mapper.resolve("Beta_vulgaris_s1").species == "Beta_vulgaris"
        mapper = NameMapper(delim="_", delim_idxs=[0, 1], delim_join="-")
        hit = mapper.resolve("Beta_vulgaris_s1")
        assert (hit.species, hit.method) == ("Beta-vulgaris", "delim")

    def test_regex_rule_with_specimen_group(self):
        mapper = NameMapper(vocabulary=SPECIES, regex=r"^(?P<species>[A-Za-z]+_[a-z]+)_(?P<specimen>s\d+)$")
        hit = mapper.resolve("Beta_vulgaris_s12")
        assert (hit.species, hit.specimen, hit.specimen_key) == ("Beta_vulgaris", "s12", "s12")
        # a name that the regex does not match falls through to the vocabulary
        assert mapper.resolve("Beta_vulgaris").species == "Beta_vulgaris"

    def test_specimen_key_defaults_to_the_leftover(self):
        hit = NameMapper(vocabulary=SPECIES).resolve("AmaranthusPumilus-Sample1")
        assert hit.specimen is None and hit.specimen_key == "sample1"

    def test_without_a_vocabulary_the_key_or_name_is_the_species(self):
        assert NameMapper().resolve("Beta_vulgaris_s1").species == "Beta_vulgaris_s1"
        assert NameMapper().resolve("Beta_vulgaris_s1").method == "self"
        assert NameMapper(delim="_", delim_idxs=[0, 1]).resolve("Beta_vulgaris_s1").species == "Beta_vulgaris"

    def test_delimiter_and_regex_cannot_be_combined(self):
        with pytest.raises(TwigError, match="either a delimiter"):
            NameMapper(delim="_", regex="(?P<species>.+)")

    def test_regex_needs_a_species_group(self):
        with pytest.raises(TwigError, match="named group"):
            NameMapper(regex=r"\w+")
        with pytest.raises(TwigError, match="invalid name regex"):
            NameMapper(regex="(")

    def test_from_args(self, tmp_path):
        path = tmp_path / "map.tsv"
        path.write_text("x1 Beta_vulgaris\n")
        args = SimpleNamespace(imap=path, name_regex=None, delim=None, delim_idxs=[0], delim_join=None,
                               normalize=None, accept_inferred=True)
        mapper = NameMapper.from_args(args, vocabulary=SPECIES)
        assert mapper.resolve("x1").species == "Beta_vulgaris"
        assert mapper.resolve("A_hybridus_s1").status == "ok"


###########################################################################
# writers
###########################################################################

class TestWriters:
    @pytest.fixture
    def hits(self):
        mapper = NameMapper(vocabulary=SPECIES)
        names = ["Beta_vulgaris", "AmaranthusPumilus-Sample1", "A_hybridus_s1", "A_pumilus_s1", "Zea_mays"]
        return [mapper.resolve(n) for n in names]

    def test_table_columns_and_empty_cells(self, hits):
        out = io.StringIO()
        write_hit_table(hits, out, counts={"Beta_vulgaris": 3})
        lines = out.getvalue().splitlines()
        assert lines[0] == "name\tspecies\tspecimen\tmethod\tstatus\tleftover\tcandidates\tn_files"
        first = lines[1].split("\t")
        assert first == ["Beta_vulgaris", "Beta_vulgaris", "-", "exact", "ok", "-", "-", "3"]
        assert lines[2].split("\t")[5] == "sample1"
        assert "Amaranthus_pumilus(" in lines[4].split("\t")[6]    # candidates of the ambiguous name

    def test_frozen_map_round_trips(self, hits, tmp_path):
        path = tmp_path / "frozen.tsv"
        write_frozen_map(hits, path)
        text = path.read_text()
        assert "# AMBIGUOUS A_pumilus_s1" in text and "# UNMAPPED Zea_mays" in text
        assert "A_hybridus_s1 Amaranthus_hybridus  # inferred" in text
        rows = parse_map_file(path)
        assert {r.sample: r.species for r in rows} == {
            "Beta_vulgaris": "Beta_vulgaris",
            "AmaranthusPumilus-Sample1": "Amaranthus_pumilus",
            "A_hybridus_s1": "Amaranthus_hybridus",
        }
        # a reloaded map resolves the inferred name without any abbreviation rule
        assert NameMapper(map_rows=rows, vocabulary=SPECIES).resolve("A_hybridus_s1").status == "ok"

    def test_specimen_keys_survive_a_frozen_map(self, tmp_path):
        # map rows used to give no leftover, so every sample of a species had the key '' after a reload
        mapper = NameMapper(vocabulary=SPECIES, accept_inferred=True)
        names = ["AmaranthusPumilus-Sample1", "AmaranthusPumilus-Sample2", "A_hybridus_s1", "s012_Beta_vulgaris", "Beta_vulgaris"]
        before = [mapper.resolve(n) for n in names]
        assert [h.specimen_key for h in before] == ["sample1", "sample2", "s1", "s012", ""]
        path = tmp_path / "frozen.tsv"
        write_frozen_map(before, path)
        reloaded = NameMapper(map_rows=parse_map_file(path), vocabulary=SPECIES)
        after = [reloaded.resolve(n) for n in names]
        assert {h.method for h in after} == {"map"}
        assert [h.specimen_key for h in after] == [h.specimen_key for h in before]
        assert [h.leftover for h in after] == [h.leftover for h in before]

    def test_the_leftover_of_a_map_row_whose_species_is_not_in_the_name(self, tmp_path):
        path = tmp_path / "map.tsv"
        path.write_text("x1 Beta_vulgaris\n~^Ahyb[_-] Amaranthus_hybridus\nA_pumilus_s3 Amaranthus_pumilus\n")
        mapper = NameMapper(map_rows=parse_map_file(path), vocabulary=SPECIES)
        assert mapper.resolve("x1").leftover == ("x1",)
        assert mapper.resolve("Ahyb-7").leftover == ("ahyb", "7")
        # the abbreviated genus is part of the species, as it is for a match by abbreviation
        assert mapper.resolve("A_pumilus_s3").leftover == ("s3",)

    def test_frozen_map_has_no_specimen_column_unless_explicit(self, tmp_path):
        mapper = NameMapper(vocabulary=SPECIES, regex=r"^(?P<species>[A-Za-z]+_[a-z]+)_(?P<specimen>s\d+)$")
        hits = [mapper.resolve("Beta_vulgaris_s1"), mapper.resolve("Amaranthus_pumilus_x1")]
        path = tmp_path / "frozen.tsv"
        write_frozen_map(hits, path)
        lines = [i for i in path.read_text().splitlines() if not i.startswith("#")]
        assert lines == ["Beta_vulgaris_s1 Beta_vulgaris s1", "Amaranthus_pumilus_x1 Amaranthus_pumilus"]

    def test_frozen_map_refuses_to_overwrite(self, hits, tmp_path):
        path = tmp_path / "frozen.tsv"
        write_frozen_map(hits, path)
        with pytest.raises(TwigError, match="use --force"):
            write_frozen_map(hits, path)
        write_frozen_map(hits, path, force=True)
