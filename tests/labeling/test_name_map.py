"""Tests for the name-map command (twig.labeling.name_map)."""

import pytest

pytest.importorskip("toytree")

from twig.labeling.name_map import needs_attention, read_units
from twig.labeling.core.names import NameHit
from twig.utils import TwigError


def parse_table(text):
    """Return the rows of a name-map table as dicts."""
    lines = text.rstrip("\n").split("\n")
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


@pytest.fixture
def inputs(names_dir):
    """The three gene trees and two alignments of the test data, and the species tree."""
    return f"-i {names_dir}/trees {names_dir}/aln -s {names_dir}/species.nwk"


class TestTable:
    def test_every_name_is_resolved_with_method_and_status(self, run_twig, inputs):
        code, out, err = run_twig(f"name-map {inputs}")
        assert code == 0
        rows = {r["name"]: r for r in parse_table(out)}
        assert len(rows) == 13
        assert rows["Beta_vulgaris"]["method"] == "exact"
        assert rows["AmaranthusPumilus-Sample1"]["species"] == "Amaranthus_pumilus"
        assert rows["AmaranthusPumilus-Sample1"]["method"] == "token"
        assert rows["AmaranthusPumilus-Sample1"]["leftover"] == "sample1"
        assert rows["s012_Amaranthus_hybridus"]["species"] == "Amaranthus_hybridus"
        assert rows["A_hybridus_s1"]["status"] == "inferred" and rows["A_hybridus_s1"]["method"] == "abbrev"
        assert rows["A_pumilus_sample1"]["status"] == "ambiguous" and rows["A_pumilus_sample1"]["species"] == "-"
        assert "Atriplex_pumilus(" in rows["A_pumilus_sample1"]["candidates"]
        assert rows["Zea_mays_s1"]["status"] == "unmapped"
        assert len(rows["Zea_mays_s1"]["candidates"].split(";")) == 3

    def test_number_of_files_in_which_a_name_occurs(self, run_twig, inputs):
        rows = {r["name"]: r for r in parse_table(run_twig(f"name-map {inputs}")[1])}
        assert rows["Spinacia_oleracea"]["n_files"] == "5"      # three trees and two alignments
        assert rows["Amaranthus_pumilus"]["n_files"] == "2"
        assert rows["Zea_mays_s1"]["n_files"] == "1"

    def test_show_problems_prints_only_names_that_need_attention(self, run_twig, inputs):
        rows = parse_table(run_twig(f"name-map {inputs} --show problems")[1])
        assert [r["name"] for r in rows] == ["Beta_vulgaris_maritima_s1", "A_pumilus_sample1", "A_hybridus_s1", "Zea_mays_s1"]

    def test_summary_goes_to_stderr(self, run_twig, inputs):
        _, out, err = run_twig(f"name-map {inputs}")
        assert "CMD:" in err and "INFO" in err
        assert "resolved 13 unique names from 5 files (3 trees, 2 FASTA files)" in err
        assert "status ambiguous 1, inferred 1, ok 10, unmapped 1" in err
        assert "1 name is ambiguous, matching more than one species" in err
        assert "1 name is unmapped, matching no species" in err
        assert "1 name is inferred, matched only by genus abbreviation" in err
        assert "1 word is left after the species name in the names of Beta_vulgaris. maritima" in err
        assert "Amaranthus_pumilus has up to 2 samples in one tree or FASTA file" in err
        assert "the species vocabulary has 5 species" in err
        assert "\tstatus\t" not in err      # the table itself is only written to stdout

    def test_species_without_samples_are_reported(self, run_twig, names_dir):
        sp = f"-s {names_dir}/species.nwk"
        # in OG0003 only the ambiguous A_pumilus_sample1 could have been Amaranthus_pumilus
        _, _, err = run_twig(f"name-map -i {names_dir}/trees/OG0003.treefile {sp}")
        assert "1 species has no samples. Amaranthus_pumilus" in err
        _, _, err = run_twig(f"name-map -i {names_dir}/trees/OG0001.treefile {sp}")
        assert "no samples" not in err

    def test_counts_in_the_summary_have_the_right_number(self, run_twig, names_dir, tmp_path):
        sp = tmp_path / "sp.nwk"
        sp.write_text("(Zea_mays,(Oryza_sativa,(Amaranthus_pumilus,Beta_vulgaris)));\n")
        tree = tmp_path / "OG1.nwk"
        tree.write_text("(Amaranthus_pumilus,Beta_vulgaris_maritima_s1,Beta_vulgaris_x_s2);\n")
        code, _, err = run_twig(f"name-map -i {tree} -s {sp} --strict")
        assert "2 species have no samples. Zea_mays, Oryza_sativa" in err and "speciess" not in err
        assert "resolved 3 unique names from 1 file (1 tree, 0 FASTA files)" in err
        assert "2 words are left after the species name in the names of Beta_vulgaris. maritima, x" in err
        assert code == 1 and "2 names need attention" in err
        tree.write_text("(Amaranthus_pumilus,Beta_vulgaris_maritima_s1);\n")
        code, _, err = run_twig(f"name-map -i {tree} -s {sp} --strict")
        assert code == 1 and "1 name needs attention" in err


class TestOutputs:
    def test_table_to_a_file(self, run_twig, inputs, tmp_path):
        out_path = tmp_path / "table.tsv"
        code, out, err = run_twig(f"name-map {inputs} -o {out_path}")
        assert code == 0 and out == ""
        assert len(parse_table(out_path.read_text())) == 13
        assert f"wrote the table of 13 names to {out_path}" in err

    def test_existing_output_needs_force(self, run_twig, inputs, tmp_path):
        out_path = tmp_path / "table.tsv"
        out_path.write_text("old")
        code, _, err = run_twig(f"name-map {inputs} -o {out_path}")
        assert code == 1 and "use --force" in err and out_path.read_text() == "old"
        code, _, _ = run_twig(f"name-map {inputs} -o {out_path} -f")
        assert code == 0 and out_path.read_text().startswith("name\tspecies")

    def test_an_output_that_is_a_directory_is_a_clear_error(self, run_twig, inputs, tmp_path):
        folder = tmp_path / "names"
        folder.mkdir()
        for option in ("-o", "--freeze"):
            code, _, err = run_twig(f"name-map {inputs} {option} {folder} -f")
            assert code == 1 and "is a directory and not a regular file" in err and "Traceback" not in err
        assert folder.is_dir() and list(folder.iterdir()) == []

    def test_an_output_that_is_an_input_is_refused_even_with_force(self, run_twig, names_dir, tmp_path):
        tree = tmp_path / "OG1.treefile"
        tree.write_text("(Beta_vulgaris,Spinacia_oleracea);\n")
        species = tmp_path / "species.nwk"
        species.write_text("(Beta_vulgaris,Spinacia_oleracea);\n")
        mapping = tmp_path / "map.tsv"
        mapping.write_text("Beta_vulgaris\n")
        for option, path in (("-o", tree), ("--freeze", species), ("--freeze", mapping), ("-o", mapping)):
            before = path.read_text()
            code, _, err = run_twig(f"name-map -i {tmp_path} -s {species} -I {mapping} {option} {path} -f")
            assert code == 1 and "is also an input file, so it is not overwritten" in err, (option, path)
            assert path.read_text() == before

    def test_out_and_freeze_must_differ(self, run_twig, inputs, tmp_path):
        path = tmp_path / "same.tsv"
        code, _, err = run_twig(f"name-map {inputs} -o {path} --freeze {path}")
        assert code == 1 and "must be different files" in err

    def test_freeze_writes_a_map_that_other_runs_can_use(self, run_twig, inputs, tmp_path, names_dir):
        frozen = tmp_path / "names.tsv"
        code, _, err = run_twig(f"name-map {inputs} --freeze {frozen}")
        assert code == 0 and "use it with -I" in err
        text = frozen.read_text()
        assert "A_hybridus_s1 Amaranthus_hybridus  # inferred" in text
        assert "# AMBIGUOUS A_pumilus_sample1" in text and "# UNMAPPED Zea_mays_s1" in text
        # with the frozen map the inferred name resolves through the map, and only unresolved names remain
        code, out, _ = run_twig(f"name-map {inputs} -I {frozen} --show problems")
        assert code == 0
        assert [r["name"] for r in parse_table(out)] == ["A_pumilus_sample1", "Zea_mays_s1"]
        # a frozen map is not overwritten without --force
        code, _, err = run_twig(f"name-map {inputs} --freeze {frozen}")
        assert code == 1 and "use --force" in err


    def test_a_closed_pipe_gives_no_traceback(self, names_dir):
        import subprocess, sys
        # the reader exits at once, so the child gets a broken pipe when it prints its table
        cmd = [sys.executable, "-c", "import sys; sys.argv[:] = ['twig', 'name-map', '-i', sys.argv[1], '-s', sys.argv[2]]; from twig.cli.main import main; sys.exit(main())",
               str(names_dir / "trees"), str(names_dir / "species.nwk")]
        child = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        child.stdout.close()
        _, err = child.communicate()
        assert "Traceback" not in err and "BrokenPipeError" not in err


class TestStrict:
    def test_strict_fails_but_still_writes_the_table_and_no_map(self, run_twig, inputs, tmp_path):
        table, frozen = tmp_path / "table.tsv", tmp_path / "names.tsv"
        code, _, err = run_twig(f"name-map {inputs} --strict -o {table} --freeze {frozen}")
        assert code == 1 and "4 names need attention" in err
        assert table.exists() and not frozen.exists()

    def test_strict_passes_when_everything_is_resolved(self, run_twig, names_dir):
        code, out, _ = run_twig(f"name-map -i {names_dir}/trees/OG0001.treefile {names_dir}/aln/OG0001.nt.fa -s {names_dir}/species.nwk --strict")
        assert code == 0 and len(parse_table(out)) == 5

    def test_inferred_matches_pass_strict_when_accepted(self, run_twig, names_dir, tmp_path):
        tree = tmp_path / "g.nwk"
        tree.write_text("(A_hybridus_s1,Beta_vulgaris);\n")
        sp = f"-s {names_dir}/species.nwk"
        assert run_twig(f"name-map -i {tree} {sp} --strict")[0] == 1
        assert run_twig(f"name-map -i {tree} {sp} --strict --accept-inferred")[0] == 0


class TestMappingOptions:
    def test_species_list_instead_of_a_species_tree(self, run_twig, names_dir, tmp_path):
        listing = tmp_path / "species.txt"
        listing.write_text("# species\nBeta_vulgaris\nSpinacia_oleracea   some note\n")
        code, out, err = run_twig(f"name-map -i {names_dir}/trees/OG0001.treefile -V {listing}")
        assert code == 0
        rows = {r["name"]: r["status"] for r in parse_table(out)}
        assert rows["Beta_vulgaris"] == "ok" and rows["Amaranthus_pumilus"] == "unmapped"
        assert "the species vocabulary has 2 species" in err

    def test_a_byte_order_mark_is_not_part_of_the_first_species(self, run_twig, names_dir, tmp_path):
        listing = tmp_path / "species.txt"
        listing.write_bytes("\ufeffBeta_vulgaris\nSpinacia_oleracea\n".encode("utf-8"))
        alignment = tmp_path / "OG1.fa"
        alignment.write_bytes("\ufeff>Beta_vulgaris\nACGT\n>Spinacia_oleracea\nACGT\n".encode("utf-8"))
        code, out, err = run_twig(f"name-map -i {alignment} -V {listing} --strict")
        assert code == 0, err
        rows = {r["name"]: (r["species"], r["method"]) for r in parse_table(out)}
        assert rows == {"Beta_vulgaris": ("Beta_vulgaris", "exact"), "Spinacia_oleracea": ("Spinacia_oleracea", "exact")}
        assert "\ufeff" not in out + err

    def test_a_species_list_that_is_not_text_is_a_clear_error(self, run_twig, names_dir, tmp_path):
        listing = tmp_path / "species.txt"
        listing.write_bytes(b"Amaranthus_pumilus \xe9\n")
        code, _, err = run_twig(f"name-map -i {names_dir}/trees/OG0001.treefile -V {listing}")
        assert code == 1 and "species.txt is not valid UTF-8 text" in err and "Traceback" not in err

    def test_species_tree_and_list_are_combined(self, run_twig, names_dir, tmp_path):
        listing = tmp_path / "species.txt"
        listing.write_text("Zea_mays\n")
        code, out, err = run_twig(f"name-map -i {names_dir}/trees/OG0003.treefile -s {names_dir}/species.nwk -V {listing}")
        rows = {r["name"]: r["status"] for r in parse_table(out)}
        assert rows["Zea_mays_s1"] == "ok" and "the species vocabulary has 6 species" in err

    def test_without_a_vocabulary_names_are_their_own_species_unless_parsed(self, run_twig, names_dir):
        tree = f"{names_dir}/trees/OG0002.treefile"
        code, out, err = run_twig(f"name-map -i {tree}")
        assert code == 0 and "names are not matched to species" in err
        assert {r["name"]: r["species"] for r in parse_table(out)}["Beta_vulgaris_maritima_s1"] == "Beta_vulgaris_maritima_s1"
        code, out, _ = run_twig(f"name-map -i {tree} -d _ -di 0 1")
        rows = {r["name"]: r for r in parse_table(out)}
        assert rows["Beta_vulgaris_maritima_s1"]["species"] == "Beta_vulgaris" and rows["Beta_vulgaris_maritima_s1"]["method"] == "delim"

    def test_delimiter_rule_with_a_vocabulary(self, run_twig, names_dir, tmp_path):
        tree = tmp_path / "g.nwk"
        tree.write_text("(OG1-AmaranthusPumilus-s1,OG1-BetaVulgaris-s1);\n")
        code, out, _ = run_twig(f"name-map -i {tree} -s {names_dir}/species.nwk -d - -di 1")
        rows = {r["name"]: r for r in parse_table(out)}
        assert rows["OG1-AmaranthusPumilus-s1"]["species"] == "Amaranthus_pumilus"
        assert rows["OG1-BetaVulgaris-s1"]["species"] == "Beta_vulgaris"

    def test_regex_rule_with_specimen(self, run_twig, names_dir, tmp_path):
        tree = tmp_path / "g.nwk"
        tree.write_text("(Beta_vulgaris_s1,Beta_vulgaris_s2);\n")
        code, out, _ = run_twig(f"name-map -i {tree} -s {names_dir}/species.nwk -x ^(?P<species>[A-Za-z]+_[a-z]+)_(?P<specimen>s[0-9]+)$")
        rows = parse_table(out)
        assert [(r["species"], r["specimen"]) for r in rows] == [("Beta_vulgaris", "s1"), ("Beta_vulgaris", "s2")]

    def test_map_file_with_regex_row(self, run_twig, names_dir, tmp_path):
        mapping = tmp_path / "map.tsv"
        mapping.write_text("~^Zea Beta_vulgaris\n")
        code, out, _ = run_twig(f"name-map -i {names_dir}/trees/OG0003.treefile -s {names_dir}/species.nwk -I {mapping}")
        rows = {r["name"]: r for r in parse_table(out)}
        assert rows["Zea_mays_s1"]["species"] == "Beta_vulgaris" and rows["Zea_mays_s1"]["method"] == "map-regex"

    def test_map_file_species_must_exist(self, run_twig, names_dir, tmp_path):
        mapping = tmp_path / "map.tsv"
        mapping.write_text("x Quercus_robur\n")
        code, _, err = run_twig(f"name-map -i {names_dir}/trees/OG0003.treefile -s {names_dir}/species.nwk -I {mapping}")
        assert code == 1 and "is not in the species tree" in err

    def test_normalization_can_be_switched_off(self, run_twig, names_dir, tmp_path):
        tree = tmp_path / "g.nwk"
        tree.write_text("(beta_vulgaris,Beta_vulgaris);\n")
        sp = f"-s {names_dir}/species.nwk"
        rows = {r["name"]: r["status"] for r in parse_table(run_twig(f"name-map -i {tree} {sp}")[1])}
        assert rows == {"beta_vulgaris": "ok", "Beta_vulgaris": "ok"}
        rows = {r["name"]: r["status"] for r in parse_table(run_twig(f"name-map -i {tree} {sp} -N none")[1])}
        assert rows == {"beta_vulgaris": "unmapped", "Beta_vulgaris": "ok"}


class TestInputs:
    def test_glob_directory_and_files_can_be_mixed(self, run_twig, names_dir):
        code, out, err = run_twig(f"name-map -i {names_dir}/trees/*.treefile {names_dir}/aln {names_dir}/trees/OG0001.treefile -s {names_dir}/species.nwk")
        assert code == 0 and "from 5 files" in err

    def test_multi_tree_files_count_as_one_file(self, run_twig, names_dir, tmp_path):
        trees = tmp_path / "many.nwk"
        trees.write_text("(Beta_vulgaris,Spinacia_oleracea);\n(Beta_vulgaris,Atriplex_pumilus);\n")
        code, out, err = run_twig(f"name-map -i {trees} -s {names_dir}/species.nwk")
        rows = {r["name"]: r["n_files"] for r in parse_table(out)}
        assert rows["Beta_vulgaris"] == "1" and "from 1 file (2 trees, 0 FASTA files)" in err

    def test_unknown_files_are_skipped_in_directories_but_not_when_named(self, run_twig, names_dir, tmp_path):
        folder = tmp_path / "mixed"
        folder.mkdir()
        (folder / "OG1.nwk").write_text("(Beta_vulgaris,Spinacia_oleracea);\n")
        (folder / "notes.txt").write_text("this is not a tree\n")
        code, out, err = run_twig(f"name-map -i {folder} -s {names_dir}/species.nwk")
        assert code == 0 and "skipping" in err and "notes.txt" in err
        code, _, err = run_twig(f"name-map -i {folder}/notes.txt -s {names_dir}/species.nwk")
        assert code == 1 and "neither a tree file nor a FASTA file" in err

    def test_a_phylip_alignment_in_a_directory_is_skipped(self, run_twig, names_dir, tmp_path):
        # its first line of two numbers looks like the header of a PAML tree file
        folder = tmp_path / "aln"
        folder.mkdir()
        (folder / "OG1.phy").write_text("2 8\nAmaranthus_pumilus ACGTACGT\nBeta_vulgaris ACGTACGT\n")
        (folder / "OG2.fa").write_text(">Beta_vulgaris\nACGT\n")
        code, out, err = run_twig(f"name-map -i {folder} -s {names_dir}/species.nwk")
        assert code == 0, err
        assert "skipping" in err and "OG1.phy, which is neither a tree file nor a FASTA file" in err
        assert [r["name"] for r in parse_table(out)] == ["Beta_vulgaris"]

    def test_errors(self, run_twig, names_dir, tmp_path):
        sp = f"-s {names_dir}/species.nwk"
        code, _, err = run_twig(f"name-map -i {tmp_path}/missing.nwk {sp}")
        assert code == 1 and "input file not found" in err
        bad = tmp_path / "bad.nwk"
        bad.write_text("((A,B),C;\n")
        code, _, err = run_twig(f"name-map -i {bad} {sp}")
        assert code == 1 and "bad.nwk, tree 1" in err
        code, _, err = run_twig(f"name-map -i {names_dir}/trees {tmp_path}/missing.nwk -s {tmp_path}/nospecies.nwk")
        assert code == 1

    def test_duplicate_species_tips_are_rejected(self, run_twig, names_dir, tmp_path):
        sp = tmp_path / "dup.nwk"
        sp.write_text("((A,B),(A,C));\n")
        code, _, err = run_twig(f"name-map -i {names_dir}/trees -s {sp}")
        assert code == 1 and "tip names are not unique" in err


class TestSummaryLimits:
    def test_long_lists_are_cut(self, run_twig, tmp_path):
        species = [f"Genus{i}_species" for i in range(10)]
        sp = tmp_path / "sp.nwk"
        sp.write_text("(" + ",".join(species) + ");\n")
        tree = tmp_path / "g.nwk"
        tree.write_text("(" + ",".join(f"{s},{s}_s2" for s in species) + ");\n")
        code, out, err = run_twig(f"name-map -i {tree} -s {sp}")
        assert code == 0
        assert err.count("samples in one tree or FASTA file") == 8
        assert "2 more species have several samples in a tree or FASTA file" in err


class TestHelpers:
    def test_needs_attention(self):
        ok = NameHit("a", "A", None, "exact", "ok")
        assert not needs_attention(ok)
        assert needs_attention(NameHit("a", "A", None, "abbrev", "inferred"))
        assert needs_attention(NameHit("a", "A", None, "token", "ok", leftover=("maritima", "s1")))
        assert not needs_attention(NameHit("a", "A", None, "token", "ok", leftover=("s1", "g7")))

    def test_read_units_requires_names(self, tmp_path):
        empty = tmp_path / "empty"
        empty.mkdir()
        (empty / "notes.txt").write_text("nothing\n")
        with pytest.raises(TwigError, match="no names were found"):
            read_units([empty])
