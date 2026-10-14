"""Tests for input files, loci and safe output (twig.labeling.core.batch)."""

import json
import os
from pathlib import Path

import pytest

from twig.utils import TwigError
from twig.labeling.core import batch


@pytest.fixture
def files(tmp_path):
    """A directory with trees, alignments, a hidden file and a subdirectory."""
    for name in ("OG0002.treefile", "OG0001.treefile", "OG0003.treefile"):
        (tmp_path / "trees").mkdir(exist_ok=True)
        (tmp_path / "trees" / name).write_text("(A,B);\n")
    for name in ("OG0001.nt.fa", "OG0002.nt.fa", "OG9999.nt.fa"):
        (tmp_path / "aln").mkdir(exist_ok=True)
        (tmp_path / "aln" / name).write_text(">A\nAC\n")
    (tmp_path / "trees" / ".hidden").write_text("x")
    (tmp_path / "trees" / "subdir").mkdir()
    (tmp_path / "trees" / "OG0004.treefile.tmp").write_text("x")
    return tmp_path


class TestExpand:
    def test_glob_is_sorted(self, files):
        found = batch.expand_inputs([str(files / "trees" / "*.treefile")])
        assert [p.name for p in found] == ["OG0001.treefile", "OG0002.treefile", "OG0003.treefile"]

    def test_directory_gives_regular_visible_files(self, files):
        found = batch.expand_inputs([files / "trees"])
        assert [p.name for p in found] == ["OG0001.treefile", "OG0002.treefile", "OG0003.treefile"]

    def test_a_glob_leaves_out_hidden_and_temporary_files_like_a_directory(self, files):
        # the .tmp file is a partial output of a run that was killed
        found = batch.expand_inputs([str(files / "trees" / "*")])
        assert [p.name for p in found] == ["OG0001.treefile", "OG0002.treefile", "OG0003.treefile"]
        with pytest.raises(TwigError, match="no files match"):
            batch.expand_inputs([str(files / "trees" / "*.tmp")])

    def test_files_directories_and_globs_can_be_mixed_and_deduplicated(self, files):
        items = [files / "trees" / "OG0003.treefile", str(files / "trees" / "*.treefile"), files / "aln"]
        found = batch.expand_inputs(items)
        assert [p.name for p in found] == [
            "OG0003.treefile", "OG0001.treefile", "OG0002.treefile",
            "OG0001.nt.fa", "OG0002.nt.fa", "OG9999.nt.fa",
        ]

    def test_errors(self, files):
        with pytest.raises(TwigError, match="no files match"):
            batch.expand_inputs([str(files / "trees" / "*.nope")])
        with pytest.raises(TwigError, match="input file not found"):
            batch.expand_inputs([files / "missing.nwk"])
        (files / "empty").mkdir()
        with pytest.raises(TwigError, match="contains no files"):
            batch.expand_inputs([files / "empty"])

    def test_many_files_do_not_need_a_shell(self, tmp_path):
        for i in range(2500):
            (tmp_path / f"OG{i:05d}.treefile").write_text("(A,B);")
        assert len(batch.expand_inputs([str(tmp_path / "*.treefile")])) == 2500


class TestSniff:
    @pytest.mark.parametrize("text, kind", [
        (">a\nACGT\n", "fasta"),
        ("# comment\n\n>a\nACGT\n", "fasta"),
        ("(A:1,B:2);\n", "tree"),
        ("[&R] (A,B);\n", "tree"),
        ("OG0001\t(A:1,B:2);\n", "tree"),
        ("3 2\n(A,B,C);\n(A,C,B);\n", "tree"),
        ("3 2\n\n[&R] (A,B,C);\n", "tree"),
        ("[&R]\n(A,B);\n", "tree"),
        ("[a comment] [&R]\n\nOG1\t(A,B);\n", "tree"),
        ("2 8\nAmaranthus_pumilus ACGTACGT\nBeta_vulgaris ACGTACGT\n", "unknown"),
        ("  2   8\nA  ACGTACGT\nB  ACGTACGT\n", "unknown"),
        ("3 2\n", "unknown"),
        ("#NEXUS\nbegin trees;\n", "tree"),
        ("name\tvalue\n", "unknown"),
        ("", "unknown"),
    ])
    def test_kinds(self, tmp_path, text, kind):
        path = tmp_path / "x"
        path.write_text(text)
        assert batch.sniff_kind(path) == kind

    def test_a_byte_order_mark_and_a_comment_after_the_locus_are_fine(self, tmp_path):
        # a tree file saved with a byte order mark was neither a tree file nor a FASTA file
        path = tmp_path / "x"
        path.write_bytes(b"\xef\xbb\xbf(A,B);\n")
        assert batch.sniff_kind(path) == "tree"
        path.write_text("OG0001\t[&R] (A,B);\n")
        assert batch.sniff_kind(path) == "tree"


class TestLoci:
    @pytest.mark.parametrize("name, locus", [
        ("OG0001.treefile", "OG0001"),
        ("OG0001.nt.fa", "OG0001"),
        ("OG0001", "OG0001"),
    ])
    def test_default_locus_regex(self, name, locus):
        assert batch.locus_name(name) == locus

    def test_custom_regex_and_errors(self):
        assert batch.locus_name("run1_OG0001.fa", r"(OG\d+)") == "OG0001"
        assert batch.locus_name("OG0001.fa", r"^OG\d+") == "OG0001"      # no group, so the whole match
        with pytest.raises(TwigError, match="does not match"):
            batch.locus_name("x.fa", r"^OG(\d+)")
        with pytest.raises(TwigError, match="invalid locus regex"):
            batch.locus_name("x.fa", "(")
        with pytest.raises(TwigError, match="empty name"):
            batch.locus_name(".hidden", r"^([^.]*)")

    def test_pairing_trees_and_alignments(self, files):
        trees = batch.expand_inputs([files / "trees" / "*.treefile"])
        alns = batch.expand_inputs([files / "aln"])
        pairing = batch.pair_by_locus(trees, alns)
        assert [(p[0], p[2] is not None) for p in pairing.pairs] == [("OG0001", True), ("OG0002", True), ("OG0003", False)]
        assert pairing.unpaired_trees == ("OG0003",)
        assert pairing.orphan_alignments == ("OG9999",)

    def test_trees_only(self, files):
        trees = batch.expand_inputs([files / "trees" / "*.treefile"])
        pairing = batch.pair_by_locus(trees)
        assert pairing.unpaired_trees == () and pairing.orphan_alignments == ()
        assert all(p[2] is None for p in pairing.pairs)

    def test_duplicate_locus_names_are_an_error(self, files):
        (files / "trees" / "OG0001.contree").write_text("(A,B);")
        trees = batch.expand_inputs([files / "trees" / "*"])
        with pytest.raises(TwigError, match="locus name 'OG0001' is shared by tree files"):
            batch.pair_by_locus(trees)


class TestTreeRecords:
    def test_loci_of_single_and_multi_tree_files(self, tmp_path):
        one = tmp_path / "OG0001.treefile"
        one.write_text("(A,B);\n")
        many = tmp_path / "all.trees"
        many.write_text("(A,B);\n(C,D);\n(E,F);\n")
        named = tmp_path / "set.nwk"
        named.write_text("geneX\t(A,B);\ngeneY\t(C,D);\n")
        records = list(batch.tree_records([one, many, named]))
        assert [r[0] for r in records] == ["OG0001", "all_1", "all_2", "all_3", "geneX", "geneY"]
        assert records[0][1] == "(A,B);" and records[2][1] == "(C,D);" and records[4][2] == named

    def test_multi_tree_loci_match_alignment_names_under_the_default_regex(self, tmp_path):
        many = tmp_path / "all.trees"
        many.write_text("(A,B);\n(C,D);\n")
        loci = [r[0] for r in batch.tree_records([many])]
        # the alignment of the second tree would be named all_2.nt.fa
        assert loci == [batch.locus_name("all_1.nt.fa"), batch.locus_name("all_2.nt.fa")]

    def test_errors(self, tmp_path):
        a = tmp_path / "OG1.treefile"
        a.write_text("(A,B);\n")
        b = tmp_path / "OG1.contree"
        b.write_text("(A,B);\n")
        with pytest.raises(TwigError, match="locus name 'OG1' is used by trees in"):
            list(batch.tree_records([a, b]))
        empty = tmp_path / "empty.nwk"
        empty.write_text("\n# nothing\n")
        with pytest.raises(TwigError, match="no trees found"):
            list(batch.tree_records([empty]))


class TestManifest:
    def test_rows_with_paths_relative_to_cwd_or_manifest(self, files, monkeypatch):
        manifest = files / "lists" / "manifest.tsv"
        manifest.parent.mkdir()
        manifest.write_text(
            "# locus tree alignment\n"
            "OG0001 trees/OG0001.treefile aln/OG0001.nt.fa\n"
            "OG0002 ../trees/OG0002.treefile\n"
        )
        monkeypatch.chdir(files)
        rows = batch.read_manifest(manifest)
        assert [r.locus for r in rows] == ["OG0001", "OG0002"]
        assert rows[0].alignment == Path("aln/OG0001.nt.fa") and rows[1].alignment is None
        assert rows[1].tree == manifest.parent / "../trees/OG0002.treefile"

    @pytest.mark.parametrize("text, message", [
        ("OG1\n", r"expected 2 or 3 columns"),
        ("OG1 nope.nwk\n", r"file not found: nope.nwk"),
        ("", r"no manifest rows found"),
    ])
    def test_errors(self, tmp_path, text, message):
        path = tmp_path / "manifest.tsv"
        path.write_text(text)
        with pytest.raises(TwigError, match=message):
            batch.read_manifest(path)

    def test_a_byte_order_mark_is_not_part_of_the_first_locus(self, tmp_path):
        (tmp_path / "t.nwk").write_text("(A,B);")
        path = tmp_path / "manifest.tsv"
        path.write_bytes(f"\ufeffOG1 {tmp_path / 't.nwk'}\n".encode("utf-8"))
        assert [r.locus for r in batch.read_manifest(path)] == ["OG1"]
        path.write_bytes(b"OG\xe9 t.nwk\n")
        with pytest.raises(TwigError, match="not valid UTF-8"):
            batch.read_manifest(path)

    def test_duplicate_locus(self, tmp_path):
        (tmp_path / "t.nwk").write_text("(A,B);")
        path = tmp_path / "manifest.tsv"
        path.write_text("OG1 t.nwk\nOG1 t.nwk\n")
        with pytest.raises(TwigError, match=r"manifest.tsv:2: duplicate locus 'OG1'"):
            batch.read_manifest(path)


class TestAtomicWrite:
    def test_writes_text_and_bytes_and_leaves_no_temporary_file(self, tmp_path):
        path = tmp_path / "sub" / "out.txt"
        batch.atomic_write(path, "hello\n")
        assert path.read_text() == "hello\n"
        batch.atomic_write(tmp_path / "out.bin", b"\x00\x01")
        assert (tmp_path / "out.bin").read_bytes() == b"\x00\x01"
        assert sorted(os.listdir(tmp_path)) == ["out.bin", "sub"]
        assert os.listdir(tmp_path / "sub") == ["out.txt"]

    def test_existing_output_needs_force(self, tmp_path):
        path = tmp_path / "out.txt"
        batch.atomic_write(path, "one")
        with pytest.raises(TwigError, match="use --force to overwrite"):
            batch.atomic_write(path, "two")
        assert path.read_text() == "one"
        batch.atomic_write(path, "two", force=True)
        assert path.read_text() == "two"

    def test_failed_write_keeps_the_old_file_and_removes_the_temporary_file(self, tmp_path, monkeypatch):
        path = tmp_path / "out.txt"
        batch.atomic_write(path, "old")

        def boom(*args, **kwargs):
            raise OSError("disk full")

        monkeypatch.setattr(batch.os, "replace", boom)
        with pytest.raises(TwigError, match=r"cannot write .*out.txt \(disk full\)"):
            batch.atomic_write(path, "new", force=True)
        assert path.read_text() == "old"
        assert os.listdir(tmp_path) == ["out.txt"]

    def test_a_directory_is_never_replaced(self, tmp_path):
        folder = tmp_path / "out"
        folder.mkdir()
        with pytest.raises(TwigError, match="is a directory and not a regular file, so it is never replaced, even with --force"):
            batch.atomic_write(folder, "text", force=True)
        assert folder.is_dir() and os.listdir(tmp_path) == ["out"]

    @pytest.mark.skipif(not hasattr(os, "mkfifo"), reason="needs named pipes")
    def test_a_device_or_a_link_to_one_is_never_replaced(self, tmp_path):
        # with --force a device such as /dev/null used to be renamed over and become a regular file
        pipe = tmp_path / "pipe"
        os.mkfifo(pipe)
        link = tmp_path / "null"
        link.symlink_to(os.devnull)
        for path in (pipe, link):
            with pytest.raises(TwigError, match="is a device, a pipe or a socket and not a regular file"):
                batch.atomic_write(path, "text", force=True)
        assert not pipe.is_file() and link.is_symlink() and os.readlink(link) == os.devnull
        assert sorted(os.listdir(tmp_path)) == ["null", "pipe"]

    def test_a_link_to_a_regular_file_can_be_overwritten(self, tmp_path):
        target = tmp_path / "target.txt"
        target.write_text("old")
        link = tmp_path / "link.txt"
        link.symlink_to(target)
        batch.atomic_write(link, "new", force=True)
        assert link.read_text() == "new"


class TestWriteStdout:
    def test_writes_text(self, capsys):
        batch.write_stdout("hello\n")
        assert capsys.readouterr().out == "hello\n"

    def test_a_closed_pipe_is_not_an_error(self, monkeypatch):
        class Closed:
            def write(self, text):
                raise BrokenPipeError()

            def fileno(self):
                raise ValueError("no file descriptor")

        monkeypatch.setattr(batch.sys, "stdout", Closed())
        batch.write_stdout("hello\n")      # must not raise


class TestRunRecord:
    SETTINGS = {"match": "split", "root_edges": "unlabel", "jobs": 1, "log_level": "INFO"}

    def test_record_content(self, tmp_path):
        record = batch.write_run_record(tmp_path, "twig tree-label-transfer ...", self.SETTINGS, "(A,B);")
        saved = json.loads((tmp_path / "twig-run.json").read_text())
        assert saved == record
        assert saved["command"] == "twig tree-label-transfer ..."
        assert saved["sptree_hash"] == batch.text_hash("(A,B);")
        assert saved["twig_version"]

    def test_resume_accepts_the_same_run_and_ignores_runtime_options(self, tmp_path):
        batch.write_run_record(tmp_path, "cmd", self.SETTINGS, "(A,B);")
        changed = {**self.SETTINGS, "jobs": 8, "log_level": "DEBUG"}
        assert batch.check_resume(tmp_path, changed, "(A,B);")["command"] == "cmd"

    def test_resume_refuses_different_settings_or_species_tree(self, tmp_path):
        batch.write_run_record(tmp_path, "cmd", self.SETTINGS, "(A,B);")
        with pytest.raises(TwigError, match="the settings differ"):
            batch.check_resume(tmp_path, {**self.SETTINGS, "match": "clade"}, "(A,B);")
        with pytest.raises(TwigError, match="the species tree differ"):
            batch.check_resume(tmp_path, self.SETTINGS, "((A,B),C);")
        with pytest.raises(TwigError, match="the settings and species tree differ"):
            batch.check_resume(tmp_path, {**self.SETTINGS, "match": "clade"}, "((A,B),C);")

    def test_text_hash(self):
        assert batch.text_hash(None) is None
        assert batch.text_hash("(A,B);") == batch.text_hash("(A,B);") != batch.text_hash("(A,C);")

    def test_resume_without_a_record(self, tmp_path):
        with pytest.raises(TwigError, match="cannot resume"):
            batch.check_resume(tmp_path, self.SETTINGS)
