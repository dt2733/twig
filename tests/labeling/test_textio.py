"""Tests for reading text files (twig.labeling.core.textio)."""

import pytest

from twig.utils import TwigError
from twig.labeling.core.textio import read_text, reading


class TestReadText:
    def test_the_text_of_a_file(self, tmp_path):
        path = tmp_path / "a.txt"
        path.write_text("one\ntwo\n")
        assert read_text(path, "list") == "one\ntwo\n"

    def test_a_byte_order_mark_is_dropped(self, tmp_path):
        path = tmp_path / "a.txt"
        path.write_bytes("\ufeffone\n".encode("utf-8"))
        assert read_text(path, "list") == "one\n"

    def test_a_missing_file_is_named_and_not_found(self, tmp_path):
        with pytest.raises(TwigError) as info:
            read_text(tmp_path / "missing.txt", "map file")
        assert str(info.value) == f"map file {tmp_path / 'missing.txt'} not found"

    def test_a_directory_cannot_be_read(self, tmp_path):
        with pytest.raises(TwigError) as info:
            read_text(tmp_path, "tree file")
        assert str(info.value).startswith(f"cannot read the tree file {tmp_path} (")
        assert ":" not in str(info.value).replace(str(tmp_path), "")

    def test_a_file_that_is_not_text(self, tmp_path):
        path = tmp_path / "a.txt"
        path.write_bytes(b"\xff\xfe\x00")
        with pytest.raises(TwigError) as info:
            read_text(path, "species list")
        assert str(info.value) == f"species list {path} is not valid UTF-8 text"

    def test_the_kind_of_file_is_part_of_every_message(self, tmp_path):
        for what in ("tree file", "FASTA file", "selector list", "label spec", "manifest"):
            with pytest.raises(TwigError, match=f"^{what} "):
                read_text(tmp_path / "missing.txt", what)

    def test_other_errors_pass_through_the_context_manager(self, tmp_path):
        with pytest.raises(TwigError, match="my own"):
            with reading(tmp_path / "x", "file"):
                raise TwigError("my own error")
        with pytest.raises(KeyError):
            with reading(tmp_path / "x", "file"):
                raise KeyError("x")


class TestTheSameMessageInEveryCommand:
    """A missing input file used to be a raw OSError in one command and another wording in the next."""

    def test_missing_input_files(self, run_twig, tmp_path):
        missing = tmp_path / "missing.nwk"
        results = {
            "tree-label": run_twig("tree-label", "-i", missing),
            "name-map": run_twig("name-map", "-i", missing),
        }
        for command, (code, out, err) in results.items():
            assert code == 1 and out == "", command
            assert f"{missing} not found" in err, command
            assert "Traceback" not in err and "Errno" not in err, command
        assert "tree file" in results["tree-label"][2] and "input file" in results["name-map"][2]
