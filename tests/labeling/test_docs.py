"""Tests that keep the labeling documentation honest."""

import re
import subprocess
from pathlib import Path

import pytest

DOCS = Path(__file__).parents[2] / "docs" / "labeling"
PAGES = sorted(DOCS.glob("*.md"))


def prose_lines(text):
    """Yield (line number, line) for lines outside fenced code blocks."""
    in_fence = False
    for number, line in enumerate(text.splitlines(), start=1):
        if line.lstrip().startswith("```"):
            in_fence = not in_fence
        elif not in_fence:
            yield number, line


def code_blocks(text):
    """Return the contents of fenced code blocks as a list of (language, text)."""
    return [(m.group(1), m.group(2)) for m in re.finditer(r"```(\w*)\n(.*?)```", text, re.S)]


def block_after(text, command):
    """Return the first code block that follows a bash block containing the command."""
    blocks = code_blocks(text)
    for index, (lang, body) in enumerate(blocks):
        if lang == "bash" and command in body:
            return blocks[index + 1][1]
    raise AssertionError(f"no code block follows {command!r}")


def fields(line):
    return re.split(r"\s{2,}", line.strip())


def test_pages_exist():
    assert {p.name for p in PAGES} >= {"design.md", "overview.md", "name-map.md"}


@pytest.mark.parametrize("page", PAGES, ids=lambda p: p.name)
class TestWritingRules:
    def test_no_em_dashes(self, page):
        assert "—" not in page.read_text() and "–" not in page.read_text()

    def test_no_colons_in_prose(self, page):
        bad = []
        for number, line in prose_lines(page.read_text()):
            without_code = re.sub(r"`[^`]*`", "", line)
            if ":" in without_code:
                bad.append(f"{page.name}:{number}: {line.strip()[:100]}")
        assert not bad, "\n".join(bad)

    def test_no_text_that_markdown_would_read_as_html(self, page):
        bad = []
        for number, line in prose_lines(page.read_text()):
            without_code = re.sub(r"`[^`]*`", "", line)
            if re.search(r"<[A-Za-z?/!][^>]*>", without_code):
                bad.append(f"{page.name}:{number}: {line.strip()[:100]}")
        assert not bad, "\n".join(bad)

    def test_code_blocks_are_balanced(self, page):
        assert page.read_text().count("```") % 2 == 0


class TestNameMapPage:
    @pytest.fixture
    def text(self):
        return (DOCS / "name-map.md").read_text()

    def test_help_text_is_current(self, text):
        from twig.labeling.parsers import get_parser_name_map
        documented = block_after(text, "twig name-map -h")
        # like the other pages, this one shows the help without its list of examples
        actual = get_parser_name_map().format_help().split("\nExamples\n")[0]
        assert documented.rstrip("\n") == actual.rstrip("\n")

    def test_example_table_is_current(self, text, names_dir):
        documented = block_after(text, "twig name-map -i trees aln -s species.nwk\n")
        result = subprocess.run(
            ["twig", "name-map", "-i", "trees", "aln", "-s", "species.nwk"],
            cwd=names_dir, capture_output=True, text=True, check=True,
        )
        actual = [line.split("\t") for line in result.stdout.splitlines()]
        assert [fields(line) for line in documented.splitlines()] == actual

    def test_problem_table_is_current(self, text, names_dir):
        documented = block_after(text, "species.nwk --show problems\n")
        result = subprocess.run(
            ["twig", "name-map", "-i", "trees", "aln", "-s", "species.nwk", "--show", "problems"],
            cwd=names_dir, capture_output=True, text=True, check=True,
        )
        actual = [line.split("\t") for line in result.stdout.splitlines()]
        assert [fields(line) for line in documented.splitlines()] == actual

    def test_frozen_map_is_current(self, text, names_dir, tmp_path):
        documented = block_after(text, "--freeze names.tsv")
        frozen = tmp_path / "names.tsv"
        subprocess.run(
            ["twig", "name-map", "-i", "trees", "aln", "-s", "species.nwk", "--freeze", str(frozen)],
            cwd=names_dir, capture_output=True, text=True, check=True,
        )
        assert documented == frozen.read_text()

    def test_summary_lines_are_current(self, text, names_dir):
        documented = next(body for _, body in code_blocks(text) if body.startswith("resolved 13 unique names"))
        result = subprocess.run(
            ["twig", "name-map", "-i", "trees", "aln", "-s", "species.nwk"],
            cwd=names_dir, capture_output=True, text=True, check=True,
        )
        prefix = r"^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} \| \w+\s+\| \S+\s+\| "
        messages = [re.sub(prefix, "", line) for line in result.stderr.splitlines()]
        messages = [m for m in messages if not m.startswith(("CMD:", "the species vocabulary"))]
        assert documented.splitlines() == messages
