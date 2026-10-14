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


def normalized(line):
    """Return the non-empty fields of a line of a table, whether they are separated by tabs or by spaces."""
    return [f for f in re.split(r"\t|\s{2,}", line.strip()) if f]


def command_pairs(text, prefix):
    """Return (command, output) for every bash block that starts with prefix and is followed by a text block."""
    blocks = code_blocks(text)
    return [
        (body.strip(), blocks[i + 1][1])
        for i, (lang, body) in enumerate(blocks[:-1])
        if lang == "bash" and body.startswith(prefix) and blocks[i + 1][0] == "text"
    ]


def table_rows(text, after):
    """Return the cells of the rows of the markdown table that follows a given sentence or starts with a given line."""
    rows, started = [], False
    for line in text[text.index(after):].splitlines():
        if line.startswith("|"):
            started = True
            rows.append([c.strip() for c in line.strip().strip("|").split("|")])
        elif started:
            break
    return rows[2:]


class TestTreeLabelPage:
    @pytest.fixture
    def text(self):
        return (DOCS / "tree-label.md").read_text()

    def test_help_text_is_current(self, text):
        from twig.labeling.parsers import get_parser_tree_label
        documented = block_after(text, "twig tree-label -h")
        actual = get_parser_tree_label().format_help().split("\nExamples\n")[0]
        assert documented.rstrip("\n") == actual.rstrip("\n")

    def test_every_example_gives_the_documented_output(self, text, data_dir):
        pairs = [(c, o) for c, o in command_pairs(text, "twig tree-label") if c != "twig tree-label -h"]
        assert len(pairs) >= 8
        for command, documented in pairs:
            result = subprocess.run(command, shell=True, cwd=data_dir, capture_output=True, text=True)
            assert result.returncode == 0, (command, result.stderr)
            assert [normalized(i) for i in documented.splitlines()] == [normalized(i) for i in result.stdout.splitlines()], command

    def test_the_table_of_rules_is_current(self, text, data_dir):
        pytest.importorskip("toytree")
        from twig.labeling.core import treeio
        rows = table_rows(text, "The next table shows what each rule labels")
        assert len(rows) >= 7
        for options, expected in rows:
            words = options.strip("`").split()
            command = ["twig", "tree-label", "-i", "species.named.nwk", "-L", "Test", "dioecious.txt",
                       "--rule", *words, "--internal-names", "keep", "-l", "ERROR"]
            result = subprocess.run(command, cwd=data_dir, capture_output=True, text=True, check=True)
            tree = treeio.read_tree(result.stdout.strip())
            assert {n.name for n in tree if n.label} == set(expected.split(", ")), options

    def test_the_rules_and_their_hyphy_names_are_listed(self, text):
        from twig.labeling.core.constants import RULE_NAMES
        from twig.labeling.core.rules import resolve_rule
        rows = table_rows(text, "| rule | HyPhy name | branches covered |")
        assert [r[0].strip("`") for r in rows] == list(RULE_NAMES)
        for rule, hyphy, _ in rows:
            if hyphy != "not in HyPhy":
                assert resolve_rule(hyphy) == rule.strip("`")

    def test_the_hyphy_options_and_the_skip_mode_are_explained(self, text):
        from twig.labeling.parsers import get_parser_tree_label
        rows = table_rows(text, "| HyPhy `label-tree` | twig `tree-label` |")
        assert len(rows) == 7
        known = {o for a in get_parser_tree_label()._actions for o in a.option_strings}
        for _, twig_side in rows:
            flags = re.findall(r"(?<![\w-])(-{1,2}[A-Za-z][\w-]*)", twig_side.split("`")[1])
            assert flags and set(flags) <= known, twig_side
        assert "HyPhy's Skip mode is not the same as `--leaf-nodes skip`" in text

    def test_the_root_edge_values_are_listed(self, text):
        from twig.labeling.core.constants import ROOT_EDGE_POLICIES
        rows = table_rows(text, "| value | label of the joined branch |")
        assert [r[0].strip("`") for r in rows[:-1]] == list(ROOT_EDGE_POLICIES) and rows[-1][0] == "a label name"

    def test_the_formats_are_listed(self, text):
        from twig.labeling.core.constants import FORMAT_NAMES
        rows = table_rows(text, "| format | output |")
        assert [r[0].strip("`") for r in rows] == list(FORMAT_NAMES)

    def test_the_columns_of_the_branch_table_are_listed(self, text):
        from twig.labeling.core.writers import TableWriter, WriteContext
        rows = table_rows(text, "| column | meaning |")
        listed = [r[0].strip("`") for r in rows]
        assert listed == list(TableWriter().columns(WriteContext(species_of=str)))


class TestTreeLabelExamplesOnThePage:
    def test_the_examples_of_the_help_are_on_the_page(self):
        from twig.labeling.parsers import TREE_LABEL_EXAMPLES
        text = (DOCS / "tree-label.md").read_text()
        documented = next(body for lang, body in code_blocks(text) if lang == "text" and body.startswith("Examples\n--------"))
        assert documented.rstrip("\n") == TREE_LABEL_EXAMPLES.rstrip("\n")

    def test_the_page_says_what_fitch_means(self):
        text = (DOCS / "tree-label.md").read_text()
        assert "Walter Fitch" in text and "1971" in text and "The rule `parsimony` means the same as `fitch`" in text
