"""Tests that the examples at the end of the tree-label help give the output that the help shows."""

import re
import shlex
import subprocess
from pathlib import Path

import pytest

pytest.importorskip("toytree")

from twig.cli.main import main
from twig.labeling.core.constants import RULE_NAMES
from twig.labeling.parsers import TREE_LABEL_EXAMPLES, get_parser_tree_label

# a line of the log, which the examples leave out
LOG = re.compile(r"^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} \|")


def parse_examples(text):
    """Return (command, shown output lines) for every line that starts with a dollar sign."""
    lines = text.splitlines()
    found, i = [], 0
    while i < len(lines):
        if lines[i].startswith("$ "):
            shown, j = [], i + 1
            while j < len(lines) and lines[j].strip() and not lines[j].startswith(("$ ", "# ")):
                shown.append(lines[j])
                j += 1
            found.append((lines[i][2:], shown))
            i = j
        else:
            i += 1
    return found


def visible(stdout, stderr):
    """Return the lines that a person sees, which are the output and the messages that are not log lines."""
    return stdout.splitlines() + [line for line in stderr.splitlines() if not LOG.match(line)]


def run_in_process(command, capsys):
    """Run a twig command in this process, with a redirect to a file as the shell would do it."""
    tokens = shlex.split(command)
    redirect = None
    if ">" in tokens:
        at = tokens.index(">")
        redirect, tokens = tokens[at + 1], tokens[:at]
    assert main(" ".join(tokens[1:])) == 0, command
    captured = capsys.readouterr()
    if redirect:
        Path(redirect).write_text(captured.out)
        return visible("", captured.err)
    return visible(captured.out, captured.err)


class TestExamples:
    def test_the_help_ends_with_the_examples(self):
        help_text = get_parser_tree_label().format_help()
        assert "\nExamples\n--------\n" in help_text
        assert help_text.rstrip("\n").endswith(TREE_LABEL_EXAMPLES.rstrip("\n"))

    # options that the examples leave out, because they only change how the command talks or are explained by name-map
    NOT_SHOWN = {"--help", "--force", "--log-level", "--accept-inferred", "--normalize", "--delim-idxs", "--delim-join"}

    def test_there_is_an_example_of_every_rule_and_of_every_option(self):
        commands = " ".join(c for c, _ in parse_examples(TREE_LABEL_EXAMPLES))
        for rule in RULE_NAMES:
            assert re.search(rf"--rule {rule}\b|-L \w+:{rule}\b", commands), rule
        parser = get_parser_tree_label()
        for action in parser._actions:
            options = action.option_strings
            if options and not set(options) & self.NOT_SHOWN:
                assert any(re.search(rf"(^|\s){re.escape(o)}(\s|$)", commands) for o in options), options

    def test_every_option_of_the_command_has_help(self):
        for action in get_parser_tree_label()._actions:
            assert action.help, action.option_strings

    def test_every_example_gives_the_output_that_the_help_shows(self, tmp_path, monkeypatch, capsys):
        monkeypatch.chdir(tmp_path)
        examples = parse_examples(TREE_LABEL_EXAMPLES)
        assert len(examples) >= 40
        for command, shown in examples:
            if command.startswith("twig ") and "|" not in command:
                got = run_in_process(command, capsys)
            else:
                result = subprocess.run(command, shell=True, cwd=tmp_path, capture_output=True, text=True)
                assert result.returncode == 0, (command, result.stderr)
                got = visible(result.stdout, result.stderr)
            assert got == shown, command
