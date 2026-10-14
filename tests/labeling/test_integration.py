"""Tests for how the labeling commands live inside twig."""

import importlib
import re
import subprocess
import sys
from argparse import ArgumentParser

import pytest

from twig.cli.main import DISPATCH
from twig.labeling import parsers
from twig.labeling.core import constants
from twig.labeling.core.model import Status

COMMANDS = ["name-map", "sample-reduce", "tree-label", "tree-label-transfer"]


def writing_problems(text: str) -> list[str]:
    """Return the lines of help text with a colon or a dash that the writing rules forbid.

    The rules are no em-dashes, no en-dashes and no colons in prose. Colons are
    allowed in the lines that argparse makes (usage and group titles), in command
    lines that start with a dollar sign, in Newick output and in tables.
    """
    found = []
    for line in text.splitlines():
        if line.startswith("usage: ") or (line.endswith(":") and not line.startswith(" ")):
            continue
        if line.startswith("$ ") or "\t" in line or re.match(r"^\(.*\);$", line):
            continue
        if ":" in line or "\u2014" in line or "\u2013" in line:
            found.append(line.strip())
    return found


def run_python(code: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)


class TestStartup:
    """twig builds every parser at every start, so these modules must stay light."""

    @pytest.mark.parametrize("module", [
        "twig.cli.subcommands",
        "twig.cli.main",
        "twig.labeling.parsers",
        "twig.labeling.core.constants",
    ])
    def test_importing_does_not_load_toytree(self, module):
        result = run_python(f"import sys, {module}; assert 'toytree' not in sys.modules, sorted(m for m in sys.modules if 'toytree' in m)")
        assert result.returncode == 0, result.stderr

    def test_building_every_parser_does_not_load_toytree(self):
        code = (
            "import sys\n"
            "from twig.cli.main import setup_parsers, DISPATCH\n"
            "from twig.cli import subcommands\n"
            "p = setup_parsers(); sp = p.add_subparsers(dest='subcommand')\n"
            "for m in DISPATCH: getattr(subcommands, 'get_parser_' + m.replace('-', '_'))(sp)\n"
            "assert 'toytree' not in sys.modules\n"
        )
        result = run_python(code)
        assert result.returncode == 0, result.stderr

    def test_the_package_init_is_empty(self):
        import twig.labeling
        assert twig.labeling.__file__.endswith("__init__.py")
        with open(twig.labeling.__file__) as handle:
            assert handle.read().strip() == ""


class TestRegistration:
    def test_commands_are_dispatched_to_the_labeling_package(self):
        for name in COMMANDS:
            assert DISPATCH[name] == "..labeling"

    def test_module_and_function_names_follow_the_command_names(self):
        # main imports twig.labeling.<name> and calls run_<name>
        for name in COMMANDS:
            attr = name.replace("-", "_")
            module = importlib.import_module(f"twig.labeling.{attr}")
            assert callable(getattr(module, f"run_{attr}"))
            assert callable(getattr(parsers, f"get_parser_{attr}"))

    def test_subcommands_module_exposes_the_parser_functions(self):
        from twig.cli import subcommands
        for name in COMMANDS:
            attr = f"get_parser_{name.replace('-', '_')}"
            assert getattr(subcommands, attr) is getattr(parsers, attr)

    def test_existing_commands_are_untouched(self):
        for name in ("tree-filter", "tree-rooter", "tree-skeleton", "filter-concat", "align-pre", "csubst"):
            assert name in DISPATCH and DISPATCH[name] != "..labeling"


class TestParsers:
    @pytest.mark.parametrize("name", COMMANDS)
    def test_stand_alone_parser(self, name):
        parser = getattr(parsers, f"get_parser_{name.replace('-', '_')}")()
        assert isinstance(parser, ArgumentParser) and parser.prog == name

    def test_name_map_defaults(self):
        args = parsers.get_parser_name_map().parse_args(["-i", "x.nwk"])
        assert args.delim_idxs == [0] and args.delim_join is None and args.normalize == ["accents", "camel", "sep", "case"]
        assert args.show == "all" and not args.strict and not args.accept_inferred and args.log_level == "INFO"

    def test_taxon_mapping_group_options(self):
        args = parsers.get_parser_name_map().parse_args(
            ["-i", "x.nwk", "-d", "_", "-di", "0", "1", "-dj", "-", "-N", "case", "sep", "--accept-inferred", "-I", "map.tsv"])
        assert (args.delim, args.delim_idxs, args.delim_join, args.normalize) == ("_", [0, 1], "-", ["case", "sep"])
        assert str(args.imap) == "map.tsv"

    def test_delimiter_and_regex_are_mutually_exclusive(self):
        with pytest.raises(SystemExit):
            parsers.get_parser_name_map().parse_args(["-i", "x", "-d", "_", "-x", "(?P<species>.+)"])

    def test_choices_are_checked(self):
        parser = parsers.get_parser_name_map()
        with pytest.raises(SystemExit):
            parser.parse_args(["-i", "x", "--show", "some"])
        with pytest.raises(SystemExit):
            parser.parse_args(["-i", "x", "-N", "upper"])

    def test_the_unmapped_option_is_only_added_where_it_acts(self):
        parser = ArgumentParser()
        parsers.add_taxon_mapping_args(parser)
        assert parser.parse_args([]).unmapped == "error"
        parser = ArgumentParser()
        parsers.add_taxon_mapping_args(parser, unmapped=False)
        assert not hasattr(parser.parse_args([]), "unmapped")

    def test_help_text_has_no_dashes_and_no_prose_colons(self):
        for name in COMMANDS:
            parser = getattr(parsers, f"get_parser_{name.replace('-', '_')}")()
            assert writing_problems(parser.format_help()) == [], name

    def test_the_writing_check_finds_what_it_is_meant_to_find(self):
        text = "\n".join([
            "usage: tree-label [options]",
            "Input/Output:",
            "  -i, --input path   fine words",
            "  -o, --out path     a colon: here",
            "# a comment with a colon: here",
            "# a comment with an em\u2014dash",
            "A line with an en\u2013dash",
            "$ twig tree-label -L Test:fitch dioecious.txt",
            "(A{Test},B);",
            "A\tB: tabs mean a table",
        ])
        found = writing_problems(text)
        assert [line.split(" ")[0] for line in found] == ["-o,", "#", "#", "A"]


class TestCommandLine:
    def test_twig_help_lists_the_new_commands(self, run_twig, capsys):
        with pytest.raises(SystemExit) as info:
            run_twig("-h")
        assert info.value.code == 0
        out = capsys.readouterr().out
        for name in COMMANDS:
            assert re.search(rf"^\s+{re.escape(name)}\s", out, re.M), name

    def test_twig_help_text(self, capsys):
        from twig.cli.main import main
        with pytest.raises(SystemExit):
            main("-h")
        out = capsys.readouterr().out
        for name in COMMANDS:
            assert name in out
        assert "(planned for M2)" in out and "resolve tip and sequence names" in out
        assert "label the branches of a tree from selections of tips" in out

    @pytest.mark.parametrize("name", COMMANDS)
    def test_each_command_has_help(self, name, capsys):
        from twig.cli.main import main
        with pytest.raises(SystemExit) as info:
            main(f"{name} -h")
        assert info.value.code == 0
        assert name in capsys.readouterr().out

    @pytest.mark.parametrize("name, milestone", [("tree-label-transfer", "M2"), ("sample-reduce", "M4")])
    def test_planned_commands_say_so_and_fail_cleanly(self, run_twig, name, milestone):
        code, out, err = run_twig(name)
        assert code == 1 and out == ""
        assert f"planned for milestone {milestone}" in err and "Traceback" not in err


class TestScripts:
    """Every twig run module can be run as a script, so the labeling modules can too."""

    @pytest.mark.parametrize("name", COMMANDS)
    def test_each_run_module_prints_its_help(self, name):
        module = "twig.labeling." + name.replace("-", "_")
        result = subprocess.run([sys.executable, "-m", module, "-h"], capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        assert result.stdout.startswith(f"usage: {name}")

    def test_a_user_error_is_a_message_and_not_a_traceback(self, tmp_path):
        result = subprocess.run(
            [sys.executable, "-m", "twig.labeling.tree_label", "-i", str(tmp_path / "none.nwk")],
            capture_output=True, text=True,
        )
        assert result.returncode == 1 and "ERROR" in result.stderr and "Traceback" not in result.stderr

    def test_a_planned_command_says_so_and_fails_cleanly(self):
        result = subprocess.run([sys.executable, "-m", "twig.labeling.sample_reduce"], capture_output=True, text=True)
        assert result.returncode == 1 and "planned for milestone M4" in result.stderr and "Traceback" not in result.stderr

    def test_a_script_does_the_work_like_the_twig_command(self, tmp_path):
        tree = tmp_path / "t.nwk"
        tree.write_text("((A,B),(C,D));")
        result = subprocess.run(
            [sys.executable, "-m", "twig.labeling.tree_label", "-i", str(tree), "-L", "Test", "^[AB]$", "--rule", "tips", "-l", "ERROR"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, result.stderr
        assert result.stdout == "((A{Test},B{Test}),(C,D));\n"


class TestConstants:
    def test_status_names_match_the_status_enum(self):
        assert list(constants.STATUS_NAMES) == [s.value for s in Status]

    def test_status_formats_as_its_value(self):
        assert f"{Status.LABELED}" == "labeled" and str(Status.SPECIES_UNLABELED) == "species-unlabeled"
        assert "\t".join([Status.CONFLICT, Status.MIXED]) == "conflict\tmixed"

    def test_normalize_choices_are_all_accepted(self):
        from twig.labeling.core.names import normalize_steps
        assert normalize_steps(list(constants.NORMALIZE_STEPS)) == constants.NORMALIZE_STEPS
        assert normalize_steps(["none"]) == ()
        assert set(constants.DEFAULT_NORMALIZE) <= set(constants.NORMALIZE_STEPS)

    def test_core_package_imports_submodules_lazily(self):
        code = ("import sys, twig.labeling.core as core; assert 'twig.labeling.core.names' not in sys.modules; "
                "core.names; assert 'twig.labeling.core.names' in sys.modules; "
                "import pytest\nwith pytest.raises(AttributeError): core.nothing")
        result = run_python(code)
        assert result.returncode == 0, result.stderr

    def test_the_tie_rules_have_one_definition(self):
        import inspect
        from twig.labeling.core import labelsets, parsimony, rules
        assert constants.AMBIGUOUS_MODES == constants.TIE_RULES + ("unlabel",)
        defaults = {
            "parser": parsers.get_parser_tree_label().get_default("ambiguous"),
            "config": labelsets.LabelingConfig().ambiguous,
            "rule options": rules.RuleOptions().ambiguous,
            "make_options": rules.make_options().ambiguous,
            "parsimony": inspect.signature(parsimony.sankoff_binary).parameters["tie"].default,
        }
        assert set(defaults.values()) == {constants.DEFAULT_AMBIGUOUS}, defaults
        assert parsimony.TIE_RULES == constants.TIE_RULES and rules.TIE_RULES == constants.TIE_RULES

    def test_the_constants_that_nothing_used_are_gone(self):
        assert not hasattr(constants, "HIT_METHODS") and not hasattr(constants, "HIT_STATUSES")
