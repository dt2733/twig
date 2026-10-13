#!/usr/bin/env python

"""Argument parsers for the twig tree labeling commands.

twig builds the parsers of every command each time it starts, so this module
imports only argparse, pathlib, textwrap and the dependency-free constants.
toytree is imported only by the run modules and the library.
"""

from __future__ import annotations

from argparse import ArgumentParser
from pathlib import Path
from textwrap import dedent

from twig.labeling.core.constants import (
    DEFAULT_NORMALIZE,
    DEFAULT_NUISANCE_LABEL,
    EXISTING_CHOICES,
    FORMAT_NAMES,
    INTERNAL_NAME_CHOICES,
    LEAF_NODE_CHOICES,
    NORMALIZE_CHOICES,
    PRECEDENCE_CHOICES,
    ROOT_EDGE_POLICIES,
    ROOT_STATE_CHOICES,
    RULE_NAMES,
    SELECT_ON_CHOICES,
    SHOW_CHOICES,
    UNMAPPED_POLICIES,
)

FRAME = "-" * 67


def _formatter(prog: str):
    """Return the help formatter used by the newer twig commands."""
    # imported here because twig.cli.subcommands imports this module at its end
    from twig.cli.subcommands import SingleMetavarHelpFormatter
    return SingleMetavarHelpFormatter(prog, width=220, max_help_position=56)


def _create(parent, kwargs: dict) -> ArgumentParser:
    """Return a subparser of parent, or a stand-alone parser if parent is None."""
    if parent:
        kwargs["name"] = kwargs.pop("prog")
        return parent.add_parser(**kwargs)
    kwargs.pop("help")
    return ArgumentParser(**kwargs)


def add_taxon_mapping_args(parser: ArgumentParser, unmapped: bool = True) -> None:
    """Add the group of options that map sample names to species.

    The group is the same in every labeling command. The option to choose what
    happens to unresolved names is added only by commands that act on them.
    """
    group = parser.add_argument_group("Taxon Mapping")
    group.add_argument("-I", "--imap", type=Path, metavar="path", help="map file rows <sample> <species> [specimen]; ~regex rows allowed")
    exclusive = group.add_mutually_exclusive_group()
    exclusive.add_argument("-d", "--delim", type=str, metavar="str", help="delimiter used to split sample names")
    exclusive.add_argument("-x", "--name-regex", type=str, metavar="str", help="regex with (?P<species>...) and optional (?P<specimen>...)")
    group.add_argument("-di", "--delim-idxs", type=int, metavar="int", nargs="+", default=[0], help="token indices kept after splitting [0]")
    group.add_argument("-dj", "--delim-join", type=str, metavar="str", help="join string for kept tokens [same as --delim]")
    group.add_argument("-N", "--normalize", type=str, metavar="str", nargs="+", choices=NORMALIZE_CHOICES, default=list(DEFAULT_NORMALIZE), help="normalization steps used for matching, or none [" + " ".join(DEFAULT_NORMALIZE) + "]")
    group.add_argument("--accept-inferred", action="store_true", help="accept genus-abbreviation matches")
    if unmapped:
        group.add_argument("-U", "--unmapped", type=str, metavar="str", choices=UNMAPPED_POLICIES, default="error", help="what to do with names that cannot be resolved (error, skip, drop, ignore) [%(default)s]")


def add_writer_args(parser: ArgumentParser) -> None:
    """Add the group of options that control how trees are written.

    The group is the same in every labeling command that writes trees.
    """
    group = parser.add_argument_group("Output Format")
    group.add_argument("--root-edges", type=str, metavar="str", default="unlabel", help="label of the one branch that HyPhy and PAML see in place of the two branches at the root (" + ", ".join(ROOT_EDGE_POLICIES) + ", or a label that wins) [%(default)s]")
    group.add_argument("--keep-root", action="store_true", help="write the tree rooted as it is, without joining the two branches at the root")
    group.add_argument("--internal-names", type=str, metavar="str", choices=INTERNAL_NAME_CHOICES, default="none", help="names written for internal nodes (" + ", ".join(INTERNAL_NAME_CHOICES) + ") [%(default)s]")
    group.add_argument("--no-dists", action="store_true", help="write trees without branch lengths")
    group.add_argument("--sanitize-suffix", action="store_true", help="add _2, _3 to tip names that become equal when they are made safe for HyPhy")
    group.add_argument("--names-table", type=Path, metavar="path", help="write the table of original and written tip names")


def get_parser_name_map(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for the name-map tool."""
    KWARGS = dict(
        prog="name-map",
        usage="name-map -i TREES/FASTA [...] -s SPTREE [options]",
        help="resolve tip and sequence names to species and freeze the mapping",
        formatter_class=_formatter,
        description=dedent(f"""
            {FRAME}
            | name-map resolves tip and sequence names to species names
            {FRAME}
            | Reads gene trees and/or FASTA files and resolves every tip and
            | sequence name to a species of the species tree (or of a species
            | list). Names are matched exactly, after normalization, by whole
            | tokens anywhere in the name, and by genus abbreviation. Matches
            | by abbreviation are reported as inferred and are not accepted by
            | other twig commands unless they are frozen into a map file with
            | --freeze or allowed with --accept-inferred. The table of names is
            | written to stdout and a summary to stderr.
            {FRAME}
        """),
        epilog=dedent("""
            Examples
            --------
            # check every name in a set of gene trees and alignments
            $ twig name-map -i 'trees/*.treefile' 'aln/*.nt.fa' -s sp.nwk

            # show only the names that need attention
            $ twig name-map -i 'trees/*.treefile' -s sp.nwk --show problems

            # parse the species from names such as OG0001-AmaranthusPumilus-s1
            $ twig name-map -i 'trees/*.treefile' -s sp.nwk -d '-' -di 1

            # save the resolved mapping for reuse with other twig commands
            $ twig name-map -i 'trees/*.treefile' 'aln/*.nt.fa' -s sp.nwk --freeze names.tsv

            Example map file (columns are sample, species, and optionally specimen)
            -----------------------------------------------------------------------
            A_pumilus_sample1           Amaranthus_pumilus    s1
            AmaranthusPumilus-Sample2   Amaranthus_pumilus    s2
            ~^Ahyb[_-]                  Amaranthus_hybridus
        """),
    )
    parser = _create(parser, KWARGS)

    io = parser.add_argument_group("Input/Output")
    io.add_argument("-i", "--input", type=Path, metavar="path", nargs="+", required=True, help="trees and/or FASTA files, directories or quoted globs")
    io.add_argument("-s", "--sptree", type=Path, metavar="path", help="species tree whose tips form the vocabulary")
    io.add_argument("-V", "--vocab", type=Path, metavar="path", help="file listing species names, one per line")
    io.add_argument("-o", "--out", type=Path, metavar="path", help="write the mapping table to this file [stdout]")
    io.add_argument("--freeze", type=Path, metavar="path", help="write the resolved map for reuse with -I")

    add_taxon_mapping_args(parser, unmapped=False)

    report = parser.add_argument_group("Reporting")
    report.add_argument("--show", type=str, metavar="str", choices=SHOW_CHOICES, default="all", help="rows of the table to print (all, problems) [%(default)s]")
    report.add_argument("--strict", action="store_true", help="error if any name is unmapped, ambiguous or inferred, or leaves alphabetic leftover tokens")

    runtime = parser.add_argument_group("Runtime")
    runtime.add_argument("-f", "--force", action="store_true", help="overwrite existing output files")
    runtime.add_argument("-l", "--log-level", type=str, metavar="level", default="INFO", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    return parser


def _get_planned_parser(parent, prog: str, summary: str, milestone: str) -> ArgumentParser:
    """Return a parser for a command that is designed but not available yet."""
    KWARGS = dict(
        prog=prog,
        usage=f"{prog} [options]",
        help=f"(planned for {milestone}) {summary}",
        formatter_class=_formatter,
        description=dedent(f"""
            {FRAME}
            | {prog} is planned for milestone {milestone} and is not available yet
            {FRAME}
            | {summary}
            | The design is in docs/labeling/design.md.
            {FRAME}
        """),
    )
    parser = _create(parent, KWARGS)
    parser.add_argument("-l", "--log-level", type=str, metavar="level", default="INFO", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    return parser


def get_parser_tree_label(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for the tree-label tool."""
    KWARGS = dict(
        prog="tree-label",
        usage="tree-label -i TREE [-L NAME SELECTOR...] [options]",
        help="label the branches of a tree from selections of tips",
        formatter_class=_formatter,
        description=dedent(f"""
            {FRAME}
            | tree-label labels the branches of a tree from selections of tips
            {FRAME}
            | Tips are selected by name (list files or regular expressions), or
            | by species with --select-on species, and a rule decides which
            | branches get each label. The rules are those of HyPhy label-tree
            | (tips, all, crown, some) and real two-state parsimony (fitch), plus
            | clade and stem. Stems, the branches on which a label begins, and
            | losses, the branches on which parsimony puts the loss of a trait,
            | can be left unlabeled or given a label of their own. Species that
            | belong in neither the test nor the reference set, for any reason,
            | can be set aside with --nuisance. Labels are written in braces
            | after the name, as HyPhy reads them, and --show lists the labeled
            | branches so the result can be checked by eye.
            {FRAME}
        """),
        epilog=dedent("""
            Examples
            --------
            # label the dioecious species and the branches where parsimony puts
            # the trait, leave the stems unlabeled, label the rest Reference
            $ twig tree-label -i sp.nwk -L Test dioecious.txt --rule fitch --stems unlabel --background Reference --show > sp.labeled.nwk

            # the same, with reversions found by parsimony labeled Rev
            $ twig tree-label -i sp.nwk -L Test dioecious.txt --rule fitch --stems unlabel --losses Rev --background Reference > sp.labeled.nwk

            # set species C aside as nuisance, whatever the reason, so that the
            # history of the trait is inferred without it
            $ twig tree-label -i sp.nwk -L Test dioecious.txt --rule fitch --nuisance C --background Reference > sp.labeled.nwk

            # label tips by regular expression as HyPhy label-tree does, and write a HyPhy tree
            $ twig tree-label -i gene.nwk -L Test '^Amaranthus_(hybridus|powellii)' -F hyphy -o gene.hyphy.nwk

            # two labels, each with its own rule
            $ twig tree-label -i sp.nwk -L Test:fitch dioecious.txt -L Out:tips '^O$' --background Reference

            # one table row per branch
            $ twig tree-label -i sp.nwk -L Test dioecious.txt -F table

            Example list file for a selection (one name per line, ~ starts a regular expression)
            ------------------------------------------------------------------------------------
            A
            B
            ~^D
        """),
    )
    parser = _create(parser, KWARGS)

    io = parser.add_argument_group("Input/Output")
    io.add_argument("-i", "--input", type=Path, metavar="path", required=True, help="file with one or more Newick trees, labeled or not")
    io.add_argument("-o", "--out", type=Path, metavar="path", help="write the output to this file [stdout]")
    io.add_argument("-F", "--format", type=str, metavar="str", choices=FORMAT_NAMES, default="braces", help="output format (" + ", ".join(FORMAT_NAMES) + ") [%(default)s]")
    io.add_argument("-B", "--branch-table", type=Path, metavar="path", help="also write a table with one row per branch")
    io.add_argument("--show", action="store_true", help="print the labeled branches, clade by clade, to stderr")
    io.add_argument("--paml-codes", type=str, metavar="str", nargs="+", help="names for PAML label numbers, such as 1=Test 2=Reference")

    select = parser.add_argument_group("Selection")
    select.add_argument("-L", "--label", type=str, metavar="str", nargs="+", action="append", help="label name followed by selectors, which are list files or regular expressions; repeat for more labels, and see the examples to give a label its own rule")
    select.add_argument("--spec", type=Path, metavar="path", help="tab-separated file with the columns label, rule, selectors and options, one label per line")
    select.add_argument("--nuisance", type=str, metavar="str", nargs="+", help="selectors for tips to set aside for any reason; they get the nuisance label, no other label overwrites it, and parsimony treats their state as unknown")
    select.add_argument("--nuisance-label", type=str, metavar="str", default=DEFAULT_NUISANCE_LABEL, help="label given to nuisance tips [%(default)s]")
    select.add_argument("--rule", "--internal-nodes", dest="rule", type=str, metavar="str", default="all", help="branches covered by a label (" + ", ".join(RULE_NAMES) + "), or a HyPhy name such as 'All descendants' [%(default)s]")
    select.add_argument("--select-on", type=str, metavar="str", choices=SELECT_ON_CHOICES, default="name", help="match selectors to tip names or to the species of tips (" + ", ".join(SELECT_ON_CHOICES) + ") [%(default)s]")
    select.add_argument("--exact", action="store_true", help="treat selectors given on the command line as exact names")
    select.add_argument("--invert", action="store_true", help="label the tips that the selectors do not match")
    select.add_argument("--leaf-nodes", type=str, metavar="str", choices=LEAF_NODE_CHOICES, default="label", help="label tips, or skip them after the rule has run (" + ", ".join(LEAF_NODE_CHOICES) + ") [%(default)s]")
    select.add_argument("--stems", type=str, metavar="str", default="keep", help="keep stems, unlabel them, or give them a label [%(default)s]")
    select.add_argument("--losses", type=str, metavar="str", default="keep", help="keep branches where parsimony puts a loss, unlabel them, or give them a label [%(default)s]")
    select.add_argument("--background", type=str, metavar="str", help="label given to every branch that is still free")
    select.add_argument("--precedence", type=str, metavar="str", choices=PRECEDENCE_CHOICES, default="specific", help="which label wins a branch claimed twice (" + ", ".join(PRECEDENCE_CHOICES) + ") [%(default)s]")
    select.add_argument("--existing", type=str, metavar="str", choices=EXISTING_CHOICES, default="keep", help="keep the labels already in the tree, or clear them first (" + ", ".join(EXISTING_CHOICES) + ") [%(default)s]")
    select.add_argument("-r", "--outgroups", type=str, metavar="str", nargs="+", help="root the tree on the first of these outgroups that fits, tips (or species) joined by commas for a clade")

    pars = parser.add_argument_group("Parsimony")
    pars.add_argument("--ambiguous", type=str, metavar="str", default="acctran", help="branches that differ between equally good histories follow acctran or deltran, stay unlabeled (unlabel), or get a label [%(default)s]")
    pars.add_argument("--gain-cost", type=str, metavar="float", default="1", help="cost of gaining the trait [%(default)s]")
    pars.add_argument("--loss-cost", type=str, metavar="float", default="1", help="cost of losing the trait [%(default)s]")
    pars.add_argument("--unknown", type=str, metavar="str", nargs="+", help="selectors for tips whose state is unknown")
    pars.add_argument("--root-state", type=str, metavar="str", choices=ROOT_STATE_CHOICES, default="free", help="state of the root (" + ", ".join(ROOT_STATE_CHOICES) + ") [%(default)s]")

    add_taxon_mapping_args(parser, unmapped=False)
    add_writer_args(parser)

    runtime = parser.add_argument_group("Runtime")
    runtime.add_argument("-f", "--force", action="store_true", help="overwrite existing output files")
    runtime.add_argument("-l", "--log-level", type=str, metavar="level", default="INFO", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    return parser


def get_parser_tree_label_transfer(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for the tree-label-transfer tool (planned)."""
    return _get_planned_parser(
        parser, "tree-label-transfer", "copy the labels of a species tree onto gene trees by branch matching", "M2",
    )


def get_parser_sample_reduce(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for the sample-reduce tool (planned)."""
    return _get_planned_parser(
        parser, "sample-reduce", "reduce gene trees and alignments to one specimen per species", "M4",
    )
