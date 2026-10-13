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
    NORMALIZE_CHOICES,
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
    """Return a parser for the tree-label tool (planned)."""
    return _get_planned_parser(
        parser, "tree-label", "label the branches of a tree from tip selections", "M1",
    )


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
