#!/usr/bin/env python

"""Command-line interface for twig.

Examples
--------
twig -h
twig [subcommand]
"""

from typing import Optional
import os
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from textwrap import dedent
from loguru import logger
import importlib
from . import subcommands
from twig import __version__ as VERSION
from twig.utils import TwigError


THREAD_ENV_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "BLIS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)

DISPATCH = {
    "align-pre": "..seqlab",
    "align-cds": "..seqlab",
    "align-post": "..seqlab",
    "consensus-fasta": "..seqlab",
    "csubst": "..seqlab",
    "diamond-bl": "..seqlab",
    "diamond-pw": "..seqlab",
    "filter-concat": "..seqlab",
    "format-fasta": "..synteny",
    "format-gff": "..synteny",
    "genome-table": "..synteny",
    # "hyphy-relax": "..seqlab",
    # "isoform-prune": "..seqlab",
    # "macse-trim": "..seqlab",
    # "macse-align": "..seqlab",
    # "macse-refine": "..seqlab",
    "name-map": "..labeling",
    "partition-cds": "..seqlab",
    "sample-reduce": "..labeling",
    "tree-filter": "..treelab",
    "tree-label": "..labeling",
    "tree-label-transfer": "..labeling",
    "tree-rooter": "..treelab",
    "tree-skeleton": "..treelab",
}


def setup_parsers() -> ArgumentParser:
    """Setup and return an ArgumentParser w/ subcommands."""
    parser = ArgumentParser(
        "twig",
        usage="twig [subcommand] --help",
        formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),
        description=dedent("""
            -----------------------------------------------------
            |  %(prog)s: tree-based workflows for integrative genomics |
            -----------------------------------------------------
            """),
        epilog=dedent(r"""
            Tutorials
            ---------
            # genome file formatting
            https://eaton-lab.org/twig/formatting

            # gene orthogroup analysis
            https://eaton-lab.org/twig/orthology

            # gene synteny analysis
            https://eaton-lab.org/twig/synteny

            # sequence rate/convergence/etc analysis
            https://eaton-lab.org/twig/evolution
        """)
    )
    parser.add_argument("-v", "--version    ", action='version', version=f"twig {VERSION}")
    return parser


def main(cmd: Optional[str] = None) -> int:
    """Command line tool.

    """
    # in cli pin numpy to single-threading
    for k in THREAD_ENV_VARS:
        os.environ.setdefault(k, '1')

    # load parser and attach subparsers
    parser = setup_parsers()
    subparsers = parser.add_subparsers(
        prog="%(prog)s", required=True,
        title="subcommands",
        dest="subcommand",
        metavar="--------------",
        help="-----------------------------------------------------",
    )
    for method in DISPATCH:
        meth = method.replace("-", "_")
        subparser = getattr(subcommands, f"get_parser_{meth}")
        subparser(subparsers)

    # parse args
    args = parser.parse_args(cmd.split() if cmd else None)

    # run subcommand
    if args.subcommand in DISPATCH:
        meth = args.subcommand.replace("-", "_")
        path = DISPATCH[args.subcommand]
        mod = importlib.import_module(f"{path}.{meth}", package=__package__)
        run_func = getattr(mod, f"run_{meth}")
        try:
            run_func(args)
            return 0
        except KeyboardInterrupt:
            logger.warning("interrupted by user")
        except TwigError as exc:
            logger.error(exc)
        except Exception as exc:
            logger.error(exc)
            raise
        return 1

    # unreachable
    parser.print_help()
    return 1


if __name__ == "__main__":
    main()
