#!/usr/bin/env python

"""Run the tree-label-transfer subcommand, which is planned for milestone M2.

The command is registered so that it appears in `twig -h`, but it is not
available yet. See docs/labeling/design.md for the design.
"""

import sys

from loguru import logger

from twig.utils import TwigError, set_log_level


def run_tree_label_transfer(args) -> None:
    """Run tree-label-transfer subcommand."""
    set_log_level(args.log_level)
    raise TwigError("tree-label-transfer is planned for milestone M2 and is not available yet, see docs/labeling/design.md")


def main():
    """Run tree-label-transfer as a script, like the other twig run modules."""
    from .parsers import get_parser_tree_label_transfer
    args = get_parser_tree_label_transfer().parse_args()
    run_tree_label_transfer(args)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.warning("interrupted by user")
        sys.exit(1)
    except TwigError as exc:
        logger.error(exc)
        sys.exit(1)
    except Exception as exc:
        logger.error(exc)
        raise
