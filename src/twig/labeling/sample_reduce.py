#!/usr/bin/env python

"""Run the sample-reduce subcommand, which is planned for milestone M4.

The command is registered so that it appears in `twig -h`, but it is not
available yet. See docs/labeling/design.md for the design.
"""

import sys

from loguru import logger

from twig.utils import TwigError, set_log_level


def run_sample_reduce(args) -> None:
    """Run sample-reduce subcommand."""
    set_log_level(args.log_level)
    raise TwigError("sample-reduce is planned for milestone M4 and is not available yet, see docs/labeling/design.md")


def main():
    """Run sample-reduce as a script, like the other twig run modules."""
    from .parsers import get_parser_sample_reduce
    args = get_parser_sample_reduce().parse_args()
    run_sample_reduce(args)


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
