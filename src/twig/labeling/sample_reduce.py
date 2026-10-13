#!/usr/bin/env python

"""Run the sample-reduce subcommand, which is planned for milestone M4.

The command is registered so that it appears in `twig -h`, but it is not
available yet. See docs/labeling/design.md for the design.
"""

from twig.utils import TwigError, set_log_level


def run_sample_reduce(args) -> None:
    """Run sample-reduce subcommand."""
    set_log_level(args.log_level)
    raise TwigError("sample-reduce is planned for milestone M4 and is not available yet, see docs/labeling/design.md")
