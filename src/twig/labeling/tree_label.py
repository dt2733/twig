#!/usr/bin/env python

"""Run the tree-label subcommand, which is planned for milestone M1.

The command is registered so that it appears in `twig -h`, but it is not
available yet. See docs/labeling/design.md for the design.
"""

from twig.utils import TwigError, set_log_level


def run_tree_label(args) -> None:
    """Run tree-label subcommand."""
    set_log_level(args.log_level)
    raise TwigError("tree-label is planned for milestone M1 and is not available yet, see docs/labeling/design.md")
