#!/usr/bin/env python

"""Run the tree-label-transfer subcommand, which is planned for milestone M2.

The command is registered so that it appears in `twig -h`, but it is not
available yet. See docs/labeling/design.md for the design.
"""

from twig.utils import TwigError, set_log_level


def run_tree_label_transfer(args) -> None:
    """Run tree-label-transfer subcommand."""
    set_log_level(args.log_level)
    raise TwigError("tree-label-transfer is planned for milestone M2 and is not available yet, see docs/labeling/design.md")
