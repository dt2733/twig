#!/usr/bin/env python

"""Small data types shared by the labeling modules.

This module does not import toytree.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class Status(str, Enum):
    """The fate of a gene tree branch when labels are transferred to it."""
    LABELED = "labeled"
    SPECIES_UNLABELED = "species-unlabeled"
    UNMATCHED = "unmatched"
    COMPATIBLE = "compatible"
    CONFLICT = "conflict"
    MIXED = "mixed"
    PARALOG = "paralog"
    UNMAPPED = "unmapped"
    MASKED = "masked"

    # tables and f-strings write "labeled" and not "Status.LABELED"
    __str__ = str.__str__
    __format__ = str.__format__


@dataclass(frozen=True)
class RootMerge:
    """How the two branches at the root become one branch when unrooting.

    HyPhy and PAML read every tree as unrooted, so they see the two branches
    at the root of a rooted tree as a single branch. `dissolved` is the idx of
    the root child whose node is removed (its children move up to the root),
    `kept` is the idx of the other root child, which receives the summed
    branch length and the merged `label`. `labels` are the labels of the two
    original branches, and `conflict` says whether they disagreed.
    """
    dissolved: int
    kept: int
    labels: tuple[str | None, str | None]
    label: str | None
    conflict: bool
