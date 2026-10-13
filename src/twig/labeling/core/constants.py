#!/usr/bin/env python

"""Names shared by the labeling parsers and the labeling library.

This module has no dependencies, so the argument parsers can import it at
startup without loading toytree or any other heavy package.
"""

from __future__ import annotations

# node feature that stores the label of the branch above a node
LABEL = "label"

# default regular expression that turns a file name into a locus name
DEFAULT_LOCUS_REGEX = r"^([^.]+)"

# normalization steps used when comparing names (applied in this order)
NORMALIZE_STEPS = ("accents", "camel", "digits", "sep", "case")
DEFAULT_NORMALIZE = ("accents", "camel", "sep", "case")
NORMALIZE_CHOICES = NORMALIZE_STEPS + ("none",)

# what to do with names that cannot be resolved to a species
UNMAPPED_POLICIES = ("error", "skip", "drop", "ignore")

# rows printed by name-map
SHOW_CHOICES = ("all", "problems")

# how a name was resolved, and whether the result can be used
HIT_METHODS = (
    "map",
    "map-regex",
    "delim",
    "regex",
    "exact",
    "canon",
    "token",
    "abbrev",
    "self",
    "none",
)
HIT_STATUSES = ("ok", "inferred", "ambiguous", "unmapped")

# statuses of gene tree branches after labels are transferred
STATUS_NAMES = (
    "labeled",
    "species-unlabeled",
    "unmatched",
    "compatible",
    "conflict",
    "mixed",
    "paralog",
    "unmapped",
    "masked",
)

# syntaxes recognized when reading labels from a tree string
LABEL_SYNTAXES = ("auto", "braces", "paml", "nhx", "none")
