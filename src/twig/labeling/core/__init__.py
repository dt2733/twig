#!/usr/bin/env python

"""Library code for twig tree labeling.

Submodules are imported on first use, so importing this package does not load
toytree. Use the submodules directly, for example
`from twig.labeling.core.names import NameMapper`.
"""

from __future__ import annotations

import importlib

_SUBMODULES = (
    "batch",
    "constants",
    "fasta",
    "model",
    "names",
    "newick",
    "sanitize",
    "selectors",
    "treeio",
)

__all__ = list(_SUBMODULES)


def __getattr__(name: str):
    """Import a submodule the first time it is accessed as an attribute."""
    if name in _SUBMODULES:
        return importlib.import_module(f"{__name__}.{name}")
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
