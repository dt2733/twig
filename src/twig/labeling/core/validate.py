#!/usr/bin/env python

"""Problems found while checking trees for the programs that read them.

This module does not import toytree.
"""

from __future__ import annotations

from dataclasses import dataclass

ERROR = "error"
WARNING = "warning"


@dataclass(frozen=True)
class Issue:
    """A problem with a tree.

    An error means the tree cannot be used by the program it was written for,
    and a warning means it can be used but something should be looked at. The
    code is a short fixed word and the locus names the tree or gene.
    """
    level: str
    code: str
    message: str
    locus: str | None = None

    def __str__(self) -> str:
        where = f"{self.locus}: " if self.locus else ""
        return f"{where}{self.message}"


def errors(issues) -> list[Issue]:
    """Return the issues that are errors."""
    return [i for i in issues if i.level == ERROR]


def warnings(issues) -> list[Issue]:
    """Return the issues that are warnings."""
    return [i for i in issues if i.level == WARNING]
