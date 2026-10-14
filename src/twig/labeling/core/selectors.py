#!/usr/bin/env python

"""Select tips by name from a list file or a regular expression.

A selector value that is an existing file is read as a list with one name per
line, where lines that start with `~` are regular expressions and lines that
start with `#` are comments. Any other value is a regular expression that is
matched with `re.search`, or an exact name if `exact` is set, so that `A` does
not also select `MA`.

This module does not import toytree.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Sequence
from pathlib import Path

from twig.utils import TwigError


def _compile(pattern: str, where: str) -> re.Pattern:
    try:
        return re.compile(pattern)
    except re.error as exc:
        raise TwigError(f"{where}: invalid regular expression {pattern!r} ({exc})") from exc


def read_name_list(path: str | Path) -> tuple[list[str], list[str]]:
    """Return the exact names and the regex patterns listed in a file."""
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8-sig")       # utf-8-sig drops a byte order mark
    except OSError as exc:
        raise TwigError(f"cannot read selector list {path}: {exc}") from exc
    except UnicodeDecodeError as exc:
        raise TwigError(f"selector list {path} is not valid UTF-8 text") from exc
    names, patterns = [], []
    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("~"):
            pattern = line[1:]
            _compile(pattern, f"{path}:{lineno}")
            patterns.append(pattern)
        else:
            names.append(line)
    return names, patterns


class Selector:
    """A set of tip names chosen by lists, exact names and regular expressions.

    A name is selected if it equals a listed name or matches any pattern.

    Parameters
    ----------
    values: Sequence[str]
        Existing files are read as lists. Other values are regular expressions,
        or exact names if `exact` is True.
    exact: bool
        Treat values that are not files as exact names.
    """

    def __init__(self, values: Sequence[str | Path] = (), exact: bool = False):
        self.values = tuple(str(i) for i in values)
        self.exact = exact
        self.names: dict[str, str] = {}
        self.patterns: dict[str, re.Pattern] = {}
        for value in self.values:
            if Path(value).is_file():
                listed, patterns = read_name_list(value)
                for name in listed:
                    self.names.setdefault(name, value)
                for pattern in patterns:
                    self.patterns.setdefault(pattern, _compile(pattern, value))
            elif exact:
                self.names.setdefault(value, "command line")
            else:
                self.patterns.setdefault(value, _compile(value, "command line"))

    def __bool__(self) -> bool:
        return bool(self.names or self.patterns)

    def matches(self, name: str) -> bool:
        """Return True if name is selected."""
        if name in self.names:
            return True
        return any(p.search(name) for p in self.patterns.values())

    def select(self, names: Iterable[str]) -> list[str]:
        """Return the selected names, in the order given."""
        return [name for name in names if self.matches(name)]

    def unused(self, names: Iterable[str]) -> list[str]:
        """Return listed names and patterns that select none of the names.

        This is useful for warning about typos in a list file.
        """
        names = list(names)
        present = set(names)
        missing = [i for i in self.names if i not in present]
        for pattern, compiled in self.patterns.items():
            if not any(compiled.search(name) for name in names):
                missing.append("~" + pattern)
        return missing

    def describe(self) -> str:
        """Return a short description for log messages."""
        return f"{len(self.names)} names and {len(self.patterns)} patterns from {', '.join(self.values) or 'nothing'}"
