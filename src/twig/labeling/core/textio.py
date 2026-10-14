#!/usr/bin/env python

"""Read text files for the labeling commands.

Every input file is read in one place, so that a file that is missing, a file
that cannot be read and a file that is not text give the same kind of message
in every command, and so that the byte order mark that some editors add is
ignored.

This module does not import toytree.
"""

from __future__ import annotations

from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

from twig.utils import TwigError


@contextmanager
def reading(path: str | Path, what: str = "file") -> Iterator[None]:
    """Turn the errors of reading a file into a TwigError that names the file.

    `what` says what kind of file it is, such as "tree file" or "map file".
    """
    try:
        yield
    except FileNotFoundError:
        raise TwigError(f"{what} {path} not found") from None
    except OSError as exc:
        raise TwigError(f"cannot read the {what} {path} ({exc.strerror or exc})") from exc
    except UnicodeDecodeError as exc:
        raise TwigError(f"{what} {path} is not valid UTF-8 text") from exc


def read_text(path: str | Path, what: str = "file") -> str:
    """Return the text of a UTF-8 file without a byte order mark."""
    path = Path(path)
    with reading(path, what):
        return path.read_text(encoding="utf-8-sig")
