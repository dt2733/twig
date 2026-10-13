"""Shared fixtures for the tree labeling tests."""

from pathlib import Path

import pytest

DATA = Path(__file__).parent / "data"

# the species tree and the labeled species tree used throughout the design
SPECIES = "(O,((M1,(A,(B,C))),(M2,(D,E))));"
LABELED = (
    "(O{Reference},((M1{Reference},(A{Test},(B{Test},C{Rev}){Test})){Reference},"
    "(M2{Reference},(D{Test},E{Test})){Reference}){Reference});"
)


@pytest.fixture
def data_dir() -> Path:
    """The directory of test data files."""
    return DATA


@pytest.fixture
def names_dir() -> Path:
    """The directory of data for name mapping, with trees, alignments and a species tree."""
    return DATA / "names"


@pytest.fixture
def run_twig(capsys):
    """Return a function that runs `twig` with arguments and returns (code, stdout, stderr).

    Arguments are joined with spaces, as `twig.cli.main.main` splits its
    command on whitespace, so paths must not contain spaces.
    """
    from twig.cli.main import main

    def run(*args) -> tuple[int, str, str]:
        code = main(" ".join(str(i) for i in args))
        captured = capsys.readouterr()
        return code, captured.out, captured.err

    return run
