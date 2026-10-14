"""Tests that error and log messages follow the writing rules, which allow no colons in prose."""

import ast
from pathlib import Path

import pytest

SOURCE = Path(__file__).parents[2] / "src" / "twig" / "labeling"

# a name, a path with a line number, or a label that starts a message may be followed by a colon
PREFIXES = ("{}:{}: ", "{}: ", "{}, tree {}: ", "label {}: ", "CMD: ")


def text_of(node) -> list[str]:
    """Return the text of a string, an f-string or a sum of them, with {} for each placeholder."""
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return [node.value]
    if isinstance(node, ast.JoinedStr):
        return ["".join(v.value if isinstance(v, ast.Constant) else "{}" for v in node.values)]
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
        return text_of(node.left) + text_of(node.right)
    return []


def messages(tree):
    """Yield (line number, text) for every message that a module raises, logs or collects."""
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call) or not node.args:
            continue
        func = node.func
        name = func.attr if isinstance(func, ast.Attribute) else getattr(func, "id", "")
        owner = func.value.id if isinstance(func, ast.Attribute) and isinstance(func.value, ast.Name) else ""
        if name in ("TwigError", "ArgumentTypeError") or (owner == "logger" and name in ("warning", "info", "error", "debug")):
            candidates = node.args[:1]
        elif name == "Issue":
            candidates = node.args[2:3]
        elif name == "append" and owner in ("warnings", "notes"):
            candidates = node.args[:1]
        else:
            continue
        for arg in candidates:
            for text in text_of(arg):
                yield node.lineno, text


def colons_in_prose(text: str) -> bool:
    for prefix in PREFIXES:
        if text.startswith(prefix):
            text = text[len(prefix):]
            break
    return ":" in text


@pytest.mark.parametrize("path", sorted(SOURCE.rglob("*.py")), ids=lambda p: str(p.relative_to(SOURCE)))
def test_messages_have_no_colons_in_prose(path):
    bad = [f"{path.name}:{lineno} {text!r}" for lineno, text in messages(ast.parse(path.read_text())) if colons_in_prose(text)]
    assert not bad, "\n".join(bad)


def test_the_check_finds_a_colon_and_allows_a_location():
    tree = ast.parse(
        "raise TwigError(f'cannot read {path}: {exc}')\n"
        "raise TwigError(f'{path}:{lineno}: empty header')\n"
        "raise TwigError('tree file {} not found'.format(path))\n"
        "logger.warning(f'label {name}: no tips')\n"
        "logger.info('CMD: ' + ' '.join(argv))\n"
        "warnings.append('use this: or that')\n"
    )
    found = [text for _, text in messages(tree) if colons_in_prose(text)]
    assert found == ["cannot read {}: {}", "use this: or that"]
