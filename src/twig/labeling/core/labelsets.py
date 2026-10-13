#!/usr/bin/env python

"""Label a tree from one or more label sets.

A label set is a label name, a rule, and the selectors that choose its tips,
for example `Test:fitch dioecious.txt`. All the sets are applied to a tree and
their claims on branches are settled by a precedence, and then three switches
refine the result. `stems` can unlabel or rename the branches on which a trait
arose, `losses` can unlabel or rename reversions found by parsimony, and
`background` labels every branch that is still free.

A set can also be marked as nuisance. These are species that belong neither to
the test nor to the reference set, for any reason, such as a reversion, an
uncertain trait state or poor data. A nuisance label is never overwritten by
another set, and the tips of a nuisance set count as unknown in parsimony, so
the history of a trait is inferred from the other species. How a nuisance label
is treated by a program (set aside, removed, or counted as background) is chosen
later, when the labels are copied to gene trees.

Branches are identified by the idx of the node below them. The tree is never
changed here, the result says which label every branch gets.

This module does not import toytree.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path

from twig.utils import TwigError
from twig.labeling.core.constants import DEFAULT_NUISANCE_LABEL, PRECEDENCE_CHOICES
from twig.labeling.core.parsimony import to_fraction
from twig.labeling.core.rules import RULES, RuleInput, make_options, resolve_rule
from twig.labeling.core.sanitize import check_label_names
from twig.labeling.core.selectors import Selector

NUISANCE_RULE = "all"
SPEC_OPTIONS = ("ambiguous", "gain-cost", "loss-cost")


@dataclass(frozen=True)
class LabelSet:
    """A label name, the rule that labels branches, and the selectors of its tips.

    `values` are selector values, which are list files or regular expressions.
    The options override the defaults of the parsimony rule for this set.
    """
    name: str
    rule: str
    values: tuple[str, ...]
    ambiguous: str | None = None
    gain_cost: Fraction | None = None
    loss_cost: Fraction | None = None
    nuisance: bool = False


@dataclass(frozen=True)
class LabelingConfig:
    """Settings that apply to all label sets of a run."""
    select_on: str = "name"
    exact: bool = False
    invert: bool = False
    leaf_nodes: str = "label"
    stems: str = "keep"
    losses: str = "keep"
    background: str | None = None
    precedence: str = "specific"
    existing: str = "keep"
    unknown: tuple[str, ...] = ()
    ambiguous: str = "acctran"
    gain_cost: Fraction = Fraction(1)
    loss_cost: Fraction = Fraction(1)
    root_state: str = "free"


@dataclass
class SetReport:
    """What happened to one label set on one tree."""
    name: str
    rule: str
    nuisance: bool
    selected: int = 0
    claimed: int = 0
    won: int = 0
    cost: Fraction | None = None
    unused: list[str] = field(default_factory=list)


@dataclass
class LabelingResult:
    """The labels that a tree gets, and how they came about.

    `labels` has an entry for every branch, which is None if it is unlabeled.
    `source` says where a label came from (a set as `name:rule`, existing,
    stem, loss, ambiguous or background), and `kind` is the structural role of
    the branch, which is tip, stem, crown or internal. Stems are found before
    the stems switch is applied, so an unlabeled stem still has the kind stem.
    """
    labels: dict[int, str | None]
    source: dict[int, str]
    kind: dict[int, str]
    reports: list[SetReport]
    warnings: list[str] = field(default_factory=list)
    counts: dict[str, int] = field(default_factory=dict)


@dataclass(frozen=True)
class _Claim:
    label: str
    rank: tuple
    source: str
    set_index: int | None


###########################################################################
# reading label sets
###########################################################################

def _split_name_rule(text: str, default_rule: str) -> tuple[str, str]:
    name, _, rule = text.partition(":")
    return name, resolve_rule(rule) if rule else default_rule


def parse_label_args(items: Iterable[Sequence[str]], default_rule: str = "all") -> list[LabelSet]:
    """Return label sets from repeated `-L NAME[:RULE] SELECTOR...` arguments."""
    sets = []
    for item in items:
        if len(item) < 2:
            raise TwigError(f"a label needs a name and at least one selector, for example -L Test dioecious.txt (got {' '.join(item)!r})")
        name, rule = _split_name_rule(item[0], default_rule)
        sets.append(LabelSet(name, rule, tuple(item[1:])))
    return sets


def read_spec(path: str | Path, default_rule: str = "all") -> list[LabelSet]:
    """Return label sets from a tab-separated spec file.

    Each line has the columns label, rule, selectors and options, and the last
    column may be left out. The rule can be empty or `-` to use the default.
    The selectors are separated by whitespace. The options are separated by
    whitespace and are ambiguous=, gain-cost= and loss-cost=. A # starts a
    comment line.
    """
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise TwigError(f"cannot read label spec {path}: {exc}") from exc
    sets = []
    for lineno, raw in enumerate(text.splitlines(), start=1):
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        where = f"{path}:{lineno}"
        columns = raw.rstrip("\n").split("\t")
        if len(columns) < 3 or len(columns) > 4:
            raise TwigError(f"{where}: expected 3 or 4 tab-separated columns (label, rule, selectors, options), found {len(columns)}")
        name, rule, values = columns[0].strip(), columns[1].strip(), tuple(columns[2].split())
        if not name or not values:
            raise TwigError(f"{where}: the label and the selectors cannot be empty")
        options: dict[str, str] = {}
        for token in (columns[3].split() if len(columns) == 4 else []):
            key, equals, value = token.partition("=")
            if not equals or key not in SPEC_OPTIONS:
                raise TwigError(f"{where}: unknown option {token!r}, use {', '.join(i + '=' for i in SPEC_OPTIONS)}")
            options[key] = value
        try:
            resolved = resolve_rule(rule) if rule not in ("", "-") else default_rule
        except TwigError as exc:
            raise TwigError(f"{where}: {exc}") from None
        sets.append(LabelSet(
            name, resolved, values,
            ambiguous=options.get("ambiguous"),
            gain_cost=to_fraction(options["gain-cost"]) if "gain-cost" in options else None,
            loss_cost=to_fraction(options["loss-cost"]) if "loss-cost" in options else None,
        ))
    if not sets:
        raise TwigError(f"{path}: no label sets found")
    return sets


def nuisance_set(values: Sequence[str], label: str = DEFAULT_NUISANCE_LABEL) -> LabelSet:
    """Return the label set for nuisance species.

    The rule is `all`, so a clade whose species are all nuisance is labeled as
    a whole, and the tips of the set are unknown to parsimony.
    """
    return LabelSet(label, NUISANCE_RULE, tuple(values), nuisance=True)


###########################################################################
# selecting tips
###########################################################################

def select_tips(
    tree,
    values: Sequence[str],
    exact: bool = False,
    select_on: str = "name",
    species_of: Callable[[str], str] | None = None,
) -> tuple[frozenset[int], list[str]]:
    """Return the idx of the tips chosen by selector values, and unused entries.

    With select_on species the selectors are matched to the species of each tip
    instead of to its name, so one selection works on trees with several samples
    per species. The unused entries are listed names and patterns that matched
    no tip, which helps to find typos.
    """
    if select_on == "species" and species_of is None:
        raise TwigError("selecting on species needs a way to map names to species")
    selector = Selector(values, exact=exact)
    keys = {
        node.idx: (species_of(node.name) if select_on == "species" else node.name)
        for node in tree if node.is_leaf()
    }
    chosen = frozenset(idx for idx, key in keys.items() if selector.matches(key))
    return chosen, selector.unused(set(keys.values()))


###########################################################################
# applying label sets
###########################################################################

def _check_names(sets: Sequence[LabelSet], config: LabelingConfig) -> None:
    names = [s.name for s in sets]
    names += [i for i in (config.background,) if i]
    for switch in (config.stems, config.losses):
        if switch not in ("keep", "unlabel"):
            names.append(switch)
    if config.ambiguous not in ("acctran", "deltran", "unlabel"):
        names.append(config.ambiguous)
    names += [s.ambiguous for s in sets if s.ambiguous and s.ambiguous not in ("acctran", "deltran", "unlabel")]
    check_label_names(dict.fromkeys(names))


def check_settings(sets: Sequence[LabelSet], config: LabelingConfig) -> None:
    """Raise a TwigError if the label sets and settings cannot be applied to any tree.

    This needs no tree, so a command can call it before it reads anything.
    """
    if config.precedence not in PRECEDENCE_CHOICES:
        raise TwigError(f"unknown precedence {config.precedence!r}, choose from {list(PRECEDENCE_CHOICES)}")
    _check_names(sets, config)
    for item in sets:
        if item.rule == "tips" and config.leaf_nodes == "skip" and not item.nuisance:
            raise TwigError(f"the label {item.name!r} uses the tips rule, which labels nothing when tip labels are skipped")


def apply_label_sets(
    tree,
    sets: Sequence[LabelSet],
    config: LabelingConfig = LabelingConfig(),
    species_of: Callable[[str], str] | None = None,
) -> LabelingResult:
    """Return the labels that a tree gets from label sets.

    The steps are as follows. Each set selects tips and its rule labels
    branches. Where several sets (or the labels already in the tree) claim the
    same branch the precedence decides, and a nuisance label always wins.
    Tip labels can be skipped, stems and losses are unlabeled or renamed, and the
    background label goes to every branch that is still free, except stems,
    losses and ambiguous branches that were deliberately left unlabeled.
    """
    check_settings(sets, config)
    warnings: list[str] = []
    root = tree.treenode
    nonroot = [node for node in tree if node.up is not None]
    leaves = [node for node in tree if node.is_leaf()]

    # tips chosen by each set, and the tips whose state is unknown to parsimony
    selections, groups, reports = [], [], []
    for item in sets:
        report = SetReport(item.name, item.rule, item.nuisance)
        chosen, unused = select_tips(tree, item.values, config.exact, config.select_on, species_of)
        pieces = [select_tips(tree, [v], config.exact, config.select_on, species_of)[0] for v in item.values]
        if config.invert and not item.nuisance:
            everything = frozenset(node.idx for node in leaves)
            chosen = everything - chosen
            pieces = [chosen]
        report.selected, report.unused = len(chosen), unused
        if unused:
            warnings.append(f"label {item.name}: {len(unused)} listed names or patterns match no tip ({', '.join(unused[:5])}{' ...' if len(unused) > 5 else ''})")
        if not chosen:
            warnings.append(f"label {item.name}: the selectors choose no tips")
        selections.append(chosen)
        groups.append(tuple(pieces))
        reports.append(report)

    unknown, _ = select_tips(tree, config.unknown, config.exact, config.select_on, species_of) if config.unknown else (frozenset(), [])
    unknown = frozenset(unknown) | frozenset(i for item, sel in zip(sets, selections) if item.nuisance for i in sel)

    # every claim on a branch, with a rank where lower ranks win
    claims: dict[int, list[_Claim]] = defaultdict(list)
    if config.existing == "keep":
        for node in nonroot:
            label = getattr(node, "label", None)
            if label:
                claims[node.idx].append(_Claim(label, (1, 0), "existing", None))
    losses: set[int] = set()
    ambiguous: set[int] = set()
    for number, (item, chosen) in enumerate(zip(sets, selections)):
        options = make_options(
            gain=item.gain_cost if item.gain_cost is not None else config.gain_cost,
            loss=item.loss_cost if item.loss_cost is not None else config.loss_cost,
            ambiguous=item.ambiguous or config.ambiguous,
            root_state=config.root_state,
        )
        if item.rule == "fitch" and options.gain_cost != options.loss_cost and len(root.children) > 2:
            warnings.append(
                f"label {item.name}: the costs of gain and loss differ and the tree is unrooted, so the result depends on the position of the root"
            )
        out = RULES[item.rule](RuleInput(tree, item.name, chosen, unknown, groups[number], options))
        reports[number].claimed = len(out.labels)
        reports[number].cost = out.cost
        losses |= out.losses
        if options.ambiguous == "unlabel":
            ambiguous |= out.ambiguous
        for idx, label in out.labels.items():
            if item.nuisance:
                rank = (0, number)
            elif config.precedence == "specific":
                rank = (2 if idx in out.direct else 3, number)
            elif config.precedence == "last":
                rank = (2, -number)
            else:
                rank = (2, number)
            claims[idx].append(_Claim(label, rank, f"{item.name}:{item.rule}", number))

    # settle the claims
    labels: dict[int, str | None] = {node.idx: None for node in nonroot}
    source: dict[int, str] = {}
    winner_set: dict[int, int | None] = {}
    set_won: Counter = Counter()
    set_lost: Counter = Counter()      # branches lost to another regular set or to labels already in the tree
    for idx, found in claims.items():
        if config.precedence == "error":
            others = {c.label for c in found if c.rank[0] >= 2}
            if len(others) > 1:
                raise TwigError(f"branch {idx} is labeled {sorted(others)} by different sets, which --precedence error forbids")
        winner = min(found, key=lambda c: c.rank)
        labels[idx], source[idx], winner_set[idx] = winner.label, winner.source, winner.set_index
        if winner.set_index is not None:
            set_won[winner.set_index] += 1
        # a nuisance label is meant to win, so losing to it is not worth a warning
        if winner.rank[0] != 0:
            for claim in found:
                if claim.set_index is not None and claim.set_index != winner.set_index:
                    set_lost[claim.set_index] += 1
    for number, report in enumerate(reports):
        report.won = set_won[number]
        if set_lost[number]:
            warnings.append(
                f"label {report.name}: {set_lost[number]} of {report.claimed} branches went to another label because of --precedence {config.precedence}"
            )

    # branches labeled by a regular set, which means not by a nuisance set and not already in the tree
    def from_regular_set(idx: int) -> bool:
        number = winner_set.get(idx)
        return number is not None and not sets[number].nuisance

    # tips are not labeled if the user asked to skip them, except nuisance tips
    if config.leaf_nodes == "skip":
        for node in leaves:
            if node.up is not None and labels.get(node.idx) and from_regular_set(node.idx):
                labels[node.idx], source[node.idx] = None, ""

    # stems are labeled branches of the sets that do not continue a labeled parent
    regular = {idx for idx in source if labels.get(idx) and from_regular_set(idx)}
    stems = {
        idx for idx in regular
        if not tree[idx].is_leaf() and (tree[idx].up is root or labels.get(tree[idx].up.idx) != labels[idx])
    }

    kind: dict[int, str] = {}
    for node in nonroot:
        if node.is_leaf():
            kind[node.idx] = "tip"
        elif node.idx in stems:
            kind[node.idx] = "stem"
        elif labels.get(node.idx) and node.up is not root and labels.get(node.up.idx) == labels[node.idx]:
            kind[node.idx] = "crown"
        else:
            kind[node.idx] = "internal"

    held: set[int] = set()        # branches that the background label must not take
    counts: dict[str, int] = {"stems": len(stems)}
    if config.stems == "unlabel":
        for idx in stems:
            labels[idx], source[idx] = None, ""
        held |= stems
    elif config.stems != "keep":
        for idx in stems:
            labels[idx], source[idx] = config.stems, "stem"

    counts["losses"] = len(losses)
    for idx in losses:
        if config.losses == "unlabel":
            held.add(idx)
        elif config.losses != "keep" and labels.get(idx) is None:
            labels[idx], source[idx] = config.losses, "loss"

    counts["ambiguous"] = len(ambiguous)
    held |= ambiguous

    counts["background"] = 0
    if config.background:
        for node in nonroot:
            if labels[node.idx] is None and node.idx not in held:
                labels[node.idx], source[node.idx] = config.background, "background"
                counts["background"] += 1

    return LabelingResult(labels, source, kind, reports, warnings, counts)
