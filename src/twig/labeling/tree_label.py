#!/usr/bin/env python

"""Label the branches of trees from selections of tips.

Tips are selected by name or species, a rule turns each selection into labeled
branches, and the trees are written in one of several formats. The usual use is
to label the species tree once, and the labeled tree is then the input of
`tree-label-transfer`.

>>> twig tree-label -i sp.nwk -L Test dioecious.txt --rule fitch --stems unlabel --background Reference --show
>>> twig tree-label -i sp.nwk -L Test dioecious.txt --rule fitch --nuisance C --background Reference
>>> twig tree-label -i gene.nwk -L Test '^Amaranthus_hybridus' -F hyphy -o gene.hyphy.nwk
"""

from __future__ import annotations

import sys
from collections import Counter
from collections.abc import Callable
from contextlib import contextmanager
from dataclasses import replace
from pathlib import Path

from loguru import logger

from twig.utils import TwigError, set_log_level
from twig.labeling.core import batch, display, rooting, treeio, writers
from twig.labeling.core.constants import (
    DEFAULT_NORMALIZE,
    EXISTING_CHOICES,
    LEAF_NODE_CHOICES,
    PRECEDENCE_CHOICES,
    SELECT_ON_CHOICES,
)
from twig.labeling.core.display import count, shown
from twig.labeling.core.labelsets import (
    LabelingConfig,
    LabelingResult,
    LabelSet,
    apply_label_sets,
    check_settings,
    nuisance_set,
    parse_label_args,
    read_spec,
    unused_settings,
)
from twig.labeling.core.names import NameHit, NameMapper, read_vocab_list
from twig.labeling.core.parsimony import to_fraction
from twig.labeling.core.rules import resolve_rule
from twig.labeling.core.sanitize import NameTable, label_name_problem
from twig.labeling.core.validate import errors, warnings


# how many tip names a warning lists
SHOWN = 8

# the Taxon Mapping options that find species, with the words used for them in messages
MAPPING_RULES = (
    ("imap", "-I", "the map file"),
    ("delim", "-d", "the delimiter"),
    ("name_regex", "-x", "the regex"),
    ("vocab", "-V", "the species list"),
)


def _joined(items: list[str], word: str = "and") -> str:
    """Return items joined by commas and a last word, such as 'a, b and c'."""
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + f" {word} {items[-1]}"


@contextmanager
def at_tree(locus: str):
    """Add the name of the tree to the message of a TwigError raised inside the block."""
    try:
        yield
    except TwigError as exc:
        raise TwigError(f"{locus}: {exc}") from None


def read_paml_codes(items: list[str] | None) -> dict[str, str] | None:
    """Return the names of PAML label numbers from entries such as 1=Test.

    The names are label names, so they are checked like the other label names
    before any tree is read.
    """
    codes = {}
    for item in items or ():
        number, equals, name = item.partition("=")
        if not equals or not number.strip().isdigit() or not name.strip():
            raise TwigError(f"--paml-codes expects entries such as 1=Test, got {item!r}")
        codes[str(int(number))] = name.strip()
    for number, name in codes.items():
        problem = label_name_problem(name)
        if problem:
            raise TwigError(f"--paml-codes gives {number} the name {name!r}, which {problem}")
    return codes or None


def input_label_warnings(tree, syntax: str | None, seen: set[str]) -> list[str]:
    """Return warnings for the labels read from a tree that are not valid label names.

    Such labels are kept and written as they are, but HyPhy cannot use them, so
    each one is named once per run. PAML numbers without a name are named in
    one warning, with a hint to name them with --paml-codes. `seen` holds the
    labels already named.
    """
    bad = {}
    for label in sorted(set(treeio.labeled_branches(tree).values()) - seen):
        problem = label_name_problem(label)
        if problem:
            bad[label] = problem
    seen.update(bad)
    if syntax == "paml" and bad:
        numbers = sorted(bad, key=lambda i: (len(i), i))
        if len(numbers) == 1:
            what, it = f"the PAML label number {numbers[0]} has no name", "it"
        else:
            what, it = f"the PAML label numbers {', '.join(numbers[:-1])} and {numbers[-1]} have no names", "them"
        return [
            f"{what}, and HyPhy cannot use a number as a label. Name {it} with --paml-codes, such as --paml-codes {numbers[0]}=Test"
        ]
    return [f"the label {label!r} read from the tree {problem}, so HyPhy cannot use it" for label, problem in bad.items()]


def mapping_warnings(args, select_on: str) -> list[str]:
    """Return warnings for the Taxon Mapping options that have no effect.

    The options act only with --select-on species. -N and --accept-inferred
    act only when names are matched to a species list given with -V, and -di
    and -dj only with -d.
    """
    rules = [flag for attr, flag, _ in MAPPING_RULES if getattr(args, attr)]
    changed = {
        "-di": list(args.delim_idxs) != [0],
        "-dj": args.delim_join is not None,
        "-N": set(args.normalize) != set(DEFAULT_NORMALIZE),
        "--accept-inferred": args.accept_inferred,
        "-U": args.unmapped != "error",
    }
    changed = [flag for flag, given in changed.items() if given]
    found = []
    if select_on != "species":
        found.append((rules + changed, "without --select-on species"))
    else:
        if not args.vocab:
            found.append(([i for i in changed if i in ("-N", "--accept-inferred")], "without a species list given with -V"))
        if not args.delim:
            found.append(([i for i in changed if i in ("-di", "-dj")], "without -d"))
    return [
        f"{_joined(flags)} {'has' if len(flags) == 1 else 'have'} no effect {why}"
        for flags, why in found if flags
    ]


def unresolved_reason(hit: NameHit, rules: list[str]) -> str | None:
    """Return why a name has no species that tree-label can use, or None if it has one.

    `rules` lists the mapping rules that were given, such as 'the map file'.
    Without any, every name is its own species.
    """
    if hit.status == "ambiguous":
        return f"it matches {_joined([c[0] for c in hit.candidates])} equally well"
    if hit.status == "inferred":
        return f"only an abbreviated genus matches it to {hit.species}, which --accept-inferred accepts"
    if hit.status == "unmapped":
        near = ", ".join(c[0] for c in hit.candidates)
        return "it matches no species of the species list" + (f", the closest are {near}" if near else "")
    if hit.method == "self" and rules:
        return f"it is not resolved by {_joined(rules, 'or')}"
    return None


def species_lookup(args, mapper: NameMapper, names: list[str]) -> Callable[[str], str | None]:
    """Return the function that gives the species of a tip name, after checking every name.

    A name that the mapping options do not resolve to one species stops the
    run with a full list under -U error. Under -U ignore it is named in a
    warning and its species is None, so that no selector chooses the tip.
    """
    rules = [words for attr, _, words in MAPPING_RULES if getattr(args, attr)]
    species: dict[str, str | None] = {}
    problems = []
    for name in dict.fromkeys(names):
        hit = mapper.resolve(name)
        reason = unresolved_reason(hit, rules)
        species[name] = None if reason else hit.species
        if reason:
            problems.append(f"{name} ({reason})")
    it = "it" if len(problems) == 1 else "them"
    if problems and args.unmapped == "error":
        raise TwigError(
            f"{count(len(problems), 'tip name')} cannot be resolved to a species, so -U error stops the run\n"
            + "\n".join(f"  {i}" for i in problems)
            + f"\nAdd {it} to a map file given with -I, check the names with twig name-map, or give -U ignore to leave {it} without a species"
        )
    if problems:
        logger.warning(
            f"{count(len(problems), 'tip name')} cannot be resolved to a species, so -U ignore leaves {it} "
            f"without one and no selector chooses {it}. {shown(problems, SHOWN)}"
        )
    return species.__getitem__


def build_label_sets(args) -> list[LabelSet]:
    """Return the label sets from a spec file, from -L, and from --nuisance, in that order."""
    default_rule = resolve_rule(args.rule)
    sets: list[LabelSet] = []
    if args.spec:
        sets += read_spec(args.spec, default_rule)
    if args.label:
        sets += parse_label_args(args.label, default_rule)
    if args.nuisance:
        sets.append(nuisance_set(args.nuisance, args.nuisance_label))
    return sets


def build_config(args) -> LabelingConfig:
    """Return the settings that apply to all label sets from the parsed arguments."""
    for name, value, choices in (
        ("--select-on", args.select_on, SELECT_ON_CHOICES),
        ("--leaf-nodes", args.leaf_nodes, LEAF_NODE_CHOICES),
        ("--precedence", args.precedence, PRECEDENCE_CHOICES),
        ("--existing", args.existing, EXISTING_CHOICES),
    ):
        if value not in choices:
            raise TwigError(f"unknown value {value!r} for {name}, choose from {list(choices)}")
    return LabelingConfig(
        select_on=args.select_on,
        exact=args.exact,
        invert=args.invert,
        leaf_nodes=args.leaf_nodes,
        stems=args.stems,
        losses=args.losses,
        background=args.background,
        precedence=args.precedence,
        existing=args.existing,
        unknown=tuple(args.unknown or ()),
        ambiguous=args.ambiguous,
        gain_cost=to_fraction(args.gain_cost),
        loss_cost=to_fraction(args.loss_cost),
        root_state=args.root_state,
    )


def check_outputs(args) -> None:
    """Raise a TwigError if an output exists, is listed twice, or is the input."""
    paths = [p for p in (args.out, args.branch_table, args.names_table) if p]
    for path in paths:
        batch.check_writable(path, args.force)
    resolved = [Path(p).resolve() for p in paths]
    if len(set(resolved)) != len(resolved):
        raise TwigError("--out, --branch-table and --names-table must be different files")
    if Path(args.input).resolve() in resolved:
        raise TwigError("an output file is the same as the input file")


def describe(result: LabelingResult, locus: str) -> list[str]:
    """Return log lines that say what each label set did to a tree."""
    lines = []
    for report in result.reports:
        kind = f"{report.rule}, nuisance" if report.nuisance else report.rule
        line = (
            f"{locus}: label {report.name} ({kind}) selected {count(report.selected, 'tip')}, "
            f"covers {count(report.claimed, 'branch')} and keeps {report.kept}"
        )
        if report.cost is not None:
            line += f", parsimony cost {report.cost}"
        lines.append(line)
    return lines


def summarize(labels: dict[int, str | None], result: LabelingResult, locus: str) -> str:
    """Return one log line that counts the labeled branches of a tree."""
    counts = Counter(i for i in labels.values() if i)
    unlabeled = sum(1 for i in labels.values() if not i)
    parts = ", ".join(f"{label} {k}" for label, k in counts.most_common()) or "none"
    line = f"{locus}: {count(len(labels), 'branch')}, labeled {parts}, unlabeled {unlabeled}"
    stems = result.counts.get("stems", 0)
    if stems:
        line += f", {count(stems, 'stem')}"
    return line


def run_tree_label(args) -> None:
    """Run tree-label subcommand."""
    set_log_level(args.log_level)
    logger.info("CMD: " + " ".join(sys.argv))

    # settings are checked before any tree is read
    sets = build_label_sets(args)
    config = build_config(args)
    check_settings(sets, config)
    for message in unused_settings(sets, config):
        logger.warning(message)
    writer = writers.get_writer(args.format)
    rooting.check_root_edge_policy(args.root_edges)
    check_outputs(args)
    paml_codes = read_paml_codes(args.paml_codes)
    table_writer = writers.get_writer("table")

    for message in mapping_warnings(args, config.select_on):
        logger.warning(message)
    mapper = None
    if config.select_on == "species":
        vocabulary = read_vocab_list(args.vocab) if args.vocab else None
        mapper = NameMapper.from_args(args, vocabulary=vocabulary)
        if not any(getattr(args, attr) for attr, _, _ in MAPPING_RULES):
            logger.warning("--select-on species was given without a map file, delimiter, regex or species list, so every tip is its own species")

    records = list(batch.tree_records([args.input]))
    logger.info(f"read {count(len(records), 'tree')} from {args.input}")

    # every tip name is checked before any tree is labeled, and a tip without a species is shown by its name
    species_of = species_name = None
    if mapper is not None:
        tip_names = []
        for locus, newick, _ in records:
            with at_tree(locus):
                tip_names += treeio.tip_names(newick)
        species_of = species_lookup(args, mapper, tip_names)
        species_name = lambda name: species_of(name) or name  # noqa: E731

    if not sets:
        what = "are cleared" if config.existing == "clear" else "are kept"
        logger.info(f"no label sets were given, so the labels already in the tree {what}")

    # one table of written names for all trees, so that a sample has the same name everywhere
    names = None
    if writer.sanitizer is not None:
        names = NameTable(writer.sanitizer, suffix=args.sanitize_suffix)
        for _, newick, _ in records:
            names.update(treeio.tip_names(newick))
        names.build()
        changed = names.changed()
        if changed:
            logger.info(f"{count(len(changed), 'tip name')} changed to be safe for {writer.sanitizer}" + (
                f", see {args.names_table}" if args.names_table else ", use --names-table to save the table"))
    elif args.names_table:
        logger.warning("--names-table is only used with -F hyphy")

    pieces: list[str] = []
    table_pieces: list[str] = []
    seen_labels: set[str] = set()
    bad_input_labels: set[str] = set()
    for number, (locus, newick, path) in enumerate(records):
        parsed = treeio.parse_newick(newick, "auto", paml_codes, source=f"{path} ({locus})")
        tree = parsed.tree
        if parsed.root_label:
            logger.warning(f"{locus}: ignoring the label {parsed.root_label!r} on the root, which has no branch")
        if config.existing == "keep":
            for message in input_label_warnings(tree, parsed.syntax, bad_input_labels):
                logger.warning(f"{locus}: {message}")

        outgroup_side = None
        if args.outgroups:
            try:
                with at_tree(locus):
                    rooted = rooting.reroot(tree, args.outgroups, args.root_edges, species_of)
            except treeio.toytree_errors() as exc:
                raise TwigError(f"{locus}: cannot reroot the tree ({exc})") from exc
            tree, outgroup_side = rooted.tree, rooted.outgroup_side
            logger.info(f"{locus}: {rooted.description}")
            for note in rooted.notes:
                logger.warning(f"{locus}: {note}")

        with at_tree(locus):
            result = apply_label_sets(tree, sets, config, species_of)
        for line in describe(result, locus):
            logger.info(line)
        for message in result.warnings:
            logger.warning(f"{locus}: {message}")
        labeled = treeio.set_labels(tree, result.labels, inplace=True)
        seen_labels.update(i for i in result.labels.values() if i)
        if not any(result.labels.values()):
            logger.warning(f"{locus}: no branch is labeled")
        logger.info(summarize(result.labels, result, locus))

        ctx = writers.WriteContext(
            locus=locus,
            dists=not args.no_dists,
            root_edges=args.root_edges,
            keep_root=args.keep_root,
            internal_names=args.internal_names,
            names=names,
            outgroup_side=outgroup_side,
            kinds=result.kind,
            sources=result.source,
            species_of=species_name,
            header=number == 0,
        )
        issues = writer.check(labeled, ctx)
        for issue in warnings(issues):
            logger.warning(str(issue))
        if errors(issues):
            raise TwigError("cannot write the tree\n" + "\n".join(f"  {i}" for i in errors(issues)))

        with at_tree(locus):
            pieces.append(writer.render(labeled, ctx))
            plan = writer.plan(labeled, ctx) if isinstance(writer, writers.TreeWriter) else None
            if args.branch_table:
                table_pieces.append(table_writer.render(labeled, replace(ctx, plan=plan)))
        if args.show:
            # the listing describes the tree as it is written, like the table
            heading = [f"tree {locus}"] if len(records) > 1 else []
            sys.stderr.write("\n".join(heading + display.branch_listing(
                labeled, result.labels, result.kind, result.source, species_name, merge=plan.merge if plan else None)) + "\n")

    missing = rooting.missing_policy_label(args.root_edges, seen_labels)
    if missing:
        logger.warning(missing)

    text = "\n".join(pieces) + "\n"
    if args.branch_table:
        batch.atomic_write(args.branch_table, "\n".join(table_pieces) + "\n", force=args.force)
        logger.info(f"wrote the table of branches to {args.branch_table}")
    if args.names_table and names is not None:
        names.write_tsv(args.names_table, force=args.force)
        logger.info(f"wrote the table of tip names to {args.names_table}")
    if args.out:
        batch.atomic_write(args.out, text, force=args.force)
        logger.info(f"wrote {count(len(pieces), 'tree')} to {args.out}")
    else:
        batch.write_stdout(text)


def main():
    """Run tree-label as a script, like the other twig run modules."""
    from .parsers import get_parser_tree_label
    args = get_parser_tree_label().parse_args()
    run_tree_label(args)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.warning("interrupted by user")
        sys.exit(1)
    except TwigError as exc:
        logger.error(exc)
        sys.exit(1)
    except Exception as exc:
        logger.error(exc)
        raise
