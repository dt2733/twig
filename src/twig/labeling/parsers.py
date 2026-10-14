#!/usr/bin/env python

"""Argument parsers for the twig tree labeling commands.

twig builds the parsers of every command each time it starts, so this module
imports only argparse, pathlib, textwrap and the dependency-free constants.
toytree is imported only by the run modules and the library.
"""

from __future__ import annotations

from argparse import ArgumentParser
from pathlib import Path
from textwrap import dedent

from twig.labeling.core.constants import (
    DEFAULT_NORMALIZE,
    DEFAULT_NUISANCE_LABEL,
    EXISTING_CHOICES,
    FORMAT_NAMES,
    INTERNAL_NAME_CHOICES,
    LEAF_NODE_CHOICES,
    NORMALIZE_CHOICES,
    PRECEDENCE_CHOICES,
    ROOT_EDGE_POLICIES,
    ROOT_STATE_CHOICES,
    RULE_NAMES,
    SELECT_ON_CHOICES,
    SHOW_CHOICES,
    UNMAPPED_POLICIES,
)

FRAME = "-" * 67


def _formatter(prog: str):
    """Return the help formatter used by the newer twig commands."""
    # imported here because twig.cli.subcommands imports this module at its end
    from twig.cli.subcommands import SingleMetavarHelpFormatter
    return SingleMetavarHelpFormatter(prog, width=220, max_help_position=56)


def _create(parent, kwargs: dict) -> ArgumentParser:
    """Return a subparser of parent, or a stand-alone parser if parent is None."""
    if parent:
        kwargs["name"] = kwargs.pop("prog")
        return parent.add_parser(**kwargs)
    kwargs.pop("help")
    return ArgumentParser(**kwargs)


def add_taxon_mapping_args(parser: ArgumentParser, unmapped: bool = True) -> None:
    """Add the group of options that map sample names to species.

    The group is the same in every labeling command. The option to choose what
    happens to unresolved names is added only by commands that act on them.
    """
    group = parser.add_argument_group("Taxon Mapping")
    group.add_argument("-I", "--imap", type=Path, metavar="path", help="map file rows <sample> <species> [specimen]; ~regex rows allowed")
    exclusive = group.add_mutually_exclusive_group()
    exclusive.add_argument("-d", "--delim", type=str, metavar="str", help="delimiter used to split sample names")
    exclusive.add_argument("-x", "--name-regex", type=str, metavar="str", help="regex with (?P<species>...) and optional (?P<specimen>...)")
    group.add_argument("-di", "--delim-idxs", type=int, metavar="int", nargs="+", default=[0], help="token indices kept after splitting [0]")
    group.add_argument("-dj", "--delim-join", type=str, metavar="str", help="join string for kept tokens [same as --delim]")
    group.add_argument("-N", "--normalize", type=str, metavar="str", nargs="+", choices=NORMALIZE_CHOICES, default=list(DEFAULT_NORMALIZE), help="normalization steps used for matching, or none [" + " ".join(DEFAULT_NORMALIZE) + "]")
    group.add_argument("--accept-inferred", action="store_true", help="accept genus-abbreviation matches")
    if unmapped:
        group.add_argument("-U", "--unmapped", type=str, metavar="str", choices=UNMAPPED_POLICIES, default="error", help="what to do with names that cannot be resolved (error, skip, drop, ignore) [%(default)s]")


def add_writer_args(parser: ArgumentParser) -> None:
    """Add the group of options that control how trees are written.

    The group is the same in every labeling command that writes trees.
    """
    group = parser.add_argument_group("Output Format")
    group.add_argument("--root-edges", type=str, metavar="str", default="unlabel", help="HyPhy and PAML read a rooted tree as unrooted and see the two branches at its root as one branch. This sets the label of that branch (" + ", ".join(ROOT_EDGE_POLICIES) + "). unlabel gives no label if the two branches differ, labeled uses the only label there is, outgroup and ingroup use the label of that side, error stops, and a label name wins if either branch has it [%(default)s]")
    group.add_argument("--keep-root", action="store_true", help="write the tree rooted as it is, without joining the two branches at the root")
    group.add_argument("--internal-names", type=str, metavar="str", choices=INTERNAL_NAME_CHOICES, default="none", help="names written after internal nodes (" + ", ".join(INTERNAL_NAME_CHOICES) + "). none writes labels only, auto writes twn1, twn2 and so on, support writes the support values, and keep writes the names the tree had [%(default)s]")
    group.add_argument("--no-dists", action="store_true", help="write trees without branch lengths")
    group.add_argument("--sanitize-suffix", action="store_true", help="for HyPhy, add _2, _3 to tip names that become equal when they are made safe")
    group.add_argument("--names-table", type=Path, metavar="path", help="for HyPhy, write the table of original and written tip names")


def get_parser_name_map(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for the name-map tool."""
    KWARGS = dict(
        prog="name-map",
        usage="name-map -i TREES/FASTA [...] -s SPTREE [options]",
        help="resolve tip and sequence names to species and freeze the mapping",
        formatter_class=_formatter,
        description=dedent(f"""
            {FRAME}
            | name-map resolves tip and sequence names to species names
            {FRAME}
            | Reads gene trees and/or FASTA files and resolves every tip and
            | sequence name to a species of the species tree (or of a species
            | list). Names are matched exactly, after normalization, by whole
            | tokens anywhere in the name, and by genus abbreviation. Matches
            | by abbreviation are reported as inferred and are not accepted by
            | other twig commands unless they are frozen into a map file with
            | --freeze or allowed with --accept-inferred. The table of names is
            | written to stdout and a summary to stderr.
            {FRAME}
        """),
        epilog=dedent("""
            Examples
            --------
            # check every name in a set of gene trees and alignments
            $ twig name-map -i 'trees/*.treefile' 'aln/*.nt.fa' -s sp.nwk

            # show only the names that need attention
            $ twig name-map -i 'trees/*.treefile' -s sp.nwk --show problems

            # parse the species from names such as OG0001-AmaranthusPumilus-s1
            $ twig name-map -i 'trees/*.treefile' -s sp.nwk -d '-' -di 1

            # save the resolved mapping for reuse with other twig commands
            $ twig name-map -i 'trees/*.treefile' 'aln/*.nt.fa' -s sp.nwk --freeze names.tsv

            Example map file (columns are sample, species, and optionally specimen)
            -----------------------------------------------------------------------
            A_pumilus_sample1           Amaranthus_pumilus    s1
            AmaranthusPumilus-Sample2   Amaranthus_pumilus    s2
            ~^Ahyb[_-]                  Amaranthus_hybridus
        """),
    )
    parser = _create(parser, KWARGS)

    io = parser.add_argument_group("Input/Output")
    io.add_argument("-i", "--input", type=Path, metavar="path", nargs="+", required=True, help="trees and/or FASTA files, directories or quoted globs")
    io.add_argument("-s", "--sptree", type=Path, metavar="path", help="species tree whose tips form the vocabulary")
    io.add_argument("-V", "--vocab", type=Path, metavar="path", help="file listing species names, one per line")
    io.add_argument("-o", "--out", type=Path, metavar="path", help="write the mapping table to this file [stdout]")
    io.add_argument("--freeze", type=Path, metavar="path", help="write the resolved map for reuse with -I")

    add_taxon_mapping_args(parser, unmapped=False)

    report = parser.add_argument_group("Reporting")
    report.add_argument("--show", type=str, metavar="str", choices=SHOW_CHOICES, default="all", help="rows of the table to print (all, problems) [%(default)s]")
    report.add_argument("--strict", action="store_true", help="error if any name is unmapped, ambiguous or inferred, or leaves alphabetic leftover tokens")

    runtime = parser.add_argument_group("Runtime")
    runtime.add_argument("-f", "--force", action="store_true", help="overwrite existing output files")
    runtime.add_argument("-l", "--log-level", type=str, metavar="level", default="INFO", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    return parser


def _get_planned_parser(parent, prog: str, summary: str, milestone: str) -> ArgumentParser:
    """Return a parser for a command that is designed but not available yet."""
    KWARGS = dict(
        prog=prog,
        usage=f"{prog} [options]",
        help=f"(planned for {milestone}) {summary}",
        formatter_class=_formatter,
        description=dedent(f"""
            {FRAME}
            | {prog} is planned for milestone {milestone} and is not available yet
            {FRAME}
            | {summary}
            | The design is in docs/labeling/design.md.
            {FRAME}
        """),
    )
    parser = _create(parent, KWARGS)
    parser.add_argument("-l", "--log-level", type=str, metavar="level", default="INFO", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    return parser


TREE_LABEL_EXAMPLES = r"""Examples
--------
The examples use a small tree and a list of the species that have a trait. To
follow along, make them first. The names on the internal branches (BC, ABC, X,
DE, Y and Z) are only there so that the results can be read. A name sits on the
branch above its group, so BC is the branch above B and C, ABC the branch above
A, B and C, and Z the branch above everything except O.

$ echo '(O,((M1,(A,(B,C)BC)ABC)X,(M2,(D,E)DE)Y)Z);' > sp.nwk
$ printf 'A\nB\nD\nE\n' > dioecious.txt

A label is written in braces after the name of its branch, so A{Test} is a tip
and )DE{Test} an internal branch. The option --internal-names keep writes the
names of the internal branches, which are left out otherwise.

Which branches get the label (--rule)
-------------------------------------
The chosen species are A, B, D and E and the label is Test. The rule decides how
far the label spreads from the chosen tips.

# tips labels only the chosen tips
$ twig tree-label -i sp.nwk -L Test dioecious.txt --rule tips --internal-names keep
(O,((M1,(A{Test},(B{Test},C)BC)ABC)X,(M2,(D{Test},E{Test})DE)Y)Z);

# all labels a branch too when every tip below it is chosen, which is DE (the default rule)
$ twig tree-label -i sp.nwk -L Test dioecious.txt --rule all --internal-names keep
(O,((M1,(A{Test},(B{Test},C)BC)ABC)X,(M2,(D{Test},E{Test})DE{Test})Y)Z);

# crown is like all but leaves out the branch at the base of each group, so DE is not labeled
$ twig tree-label -i sp.nwk -L Test dioecious.txt --rule crown --internal-names keep
(O,((M1,(A{Test},(B{Test},C)BC)ABC)X,(M2,(D{Test},E{Test})DE)Y)Z);

# some labels every branch that has at least one chosen tip below it
$ twig tree-label -i sp.nwk -L Test dioecious.txt --rule some --internal-names keep
(O,((M1,(A{Test},(B{Test},C)BC{Test})ABC{Test})X{Test},(M2,(D{Test},E{Test})DE{Test})Y{Test})Z{Test});

# clade labels the branch above the common ancestor of the chosen tips and everything below it
# (the rule goes after the label name, joined to it with a colon)
$ twig tree-label -i sp.nwk -L Test:clade '^[AB]$' --internal-names keep
(O,((M1,(A{Test},(B{Test},C{Test})BC{Test})ABC{Test})X,(M2,(D,E)DE)Y)Z);

# stem labels only the branch above the common ancestor
$ twig tree-label -i sp.nwk -L Test:stem '^[AB]$' --internal-names keep
(O,((M1,(A,(B,C)BC)ABC{Test})X,(M2,(D,E)DE)Y)Z);

What fitch means
----------------
Parsimony means choosing the explanation that needs the fewest changes, and
Fitch is Walter Fitch, who described this method in 1971. Given the species that
have the trait, fitch finds the fewest gains and losses of the trait on the tree
and labels the branches where the trait is present. The rule parsimony means the
same. Here the trait arose once on ABC, was lost on C and arose once on DE, which
is 3 changes.

# fitch labels the branches where parsimony says the trait is present
$ twig tree-label -i sp.nwk -L Test dioecious.txt --rule fitch --internal-names keep
(O,((M1,(A{Test},(B{Test},C)BC{Test})ABC{Test})X,(M2,(D{Test},E{Test})DE{Test})Y)Z);

Two histories can be equally good. Gaining the trait separately on A and on B
also costs 3 changes. The option --ambiguous says which history to use.

# acctran (the default) puts changes nearer the root, and deltran nearer the tips
$ twig tree-label -i sp.nwk -L Test dioecious.txt --rule fitch --ambiguous deltran --internal-names keep
(O,((M1,(A{Test},(B{Test},C)BC)ABC)X,(M2,(D{Test},E{Test})DE{Test})Y)Z);

# unlabel leaves the branches that differ without a label, and a name gives them a label of their own
$ twig tree-label -i sp.nwk -L Test dioecious.txt --rule fitch --ambiguous Unsure --internal-names keep
(O,((M1,(A{Test},(B{Test},C)BC{Unsure})ABC{Unsure})X,(M2,(D{Test},E{Test})DE{Test})Y)Z);

# making a loss cost more than a gain (--gain-cost and --loss-cost, both 1) avoids the loss on C
$ twig tree-label -i sp.nwk -L Test dioecious.txt --rule fitch --loss-cost 3 --internal-names keep
(O,((M1,(A{Test},(B{Test},C)BC)ABC)X,(M2,(D{Test},E{Test})DE{Test})Y)Z);

# making gains cost more (--gain-cost) and fixing the trait as absent at the root (--root-state 0)
# puts one gain at the base of the whole group (Z), with losses on M1, M2 and C
$ twig tree-label -i sp.nwk -L Test dioecious.txt --rule fitch --gain-cost 3 --root-state 0 --internal-names keep
(O,((M1,(A{Test},(B{Test},C)BC{Test})ABC{Test})X{Test},(M2,(D{Test},E{Test})DE{Test})Y{Test})Z{Test});

# --unknown marks tips whose state is not known, so parsimony may give them either state
$ twig tree-label -i sp.nwk -L Test dioecious.txt --rule fitch --unknown C --internal-names keep
(O,((M1,(A{Test},(B{Test},C{Test})BC{Test})ABC{Test})X,(M2,(D{Test},E{Test})DE{Test})Y)Z);

Stems, losses and the background
--------------------------------
A stem is the branch on which the trait arose, here ABC and DE. The change
happened somewhere along it, so a stem is often best left out of both groups. A
loss is a branch where the trait is lost again, here C. These switches work with
every rule.

# --stems unlabel leaves the stems without a label
$ twig tree-label -i sp.nwk -L Test dioecious.txt --rule fitch --stems unlabel --internal-names keep
(O,((M1,(A{Test},(B{Test},C)BC{Test})ABC)X,(M2,(D{Test},E{Test})DE)Y)Z);

# a name gives the stems a label of their own, so a later step can treat them as Test or leave them out
$ twig tree-label -i sp.nwk -L Test dioecious.txt --rule fitch --stems Stem --internal-names keep
(O,((M1,(A{Test},(B{Test},C)BC{Test})ABC{Stem})X,(M2,(D{Test},E{Test})DE{Stem})Y)Z);

# with crown the stem is the branch that the rule leaves out at the base of each group, here DE
$ twig tree-label -i sp.nwk -L Test dioecious.txt --rule crown --stems Stem --internal-names keep
(O,((M1,(A{Test},(B{Test},C)BC)ABC)X,(M2,(D{Test},E{Test})DE{Stem})Y)Z);

# --losses gives the branches where the trait is lost a label of their own
$ twig tree-label -i sp.nwk -L Test dioecious.txt --rule fitch --losses Rev --internal-names keep
(O,((M1,(A{Test},(B{Test},C{Rev})BC{Test})ABC{Test})X,(M2,(D{Test},E{Test})DE{Test})Y)Z);

# --background labels every branch that is still free, which makes the reference branches explicit
$ twig tree-label -i sp.nwk -L Test dioecious.txt --rule fitch --stems unlabel --losses Rev --background Reference --internal-names keep
(O{Reference},((M1{Reference},(A{Test},(B{Test},C{Rev})BC{Test})ABC)X{Reference},(M2{Reference},(D{Test},E{Test})DE)Y{Reference})Z{Reference});

# --leaf-nodes skip removes the labels of the tips after the rule has run
$ twig tree-label -i sp.nwk -L Test dioecious.txt --rule all --leaf-nodes skip --internal-names keep
(O,((M1,(A,(B,C)BC)ABC)X,(M2,(D,E)DE{Test})Y)Z);

Nuisance species
----------------
Some species belong in neither group, for any reason, such as a reversion, an
uncertain trait or poor data. A nuisance species gets a label of its own that no
other label overwrites, and parsimony ignores its state, so the history is worked
out from the other species.

# set C aside, which also changes the history (C is no longer a loss)
$ twig tree-label -i sp.nwk -L Test dioecious.txt --rule fitch --nuisance C --internal-names keep
(O,((M1,(A{Test},(B{Test},C{Nuisance})BC{Test})ABC{Test})X,(M2,(D{Test},E{Test})DE{Test})Y)Z);

# the name of the label can be chosen
$ twig tree-label -i sp.nwk -L Test dioecious.txt --rule fitch --nuisance C --nuisance-label Rev --stems unlabel --background Reference --internal-names keep
(O{Reference},((M1{Reference},(A{Test},(B{Test},C{Rev})BC{Test})ABC)X{Reference},(M2{Reference},(D{Test},E{Test})DE)Y{Reference})Z{Reference});

Choosing the tips
-----------------
A selector is a list file with one name per line (a line that starts with ~ is a
regular expression), or a regular expression that is matched anywhere in the name.

# the regular expression M matches M1 and M2
$ twig tree-label -i sp.nwk -L Mono M --rule tips
(O,((M1{Mono},(A,(B,C))),(M2{Mono},(D,E))));

# --exact means exactly this name, so M1 matches M1 only
$ twig tree-label -i sp.nwk -L Mono M1 --rule tips --exact
(O,((M1{Mono},(A,(B,C))),(M2,(D,E))));

# --invert labels the tips that are not chosen
$ twig tree-label -i sp.nwk -L Mono dioecious.txt --rule tips --invert
(O{Mono},((M1{Mono},(A,(B,C{Mono}))),(M2{Mono},(D,E))));

When the tips are samples of species, --select-on species matches the species
instead of the names. The species come from the Taxon Mapping options, here -d _
takes the part of each name before the first underscore.

$ echo '((A_1,A_2),(B_1,C_1));' > gene.nwk
$ twig tree-label -i gene.nwk -L Test A --select-on species -d _ --rule tips
((A_1{Test},A_2{Test}),(B_1,C_1));

# -x finds the species with a regular expression that has a group named species
$ twig tree-label -i gene.nwk -L Test A --select-on species -x '^(?P<species>[^_]+)_' --rule tips
((A_1{Test},A_2{Test}),(B_1,C_1));

# -I reads a map file with a sample and its species on each line
$ printf 'A_1 A\nA_2 A\nB_1 B\nC_1 C\n' > map.txt
$ twig tree-label -i gene.nwk -L Test A --select-on species -I map.txt --rule tips
((A_1{Test},A_2{Test}),(B_1,C_1));

Several labels
--------------
# repeat -L for more labels, each with its own rule
$ twig tree-label -i sp.nwk -L Out:tips '^O$' -L Test:fitch dioecious.txt --internal-names keep
(O{Out},((M1,(A{Test},(B{Test},C)BC{Test})ABC{Test})X,(M2,(D{Test},E{Test})DE{Test})Y)Z);

# two origins of a trait can have a label each
$ twig tree-label -i sp.nwk -L Test1:fitch '^[AB]$' -L Test2:fitch '^[DE]$' --stems unlabel --background Reference --internal-names keep
(O{Reference},((M1{Reference},(A{Test1},(B{Test1},C{Reference})BC{Test1})ABC)X{Reference},(M2{Reference},(D{Test2},E{Test2})DE)Y{Reference})Z{Reference});

# or list the labels in a file with the columns label, rule and selectors
$ printf 'Out\ttips\t^O$\nTest\tfitch\tdioecious.txt\n' > labels.tsv
$ twig tree-label -i sp.nwk --spec labels.tsv --internal-names keep
(O{Out},((M1,(A{Test},(B{Test},C)BC{Test})ABC{Test})X,(M2,(D{Test},E{Test})DE{Test})Y)Z);

When two labels claim the same branch, --precedence decides. Here Test spreads
over the clade above A and B, which includes C, and Rev is chosen for C directly.
The default specific lets a label chosen directly for a tip beat a label that
spreads from a clade.

# C keeps Rev, although Test covers it
$ twig tree-label -i sp.nwk -L Test:clade '^[AB]$' -L Rev:tips '^C$' --internal-names keep
(O,((M1,(A{Test},(B{Test},C{Rev})BC{Test})ABC{Test})X,(M2,(D,E)DE)Y)Z);

# first lets the label that was given first win, which is Test
$ twig tree-label -i sp.nwk -L Test:clade '^[AB]$' -L Rev:tips '^C$' --precedence first --internal-names keep
(O,((M1,(A{Test},(B{Test},C{Test})BC{Test})ABC{Test})X,(M2,(D,E)DE)Y)Z);

Labels already in a tree
------------------------
# labels that are already in the tree are kept, and win over new labels
$ echo '((A{Old},B),(C,D),E);' > old.nwk
$ twig tree-label -i old.nwk -L Test '^[AB]$' --rule tips
((A{Old},B{Test}),(C,D),E);

# --existing clear removes them first
$ twig tree-label -i old.nwk -L Test '^[AB]$' --rule tips --existing clear
((A{Test},B{Test}),(C,D),E);

# PAML labels such as #1 are read too, and --paml-codes gives the numbers names
$ echo '(A #1,B,(C,D) #2);' > paml.nwk
$ twig tree-label -i paml.nwk --paml-codes 1=Test 2=Reference
(A{Test},B,(C,D){Reference});

# -F newick writes the tree without labels
$ twig tree-label -i old.nwk -F newick
((A,B),(C,D),E);

# --stems acts on the labels already in a tree too, which makes the stem variant of a tree labeled elsewhere
$ echo '(O,((A{Test},B{Test}){Test},C));' > test.nwk
$ twig tree-label -i test.nwk --stems Stem
(O,((A{Test},B{Test}){Stem},C));

Rooting and writing
-------------------
# -r roots the tree on an outgroup first (tips joined by commas make a clade)
$ echo '((M1,(A,(B,C))),(M2,(D,E)),O);' > unrooted.nwk
$ twig tree-label -i unrooted.nwk -r O -L Test dioecious.txt --rule fitch
(O,((M1,(A{Test},(B{Test},C){Test}){Test}),(M2,(D{Test},E{Test}){Test})));

# -F hyphy writes a tree for HyPhy, which reads trees as unrooted. The two branches
# at the root are written as one branch, and --root-edges chooses its label
$ echo '(A{X}:1,(B{Y}:1,C{Y}:1){Z}:1);' > roots.nwk
$ twig tree-label -i roots.nwk -F hyphy
(A:2,B{Y}:1,C{Y}:1);
$ twig tree-label -i roots.nwk -F hyphy --root-edges Z
(A{Z}:2,B{Y}:1,C{Y}:1);

# --show lists the labeled branches on the screen, clade by clade, and the tree goes to the file
$ echo '((A,B),(C,D));' > small.nwk
$ twig tree-label -i small.nwk -L Test '^[AB]$' --show > small.labeled.nwk
labeled branches, clade by clade (X+Y is the branch above the common ancestor of X and Y)
  Test  stem  A+B
  Test  tip   A
  Test  tip   B
unlabeled branches
  -  internal  C+D
  -  tip       C
  -  tip       D
summary  3 Test, 3 unlabeled

# -o writes the tree to a file, and -B also writes a table with one row per branch
$ twig tree-label -i small.nwk -L Test '^[AB]$' -o small.out.nwk -B small.tsv
$ cut -f 3,5,7,8 small.tsv
node	n_tips	label	kind
A	1	Test	tip
B	1	Test	tip
C	1		tip
D	1		tip
S4	2	Test	stem
S5	2		internal

# --internal-names auto names the internal nodes twn1, twn2 and so on, so that the
# results of HyPhy for each branch can be joined to the table
$ twig tree-label -i small.nwk -L Test '^[AB]$' --internal-names auto
((A{Test},B{Test})twn1{Test},(C,D)twn2);

# --no-dists writes the tree without branch lengths
$ echo '((A:0.1,B:0.2):0.3,C:0.4);' > lengths.nwk
$ twig tree-label -i lengths.nwk -L Test '^[AB]$' --no-dists
((A{Test},B{Test}){Test},C);

# for HyPhy, tip names are made safe (here a dash and a dot become underscores)
$ echo '(A-1,A.2,(B,C));' > dashes.nwk
$ twig tree-label -i dashes.nwk -L Test A -F hyphy --rule tips --names-table names.tsv
(A_1{Test},A_2{Test},(B,C));

# names that become equal when they are made safe stop the run, unless --sanitize-suffix adds _2, _3
$ echo '(A-1,A.1,(B,C));' > same.nwk
$ twig tree-label -i same.nwk -F hyphy --sanitize-suffix
(A_1,A_1_2,(B,C));

# --keep-root writes a rooted tree as it is, instead of joining the two branches at the root
$ twig tree-label -i roots.nwk -F hyphy --keep-root
(A{X}:1,(B{Y}:1,C{Y}:1){Z}:1);"""


def get_parser_tree_label(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for the tree-label tool."""
    KWARGS = dict(
        prog="tree-label",
        usage="tree-label -i TREE [-L NAME SELECTOR...] [options]",
        help="label the branches of a tree from selections of tips",
        formatter_class=_formatter,
        description=dedent(f"""
            {FRAME}
            | tree-label puts labels on the branches of a tree
            {FRAME}
            | Programs such as HyPhy compare groups of branches, for example
            | the branches of species that have a trait (Test) with those of
            | species that do not (Reference). tree-label decides which
            | branches belong to which group, and writes the tree with the
            | label in braces after each branch, such as A{{Test}}. You choose
            | the tips with list files, regular expressions or species, and
            |   a rule spreads the label from those tips to branches
            |   (tips, all, crown, some), or from their common ancestor
            |   (clade, stem),
            |   parsimony works out where a trait arose and was lost (fitch),
            |   several labels can be given at once or read from a file, and
            |   labels that the tree already has can be kept, changed or
            |   converted, for example from PAML.
            | Stems, losses, the background and nuisance species have options
            | of their own. The examples below show each of them, and --show
            | lists the result so that it can be checked by eye.
            {FRAME}
        """),
        epilog=TREE_LABEL_EXAMPLES,
    )
    parser = _create(parser, KWARGS)

    io = parser.add_argument_group("Input/Output")
    io.add_argument("-i", "--input", type=Path, metavar="path", required=True, help="file with one or more Newick trees, labeled or not")
    io.add_argument("-o", "--out", type=Path, metavar="path", help="write the output to this file [stdout]")
    io.add_argument("-F", "--format", type=str, metavar="str", choices=FORMAT_NAMES, default="braces", help="what to write (" + ", ".join(FORMAT_NAMES) + "). braces is the tree with labels in braces after the names, hyphy is the same made ready for HyPhy (safe tip names, written unrooted), newick is the plain tree without labels, and table has one row per branch [%(default)s]")
    io.add_argument("-B", "--branch-table", type=Path, metavar="path", help="also write a table with one row per branch, with its label, its role, the tips below it and where the label came from")
    io.add_argument("--show", action="store_true", help="print the labeled branches, clade by clade, to the screen so that the result can be checked by eye")
    io.add_argument("--paml-codes", type=str, metavar="str", nargs="+", help="names for the PAML label numbers in the input, such as 1=Test 2=Reference")

    select = parser.add_argument_group("Selection")
    select.add_argument("-L", "--label", type=str, metavar="str", nargs="+", action="append", help="a label name followed by the tips it goes to, which are list files or regular expressions. Repeat -L for more labels. A rule can follow the name, joined to it with a colon, to give that label its own rule")
    select.add_argument("--spec", type=Path, metavar="path", help="read the labels from a tab-separated file with the columns label, rule, selectors and options, one label per line, instead of using -L")
    select.add_argument("--nuisance", type=str, metavar="str", nargs="+", help="species to set aside for any reason, such as a reversion, an uncertain trait or poor data. They get a label of their own that no other label overwrites, and parsimony treats their state as unknown")
    select.add_argument("--nuisance-label", type=str, metavar="str", default=DEFAULT_NUISANCE_LABEL, help="the label given to nuisance tips [%(default)s]")
    select.add_argument("--rule", "--internal-nodes", dest="rule", type=str, metavar="str", default="all", help="how far a label spreads from the chosen tips (" + ", ".join(RULE_NAMES) + "). tips is only the chosen tips, all adds the branches whose tips are all chosen, crown is all without the branch at the base of each group, some is every branch with a chosen tip below it, fitch (also called parsimony) is the branches where the trait is inferred to be present, and clade and stem start at the common ancestor of the chosen tips and label its branch and everything below it, or only its branch. HyPhy names such as 'All descendants' work too [%(default)s]")
    select.add_argument("--select-on", type=str, metavar="str", choices=SELECT_ON_CHOICES, default="name", help="match the chosen names to the names of the tips, or to the species of the tips when several samples belong to one species (" + ", ".join(SELECT_ON_CHOICES) + ") [%(default)s]")
    select.add_argument("--exact", action="store_true", help="treat the names given on the command line as exact names, so that A does not also match MA")
    select.add_argument("--invert", action="store_true", help="label the tips that are not chosen")
    select.add_argument("--leaf-nodes", type=str, metavar="str", choices=LEAF_NODE_CHOICES, default="label", help="label the tips as usual, or remove the labels of tips after the rule has run (" + ", ".join(LEAF_NODE_CHOICES) + ") [%(default)s]")
    select.add_argument("--stems", type=str, metavar="str", default="keep", help="a stem is the branch on which a trait arose, where its label begins (with fitch a branch with a gain, with crown the branch the rule leaves out at the base of each group). Choose keep, unlabel, or a label name, which gives stems a label of their own. Labels that are already in the tree count too [%(default)s]")
    select.add_argument("--losses", type=str, metavar="str", default="keep", help="a loss is a branch where parsimony puts the loss of the trait again. Choose keep, unlabel, or a label name, which gives losses a label of their own [%(default)s]")
    select.add_argument("--background", type=str, metavar="str", help="the label given to every branch that is still free after everything else")
    select.add_argument("--precedence", type=str, metavar="str", choices=PRECEDENCE_CHOICES, default="specific", help="which label wins a branch that two labels claim (" + ", ".join(PRECEDENCE_CHOICES) + "). specific lets a label chosen directly for a tip beat one that spreads from a clade, first and last let the earlier or later label win, and error stops [%(default)s]")
    select.add_argument("--existing", type=str, metavar="str", choices=EXISTING_CHOICES, default="keep", help="keep the labels that are already in the tree, and let them win over new labels, or clear them first (" + ", ".join(EXISTING_CHOICES) + ") [%(default)s]")
    select.add_argument("-r", "--outgroups", type=str, metavar="str", nargs="+", help="root the tree on an outgroup first. Give tips (or species), joined by commas for a clade, and the first outgroup that forms one side of a branch is used")

    pars = parser.add_argument_group("Parsimony")
    pars.add_argument("--ambiguous", type=str, metavar="str", default="acctran", help="two histories can be equally good. acctran puts the changes nearer the root, deltran nearer the tips, unlabel leaves the branches that differ without a label, and a label name gives them a label of their own [%(default)s]")
    pars.add_argument("--gain-cost", type=str, metavar="float", default="1", help="how costly it is to gain the trait, compared with losing it. A higher cost makes parsimony avoid gains [%(default)s]")
    pars.add_argument("--loss-cost", type=str, metavar="float", default="1", help="how costly it is to lose the trait, compared with gaining it. A higher cost makes parsimony avoid losses [%(default)s]")
    pars.add_argument("--unknown", type=str, metavar="str", nargs="+", help="tips whose state is not known, so that parsimony may give them either state")
    pars.add_argument("--root-state", type=str, metavar="str", choices=ROOT_STATE_CHOICES, default="free", help="force the root to lack the trait (0) or to have it (1), or let parsimony choose (" + ", ".join(ROOT_STATE_CHOICES) + ") [%(default)s]")

    add_taxon_mapping_args(parser, unmapped=False)
    add_writer_args(parser)

    runtime = parser.add_argument_group("Runtime")
    runtime.add_argument("-f", "--force", action="store_true", help="overwrite existing output files")
    runtime.add_argument("-l", "--log-level", type=str, metavar="level", default="INFO", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    return parser


def get_parser_tree_label_transfer(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for the tree-label-transfer tool (planned)."""
    return _get_planned_parser(
        parser, "tree-label-transfer", "copy the labels of a species tree onto gene trees by branch matching", "M2",
    )


def get_parser_sample_reduce(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for the sample-reduce tool (planned)."""
    return _get_planned_parser(
        parser, "sample-reduce", "reduce gene trees and alignments to one specimen per species", "M4",
    )
