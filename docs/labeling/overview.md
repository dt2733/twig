
# Tree labeling

Programs that test for changes in natural selection, such as HyPhy's RELAX, BUSTED-PH, BUSTED and aBSREL, read a tree in which every branch can carry a label. The label puts the branch in a set, for example the test set or the reference set, and the method compares the sets. The tree labeling tools in twig build these labeled trees for thousands of genes at once, starting from one species tree that you label by hand.

This page explains the ideas that the tools share. The complete design, including the reasons for each choice, is in [design.md](design.md).

## Commands and their status

The commands are developed in milestones, and a command that is planned but not available yet prints a message that says so.

| command | purpose | status |
|---|---|---|
| `twig name-map` | resolve every tip and sequence name to a species, print the table for checking, and freeze it for reuse | available, see [name-map](name-map.md) |
| `twig tree-label` | label one tree, usually the species tree, from tip selections | available, see [tree-label](tree-label.md) |
| `twig tree-label-transfer` | copy the labels of the species tree onto each gene tree and write the trees with their alignments | planned for milestone M2 |
| `twig sample-reduce` | reduce gene trees and alignments to one specimen per species | planned for milestone M4 |

## Requirements

The tree labeling code needs Python 3.10 or newer, as does the rest of the twig command line, and toytree 3.0.10 or newer, which is installed with twig. Its tests pass with toytree 3.0.10, 3.0.11 and the development version. HyPhy is not needed to create labeled trees, and twig does not run it.

## The problem

You know the history of a trait, such as dioecy, on the species tree, and you want to ask whether selection changed along the branches where the trait changed. Gene trees often differ from the species tree, through incomplete lineage sorting or estimation error. Forcing every gene onto the species tree puts substitutions on the wrong branches and inflates false positives. Running each gene on its own tree avoids that, but then every gene tree needs its own labels.

## Species tree branch matching

You label the species tree once. For each gene, a branch of the gene tree takes the label of the species tree branch that has exactly the same set of species below it. A gene branch with no counterpart, because the gene tree disagrees with the species tree there, stays unlabeled and drops out of the comparison. Tips always match, because a tip has exactly one species.

A stem, which is the branch that leads to a clade with the trait, is therefore labeled only in genes whose tree recovers that clade.

## Labels live on branches

A label belongs to a branch, and twig stores it on the node below the branch. The root has no branch, so it never carries a label. Labels follow their branches when a tree is rerooted. Trees can be read with HyPhy braces such as `A{Test}`, with PAML marks such as `A #1` and clade labels such as `$1`, and with the `label` entry of NHX and FigTree annotations. Trees are written in the HyPhy form, with the label after the name and before the branch length.

```
(O{Reference},((M1{Reference},(A{Test},(B{Test},C{Rev}){Test})){Reference},(M2{Reference},(D{Test},E{Test})){Reference}){Reference});
```

## Stems and the point of change

The change from one state to another happened somewhere along a stem, at an unknown point. A stem can be left unlabeled, which tests only the genes of the clade that keeps the trait, or it can be labeled, which also includes the change itself. Running both and comparing them separates genes that are associated with the origin of the trait from genes that are associated with its maintenance. The tools are designed so that one labeled species tree supports both analyses.

## Nuisance species

A nuisance species is one that you do not want in the test set or in the reference set. A reversion to the ancestral state is one reason, and there are many others, such as an uncertain trait state, a hybrid, or poor sequence quality. Most people cannot know in advance that a species is a reversion, so twig does not depend on detecting reversions. You mark any species or clade as nuisance with `--nuisance` in [tree-label](tree-label.md), whether or not parsimony finds a loss there. Nuisance species get a label of their own, no other label overwrites it, and parsimony treats their state as unknown, so the history of the trait is inferred from the other species.

How nuisance species are handled by a program is a choice made for each analysis, because it depends on the program. RELAX can set them aside as unclassified branches that are fitted but left out of the comparison. BUSTED-PH has no such set, so a nuisance species must either be removed from the tree and the alignment or accepted as background. That choice is made when the labels are copied to gene trees, which `tree-label-transfer` will do. The choice of label name is yours, and the label `Rev` in the tree above is only an example.

## What an unlabeled branch means

Leaving a branch unlabeled is not the same in every program, and this decides many defaults in the tools.

| program | test set | reference set | unlabeled branches |
|---|---|---|---|
| RELAX, Classic mode with `--test` and `--reference` | test label | reference label | unclassified, fitted with their own distribution and left out of the comparison |
| RELAX with only `--test` | test label | every other branch | none, everything else silently becomes reference |
| BUSTED and BUSTED-PH | the label given to `--branches` | every other branch | none, everything else is background |
| aBSREL | the label given to `--branches`, each branch tested separately | none | modeled but not tested |
| PAML codeml | `#1` and further codes | `#0` | `#0` |

For every program, twig will report how many branches, and how much branch length, end up in the background without being reference by choice.

## Branch statuses

When labels are transferred to a gene tree, every branch gets a status that says how it was labeled. The statuses appear in the tables that `tree-label-transfer` will write.

| status | meaning |
|---|---|
| `labeled` | matched a species tree branch that carries a label |
| `species-unlabeled` | matched a species tree branch that was deliberately left unlabeled, such as a stem |
| `unmatched` | its set of species does not exist in the species tree, which means the gene tree disagrees there |
| `compatible` | it resolves a polytomy of the species tree without contradicting it |
| `conflict` | it matched a species tree branch that was formed by merging branches whose labels disagree |
| `mixed` | samples of one species lie on both sides of it |
| `paralog` | it lies at or below a duplication node and was masked |
| `unmapped` | its tip could not be assigned to a species and the run was told to keep it |
| `masked` | a mask changed its label |

## Samples, specimens and species

A sample is a tip of a gene tree or a row of an alignment, under its raw name. A specimen is one individual. A species is a tip of the species tree, and name mapping always resolves samples to that vocabulary. Gene trees can hold several specimens of a species and several copies of a gene per specimen. The names of the same species are often written in several ways, such as `A_pumilus_sample1` and `AmaranthusPumilus-Sample1`, which is why the first tool in the workflow is [name-map](name-map.md).

## Why low-support branches are not collapsed

The tools never collapse branches with low support into polytomies. A collapsed branch disappears together with its label, labeled stems are often the short and weakly supported branches that collapsing would remove, HyPhy accounts for uncertainty through its own estimates, and PAML accepts at most three branches at any node. If you want to collapse branches anyway, run `twig tree-filter` with its support option first.

## Files and tests

All of the code is in `src/twig/labeling`, the tests and their example data are in `tests/labeling`, and these pages are in `docs/labeling`. The tests run with `pytest tests/labeling` after `pip install pytest`.
