# tree-label

The `twig tree-label` tool labels the branches of a tree from selections of tips. You usually run it once, on your species tree, and the labeled species tree is then the input of `tree-label-transfer`, which copies the labels onto every gene tree. It can also label a single gene tree directly. It writes the result in the form HyPhy reads, as a plain tree, or as a table with one row per branch.

You choose tips by name, with list files or regular expressions, or by species. A rule then decides which branches get the label. The rules include those of the HyPhy `label-tree` script and real two-state parsimony, which the HyPhy option of that name does not do. Three switches refine the result. Stems can be left unlabeled or given their own label, losses found by parsimony can be left unlabeled or given their own label, and every branch that is still free can get a background label. Species that belong in neither the test set nor the reference set, for any reason, can be set aside as nuisance species.

The example files used on this page are in `tests/labeling/data`. The species tree `species.nwk` is `(O,((M1,(A,(B,C))),(M2,(D,E))));`, where O is the outgroup, M1 and M2 are monoecious, A, B, D and E are dioecious, and C is the species that lost dioecy again. The file `dioecious.txt` lists A, B, D and E.

## The tree-label command

```bash
twig tree-label -h
```

```text
usage: tree-label -i TREE [-L NAME SELECTOR...] [options]

-------------------------------------------------------------------
| tree-label labels the branches of a tree from selections of tips
-------------------------------------------------------------------
| Tips are selected by name (list files or regular expressions), or
| by species with --select-on species, and a rule decides which
| branches get each label. The rules are those of HyPhy label-tree
| (tips, all, crown, some) and real two-state parsimony (fitch), plus
| clade and stem. Stems, the branches on which a label begins, and
| losses, the branches on which parsimony puts the loss of a trait,
| can be left unlabeled or given a label of their own. Species that
| belong in neither the test nor the reference set, for any reason,
| can be set aside with --nuisance. Labels are written in braces
| after the name, as HyPhy reads them, and --show lists the labeled
| branches so the result can be checked by eye.
-------------------------------------------------------------------

options:
  -h, --help                       show this help message and exit

Input/Output:
  -i, --input path                 file with one or more Newick trees, labeled or not
  -o, --out path                   write the output to this file [stdout]
  -F, --format str                 output format (braces, hyphy, newick, table) [braces]
  -B, --branch-table path          also write a table with one row per branch
  --show                           print the labeled branches, clade by clade, to stderr
  --paml-codes str [str ...]       names for PAML label numbers, such as 1=Test 2=Reference

Selection:
  -L, --label str [str ...]        label name followed by selectors, which are list files or regular expressions; repeat for more labels, and see the examples to give a label its own rule
  --spec path                      tab-separated file with the columns label, rule, selectors and options, one label per line
  --nuisance str [str ...]         selectors for tips to set aside for any reason; they get the nuisance label, no other label overwrites it, and parsimony treats their state as unknown
  --nuisance-label str             label given to nuisance tips [Nuisance]
  --rule, --internal-nodes str     branches covered by a label (tips, all, crown, some, fitch, clade, stem), or a HyPhy name such as 'All descendants' [all]
  --select-on str                  match selectors to tip names or to the species of tips (name, species) [name]
  --exact                          treat selectors given on the command line as exact names
  --invert                         label the tips that the selectors do not match
  --leaf-nodes str                 label tips, or skip them after the rule has run (label, skip) [label]
  --stems str                      keep stems, unlabel them, or give them a label [keep]
  --losses str                     keep branches where parsimony puts a loss, unlabel them, or give them a label [keep]
  --background str                 label given to every branch that is still free
  --precedence str                 which label wins a branch claimed twice (specific, first, last, error) [specific]
  --existing str                   keep the labels already in the tree, or clear them first (keep, clear) [keep]
  -r, --outgroups str [str ...]    root the tree on the first of these outgroups that fits, tips (or species) joined by commas for a clade

Parsimony:
  --ambiguous str                  branches that differ between equally good histories follow acctran or deltran, stay unlabeled (unlabel), or get a label [acctran]
  --gain-cost float                cost of gaining the trait [1]
  --loss-cost float                cost of losing the trait [1]
  --unknown str [str ...]          selectors for tips whose state is unknown
  --root-state str                 state of the root (free, 0, 1) [free]

Taxon Mapping:
  -I, --imap path                  map file rows <sample> <species> [specimen]; ~regex rows allowed
  -d, --delim str                  delimiter used to split sample names
  -x, --name-regex str             regex with (?P<species>...) and optional (?P<specimen>...)
  -di, --delim-idxs int [int ...]  token indices kept after splitting [0]
  -dj, --delim-join str            join string for kept tokens [same as --delim]
  -N, --normalize str [str ...]    normalization steps used for matching, or none [accents camel sep case]
  --accept-inferred                accept genus-abbreviation matches

Output Format:
  --root-edges str                 label of the one branch that HyPhy and PAML see in place of the two branches at the root (unlabel, labeled, outgroup, ingroup, error, or a label that wins) [unlabel]
  --keep-root                      write the tree rooted as it is, without joining the two branches at the root
  --internal-names str             names written for internal nodes (none, auto, support, keep) [none]
  --no-dists                       write trees without branch lengths
  --sanitize-suffix                add _2, _3 to tip names that become equal when they are made safe for HyPhy
  --names-table path               write the table of original and written tip names

Runtime:
  -f, --force                      overwrite existing output files
  -l, --log-level level            stderr logging level (DEBUG, [INFO], WARNING, ERROR)
```

## Selecting tips

A selector is either a list file or a regular expression. A value that names an existing file is read as a list with one name per line, where a line that starts with `~` is a regular expression and a line that starts with `#` is a comment. Any other value is a regular expression that is matched with `re.search`, so `A` matches every name that contains an A. The option `--exact` treats the values given on the command line as exact names, and `--invert` labels the tips that the selectors do not match.

With `--select-on species` the selectors are matched to the species of each tip, which makes one selection work on gene trees that have several samples per species. The species come from the Taxon Mapping options, which are the same as in [name-map](name-map.md). For example `-d _` takes the part of each name before the first underscore. Without a map file, delimiter or regex, every tip is its own species, and the command says so.

Every label reports how many tips it selected, how many branches it covers, and how many it keeps after precedence, so a label that another label overrode never goes unnoticed. Names or patterns that match no tip are reported too, which helps to find typos.

## Several labels

Repeat `-L` for more labels. Each label can have its own rule, which is written after the label name with a colon, for example `-L Test:fitch dioecious.txt -L Out:tips '^O$'`. The labels can also be listed in a spec file given with `--spec`. It has one tab-separated line for each label with the columns label, rule, selectors and options, where the rule can be left empty and the options can be `ambiguous=`, `gain-cost=` and `loss-cost=`. A line that starts with `#` is a comment.

When two labels claim the same branch, `--precedence` decides. The default `specific` lets a label that a tip received by direct selection beat a label that a branch inherits from the clade around it, so a species that you select for its own label is never overwritten by the label of its clade. The value `first` lets earlier labels win, as HyPhy does, `last` lets later labels win, and `error` stops the run.

Labels that are already in the tree are kept by default and have priority over new labels, as in HyPhy. The option `--existing clear` removes them first. If you give no label at all, the labels in the tree are kept, so `tree-label` can also be used only to convert a labeled tree to another format.

## Rules

The option `--rule`, also called `--internal-nodes`, decides which branches a label covers. The HyPhy names are accepted in any case, so `--internal-nodes "All descendants"` works.

| rule | HyPhy name | branches covered |
|---|---|---|
| `tips` | None | the selected tips only |
| `all` | All descendants | every branch all of whose descendant tips are selected |
| `crown` | All descendants, no MRCA | the result of `all` without the internal branches whose parent branch is not labeled, so isolated selected tips stay labeled |
| `some` | Some descendants | every branch with at least one selected tip below it |
| `fitch` | Parsimony | the branches on which two-state parsimony places the trait, see below |
| `clade` | not in HyPhy | for each selector, the branch above the common ancestor of its tips and everything below it |
| `stem` | not in HyPhy | for each selector, only the branch above the common ancestor of its tips |

The option `--leaf-nodes skip` removes the labels of tips after the rule has run. With `fitch` the tips still provide the data. The rule `tips` together with `skip` is an error because it would label nothing. The rules `clade` and `stem` are errors when the common ancestor is the root, because the root has no branch.

The next table shows what each rule labels on the species tree when the selection is A, B, D and E. The internal branches are named in `species.named.nwk` as BC, ABC, X (M1 and its sister), DE, Y (M2 and its sister) and Z (the whole ingroup).

| rule and options | labeled branches |
|---|---|
| `tips` | A, B, D, E |
| `all` | A, B, D, E, DE |
| `crown` | A, B, D, E |
| `some` | A, B, D, E, BC, ABC, X, DE, Y, Z |
| `fitch` | A, B, BC, ABC, D, E, DE |
| `fitch --stems unlabel` | A, B, BC, D, E |
| `fitch --ambiguous deltran` | A, B, D, E, DE |

For example, the `fitch` row is produced by this command. The option `--internal-names keep` writes the names of the internal branches so that they can be read in the output.

```bash
twig tree-label -i species.named.nwk -L Test dioecious.txt --rule fitch --internal-names keep
```

```text
(O,((M1,(A{Test},(B{Test},C)BC{Test})ABC{Test})X,(M2,(D{Test},E{Test})DE{Test})Y)Z);
```

## Parsimony

The `fitch` rule reconstructs a trait that is either absent or present, using parsimony, and labels every branch on which the trait is present. Selected tips have the trait, tips that are not selected do not, and tips given with `--unknown` cost nothing in either state. Gaining the trait costs `--gain-cost` and losing it costs `--loss-cost`, both 1 by default. The costs are kept as exact fractions, so equal totals are recognized without rounding errors. The state of the root is free by default, and `--root-state 0` or `--root-state 1` fixes it. If the root has three or more children and the two costs differ, the result depends on where the tree is rooted, and the command warns.

Often two histories are equally good. On the species tree the trait can have been gained once on the branch above A, B and C and then lost on C, or it can have been gained separately on A and on B. Both histories cost 3 in total, which includes one gain on the branch above D and E. The option `--ambiguous` says what to do with the branches that differ between the equally good histories. The default `acctran` moves changes toward the root, which gives the first history. The value `deltran` delays changes toward the tips, which gives the second history. The value `unlabel` leaves those branches unlabeled, and a label name gives them a label of their own so that they can be handled separately later. The same tree under `deltran` has no label on the branches above A, B and C.

```bash
twig tree-label -i species.nwk -L Test dioecious.txt --rule fitch --ambiguous deltran
```

```text
(O,((M1,(A{Test},(B{Test},C))),(M2,(D{Test},E{Test}){Test})));
```

## Stems, losses and background

These three switches work with every rule and turn the result of a rule into the labeling that your question needs.

A stem is a labeled internal branch whose parent branch does not carry the same label. It is the branch on which the label begins, and the change happened somewhere along its length. The option `--stems unlabel` leaves stems out of both the test and the reference set, and a name such as `--stems Stem` gives them a class of their own. With a class of their own, one labeled species tree supports both an analysis of the crown only and an analysis of stem plus crown, because the choice can be made later for each run.

A loss is a branch with the trait absent whose parent has it present, which is a reversion to the ancestral state. Only parsimony finds losses. The option `--losses unlabel` leaves them unlabeled, and `--losses Rev` gives them the label Rev. You do not need to know which species are reversions, because any species can be set aside as a nuisance species, as described in the next section.

The option `--background Reference` gives the label Reference to every branch that no label covers, except the stems, losses and ambiguous branches that you deliberately left unlabeled. This makes the reference branches explicit.

## Nuisance species

A nuisance species is one that you do not want in the test set or in the reference set. A reversion is one reason, and there are many others, such as an uncertain trait state, a hybrid, or poor sequence quality. You mark species with `--nuisance`, using the same selectors as for any label, and you choose the label with `--nuisance-label`, which is Nuisance by default.

A nuisance label has three effects. No other label overwrites it, whatever the precedence and whether or not the tree already had a label there. It is not touched by `--stems`, `--leaf-nodes skip` or `--background`. And the tips are treated as unknown by parsimony, so the history of the trait is inferred from the other species without them.

```bash
twig tree-label -i species.nwk -L Test dioecious.txt --rule fitch --stems unlabel --nuisance C --background Reference
```

```text
(O{Reference},((M1{Reference},(A{Test},(B{Test},C{Nuisance}){Test})){Reference},(M2{Reference},(D{Test},E{Test})){Reference}){Reference});
```

Here C is set aside. Parsimony no longer has to explain it, so the cost of the history falls from 3 to 2, and C keeps its own label. Take care that unknown states can change the history elsewhere. If M1 and M2 were set aside as well, for example, a single gain above the whole ingroup would become the best history and the branches above M1, M2 and the ingroup would be labeled Test.

What a program does with the nuisance label is decided later, for each analysis, when the labels are copied to gene trees. RELAX can set such branches aside, and BUSTED-PH has no such set, so the species must either be removed from the trees and alignments or be accepted as background. That choice belongs to `tree-label-transfer`.

## Rooting and the two branches at the root

HyPhy and PAML read every tree as unrooted. A rooted tree has two branches at its root, and the programs see them as one branch. HyPhy keeps only the label of one of them, and which one depends on the order of the two children in the Newick string, so the result can change silently. For this reason twig always joins the two branches itself, with a rule that you choose with `--root-edges`, whenever it writes the `hyphy` format.

| value | label of the joined branch |
|---|---|
| `unlabel` | the shared label if the two branches agree, otherwise no label and a warning |
| `labeled` | the label of the only labeled branch, and no label if both are labeled differently |
| `outgroup` | the label of the branch on the side of the outgroup |
| `ingroup` | the label of the branch on the side of the ingroup |
| `error` | none, the run stops if the labels differ |
| a label name | that label if either branch carries it |

On the species tree of this page the two branches at the root are O and the whole ingroup, and in the labeled tree both are Reference, so the joined branch is Reference. In the next tree they differ, with X on O and Z on the ingroup.

```bash
twig tree-label -i root-edges.nwk -F hyphy --root-edges unlabel
```

```text
(O:2,(M1{R}:1,(A{T}:1,(B{T}:1,C{V}:1){T}:1){S}:1){R}:1,(M2{R}:1,(D{T}:1,E{T}:1){W}:1){R}:1);
```

```bash
twig tree-label -i root-edges.nwk -F hyphy --root-edges Z
```

```text
(O{Z}:2,(M1{R}:1,(A{T}:1,(B{T}:1,C{V}:1){T}:1){S}:1){R}:1,(M2{R}:1,(D{T}:1,E{T}:1){W}:1){R}:1);
```

The joined branch has the summed length of the two branches. For the `outgroup` and `ingroup` values the side of the outgroup is the one that contains the outgroup that `-r` used for rooting, and without `-r` it is the smaller side. If the two sides are equally large, twig asks you to give `-r`. The option `--keep-root` skips the joining and writes the tree as it is, which is useful for programs that read rooted trees.

The option `-r` roots the tree on an outgroup before any label is applied. Give tips or species, joined by commas for a clade, and the first outgroup that forms one side of a branch is used. Labels stay on their branches when the tree is rerooted, and the two old root branches are joined by the same `--root-edges` value. An unrooted tree can be rooted in this way.

```bash
twig tree-label -i species.unrooted.nwk -r O -L Test dioecious.txt --rule fitch -F hyphy
```

```text
(O,(M1,(A{Test},(B{Test},C){Test}){Test}),(M2,(D{Test},E{Test}){Test}));
```

## Writing the tree

The option `-F` chooses the output format.

| format | output |
|---|---|
| `braces` | the labels in braces after the name, as HyPhy reads them, with the tree rooted as it was given |
| `hyphy` | the same, with tip names made safe for HyPhy and the tree written unrooted as described above |
| `newick` | a plain tree with the labels removed |
| `table` | one row per branch, described below |

Branch lengths are written when the input had them, and `--no-dists` leaves them out. A tree that was read without branch lengths is written without them. Internal node names are not written by default. The option `--internal-names` can write the support values (`support`), the names the tree had (`keep`), or new names `twn1`, `twn2` and so on in postorder (`auto`). The new names let the per-branch results of HyPhy be joined to the branch table.

For the `hyphy` format twig makes every tip name safe in the way HyPhy would. Characters other than letters, digits and underscores become underscores, and a name that HyPhy reserves gets an underscore added. Names that become equal, without regard to case, stop the run with a list, or get `_2`, `_3` and so on when `--sanitize-suffix` is given. The option `--names-table` saves the table of original and written names. The names are made safe once for all the trees of the input, so a sample is written with the same name in every tree. Label names must start with a letter, contain only letters, digits and underscores, and differ from the HyPhy branch sets All, Internal, Leaves, Unlabeled and fg and from the words that twig uses for mask actions and branch statuses.

```bash
twig tree-label -i gene.nwk -L Test dioecious.txt --select-on species -d _ --rule all -F hyphy --internal-names auto
```

```text
(O_1:0.18,(M1_1:0.2,(A_1{Test}:0.05,A_2{Test}:0.04,(B_1{Test}:0.1,C_1:0.12)twn1:0.02)twn2:0.03)twn3:0.04,(M2_1:0.2,(D_1{Test}:0.1,E_1{Test}:0.11)twn4{Test}:0.03)twn5:0.05);
```

In this gene tree A has two samples. The selector A matched the species of each tip, so both samples are labeled, and the branch above O has the summed length 0.12 and 0.06.

## The branch table

The option `-B` writes a table with one row per branch next to the tree, and `-F table` writes it instead of the tree. The root has no branch and has no row. The table describes the tree as it is written, so with `-F hyphy` the dissolved root branch is marked.

| column | meaning |
|---|---|
| `locus` | the name of the tree, which is the file name before its first dot, followed by a number if the file holds several trees |
| `idx` | the number of the node in toytree |
| `node` | the tip name, the name of the internal node if it is unique, or S followed by the idx |
| `is_tip` | yes or no |
| `n_tips` | the number of tips below the branch |
| `tips` | the names of the tips below the branch, cut after 200 |
| `species` | the species below the branch, only with `--select-on species` |
| `label` | the label as written |
| `kind` | tip, stem, crown (a labeled branch with the same label as its parent) or internal |
| `source` | where the label came from, which is the label and rule such as `Test:fitch`, or existing, stem, loss, background, or nothing if unlabeled |
| `dist` | the branch length as written |
| `support` | the support value of the node, if it has one |
| `written_name` | the name written for the node |
| `written` | no for the root child that is dissolved when the root branches are joined |
| `merged_into` | the node that receives the dissolved branch |
| `note` | root-merge for the two branches that were joined |

The kind of a branch is found before `--stems` is applied, so a stem that was left unlabeled still has the kind stem. The next command shows some of the columns for the species tree with unlabeled stems.

```bash
twig tree-label -i species.nwk -L Test dioecious.txt --rule fitch --stems unlabel -F table | cut -f 3,5,7-9
```

```text
node  n_tips  label  kind      source
O     1              tip
M1    1              tip
A     1       Test   tip       Test:fitch
B     1       Test   tip       Test:fitch
C     1              tip
M2    1              tip
D     1       Test   tip       Test:fitch
E     1       Test   tip       Test:fitch
S8    2       Test   crown     Test:fitch
S9    3              stem
S10   4              internal
S11   2              stem
S12   3              internal
S13   7              internal
```

## Checking the result

The option `--show` prints one line for every labeled branch to stderr, in the order of the tree, followed by the unlabeled branches. Each line gives the label, the role of the branch and the tips below it, where X+Y means the branch above the common ancestor of X and Y. The roles are tip, stem, loss, crown and internal. The stems that were left without a label are listed as well, so that they are not missed.

## Step 1 of an analysis

This command labels the species tree of this page for a RELAX or BUSTED-PH analysis. Parsimony finds where dioecy arose and was lost, the stems are left unlabeled, the loss on C gets the label Rev, and everything else is Reference. The listing is printed to stderr, and the labeled tree is written to stdout. Here both are shown, and the option `-l ERROR` leaves out the log messages.

```bash
twig tree-label -i species.nwk -L Test dioecious.txt --rule fitch --stems unlabel --losses Rev --background Reference --show -l ERROR 2>&1
```

```text
labeled branches, clade by clade (X+Y is the branch above the common ancestor of X and Y)
  Reference  tip       O
  Reference  internal  M1+A+B+C+M2+D+E
  Reference  internal  M1+A+B+C
  Reference  tip       M1
  Test       tip       A
  Test       crown     B+C
  Test       tip       B
  Rev        loss      C
  Reference  internal  M2+D+E
  Reference  tip       M2
  Test       tip       D
  Test       tip       E
unlabeled branches
  -  stem  A+B+C
  -  stem  D+E
summary  6 Reference, 5 Test, 1 Rev, 2 unlabeled
(O{Reference},((M1{Reference},(A{Test},(B{Test},C{Rev}){Test})){Reference},(M2{Reference},(D{Test},E{Test})){Reference}){Reference});
```
