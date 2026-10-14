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
| tree-label puts labels on the branches of a tree
-------------------------------------------------------------------
| Programs such as HyPhy compare groups of branches, for example
| the branches of species that have a trait (Test) with those of
| species that do not (Reference). tree-label decides which
| branches belong to which group, and writes the tree with the
| label in braces after each branch, such as A{Test}. You choose
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
-------------------------------------------------------------------

options:
  -h, --help                       show this help message and exit

Input/Output:
  -i, --input path                 file with one or more Newick trees, labeled or not
  -o, --out path                   write the output to this file [stdout]
  -F, --format str                 what to write (braces, hyphy, newick, table). braces is the tree with labels in braces after the names, hyphy is the same made ready for HyPhy (safe tip names, written unrooted), newick
                                   is the plain tree without labels, and table has one row per branch [braces]
  -B, --branch-table path          also write a table with one row per branch, with its label, its role, the tips below it and where the label came from
  --show                           print the labeled branches, clade by clade, to the screen so that the result can be checked by eye
  --paml-codes str [str ...]       names for the PAML label numbers in the input, such as 1=Test 2=Reference

Selection:
  -L, --label str [str ...]        a label name followed by the tips it goes to, which are list files or regular expressions. Repeat -L for more labels. A rule can follow the name, joined to it with a colon, to give that
                                   label its own rule
  --spec path                      read the labels from a tab-separated file with the columns label, rule, selectors and options, one label per line, instead of using -L
  --nuisance str [str ...]         species to set aside for any reason, such as a reversion, an uncertain trait or poor data. They get a label of their own that no other label overwrites, and parsimony treats their state
                                   as unknown
  --nuisance-label str             the label given to nuisance tips [Nuisance]
  --rule, --internal-nodes str     how far a label spreads from the chosen tips (tips, all, crown, some, fitch, clade, stem). tips is only the chosen tips, all adds the branches whose tips are all chosen, crown is all
                                   without the branch at the base of each group, some is every branch with a chosen tip below it, fitch (also called parsimony) is the branches where the trait is inferred to be present,
                                   and clade and stem start at the common ancestor of the chosen tips and label its branch and everything below it, or only its branch. HyPhy names such as 'All descendants' work too [all]
  --select-on str                  match the chosen names to the names of the tips, or to the species of the tips when several samples belong to one species (name, species) [name]
  --exact                          treat the names given on the command line as exact names, so that A does not also match MA
  --invert                         label the tips that are not chosen
  --leaf-nodes str                 label the tips as usual, or remove the labels of tips after the rule has run (label, skip) [label]
  --stems str                      a stem is the branch on which a trait arose, where its label begins (with fitch a branch with a gain, with crown the branch the rule leaves out at the base of each group). Choose keep,
                                   unlabel, or a label name, which gives stems a label of their own. Labels that are already in the tree count too [keep]
  --losses str                     a loss is a branch where parsimony puts the loss of the trait again. Choose keep, unlabel, or a label name, which gives losses a label of their own [keep]
  --background str                 the label given to every branch that is still free after everything else
  --precedence str                 which label wins a branch that two labels claim (specific, first, last, error). specific lets a label chosen directly for a tip beat one that spreads from a clade, first and last let
                                   the earlier or later label win, and error stops [specific]
  --existing str                   keep the labels that are already in the tree, and let them win over new labels, or clear them first (keep, clear) [keep]
  -r, --outgroups str [str ...]    root the tree on an outgroup first, which every rule except tips needs when the tree is unrooted. Give tips (or species), joined by commas for a clade, and the first outgroup that forms
                                   one side of a branch is used

Parsimony:
  --ambiguous str                  two histories can be equally good. acctran puts the changes nearer the root, deltran nearer the tips, unlabel leaves the branches that differ without a label, and a label name gives
                                   them a label of their own [acctran]
  --gain-cost float                how costly it is to gain the trait, compared with losing it. A higher cost makes parsimony avoid gains [1]
  --loss-cost float                how costly it is to lose the trait, compared with gaining it. A higher cost makes parsimony avoid losses [1]
  --unknown str [str ...]          tips whose state is not known, so that parsimony may give them either state
  --root-state str                 force the root to lack the trait (0) or to have it (1), or let parsimony choose (free, 0, 1) [free]

Taxon Mapping:
  -I, --imap path                  map file rows <sample> <species> [specimen]; ~regex rows allowed
  -d, --delim str                  delimiter used to split sample names
  -x, --name-regex str             regex with (?P<species>...) and optional (?P<specimen>...)
  -di, --delim-idxs int [int ...]  token indices kept after splitting [0]
  -dj, --delim-join str            join string for kept tokens [same as --delim]
  -N, --normalize str [str ...]    normalization steps used for matching, or none [accents camel sep case]
  --accept-inferred                accept genus-abbreviation matches

Output Format:
  --root-edges str                 HyPhy and PAML read a rooted tree as unrooted and see the two branches at its root as one branch. This sets the label of that branch (unlabel, labeled, outgroup, ingroup, error).
                                   unlabel gives no label if the two branches differ, labeled uses the only label there is, outgroup and ingroup use the label of that side, error stops, and a label name wins if either
                                   branch has it [unlabel]
  --keep-root                      write the tree rooted as it is, without joining the two branches at the root
  --internal-names str             names written after internal nodes (none, auto, support, keep). none writes labels only, auto writes twn1, twn2 and so on, support writes the support values, and keep writes the names
                                   the tree had [none]
  --no-dists                       write trees without branch lengths
  --sanitize-suffix                for HyPhy, add _2, _3 to tip names that become equal when they are made safe
  --names-table path               for HyPhy, write the table of original and written tip names

Runtime:
  -f, --force                      overwrite existing output files
  -l, --log-level level            stderr logging level (DEBUG, [INFO], WARNING, ERROR)
```

## Examples from the help

The help ends with a list of examples, which are repeated here with the output that each command prints. They show every rule and every option on one small tree, and the tests run each command and compare its output, so the list cannot drift from the program. To follow along, make the two files with the first two commands.

```text
Examples
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
(A{X}:1,(B{Y}:1,C{Y}:1){Z}:1);
```

## Selecting tips

A selector is either a list file or a regular expression. A value that names an existing file is read as a list with one name per line, where a line that starts with `~` is a regular expression and a line that starts with `#` is a comment. Any other value is a regular expression that is matched with `re.search`, so `A` matches every name that contains an A. The option `--exact` treats the values given on the command line as exact names, and `--invert` labels the tips that the selectors do not match.

With `--select-on species` the selectors are matched to the species of each tip, which makes one selection work on gene trees that have several samples per species. The species come from the Taxon Mapping options, which are the same as in [name-map](name-map.md). For example `-d _` takes the part of each name before the first underscore. Without a map file, delimiter or regex, every tip is its own species, and the command says so.

Every label reports how many tips it selected, how many branches it covers, and how many keep its label in the end, after precedence, `--leaf-nodes skip` and `--stems`. A label that loses branches to another label, or to a label that the tree already had, says how many and which option decided it, so a label that was overridden never goes unnoticed. Names or patterns that match no tip are reported too, which helps to find typos.

## Several labels

Repeat `-L` for more labels. Each label can have its own rule, which is written after the label name with a colon, for example `-L Test:fitch dioecious.txt -L Out:tips '^O$'`. The labels can also be listed in a spec file given with `--spec`. It has one tab-separated line for each label with the columns label, rule, selectors and options, where the rule can be left empty and the options can be `ambiguous=`, `gain-cost=` and `loss-cost=`. A line that starts with `#` is a comment.

When two labels claim the same branch, `--precedence` decides. The default `specific` lets a label that a tip received by direct selection beat a label that a branch inherits from the clade around it, so a species that you select for its own label is never overwritten by the label of its clade. The value `first` lets earlier labels win, as HyPhy does, `last` lets later labels win, and `error` stops the run.

Labels that are already in the tree are kept by default and have priority over new labels, as in HyPhy. The option `--existing clear` removes them first. If you give no label at all, the labels in the tree are kept, so `tree-label` can also be used only to convert a labeled tree to another format, or to apply `--stems` and `--background` to it.

## Rules

The option `--rule`, also called `--internal-nodes`, decides which branches a label covers. The HyPhy names are accepted in any case, so `--internal-nodes "All descendants"` works.

| rule | HyPhy name | branches covered |
|---|---|---|
| `tips` | None | the selected tips only |
| `all` | All descendants | every branch all of whose descendant tips are selected |
| `crown` | All descendants, no MRCA | the result of `all` without the internal branches whose parent branch is not labeled, so isolated selected tips stay labeled |
| `some` | Some descendants | every branch with at least one selected tip below it |
| `fitch` | Parsimony | the branches on which two-state parsimony places the trait, see below. The name `parsimony` works too |
| `clade` | not in HyPhy | for each selector, the branch above the common ancestor of its tips and everything below it |
| `stem` | not in HyPhy | for each selector, only the branch above the common ancestor of its tips |

The option `--leaf-nodes skip` removes the labels of tips after the rule has run. With `fitch` the tips still provide the data. A tip whose label was removed in this way is left out of the background as well. The rule `tips` together with `skip` is an error because it would label nothing. The rules `clade` and `stem` are errors when the common ancestor is the root, because the root has no branch.

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

The name `fitch` needs a word of explanation. Parsimony means choosing the explanation that needs the fewest changes, and Fitch is Walter Fitch, who described this method in 1971. The rule `parsimony` means the same as `fitch`. Given the species that have a trait, it finds the fewest gains and losses of the trait on the tree that explain them. The `fitch` rule reconstructs a trait that is either absent or present, and labels every branch on which the trait is present. Selected tips have the trait, tips that are not selected do not, and tips given with `--unknown` cost nothing in either state. Gaining the trait costs `--gain-cost` and losing it costs `--loss-cost`, both 1 by default. The costs are kept as exact fractions, so equal totals are recognized without rounding errors, and a cost may be as large as you like. The two costs cannot both be 0, because every history would then cost nothing. The state of the root is free by default, and `--root-state 0` or `--root-state 1` fixes it. If the root has three or more children, the tree is most likely unrooted, as IQ-TREE and ASTRAL write their trees, and the result depends on where its root is drawn, whatever the costs. A branch is labeled by the state of the node below it, so a basal node that falls inside the group with the trait changes the stems, the losses and the labels around it. The command warns in that case, and `-r` with the outgroup roots the tree first.

Often two histories are equally good. On the species tree the trait can have been gained once on the branch above A, B and C and then lost on C, or it can have been gained separately on A and on B. Both histories cost 3 in total, which includes one gain on the branch above D and E. The option `--ambiguous` says what to do with the branches that differ between the equally good histories. The default `acctran` moves changes toward the root, which gives the first history. The value `deltran` delays changes toward the tips, which gives the second history. The value `unlabel` leaves those branches unlabeled, and a label name gives them a label of their own so that they can be handled separately later. The same tree under `deltran` has no label on the branches above A, B and C. A tip given with `--unknown` can be one of those branches, when its state differs between the equally good histories, and it is then treated like any other of them. The modes of `--ambiguous`, `--stems` and `--losses` are written in lower case, and a value such as `Deltran` or `Keep` is refused, so that a slip of the shift key does not turn a mode into a label name.

```bash
twig tree-label -i species.nwk -L Test dioecious.txt --rule fitch --ambiguous deltran
```

```text
(O,((M1,(A{Test},(B{Test},C))),(M2,(D{Test},E{Test}){Test})));
```

## Stems, losses and background

These three switches work with every rule and turn the result of a rule into the labeling that your question needs.

A stem is the internal branch on which a trait arose, so the change happened somewhere along its length. How it is found depends on the rule. With `fitch` a stem is a branch with a gain, which means the trait is present on it and absent in its parent, so a branch below the root is a stem only when the root lacks the trait. With `crown` a stem is the branch that the rule leaves out at the base of each group, so `--stems Stem` labels it, the background never takes it, and the first crown branch below it is not a stem. With the other rules a stem is a labeled internal branch whose parent branch does not carry the same label, and a branch below the root is a stem unless every child of the root carries the same label. A tip is never a stem, and neither is a branch with the label of `--ambiguous`. When `--stems` is unlabel or a name, the labels that are already in the tree count as well, as with the rules other than `fitch` and `crown`, so a tree labeled elsewhere, such as a PAML tree, can be turned into its crown or stem variant. twig cannot tell whether such a label marks a trait or whether its stems were handled before, so it names the stems it finds among them in a warning. The labels that the same run gives to the background, to losses and to nuisance species are never taken for traits. Use this on a tree whose stems have not been handled yet, because on the output of `--stems unlabel` the branch below each old stem would become a new stem. The option `--stems unlabel` leaves stems out of both the test and the reference set, and a name such as `--stems Stem` gives them a class of their own. With a class of their own, one labeled species tree supports both an analysis of the crown only and an analysis of stem plus crown, because the choice can be made later for each run.

A loss is a branch with the trait absent whose parent has it present, which is a reversion to the ancestral state. Only parsimony finds losses, so `--losses` acts on the result of `fitch` and not on labels that are already in the tree. The option `--losses unlabel` leaves them unlabeled, and `--losses Rev` gives them the label Rev. You do not need to know which species are reversions, because any species can be set aside as a nuisance species, as described in the next section.

The option `--background Reference` gives the label Reference to every branch that no label covers, except the stems, losses, ambiguous branches and skipped tips that you deliberately left unlabeled. This makes the reference branches explicit.

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

The option `-r` roots the tree on an outgroup before any label is applied. Give tips or species, joined by commas for a clade, and the first outgroup that forms one side of a branch is used. Labels stay on their branches when the tree is rerooted, and the two old root branches are joined by the same `--root-edges` value. An unrooted tree can be rooted in this way, and it should be, because every rule except `tips` decides by the tips that lie below a branch, which depend on the root. When the root has three or more children, every label whose rule is not `tips` gives a warning that suggests `-r`.

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
| `kind` | tip, stem, crown (a labeled branch with the same label as its parent, or below a stem of its label) or internal |
| `source` | where the label came from, which is the label and rule such as `Test:fitch`, or existing, stem, loss, ambiguous, background, or nothing if unlabeled |
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

The option `--show` prints one line for every labeled branch to stderr, in the order of the tree, followed by the unlabeled branches. Each line gives the label, the role of the branch and the tips below it, where X+Y means the branch above the common ancestor of X and Y. The roles are tip, stem, loss, ambiguous, crown and internal. The stems that were left without a label are listed as well, so that they are not missed.

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
