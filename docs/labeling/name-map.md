
# name-map

The `twig name-map` tool resolves the name of every tip in a set of gene trees, and the name of every sequence in a set of FASTA files, to a species of your species tree. It prints a table so that you can check the result by eye, and it can freeze the resolved mapping into a file that the other tree labeling tools read. Run it first, because every later step depends on each tip being assigned to the right species.

Gene trees often name the same species in several ways, for example `A_pumilus_sample1` in one gene and `AmaranthusPumilus-Sample1` in another, and they can hold several individuals of a species. The tool matches names to species exactly and also by tolerant rules, but it never assigns a species by similarity. A name that matches more than one species is reported as ambiguous, and a name that matches none is reported as unmapped with the closest species listed only as a suggestion.

The example files used on this page are in `tests/labeling/data/names`, which holds three gene trees, two alignments and a species tree of five species.

## The name-map command

```bash
twig name-map -h
```

```
usage: name-map -i TREES/FASTA [...] -s SPTREE [options]

-------------------------------------------------------------------
| name-map resolves tip and sequence names to species names
-------------------------------------------------------------------
| Reads gene trees and/or FASTA files and resolves every tip and
| sequence name to a species of the species tree (or of a species
| list). Names are matched exactly, after normalization, by whole
| tokens anywhere in the name, and by genus abbreviation. Matches
| by abbreviation are reported as inferred and are not accepted by
| other twig commands unless they are frozen into a map file with
| --freeze or allowed with --accept-inferred. The table of names is
| written to stdout and a summary to stderr.
-------------------------------------------------------------------

options:
  -h, --help                       show this help message and exit

Input/Output:
  -i, --input path [path ...]      trees and/or FASTA files, directories or quoted globs
  -s, --sptree path                species tree whose tips form the vocabulary
  -V, --vocab path                 file listing species names, one per line
  -o, --out path                   write the mapping table to this file [stdout]
  --freeze path                    write the resolved map for reuse with -I

Taxon Mapping:
  -I, --imap path                  map file rows <sample> <species> [specimen]; ~regex rows allowed
  -d, --delim str                  delimiter used to split sample names
  -x, --name-regex str             regex with (?P<species>...) and optional (?P<specimen>...)
  -di, --delim-idxs int [int ...]  token indices kept after splitting [0]
  -dj, --delim-join str            join string for kept tokens [same as --delim]
  -N, --normalize str [str ...]    normalization steps used for matching, or none [accents camel sep case]
  --accept-inferred                accept genus-abbreviation matches

Reporting:
  --show str                       rows of the table to print (all, problems) [all]
  --strict                         error if any name is unmapped, ambiguous or inferred, or leaves alphabetic leftover tokens, or if a tree or FASTA file repeats a name

Runtime:
  -f, --force                      overwrite existing output files
  -l, --log-level level            stderr logging level (DEBUG, [INFO], WARNING, ERROR)
```

## Inputs

The input to `-i` is any mix of tree files, FASTA files, directories and globs. Quote globs, as in `'trees/*.treefile'`, so that twig expands them itself, which keeps thousands of files from hitting the argument limit of the shell. A file is recognized as a tree file or a FASTA file from its first lines and not from its extension. A file that is neither, such as a PHYLIP alignment, is skipped with a warning when it comes from a directory or a glob, and is an error when you name it directly. Directories and globs leave out hidden files and files whose names end in `.tmp`, which are the partial outputs of a twig run that was killed.

A tree file can hold one tree or many trees, which can span several lines. Blank lines, comment lines that start with `#` and a PAML header with the number of species and trees are skipped. Trees can carry labels, support values, branch lengths and annotations, all of which are ignored here. NEXUS files and quoted names are not supported.

The species come from the tips of the species tree given to `-s`, from a text file with one species per line given to `-V`, or from both. Labels on the species tree are ignored. Without either option the tool still lists the names, but it cannot match them to species. A byte order mark at the start of a file, which some Windows editors add, is ignored in every file that the tool reads, including map files.

## How names are matched

For each name, the first rule that gives an answer is used.

1. A row of a map file that names the sample exactly.
2. A regular expression row of a map file.
3. A species name parsed from the sample name with `-d` or `-x`, matched to the species tree by the rules below.
4. The whole sample name, matched to the species tree by the rules below.
5. Otherwise the name is unmapped.

The matching of a name to the species tree tries four steps and stops at the first that finds a species.

| method | rule | example |
|---|---|---|
| `exact` | the name is a species name | `Beta_vulgaris` |
| `canon` | the names are equal after normalization | `AmaranthusPumilus` and `Amaranthus_pumilus` |
| `token` | the words of a species name appear next to each other anywhere in the name | `s012_Amaranthus_hybridus`, `og1 amaranthus pumilus g7` |
| `abbrev` | like `token`, except that the genus may be a prefix of up to three letters | `A_hybridus_s1` for `Amaranthus_hybridus` |

When several species match in the same step, the one that covers the most words wins, so `A_hybridus_quitensis_s1` goes to `A_hybridus_quitensis` when both species exist. When several species tie, the name is ambiguous. For example `A_pumilus_sample1` is ambiguous between `Amaranthus_pumilus` and `Atriplex_pumilus`, while `AmaranthusPumilus-Sample1` still resolves to the first.

A match by abbreviation is only a guess about a genus, so it gets the status `inferred`. It is shown in the table, but other tree labeling commands do not accept it unless you freeze the mapping after reviewing it, or pass `--accept-inferred`. This way a sample of Atriplex that is missing from the species tree is never silently assigned to Amaranthus.

The words that are left over after the species part of a name are reported in the `leftover` column. They are typically a specimen identifier such as `sample1`, but they can also expose a variety or a subspecies that is being lumped into a species, such as the `maritima` in `Beta_vulgaris_maritima_s1`.

## Normalization

Normalization is used only to compare names. The names that are written to the table and to the map file are never changed. The steps are applied in the order of the table below, whatever order they are given in.

| step | effect | example |
|---|---|---|
| `accents` | accents are removed | `Ñandú` becomes `Nandu` |
| `camel` | a separator is added between a lower case and an upper case letter, and at the end of an acronym | `AmaranthusPumilus` becomes `Amaranthus_Pumilus` |
| `digits` | a separator is added between letters and digits, off by default | `Sample1` becomes `Sample_1` |
| `sep` | every run of characters that are not letters or digits becomes one `_` | `A.pumilus-s1` becomes `A_pumilus_s1` |
| `case` | upper case letters become lower case | `Amaranthus_Pumilus` becomes `amaranthus_pumilus` |

The default is `accents camel sep case`, and `-N none` turns normalization off, so that only `exact` matches can succeed.

## Reading the table

The table is written to stdout, or to the file given to `-o`, with one row for each distinct name.

```bash
cd tests/labeling/data/names
twig name-map -i trees aln -s species.nwk
```

The output is tab separated, and it is shown here with aligned columns.

```
name                       species              specimen  method  status     leftover     candidates                                                               n_files
Amaranthus_pumilus         Amaranthus_pumilus   -         exact   ok         -            -                                                                        2
Amaranthus_hybridus        Amaranthus_hybridus  -         exact   ok         -            -                                                                        2
Atriplex_pumilus           Atriplex_pumilus     -         exact   ok         -            -                                                                        3
Beta_vulgaris              Beta_vulgaris        -         exact   ok         -            -                                                                        3
Spinacia_oleracea          Spinacia_oleracea    -         exact   ok         -            -                                                                        5
AmaranthusPumilus-Sample1  Amaranthus_pumilus   -         token   ok         sample1      -                                                                        2
AmaranthusPumilus-Sample2  Amaranthus_pumilus   -         token   ok         sample2      -                                                                        2
s012_Amaranthus_hybridus   Amaranthus_hybridus  -         token   ok         s012         -                                                                        2
Atriplex_pumilus_x1        Atriplex_pumilus     -         token   ok         x1           -                                                                        2
Beta_vulgaris_maritima_s1  Beta_vulgaris        -         token   ok         maritima,s1  -                                                                        2
A_pumilus_sample1          -                    -         abbrev  ambiguous  -            Atriplex_pumilus(0.545);Amaranthus_pumilus(0.514)                        1
A_hybridus_s1              Amaranthus_hybridus  -         abbrev  inferred   s1           -                                                                        1
Zea_mays_s1                -                    -         none    unmapped   -            Beta_vulgaris(0.417);Amaranthus_pumilus(0.414);Amaranthus_hybridus(0.4)  1
```

| column | content |
|---|---|
| `name` | the raw name of a tip or a sequence |
| `species` | the species it was assigned to, or `-` if it was not assigned |
| `specimen` | the specimen, shown only when a map row or a regex group gave one explicitly |
| `method` | how the name was resolved, one of `map`, `map-regex`, `delim`, `regex`, `exact`, `canon`, `token`, `abbrev`, `self` and `none` |
| `status` | `ok`, `inferred`, `ambiguous` or `unmapped` |
| `leftover` | the words of the name that are not part of the species, separated by commas |
| `candidates` | for names that are not resolved, the species that tie or the three most similar species with their similarity, which is only a suggestion |
| `n_files` | the number of input files in which the name occurs |

The method `self` is used when there is no species tree or species list and no parsing rule, in which case a name is its own species. The option `--show problems` prints only the rows that need attention, which are the names that are not `ok` and the names with a leftover word made only of letters, such as `maritima`. Leftover words with digits, such as `sample1` or `s012`, are usually specimen identifiers and do not count. Neither do the leftover words of a name that a map row resolves, because the row says what the name means.

```bash
twig name-map -i trees aln -s species.nwk --show problems
```

```
name                       species              specimen  method  status     leftover     candidates                                                               n_files
Beta_vulgaris_maritima_s1  Beta_vulgaris        -         token   ok         maritima,s1  -                                                                        2
A_pumilus_sample1          -                    -         abbrev  ambiguous  -            Atriplex_pumilus(0.545);Amaranthus_pumilus(0.514)                        1
A_hybridus_s1              Amaranthus_hybridus  -         abbrev  inferred   s1           -                                                                        1
Zea_mays_s1                -                    -         none    unmapped   -            Beta_vulgaris(0.417);Amaranthus_pumilus(0.414);Amaranthus_hybridus(0.4)  1
```

## Reading the summary

A summary is written to stderr, which you can keep with `2> summary.log`. For the example above it reads as follows, where the date and the name of the source file that begin each message are left out.

```
resolved 13 unique names from 5 files (3 trees, 2 FASTA files)
matched by abbrev 2, exact 5, none 1, token 5
status ambiguous 1, inferred 1, ok 10, unmapped 1
1 name is ambiguous, matching more than one species. A_pumilus_sample1 (Atriplex_pumilus, Amaranthus_pumilus)
1 name is unmapped, matching no species. Zea_mays_s1 (Beta_vulgaris, Amaranthus_pumilus, Amaranthus_hybridus)
1 name is inferred, matched only by genus abbreviation (review and freeze with --freeze). A_hybridus_s1 (Amaranthus_hybridus)
1 word is left after the species name in the names of Beta_vulgaris. maritima
Amaranthus_pumilus has up to 2 samples in one tree or FASTA file, which happens in 2 of the 4 where it occurs
```

Besides the counts and the names that need attention, the summary names every tree or FASTA file in which a name occurs more than once, which means that the file is broken, and it lists species of the species tree that have no samples, the leftover words made only of letters for each species, except in names that a map row resolves, and the species that have more than one sample in a single tree or file. The last two help to find a delimiter that splits names too finely or too coarsely, and a subspecies that is lumped into its species.

## Parsing names

Names that follow a pattern can be parsed before they are matched. The option `-d` splits a name on a delimiter and keeps the pieces listed in `-di`, joined by `-dj`, which defaults to the delimiter itself. The tool matches the parsed name to the species tree first and falls back to the whole name if the parsed name does not match. A name that cannot be split as asked is not an error, it just goes straight to the whole name.

```bash
twig name-map -i genes.nwk -s species.nwk -d - -di 1
```

```
name                      species             specimen  method  status  leftover  candidates  n_files
OG1-AmaranthusPumilus-s1  Amaranthus_pumilus  -         canon   ok      og1,s1    -           1
OG1-BetaVulgaris-s1       Beta_vulgaris       -         canon   ok      og1,s1    -           1
OG1-Atriplex_pumilus-s2   Atriplex_pumilus    -         exact   ok      og1,s2    -           1
```

The option `-x` does the same with a regular expression that has a named group `species` and, optionally, a named group `specimen`. A name that the expression does not match falls through to the whole name. The options `-d` and `-x` cannot be used together.

```bash
twig name-map -i genes.nwk -s species.nwk -x '^(?P<species>[A-Za-z]+_[a-z]+)_(?P<specimen>s[0-9]+)$'
```

```
name                  species            specimen  method  status  leftover  candidates  n_files
Beta_vulgaris_s1      Beta_vulgaris      s1        exact   ok      s1        -           1
Beta_vulgaris_s2      Beta_vulgaris      s2        exact   ok      s2        -           1
Spinacia_oleracea_s1  Spinacia_oleracea  s1        exact   ok      s1        -           1
```

## Map files

When names cannot be matched by rules, or when you have decided what a name means, give the answer in a map file with `-I`. Each row has up to three columns separated by whitespace, which are the sample, the species and optionally the specimen. A `#` that starts a word begins a comment, and blank lines are ignored. A row with one column maps a name to itself, which is how the IMAP files of `tree-filter` work. A sample that starts with `~` is a regular expression that is applied with `re.search`, and these rows are tried in the order of the file after all the rows that name a sample exactly.

```
# sample                    species               specimen (optional)
A_pumilus_sample1           Amaranthus_pumilus    s1
AmaranthusPumilus-Sample2   Amaranthus_pumilus    s2
~^Ahyb[_-]                  Amaranthus_hybridus
Beta_vulgaris
```

A map row takes precedence over every other rule. A row that names a sample exactly is compared with the whole name and also with the name parsed by `-d` or `-x`, so the IMAP files of `tree-filter` work unchanged. When a species tree or list is given, the species in each row must be one of its species, or the tool stops and names the file and the line.

## Freezing a mapping

After you have reviewed the table, `--freeze` writes the resolved names as a map file, so that every later run uses exactly the mapping you reviewed.

```bash
twig name-map -i trees aln -s species.nwk --freeze names.tsv
```

```
# map file written by twig name-map
# <sample> <species> [<specimen>]
Amaranthus_pumilus Amaranthus_pumilus
Amaranthus_hybridus Amaranthus_hybridus
Atriplex_pumilus Atriplex_pumilus
Beta_vulgaris Beta_vulgaris
Spinacia_oleracea Spinacia_oleracea
AmaranthusPumilus-Sample1 Amaranthus_pumilus
AmaranthusPumilus-Sample2 Amaranthus_pumilus
s012_Amaranthus_hybridus Amaranthus_hybridus
Atriplex_pumilus_x1 Atriplex_pumilus
Beta_vulgaris_maritima_s1 Beta_vulgaris
A_hybridus_s1 Amaranthus_hybridus  # inferred by genus abbreviation, review
# AMBIGUOUS A_pumilus_sample1 Atriplex_pumilus(0.545);Amaranthus_pumilus(0.514)
# UNMAPPED Zea_mays_s1 Beta_vulgaris(0.417);Amaranthus_pumilus(0.414);Amaranthus_hybridus(0.4)
```

Names with the status `ok` or `inferred` are written. Rows that came from an abbreviation carry a comment as a reminder to review them, and the comment can stay. Ambiguous and unmapped names are listed in comments, and you resolve them by adding a row, for example `A_pumilus_sample1 Amaranthus_pumilus`. The specimen column is written only when a map row or a regex group gave a specimen, because otherwise a reloaded map would make every sample its own specimen and would silently switch off the detection of duplicate samples. A name that a map row resolves keeps its leftover words, found by matching the name to the species of the row, so a sample has the same leftover and specimen before and after freezing.

Running the same command again with the map file leaves only the names you have not yet resolved.

```bash
twig name-map -i trees aln -s species.nwk -I names.tsv --show problems
```

An existing output file is never overwritten unless you pass `-f`. Even with `-f` an output is never written over one of the input files, the species tree, the species list or the map file, and never over a path that is not a regular file, such as a directory or a device like `/dev/stdout`. The table goes to stdout when `-o` is left out.

## Strict mode

With `--strict` the tool stops with an error and exit status 1 if any name is unmapped, ambiguous or inferred, or has a leftover word made only of letters as described for `--show problems`, or if a tree or a FASTA file holds a name more than once. The table is still written, so you can see what needs attention, but no map file is frozen. This is useful in a script that should stop before a long analysis starts.

## Large data sets

The tool reads only the tip names of trees and the headers of FASTA files. It never reads sequences or builds tree objects, so thousands of loci take seconds. The table lists each distinct name once, with the number of files in which it occurs.
