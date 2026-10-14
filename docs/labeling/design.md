# Tree labeling for twig, vision and implementation plan

This document is self-contained. Part A gives any reader, human or AI, the background needed to understand why this work exists and what has already been decided. Part B fixes the boundaries of the work inside the existing twig code. Parts C onward give the design and the technical details.

---

# Part A. Background

## A.1 The repository

- **What twig is.** twig ("tree-based workflows for integrative genomics") is a Python command-line toolkit for comparative and phylogenomic analysis, written by Deren Eaton's lab (eaton-lab/twig). This work happens in a fork, `dt2733/twig`, on the branch `tree-labeling`. The fork's default branch is `main` and the package version is 0.0.5.
- **Layout.** Code lives in `src/twig/` in subpackages `cli` (the command line), `seqlab` (BLAST, alignment, trimming, concatenation), `treelab` (gene tree tools), `synteny` (genome file formatting), `hgt` (unfinished) and `utils` (logging, errors, a process pool).
- **How a command runs.** `src/twig/cli/main.py` holds a `DISPATCH` dictionary mapping each subcommand name to a subpackage, for example `"tree-filter": "..treelab"`. At startup `main` builds every subcommand's argument parser by calling `get_parser_<name>` from `src/twig/cli/subcommands.py`. When a command is chosen, it imports `twig.<subpackage>.<name>` lazily and calls `run_<name>(args)`. `TwigError` (from `twig.utils`) is the error type for user mistakes and is printed without a traceback.
- **Tree engine.** All tree work uses toytree (Deren Eaton's tree library, `toytree>=3.0.10` in `pyproject.toml`). The existing tree tools are `tree-filter`, `tree-rooter` and `tree-skeleton`. They parse sample names into taxon names with `--delim`, `--delim-idxs`, `--delim-join` and an `--imap` file, and that parsing is duplicated across several modules.
- **State of the code.** There are no tests. Documentation is a mkdocs-material site in `docs/`, deployed by a GitHub Action on pushes to `main`. Several small bugs exist in other tools (listed in the appendix). They are not part of this work.
- **Why the fork matters.** Keeping existing files untouched keeps it easy to pull future changes from eaton-lab/twig without merge conflicts. Part B is built around that.

## A.2 The scientific problem

- **The study.** The user studies the evolution of a sexual system (dioecy versus monoecy) in a plant clade. The species tree contains two independent origins of dioecy, and inside one dioecious clade a single species has reverted to monoecy. The data are thousands of single-copy orthologs, with one individual per species, and a gene tree for each.
- **The question.** Did natural selection change in association with dioecy, for example relaxed or intensified purifying selection, or episodic positive selection? The tools are HyPhy's RELAX (estimates a selection-intensity parameter k on test branches relative to reference branches), BUSTED-PH (asks whether episodic positive selection is associated with foreground branches), aBSREL (tests branches one at a time) and BUSTED.
- **Why labels are needed.** These programs take a Newick tree in which branches carry labels in curly braces, for example `A{Test}:0.1`. Each branch is assigned to a set by its label. RELAX allows a third, unclassified set of branches that is fitted but left out of the comparison. The built-in BUSTED-PH has no such set, so every branch not in the foreground counts as background. (A legacy BUSTED-PH script had a "nuisance" set, but it was dropped when the method was built into HyPhy.)
- **Gene tree discordance.** The trait's history is known on the species tree, but gene trees often differ from it, through incomplete lineage sorting or estimation error. Forcing every gene onto the species tree misplaces substitutions and inflates false positives (Mendes and Hahn). Running on each gene's own tree avoids that, but then each gene tree needs its own labels.
- **The chosen solution, species-tree branch matching.** The biologist labels the species tree once. For each gene, a gene-tree branch takes the label of the species-tree branch that has exactly the same set of species below it. Branches with no counterpart stay unlabeled and drop out of the comparison. Tips always match. A stem (the branch leading to a trait clade) is labeled only in genes whose tree recovers that clade, which reproduces the filter used in the Jaltomata study (Wu and colleagues) and the logic of the Hahn lab's hyphy-interface. No published RELAX study is known to use this rule, so twig would be the first general implementation.
- **Stems.** The transition to dioecy happened somewhere along each stem, at an unknown point. The stem can be left unlabeled (crown only), or labeled. Running both versions and comparing them separates genes associated with the origin from genes associated with the maintained state.
- **Reversions.** The reverted species should be neither test nor reference. In RELAX it can be unclassified. In BUSTED-PH, which has no unclassified set, it has to be pruned from the tree and alignment.
- **HyPhy's own labeling tool.** `hyphy label-tree` labels one tree from a list or regex of tips, with internal-node rules "None", "All descendants", "All descendants, no MRCA", "Some descendants" and "Parsimony". Reading its source showed that "Parsimony" behaves exactly like "All descendants", because unselected tips are never given a state. It also cannot transfer labels between trees.
- **Names and samples.** Gene trees in general may contain several individuals per species, with inconsistent naming such as `A_pumilus_sample1` and `AmaranthusPumilus-Sample1`, and may contain paralogs (several gene copies per individual). To be consistent with a species tree, trees can be reduced to one sample per species in three ways. (a) Prune the existing tree, so the extra specimens still informed the topology. (b) Remove the extra rows from the existing alignment and re-estimate the tree. (c) Remove the specimens, realign and re-estimate. A fourth way is to collapse conspecific samples that form a clade into one tip.

## A.3 Decisions already made by the user

These are settled and should not be reopened.

1. All tree manipulation uses toytree.
2. HyPhy is the first target, but labeling must be general enough for other programs that need labeled trees.
3. Species-tree branch matching is the main mode. HyPhy label-tree style rules are also provided.
4. Low-support branches are never collapsed. HyPhy handles uncertainty in its own estimates, and collapsing deletes branches together with their labels.
5. For each category of branch, the user can choose to leave it unlabeled, relabel it, or prune its species from the tree (and alignment).
6. Name mapping must cope with inconsistent formats and multiple samples per species.
7. All three reduction strategies, plus collapse, are offered, with several rules for which sample to keep.
8. The new work stays in its own module and touches existing twig files only where unavoidable (Part B).
9. Writing style for all prose, help text and docs is to use no em-dashes and to avoid colons in prose sentences (code, Newick and command examples may contain colons).

## A.4 Environment facts

- toytree and HyPhy are not installed in the cloud container. Tests install toytree with `pip install toytree`. HyPhy behavior is taken from its source code and documentation, and an optional test runs HyPhy only when the binary is present.
- The container runs Python 3.13. The existing `subcommands.py` already needs Python 3.10 or newer (it evaluates `ArgumentParser | None`), although `pyproject.toml` says 3.8. The new code also needs 3.10. `pyproject.toml` is not changed; the docs state the requirement.
- The research behind this plan read twig's code, toytree's source (3.0.10, 3.0.11 and master), and HyPhy's source (`label-tree.bf`, `trees.bf`, `RELAX.bf`, `BUSTED.bf`, `BUSTED-PH.bf`, `alignments.bf`). Two independent adversarial reviews of a first draft were folded in.

## A.5 Glossary

| term | meaning |
|---|---|
| branch, edge | the line between a node and its parent. A branch's label is stored on the node below it |
| label | text attached to a branch, such as `Test` or `Reference`, written `{Test}` for HyPhy |
| test, foreground | the branches whose selection regime is being tested (dioecious lineages) |
| reference, background | the branches compared against (monoecious lineages). In BUSTED-PH every non-foreground branch is background |
| unclassified, nuisance | branches fitted but excluded from the comparison (RELAX only) |
| stem | the branch leading to a trait clade, along which the transition happened |
| crown | the branches inside a trait clade |
| reversion, loss | a return to the ancestral state inside a trait clade |
| split, bipartition | the division of a tree's tips into the two sides of a branch. Splits identify branches without needing a root |
| discordance | a gene tree branch whose split does not exist in the species tree |
| sample, specimen, species | a tip or alignment row; one individual; one tip of the species tree |
| paralog, in-paralog | a gene copy from a duplication; a duplication confined to one species |

---

# Part B. Boundaries inside twig

## B.1 Principle

Everything new lives in three new directories. Existing twig code is used through imports where it is stable infrastructure, and is otherwise left alone.

```
src/twig/labeling/      all new code (commands and library)
tests/labeling/         all new tests and test data
docs/labeling/          all new documentation, including this design document
```

No new dependencies are added. The code uses the standard library, toytree and loguru, which twig already requires (numpy is used only where twig already uses it).

## B.2 Every existing file this work touches

| file | change | why it cannot be avoided |
|---|---|---|
| `src/twig/cli/main.py` | four lines added to `DISPATCH`, for example `"tree-label": "..labeling",` | `main` only runs subcommands listed there |
| `src/twig/cli/subcommands.py` | two lines appended at the very end, a comment and `from twig.labeling.parsers import get_parser_name_map, get_parser_sample_reduce, get_parser_tree_label, get_parser_tree_label_transfer  # noqa: E402,F401` | `main` looks up `get_parser_<name>` as attributes of this module. Appending at the end means the module is fully defined before the import runs, so there is no circular-import risk |
| `mkdocs.yml` (optional, only if the user agrees) | one nav block listing `docs/labeling/` pages | without it, mkdocs still builds the pages but leaves them out of the site menu |

Nothing else changes. In particular `pyproject.toml` (setuptools finds the new subpackage automatically), `README.md`, the GitHub workflow, and every file in `utils/`, `treelab/`, `seqlab/`, `synteny/` and `hgt/` stay as they are.

## B.3 Existing code imported read-only

Only stable infrastructure is imported.

| import | used for |
|---|---|
| `twig.utils.TwigError` | the error type `main` already handles |
| `twig.utils.logger_setup.set_log_level`, `setup_loguru_worker` | the same logging behavior as every other command |
| `twig.utils.parallel.run_with_pool_iter`, `safe_popen` | parallel genes, and running external commands with clean Ctrl-C handling |
| `twig.cli.subcommands.SingleMetavarHelpFormatter` | the same help layout as the newer twig commands, imported inside each `get_parser_*` function |

## B.4 Logic reimplemented inside `labeling` instead of shared

Some existing functions do part of what is needed, but changing them would affect other commands, and importing domain logic couples this work to code that may change upstream. These are rewritten inside `labeling/core/`, with a comment naming the twig function each one mirrors.

| need | existing function it mirrors | reason not to reuse |
|---|---|---|
| strict FASTA reading | `seqlab/filter_concat.parse_fasta_alignment` | five different FASTA parsers exist with different behavior, and none should change |
| imap and delimiter parsing | `filter_concat.parse_imap`, `parse_taxon`, tree-filter's imap loop | tree-filter's version has bugs, and the new mapper needs more features |
| ordered rooting clades | `treelab/tree_rooter.get_rooting_clades` | about twenty lines, and moving it would edit tree-rooter |
| presence and overlap matrix | `seqlab/align_post.present_matrix`, `overlap_matrix` | small, and align-post stays untouched |

The duplication is deliberate. A later, separate effort could offer the shared versions to the rest of twig (appendix), but that is outside this plan.

## B.5 Safety rules for living inside twig

1. **Startup cost.** `main` builds every parser on every `twig` call. So `twig/labeling/__init__.py` stays empty and `labeling/parsers.py` imports only `argparse`, `pathlib` and `textwrap`. toytree is imported only inside the four run modules and `labeling/core/`. A test asserts that importing `twig.cli.subcommands` does not load toytree.
2. **Startup failures.** A bug in one of the four `get_parser_*` functions would break every twig command. A test builds the full `twig` parser and runs `twig -h`.
3. **Flags.** argparse subparsers keep separate option namespaces, so the new commands' flags cannot collide with flags of other commands.
4. **Names.** The four command names are checked against the current `DISPATCH` and against upstream eaton-lab/twig before they are added.
5. **Removal.** Removing the feature means deleting the three directories and the six added lines.

---

# Part C. What the user gets

## C.1 Four commands

| command | purpose |
|---|---|
| `twig name-map` | resolve every tip and sequence name to a species, print the table for checking, and freeze it for reuse |
| `twig tree-label` | label one tree (usually the species tree) from tip selections with HyPhy-style rules, real parsimony, and explicit stem, loss and background handling |
| `twig tree-label-transfer` | copy the labeled species tree's labels onto each gene tree by branch matching, apply masks (unlabel, relabel, prune), validate against the target program, and write each tree with its matching alignment |
| `twig sample-reduce` | reduce gene trees and alignments to one specimen per species by pruning, row removal plus re-estimation, realignment plus re-estimation, or collapse of monophyletic conspecific samples |

## C.2 The user's analysis, end to end

The species tree from the user's figures is `(O,((M1,(A,(B,C))),(M2,(D,E))));`. O is the outgroup, M1 and M2 are monoecious, A, B, D and E are dioecious, and C is the reversion. A file `dioecious.txt` lists A, B, D and E.

```bash
# 1. label the species tree once, stems unlabeled, reversion C labeled Rev, everything else Reference
twig tree-label -i sp.nwk -L Test dioecious.txt --rule fitch --stems unlabel --losses Rev --background Reference --show > sp.labeled.nwk

# 2. confirm every gene tip resolves to a species tip
twig name-map -i 'trees/*.treefile' 'aln/*.nt.fa' -s sp.labeled.nwk --show problems

# 3. dry run, reports only
twig tree-label-transfer -s sp.labeled.nwk -i 'trees/*.treefile' -a 'aln/*.nt.fa' -P relax -O relax/ --check

# 4a. RELAX inputs (Rev joins the Unclassified set)
twig tree-label-transfer -s sp.labeled.nwk -i 'trees/*.treefile' -a 'aln/*.nt.fa' -P relax -O relax/

# 4b. BUSTED-PH inputs (no Unclassified set, so the reversion species is pruned from trees and alignments)
twig tree-label-transfer -s sp.labeled.nwk -i 'trees/*.treefile' -a 'aln/*.nt.fa' -P busted-ph -m Rev=prune -O bph/

# 5. run HyPhy with the generated command file
parallel < relax/hyphy_commands.sh
```

Step 1 writes the labeled tree from the user's first figure.

```
(O{Reference},((M1{Reference},(A{Test},(B{Test},C{Rev}){Test})){Reference},(M2{Reference},(D{Test},E{Test})){Reference}){Reference});
```

The A+B+C stem and the D+E stem have no braces, so they are unlabeled. `--show` prints every labeled branch with the species below it so the result can be checked by eye. To compare stem plus crown against crown only, the stems can instead get their own label (`--stems Stem`), and each transfer run then chooses with `-m Stem=Test` or `-m Stem=unlabel`.

## C.3 What "unlabeled" means in each program

This drives every default below.

| program | test set | reference set | unlabeled branches |
|---|---|---|---|
| RELAX, Classic mode with `--test` and `--reference` | test label | reference label | Unclassified, with their own ω distribution, outside the comparison |
| RELAX with only `--test` | test label | every other branch | none, everything else silently becomes Reference |
| BUSTED, BUSTED-PH | `--branches` label | every other branch | none, everything else is background |
| aBSREL | `--branches` label, each tested separately | none | modeled but not tested |
| PAML codeml | `#1` (and `#2` and so on) | `#0` | `#0` |

Leaving a branch unlabeled is neutral in RELAX but puts it in the background in BUSTED-PH. For every preset, twig reports the branches and branch length that land in the background without being Reference by choice (`n_uncertain_bg`, `len_uncertain_bg`).

---

# Part D. Concepts

## D.1 Labels live on branches

- **Storage.** A branch's label is the toytree node feature `label` on the node below the branch. `None` means unlabeled. The root has no branch and never carries a label.
- **Rerooting.** `label` is registered in `tree.edge_features` after every read, prune and copy, so `tree.root(..., edge_features=["label"])` moves labels correctly when the root moves.
- **Root edges.** A rooted tree has two branches at the root. HyPhy and PAML read every tree as unrooted, so they see those two branches as one. twig always merges them itself with an explicit policy (F.4), rather than letting HyPhy pick one silently.

## D.2 Branch statuses

| status | meaning |
|---|---|
| `labeled` | matched a species-tree branch that carries a label |
| `species-unlabeled` | matched a species-tree branch deliberately left unlabeled (for example a stem) |
| `unmatched` | its species split does not exist in the species tree (discordance) |
| `compatible` | it resolves a species-tree polytomy without contradicting it |
| `conflict` | it matched a species-tree branch formed by merging branches whose labels disagree |
| `mixed` | samples of one species lie on both sides of it |
| `paralog` | it lies at or below a duplication node and was masked |
| `unmapped` | its tip could not be mapped to a species and the run was told to keep it |
| `masked` | a mask changed its label (the original status is kept in the branch table) |

## D.3 Samples, specimens and species

- **Sample.** A tip of a gene tree or a row of an alignment, under its raw name.
- **Specimen.** One individual. By default it is the canonical leftover after the species part of the name is removed (`sample1` in both `A_pumilus_sample1` and `AmaranthusPumilus-Sample1`). A regex group or a map-file column can set it explicitly.
- **Species.** A tip of the species tree. Mapping always resolves to this vocabulary.

---

# Part E. Architecture

## E.1 Module layout

```
src/twig/labeling/
  __init__.py              empty, so importing the package is free
  parsers.py               get_parser_* for the four commands (argparse only, no toytree)
  name_map.py              run_name_map             (module names must equal the command names,
  tree_label.py            run_tree_label            because main imports twig.labeling.<name>
  tree_label_transfer.py   run_tree_label_transfer   and calls run_<name>)
  sample_reduce.py         run_sample_reduce
  core/                    library code, no argparse anywhere
    __init__.py            lazy re-exports
    constants.py           feature name, status names, registry name lists for help text (no toytree)
    model.py               Status, MaskRule, BranchRecord, GeneResult, config dataclasses
    names.py               NameMapper, map files, normalization, vocabulary matching (no toytree)
    fasta.py               strict FASTA reader and writer, row subsetting, empty columns, occupancy, codon checks
    selectors.py           list-file-or-regex tip selectors
    batch.py               glob and directory expansion, locus naming, pairing, manifests, multi-tree files,
                           atomic writes, resume bookkeeping
    treeio.py              label pre-tokenizer, read_tree, iter_trees, strip_labels, support handling
    newick.py              iterative Newick serializer used by every writer
    rooting.py             rooting clades, species-consistent rooting, root-edge merge planning
    rules.py               internal-node rule registry
    parsimony.py           two-state Sankoff reconstruction
    transfer.py            SpeciesIndex, restriction to present species, branch classifier
    paralogs.py            duplication detection and policies
    mask.py                mask parsing and application
    sanitize.py            HyPhy and PAML name rules, global reversible name tables, label-name checks
    writers.py             writer registry
    presets.py             program presets, requirements and command templates
    validate.py            Issue type, checks for trees, alignments and presets
    report.py              gene report, branch table, species-tree recovery table, stderr summary
    display.py             text listing of labeled branches by clade (--show)
    keep.py                keep-one selector registry
    reduce.py              reduction strategies and external command templates
```

## E.2 Conventions followed

- **Parsers.** The style of the newer twig commands (`filter-concat`, `align-pre`) is used, meaning `SingleMetavarHelpFormatter(prog, width=220, max_help_position=56)`, argument groups, `[%(default)s]` in help strings, 67-dash description frames, and an `Examples` epilog. Two helpers inside `parsers.py`, `add_taxon_mapping_args` and `add_writer_args`, keep shared groups identical across the four commands. Choice names are validated in `run_*` against the registries (E.4), because `parsers.py` must not import the library.
- **Errors.** User mistakes raise `TwigError`. Per-gene problems never stop a batch. They skip the gene and record the reason, unless `--on-fail error` is given.
- **Logging.** `run_*` calls `set_log_level(args.log_level)` and then `logger.disable("toytree")`, because toytree enables its own loguru sink on import. Each run logs a `CMD:` line and a summary to stderr.
- **Outputs.** Single-stream output goes to stdout or `-o`. Per-locus outputs go to `-O/--outdir`, and reports are then written there by default.
- **Overwriting and resuming.** An existing output raises `TwigError` unless `-f/--force`. A separate `--resume` skips loci whose outputs are complete. Outputs are written to `*.tmp` and renamed with `os.replace`, with the tree written last, so a killed run never leaves a file that looks finished. `OUTDIR/twig-run.json` records the command, a hash of the settings and of the species tree, and the twig and toytree versions, and `--resume` refuses to continue under different settings.
- **Library and command separation.** `run_*` turns `args` into frozen config dataclasses. The library never sees argparse, so it can be used from Python and tested directly.

## E.3 Core data structures (`core/model.py`)

```python
LABEL = "label"

class Status(str, Enum):
    LABELED = "labeled"; SPECIES_UNLABELED = "species-unlabeled"; UNMATCHED = "unmatched"
    COMPATIBLE = "compatible"; CONFLICT = "conflict"; MIXED = "mixed"; PARALOG = "paralog"
    UNMAPPED = "unmapped"; MASKED = "masked"
    __str__ = str.__str__          # TSVs and f-strings write "labeled", not "Status.LABELED"
    __format__ = str.__format__

@dataclass(frozen=True)
class MaskRule:
    cls: str                       # a label name or a status name
    action: str                    # keep | unlabel | prune | relabel
    target: str | None = None

@dataclass
class BranchRecord:
    locus: str; node: str; is_tip: bool; n_tips: int; tips: tuple[str, ...]
    species: tuple[str, ...]; split: str; label: str | None; status: Status
    source: str; sp_branch: str | None; origin: str | None
    dist: float; support: float | None; written: bool = True; merged_into: str | None = None; note: str = ""

@dataclass
class GeneResult:
    locus: str; status: str; reason: str
    outputs: dict[str, str]         # writer name -> rendered text (picklable for worker processes)
    alignment: dict[str, str] | None
    records: list[BranchRecord]; stats: dict[str, Any]; issues: list["Issue"]

@dataclass(frozen=True)
class TransferConfig:
    match: str = "split"; merge: str = "unanimous"; merge_tips: str = "tipward"
    overlap: str = "strict"; intraspecific: str = "inherit"
    root: str = "auto"; root_fail: str = "unrooted"; outgroups: tuple[tuple[str, ...], ...] = ()
    repeats: str = "error"; paralogs: str = "skip"; paralog_tips: str = "unlabel"
    masks: tuple[MaskRule, ...] = (); protect: tuple[str, ...] = (); root_edges: str = "unlabel"
    min_tips: int = 4; min_species: int = 3
    @classmethod
    def from_args(cls, args) -> "TransferConfig": ...
```

## E.4 Registries

Each extension point is a dictionary filled by a decorator. Adding a rule, merge policy, keep-one selector, writer or preset means writing one function plus its tests, without touching the commands.

```python
RULES: dict[str, Rule]            # core/rules.py        register_rule(name, *aliases, help="")
MERGES: dict[str, MergeFunc]      # core/transfer.py     register_merge(name)
SELECTORS: dict[str, Selector]    # core/keep.py         register_selector(name, needs=(), tolerant=True, help="")
WRITERS: dict[str, type[Writer]]  # core/writers.py      register_writer(cls)
PRESETS: dict[str, Preset]        # core/presets.py      register_preset(...)
```

`core/constants.py` holds the name lists and one-line help used by `parsers.py`. A unit test checks that those lists equal the registry keys.

## E.5 toytree rules (enforced in `treeio.py`, `newick.py` and `rooting.py`, and covered by tests)

1. **Read through the pre-tokenizer.** Label syntaxes are rewritten at the string level before toytree parses anything (F.1). This avoids three verified problems. toytree's distance parser raises on `A:0.1#1`, toytree master moves `{...}` into a `trait` feature on fully labeled trees, and a single `95{Test}` internal name breaks toytree's guess between names and support values.
2. **Query by node index.** Always pass integer `idx` values to `root()`, `mod.prune()` and `get_nodes()`. String queries also match internal names, and a name starting with `~` is treated as a regex.
3. **Keep labels on edges.** Re-register `label` as an edge feature after every read, prune or copy, and pass `edge_features=["label"]` to every `root()` call.
4. **Merge root labels before rerooting.** Rerooting an already rooted tree drops one root child's custom label, so twig plans the merge first (F.4).
5. **Pruning.** Use the tree returned by `tree.mod.prune(*idxs)`, then reset the root length to 0 and re-register the edge feature. Never use `drop_tips` followed by `remove_unary_nodes`, which leaves a unary root in 3.0.x.
6. **Never prune the labeled species tree.** toytree keeps only the tipward node's features when it merges edges. Restriction is computed directly from clade sets (G.3).
7. **No `collapse_nodes`.** Its default `min_dist=1e-6` silently collapses zero-length branches as well.
8. **No toytree `unroot()` and no toytree `write()` for labeled output.** `unroot()` leaves custom labels on the new root, where they are never written, and released `write()` cannot put `{}` text on some nodes only. twig uses its own serializer.
9. **No toytree parsimony or ancestral-state code.** In 3.0.x it is a downpass-only stub.
10. **Monophyly on possibly unrooted trees.** Use a split test (a sample set S is a clade side if some non-root node has exactly S, or exactly the complement of S, below it), because `is_monophyletic` raises on unrooted trees.
11. **Exceptions.** Catch `ToytreeError`, `NodeDataError` and `ToytreeRegexError` per gene (`from toytree.utils import ToytreeError`). Other exceptions inside a gene are logged with a traceback at DEBUG and recorded as `internal-error`, so bugs are not hidden.
12. **Version matrix.** Tests run on toytree 3.0.10, 3.0.11 and master. Code never relies on `set_node_data(edge=)`, a top-level `toytree.ToytreeError`, or `write_single_feature`.

---

# Part F. Names, reading and writing

## F.1 Name mapping (`core/names.py`, `twig name-map`)

### Map file (`-I/--imap`)

```
# sample                    species               specimen (optional)
A_pumilus_sample1           Amaranthus_pumilus    s1
AmaranthusPumilus-Sample2   Amaranthus_pumilus    s2
~^Ahyb[_-]                  Amaranthus_hybridus
Beta_vulgaris
```

- **Format.** Columns are separated by whitespace, `#` starts a comment and blank lines are ignored. A row has one, two or three columns, and a one-column row maps a name to itself, as tree-filter's IMAP does.
- **Regex rows.** A first column starting with `~` is a regular expression applied with `re.search`, the convention toytree uses. Regex rows are tried in file order after all exact rows.
- **Matching.** Exact rows are matched against the raw name and then against the extracted key, so existing tree-filter IMAP files work unchanged.
- **Errors.** Duplicate exact rows, more than three columns or an empty file raise `TwigError` with `path:line`. When a species tree or vocabulary is given, every species value must resolve to it.

### Extraction

Extraction is optional when a species tree is given, because vocabulary matching on the whole name already handles mixed formats.

- **Delimiter rule** (`-d`, `-di`, `-dj`). Same meaning as in tree-filter. In the new commands the join string defaults to the delimiter itself. A name that cannot be split as asked falls through to the next step instead of stopping the run.
- **Regex rule** (`-x/--name-regex`). A named group `species` gives the key, and an optional named group `specimen` gives the specimen. Names the regex does not match fall through.

### Normalization (`-N/--normalize`)

Normalization is used only for comparisons. Names written to output are never changed by it.

| step | effect | example |
|---|---|---|
| `accents` | Unicode NFKD, combining marks removed | `Ñandú` to `Nandu` |
| `camel` | separator at lower-to-upper and acronym boundaries | `AmaranthusPumilus` to `Amaranthus_Pumilus` |
| `digits` | separator at letter-digit boundaries (off by default) | `Sample1` to `Sample_1` |
| `sep` | every run of non-alphanumeric characters becomes `_`, trimmed at the ends | `A.pumilus-s1` to `A_pumilus_s1` |
| `case` | casefold | `Amaranthus_Pumilus` to `amaranthus_pumilus` |

The default is `accents camel sep case`, and `-N none` switches normalization off.

### Vocabulary matching

The vocabulary is the species tree's tips (`-s`) or a list file (`-V`).

1. `exact`. The extracted key equals a vocabulary entry.
2. `canon`. The canonical forms are equal.
3. `token`. The vocabulary entry's tokens appear as a contiguous run anywhere in the name's tokens, so `amaranthus pumilus sample1`, `s012 amaranthus pumilus` and `og1 amaranthus pumilus g7` all match `Amaranthus_pumilus`.
4. `abbrev`. Like `token`, except that a non-final vocabulary token (the genus) may be matched by a prefix of at most three characters, so `A_pumilus_sample1` matches `Amaranthus_pumilus`.
5. **Longest match wins.** The entry covering more tokens wins, so `A_hybridus_quitensis_s1` goes to `A_hybridus_quitensis` when both exist. Equal candidates give status `ambiguous`, so `A_pumilus` is ambiguous between `Amaranthus_pumilus` and `Atriplex_pumilus`, while `AmaranthusPumilus` still resolves.

**Safety rule.** The labeling commands accept `map`, `map-regex`, `exact`, `canon` and `token` matches automatically. An `abbrev` match is accepted only inside `name-map`, through a frozen map, or with `--accept-inferred`, so an Atriplex sample missing from the tree is never silently assigned to Amaranthus. Unresolved names list their three closest vocabulary entries (by `difflib` similarity) as suggestions, and similarity is never used to assign a species.

### Precedence and unresolved names

The first step that resolves a name wins, in the order exact map row, regex map row, extraction followed by vocabulary steps, vocabulary steps on the whole canonical name, and finally unmapped.

| `-U/--unmapped` | effect |
|---|---|
| `error` (default) | a header-only pre-scan of all inputs runs first, and any unmapped or ambiguous name stops the run with a full list |
| `skip` | genes containing such names are skipped and reported |
| `drop` | such tips and rows are removed before labeling, and each is listed |
| `ignore` | such tips stay, are left out of the species sets used for matching, and get status `unmapped` |

### API

```python
@dataclass(frozen=True)
class NameHit:
    name: str; species: str | None; specimen: str | None
    method: str; status: str              # ok | ambiguous | unmapped | inferred
    leftover: tuple[str, ...]; candidates: tuple[tuple[str, float], ...] = ()

class NameMapper:
    def __init__(self, map_rows=(), regex=None, delim=None, delim_idxs=(0,), delim_join=None,
                 vocabulary=None, normalize=DEFAULT_STEPS, accept_inferred=False): ...
    @classmethod
    def from_args(cls, args, vocabulary=None) -> "NameMapper": ...
    def resolve(self, name: str) -> NameHit           # memoized
    def resolve_all(self, names) -> dict[str, NameHit]
    def hits(self) -> list[NameHit]

def parse_map_file(path) -> list[MapRow]
def canonical(text, steps=DEFAULT_STEPS) -> str
def write_hit_table(hits, handle) -> None
def write_frozen_map(hits, path) -> None    # writes a specimen column only when one was given explicitly
```

The frozen map leaves the specimen column empty unless a regex group or map column set it. Otherwise a reloaded map would make every sample its own specimen and silently disable duplicate detection.

### `twig name-map` options

| group | flag | default | help |
|---|---|---|---|
| Input/Output | `-i, --input path+` | required | trees and/or FASTA files, directories or quoted globs |
| | `-s, --sptree path` | | species tree whose tips form the vocabulary |
| | `-V, --vocab path` | | file listing species names |
| | `-o, --out path` | stdout | mapping table TSV |
| | `--freeze path` | | write the resolved map for reuse with -I |
| Taxon Mapping | shared group | | |
| Reporting | `--show str` | `all` | rows to print (all, problems) |
| | `--strict` | off | error if any name is unmapped, ambiguous or inferred, or leaves alphabetic leftover tokens |

Output columns are `name species specimen method status leftover candidates n_files`. The stderr summary gives counts per method and status, vocabulary species with no samples, the distinct leftover tokens per species (which exposes infraspecific names lumped into a species), and the samples per species per file (which exposes a delimiter that splits too finely or too coarsely).

### Shared Taxon Mapping group (`add_taxon_mapping_args`)

| flag | default | help |
|---|---|---|
| `-I, --imap path` | | map file rows `<sample> <species> [specimen]`; ~regex rows allowed |
| `-d, --delim str` | | delimiter used to split sample names |
| `-di, --delim-idxs int+` | `[0]` | token indices kept after splitting |
| `-dj, --delim-join str` | same as delim | join string for kept tokens |
| `-x, --name-regex str` | | regex with `(?P<species>...)` and optional `(?P<specimen>...)` |
| `-N, --normalize str+` | `accents camel sep case` | normalization steps used for matching (or none) |
| `--accept-inferred` | off | accept genus-abbreviation matches |
| `-U, --unmapped str` | `error` | unresolved names (error, skip, drop, ignore) |

## F.2 Reading labeled trees (`core/treeio.py`)

```python
def iter_trees(path) -> Iterator[tuple[str | None, str]]
def read_tree(text_or_path, label_syntax="auto", paml_codes=None) -> ToyTree
def strip_labels(tree) -> ToyTree
```

- **Splitting a file into trees.** `iter_trees` reads up to each `;` outside `[...]` comments, so multi-line Newick works. It skips blank lines, `#` comment lines and a PAML `ns ntree` header, handles CRLF, accepts an optional `locus<TAB>newick` form, and gives a clear `TwigError` for NEXUS input and for quoted labels, which toytree does not support.
- **Label pre-tokenizer.** Before toytree parses the string, a small tokenizer that skips existing `[...]` comments rewrites label syntaxes.
  - `{X}` becomes `[&twiglab=X]`. This covers HyPhy and phylotree.js output (`A{Test}`, `){Test}`, `)Node3{Test}`, `)95{Test}`).
  - PAML `#k` becomes `[&twigpaml=k]` and `$k` becomes `[&twigpamlc=k]`. `$k` is expanded to every descendant after parsing, with `#` beating `$` and inner marks beating outer ones, as in codeml. Codes map to names through `--paml-codes 1=Test 2=Reference`.
  - NHX `[&&NHX:label=X]` and FigTree `[&label="X"]` are read natively by toytree, and twig strips quotes and casts the value to text.
  - A node that would end up with two `[...]` blocks raises a clear error.
- **Internal names and support.** Trees are read with `internal_labels="name"`. A fully numeric internal name then moves to `support`, and IQ-TREE's `80.5/95` (SH-aLRT and UFBoot) becomes two support features.
- **Branch lengths.** twig records whether the input had lengths. toytree fills missing lengths with 1.0, so writers omit lengths rather than inventing them.
- **Species-tree checks.** Tip names must be unique, a label on the root is ignored with a warning, and internal branch IDs are the node's name if it has one, otherwise `S{idx}`.

## F.3 Writing (`core/newick.py`, `core/writers.py`)

One iterative serializer builds every output, so token order is always `name{Label}:len` for tips and `)name{Label}:len` for internal nodes. It can write a tree unrooted in one pass by placing the dissolved root child's children inline (F.4).

```python
def to_newick(tree, *, tip_text, internal_text, label_text, dists=True, dist_fmt="%.12g", root_merge=None) -> str

class Writer:
    name: ClassVar[str]; suffix: ClassVar[str]; unrooted: ClassVar[bool]; sanitizer: ClassVar[str | None]
    def check(self, tree, ctx) -> list[Issue]: ...
    def render(self, tree, ctx) -> str: ...
```

| writer | suffix | output | rooting | milestone |
|---|---|---|---|---|
| `braces` | `.labeled.nwk` | twig's round-trip format, `name{Label}:len` | as given | M1 |
| `hyphy` | `.hyphy.nwk` | `name{Label}:len`, names sanitized for HyPhy, internal names optional | unrooted with policy | M1 |
| `table` | `.branches.tsv` | one row per branch | not applicable | M1 |
| `newick` | `.nwk` | labels stripped, for tree builders and viewers | as given | M1 |
| `paml` | `.paml.nwk` | `name #k`, `(...) #k`, no internal names or support | unrooted | M5 |
| `nhx` | `.nhx` | `[&&NHX:label=X]` | rooted | M5 |
| `csubst` | `.foreground.tsv` | foreground lineage file, one lineage per labeled origin | rooted tree alongside | M5 |
| `lists` | `.labels.tsv` | `label<TAB>species;species`, for table-driven tools such as RERconverge and PhyloAcc | not applicable | M5 |

### HyPhy specifics

- **Tip names.** Names are sanitized once across all genes during the pre-scan, so a sample gets the same written name in every gene. Every written name is unchanged by HyPhy's own identifier conversion, unique without regard to case (HyPhy matches tips to sequences case-insensitively), and not a HyPhy reserved word (for example `Beta`). Collisions stop the run with a list, or get `_2`, `_3` under `--sanitize-suffix`. The table goes to `OUTDIR/sanitized_names.tsv`, and the paired alignment uses the same names, so HyPhy never renames anything and its JSON output joins back to twig's branch table.
- **Internal names** (`--internal-names none|auto|support|keep`). `none` (default) writes only labels. `auto` names internal nodes `twn1`, `twn2` and so on in postorder, so HyPhy's per-branch results join to the branch table, at the cost of dropping support values. `support` writes support values, which HyPhy reads before auto-naming the node. `keep` sanitizes existing names.
- **Label names** must match `^[A-Za-z][A-Za-z0-9_]*$`. `All`, `Internal`, `Leaves` and `Unlabeled` in any case are rejected because they collide with HyPhy's built-in branch-set options, and twig's own action and status words are rejected too.
- **Codon alignments.** HyPhy's codon loader requires a length divisible by 3 and no stop codons. twig checks length, in-frame stops (`--genetic-code`, default standard) and the nucleotide alphabet before writing, with `--stops skip|mask|error` (default `skip`).
- **Short branches.** HyPhy's `--kill-zero-lengths Yes` default deletes internal branches it estimates at zero length, and labeled stems are the most likely to be affected. IQ-TREE's minimum length is 1e-6, so twig counts labeled internal branches at or below `--short-len 1e-6` (column `short_labeled`), and the generated HyPhy commands pass `--kill-zero-lengths Constrain`. This is a heuristic, because HyPhy decides on its own estimates.

## F.4 Root edges

HyPhy reads every tree as unrooted. Given a two-child root, it promotes the first internal root child, and the other child's branch stands for both root edges with only that child's label. The surviving label therefore depends on child order in the Newick string, which is a silent source of error. twig always unroots the tree itself before writing for HyPhy and PAML.

```python
@dataclass(frozen=True)
class RootMerge:
    dissolved: int; kept: int; labels: tuple[str | None, str | None]; label: str | None; conflict: bool

def plan_root_merge(tree, policy, priority=(), outgroup_side=None) -> RootMerge | None
```

- **Which child is dissolved.** A bifurcating internal root child when one exists, so PAML never gets more than three children at the root.
- **Policy** (`--root-edges`).

| value | merged label |
|---|---|
| `unlabel` (default) | equal labels are kept, otherwise the merged branch is unlabeled with a warning and a `root-merge` note |
| `labeled` | if exactly one of the two edges is labeled, its label is used |
| `outgroup` / `ingroup` | the label of the outgroup-side or ingroup-side edge (equal sizes go to the side containing the first `-r` outgroup) |
| `error` | the gene is skipped |
| `<LABEL>` | that label wins if either edge carries it |

- **Order of operations.** The merge is planned before validation and counting, so requirements and reports describe the branches HyPhy will actually see. The branch table marks the dissolved edge with `written=False` and `merged_into`.
- **Example.** In the user's tree the two root edges are O (Reference) and the ingroup branch (Reference). They agree, so `O{Reference}` is written with the summed length.

## F.5 Paired alignments

- **Reader.** `core/fasta.read_fasta` is strict. It rejects empty or duplicate headers and, with `aligned=True`, rows of unequal length. A header is its first whitespace-separated token.
- **What is written.** With `-a`, twig writes `OUTDIR/{locus}.fa` whose rows are exactly the written tree's tips, under the written names, in tree order.
- **Pairing tips with rows.** Tips and rows are compared through a tolerant key (exact, then IQ-TREE's renaming rule, then the mapper's canonical form), and every rename is reported. Reconciliation happens at the start of each gene, before rooting and matching, never by pruning a labeled tree afterwards. `--mismatch skip|intersect|error` (default `skip`) handles genes whose tips and rows still differ.
- **Empty columns.** `--drop-empty-columns none|site|codon` (default `none`) removes columns emptied by pruning, and `codon` keeps the reading frame.

---

# Part G. The commands in detail

## G.1 `twig tree-label`, labeling from tip selections

### Selections

- **Selectors.** A value that is an existing file is a list with one name per line, where `~` lines are regexes and `#` lines are comments. Any other value is a regex matched with `re.search`. `--exact` treats command-line values as exact names, so `A` does not also match `MA`.
- **What is matched** (`--select-on name|species`). Tip names, or the mapped species of each tip, which makes one selection work on gene trees with several samples per species.
- **Reporting.** Every set reports how many tips it selected and how many it actually labeled after precedence, so an overridden set is never silent.

### Rules (`--rule`, alias `--internal-nodes`)

| rule | HyPhy `label-tree` name | labeled branches |
|---|---|---|
| `tips` | None | the selected tips only |
| `all` (default) | All descendants | every branch all of whose descendant tips are selected |
| `crown` | All descendants, no MRCA | the `all` result minus every internal branch whose parent branch is not labeled in that result (evaluated in one pass; isolated selected tips stay labeled, matching HyPhy's code) |
| `some` | Some descendants | every non-root branch with at least one selected descendant tip |
| `fitch` | Parsimony (fixed) | branches reconstructed in the trait state by two-state parsimony |
| `clade` | (twig) | for each selector, the branch above the MRCA of its tips and everything below it |
| `stem` | (twig) | for each selector, only the branch above the MRCA of its tips |

- HyPhy names are accepted as aliases, case-insensitively, so `--internal-nodes "All descendants"` works.
- `--leaf-nodes skip` removes tip labels after the rule runs (under `fitch` the tips still provide the data). `tips` combined with `skip` is an error because it labels nothing. HyPhy's own Skip mode ignores the no-MRCA and None rules, and twig documents that difference instead of copying it.
- `clade` and `stem` raise an error when the MRCA is the root, which has no branch.

### Two-state parsimony (`core/parsimony.py`)

```python
@dataclass
class ParsimonyResult:
    states: dict[int, int]; state_sets: dict[int, frozenset[int]]; cost: Fraction
def sankoff_binary(tree, tip_states, gain=1, loss=1, root_state="free", tie="deltran") -> ParsimonyResult
```

1. **Tip states.** Selected tips are 1, unselected tips are 0, and tips given with `--unknown` cost nothing in either state. Costs are `--gain-cost` (0 to 1) and `--loss-cost` (1 to 0), stored as `Fraction` so equality tests are exact.
2. **Downpass** in idx order (children first), `d[v][s] = sum over children w of min over t of (c[s][t] + d[w][t])`. Polytomies are handled exactly.
3. **Root.** `best = min over s of d[root][s]`. With `--root-state 0|1`, `best = d[root][k]` and the other state gets infinite cost at the root.
4. **Uppass** gives the state sets of all most-parsimonious reconstructions. State s is in node v's set exactly when `d[v][s] + up[v][s] == best`.
5. **One reconstruction** is chosen in preorder, with `deltran` (default) delaying changes toward the tips and `acctran` moving them toward the root.
6. **Ambiguous nodes** (`--ambiguous deltran|acctran|unlabel|<NAME>`) follow the chosen reconstruction, stay unlabeled, or get their own label.
7. **Unrooted input.** The result depends on the root position whatever the costs, so twig warns when the root has three or more children (see H.8).

### Stems, losses and background

These switches apply to every rule and turn a reconstruction into exactly the labeling the biology needs.

- **`--stems keep|unlabel|<NAME>`.** A stem is a labeled internal branch whose parent branch does not carry the same label, the branch on which the trait arose somewhere along its length. `unlabel` leaves stems out of both test and reference. A name such as `Stem` gives them their own class, so one species tree supports both the crown-only and the stem-plus-crown analyses through masks at transfer time.
- **`--losses keep|unlabel|<NAME>`** (rule `fitch`). A loss is a branch in state 0 whose parent is in state 1, which is a reversion. `Rev` gives reversions their own class.
- **`--background <NAME>`.** Every branch not labeled by any set, not a stem and not a loss gets this label, which makes reference branches explicit for RELAX.

### Several label sets

- **Inline.** `-L/--label NAME[:RULE] SELECTOR...` can be repeated, for example `-L Test:fitch dioecious.txt -L Out:tips '^O$'`.
- **Spec file** (`--spec path`). Tab-separated `label rule selectors options`, one set per line, with options such as `ambiguous=`, `gain-cost=` and `loss-cost=`.
- **Precedence** (`--precedence specific|first|last|error`). The default `specific` lets a label given to a tip by direct selection beat a label inherited from an enclosing clade, so a reversion selected for `Rev` is never overwritten by Test. `first` reproduces HyPhy, where earlier sets win.
- **Existing labels** (`--existing keep|clear`). `keep` (default, as in HyPhy) gives labels already in the tree priority, and `clear` removes them first.
- **Rerooting** (`-r/--outgroups SPECIES...`, comma-joined for clades). twig reroots before labeling, after merging any existing root labels.

### Checking the result (`--show`)

`--show` prints one line per labeled branch to stderr, giving the label, whether it is a tip, stem, loss or crown branch, and the species below it. toytree 3.0.11's text drawing shows tip names only, so this listing is twig's own (`core/display.py`). An SVG drawing with edge colors can be added later.

### Options

| group | flag | default | help |
|---|---|---|---|
| Input/Output | `-i, --input path` | required | newick or multi-newick tree file |
| | `-o, --out path` | stdout | output tree file |
| | `-F, --format str` | `braces` | output writer (braces, hyphy, newick, table) |
| | `-B, --branch-table path` | | per-branch TSV |
| | `--show` | off | print labeled branches by clade to stderr |
| Selection | `-L, --label NAME[:RULE] SEL...` | | label tips matched by selectors (repeatable) |
| | `--spec path` | | TSV of label sets |
| | `--rule, --internal-nodes str` | `all` | tips, all, crown, some, fitch, clade, stem |
| | `--select-on str` | `name` | match selectors to tip names or mapped species |
| | `--exact`, `--invert` | off | exact names; label tips not matched |
| | `--leaf-nodes str` | `label` | label, skip |
| | `--stems`, `--losses str` | `keep` | keep, unlabel, or a label |
| | `--background str` | | label every remaining branch |
| | `--precedence str` | `specific` | specific, first, last, error |
| | `--existing str` | `keep` | keep, clear |
| | `-r, --outgroups str+` | | reroot on these tips or species first |
| Parsimony | `--ambiguous str` | `deltran` | deltran, acctran, unlabel, or a label |
| | `--gain-cost`, `--loss-cost float` | `1` | change costs |
| | `--unknown SEL+` | | tips with unknown state |
| | `--root-state str` | `free` | free, 0, 1 |
| Taxon Mapping | shared group | | used with --select-on species |
| Output Format | shared writer group | | `--root-edges`, `--keep-root`, `--internal-names`, `--no-dists` |
| Runtime | `-f, --force`, `-l, --log-level` | | |

| HyPhy `label-tree` | twig `tree-label` |
|---|---|
| `--regexp R --label L` | `-L L R` |
| `--list F --label L` | `-L L F` |
| `--internal-nodes "All descendants"` | `--rule all` (HyPhy names accepted) |
| `--leaf-nodes Skip` | `--leaf-nodes skip` |
| `--reroot X` | `-r X` |
| `--invert Yes` | `--invert` |
| `--output F` | `-o F -F hyphy` |

## G.2 `twig tree-label-transfer`, the main mode

### Per-gene pipeline

1. **Read and clean.** Parse the gene tree, strip existing labels (default), and move numeric internal names to support.
2. **Map and reconcile.** Resolve every tip to a species and specimen, apply `--unmapped`, and reconcile tips with alignment rows (F.5). All removals happen here, before any labeling.
3. **Repeats.** If any species (or specimen) occurs more than once, apply `--repeats` (G.2 repeats). Genes where every species occurs once skip all paralog logic.
4. **Species-level masks.** Apply label masks to a label view of the species tree. Species marked for pruning are removed from the gene tree and alignment.
5. **Root** only when needed, meaning for clade matching, duplication detection, rooted writers or `--root-edges outgroup|ingroup`.
6. **Restrict** the species tree to the species present, cached per species set.
7. **Classify** every gene branch.
8. **Status masks**, such as `unmatched=...`.
9. **Plan the root merge** for each requested writer (F.4).
10. **Validate** against the preset on the branches as written.
11. **Write** each format and the alignment, atomically.
12. **Report** rows for the gene report, branch table and species-tree recovery table.

### Rooting

| `--root` | behavior |
|---|---|
| `auto` (default) | root on the root split of the species tree restricted to the species present. With the outgroup present this is outgroup versus ingroup. Without it, the deepest available split is used (for example M1+ABC versus M2+DE). A polytomy there falls back to the ordered rooting clades |
| `outgroup` | root on the ordered outgroups from `-r`, or on the ordered rooting clades computed from the species tree (the tree-rooter logic, reimplemented in `core/rooting.py`) |
| `as-is` | trust the input root, for example trees rooted with `twig tree-rooter` |
| `mad` | minimal ancestor deviation, always on a copy |

The candidate side is tested as a split on sample sets before `root()` is called, which avoids toytree's error logging. `--root-fail unrooted|skip|mad` (default `unrooted`) handles failures. Split matching does not need a root, so an unrooted gene is still labeled correctly. Genes that need a root (clade matching, repeated species) and cannot be rooted are skipped with reason `root-failed`.

### Restricting the labeled species tree

The species tree is never pruned with toytree (E.5 rule 6). Restriction is computed from clade sets.

- **Restricted clades.** For each non-root species node u, `r(u)` is the set of its tips present in the gene. A node gives a restricted branch when `r(u)` is neither empty nor all present species.
- **Merged branches.** Nodes with equal `r(u)` lie on one path and become one branch, with their labels listed tipward to rootward.
- **Terminal merges** (`--merge-tips`, default `tipward`). When the merged branch ends in a tip, the tip's own label wins, so tips always keep their species' label. The other path labels go into the branch table.
- **Internal merges** (`--merge`, default `unanimous`). The shared label if all components agree (unlabeled counts as a value), otherwise unlabeled with status `conflict`. `tipward` and `rootward` are also available, and more policies can be registered later.
- **Root branch.** The two root-child branches of the restricted tree share one split, and their labels are merged with the `--root-edges` policy, so all root behavior is set in one place.

```python
@dataclass(frozen=True)
class RestrictedBranch:
    key: frozenset[str]; clade: frozenset[str]; label: str | None
    path: tuple[str, ...]; path_labels: tuple[str | None, ...]; conflict: bool; origin: str | None

class SpeciesIndex:
    def __init__(self, sptree, labels: Mapping[str, str | None]): ...
    species: frozenset[str]; branches: tuple[SpeciesBranch, ...]
    def restrict(self, present: frozenset[str], cfg) -> "RestrictedIndex": ...   # per-instance dict cache

class RestrictedIndex:
    def by_split(self, side: frozenset[str]) -> RestrictedBranch | None
    def by_clade(self, clade: frozenset[str]) -> RestrictedBranch | None
    def tip(self, species: str) -> RestrictedBranch
```

### Matching

- **Split mode** (`--match split`, default). A gene branch matches the restricted branch with the same split of present species. On every non-root branch this gives exactly the user's rule (same set of species below the branch) whenever the gene is rooted like the species tree. It differs only at the root, where HyPhy and PAML merge the two root branches anyway, and on misrooted genes, where it still labels every concordant branch correctly.
- **Clade mode** (`--match clade`). The literal rooted rule, which needs a rooted gene tree.
- **Compatible refinements.** A gene branch that resolves a species-tree polytomy is `compatible`, not `unmatched`, and is labeled when the polytomy branch and every child branch it groups carry the same label.

### The branch classifier

Overlap is computed from per-node species counts, so several samples of one species are handled correctly.

```python
total = Counter(sp[x] for x in samples if sp.get(x))
below = {}
for n in tree:                                    # idx order, children first
    below[n.idx] = (Counter({sp[n.name]: 1}) if sp.get(n.name) else Counter()) if n.is_leaf() \
                   else sum((below[c.idx] for c in n.children), Counter())

def classify(n):
    cnt = below[n.idx]; sB = frozenset(cnt)
    if not sB:
        return None, Status.UNMAPPED, "unmapped"
    if n.is_leaf():                               # tips never get an overlap check
        return from_branch(rindex.tip(only(sB)), "tip")
    if len(sB) == 1:                              # a clade of samples of one species
        return intraspecific(only(sB), cnt, total)
    if any(cnt[s] < total[s] for s in sB) and cfg.overlap == "strict":
        return None, Status.MIXED, ""
    rb = rindex.by_clade(sB) if cfg.match == "clade" else rindex.by_split(split_key(sB, P))
    if rb is None:
        return None, compatible_or_unmatched(sB), ""
    if rb.conflict:
        return None, Status.CONFLICT, ""
    return rb.label, (Status.LABELED if rb.label else Status.SPECIES_UNLABELED), "species-match"
```

`--intraspecific inherit|unlabel` decides whether branches inside a single-species clade take that species' terminal label. The default is `inherit`, but HyPhy presets default to `unlabel` and warn when any species has more than one sample, because branches within a species carry polymorphism, which inflates dN/dS. Reducing first with `sample-reduce` is the clean option.

### Repeated species and paralogs

The user's data are single-copy with one individual per species, so none of this activates for them. It keeps the tool safe on other data.

- **`--repeats error|individuals|copies`.** When a species appears more than once and no explicit specimen IDs exist, twig cannot tell individuals from gene copies, so the default stops with a message. `individuals` treats repeats as conspecific samples, and `copies` treats them as paralogs.
- **Detection.** A node is a duplication when two of its children share a species (or a specimen, when IDs exist). At a polytomy, single-species children of the same species are merged first, so `(A_1,A_2,B_1)` is not a false duplication. A duplication whose samples all belong to one species is an in-paralog inside that species' terminal branch.

| `-p/--paralogs` | effect |
|---|---|
| `skip` (default) | the gene is skipped with reason `paralog` |
| `mask` | labeling proceeds, then every branch at or below a duplication is unlabeled with status `paralog` |
| `split` | the gene is cut into maximal duplication-free subtrees (`{locus}.og1`, `.og2`), each labeled on its own |
| `error` | the run stops at the first multi-copy gene |

- **`--paralog-tips unlabel|label`** (default `unlabel`) also covers in-paralog tips, so tips in genes that keep paralogs can stay unlabeled. This avoids reading post-duplication relaxation as a trait effect.
- **Later addition.** Split subtrees lose their outgroups, so each stem merges with its sister at the root. Grafting single-copy outgroups onto each subtree (the Yang and Smith rooted-ingroup approach) is planned for a later milestone.

### Tips and deep coalescence

Every tip takes its species' terminal label. Under incomplete lineage sorting a gene's tip branch can reach back past the speciation event, so part of a labeled tip may predate the trait. This slightly dilutes signal and is standard practice. The report counts `unmatched` and `mixed` branches per gene as a measure of local discordance, and `--limit unmatched=N` drops genes that are too discordant.

### Masks (`-m/--mask CLASS=ACTION`, repeatable)

| action | meaning |
|---|---|
| `keep` | leave as is (default) |
| `unlabel` | remove the label |
| `prune` | remove species from the tree and the alignment |
| `<NAME>` | relabel to NAME (may be an existing label, so `-m Test2=Test` pools two origins) |

| class | allowed actions | applied |
|---|---|---|
| any label in the species tree | keep, unlabel, prune, relabel | on the species tree, before restriction |
| `species-unlabeled` | keep, relabel, prune (tips only) | on the species tree, before restriction |
| `nested-unlabeled` (an unlabeled tip whose parent is labeled, which captures reversions without a dedicated label) | keep, relabel, prune | on the species tree, before restriction |
| `unmatched`, `compatible`, `conflict`, `mixed`, `paralog`, `unmapped` | keep, relabel | on the gene tree, after classification |
| `unlabeled` (catch-all for unlabeled gene branches) | keep, relabel | on the gene tree, last |

- **Why label masks run on the species tree.** Masks then reach merged branches correctly. If Rev is folded into Reference before restriction, a merged branch made of Rev and Reference becomes a clean Reference branch instead of a conflict.
- **Why prune acts only on species.** A branch can only be removed by removing taxa, and an internal stem exists as long as one taxon remains on each side. So a stem can be unlabeled or relabeled, or (in BUSTED-PH) knowingly left in the background, but not pruned.
- **Protection.** `prune` never removes the rooting outgroups or species named in `--protect`. `--prune-species LIST|REGEX` prunes named species directly.
- **Two origins.** Labeling the origins `Test1` and `Test2` in the species tree allows a per-origin analysis (`-m Test2=unlabel`) or a pooled one (`-m Test1=Test -m Test2=Test`) from one tree. The branch table also has an `origin` column naming each connected labeled component.

### Presets, requirements and validation

| preset | writer | requirements (`--test-label Test`, `--ref-label Reference`) | consequence warned about |
|---|---|---|---|
| `relax` | hyphy | at least one Test and one Reference branch | other labels and unlabeled branches become Unclassified |
| `busted` | hyphy | at least one Test branch and one other | every non-Test branch is background |
| `busted-ph` | hyphy | at least one Test branch and one other | every non-Test branch is background, including stems and discordant branches, with `n_uncertain_bg` reported |
| `absrel` | hyphy | at least one Test branch | only Test branches are tested |
| `paml-bs` | paml | exactly one foreground code, at most three children per node | every other label becomes #0 |
| `csubst` | braces plus csubst | rooted output | internal labels with unlabeled tips below cannot be expressed |

- **Overrides.** `--require CLASS=N`, `--limit CLASS=N`, `--require-length CLASS=X` and `--min-species N` (default 3) add to the preset. Branch length per class (`len_<label>`) is reported, because RELAX and BUSTED-PH power depends on branch length, not branch counts.
- **Failures.** `--on-fail skip|warn|error` (default `skip`).
- **Labels outside the preset** produce a warning naming the HyPhy consequence (for example "Rev becomes Unclassified in RELAX"), not an error.
- **HyPhy commands.** With `-O`, twig writes `OUTDIR/hyphy_commands.sh` with one line per written gene, for example `hyphy relax --alignment {fa} --tree {nwk} --test Test --reference Reference --mode "Classic mode" --kill-zero-lengths Constrain --output {json}`, ready for GNU parallel or a cluster.

| validation check | when | outcome |
|---|---|---|
| unmapped or ambiguous names | pre-scan | stop with a full list (default) |
| duplicate tips or rows | per gene | skip |
| tips differ from rows after reconciliation | per gene | `--mismatch` |
| sanitized name collisions | pre-scan | stop, or suffix |
| invalid or reserved label names | start | `TwigError` |
| codon length, stops, alphabet | per gene | `--stops` |
| root-edge conflict | per gene | `--root-edges` |
| labeled short internal branches | per gene | warning and `short_labeled` column |
| requirements and limits | per gene, after the root merge | `--on-fail` |

`-c/--check` runs everything and writes only the reports.

### Reports (written into OUTDIR by default)

- **`report.tsv`** has one row per gene with the columns `locus status reason n_tips n_species max_per_species root match n_pruned n_dup n_<label>... len_<label>... n_unlabeled n_species_unlabeled n_unmatched n_compatible n_conflict n_mixed n_paralog n_uncertain_bg len_uncertain_bg stems_recovered stems_missing stems_lost_root_merge root_merge short_labeled warnings outputs`.
- **`branches.tsv`** has one row per gene branch with the columns `locus node is_tip n_tips tips species split label status source sp_branch origin path_labels dist support written merged_into note`. The `split` column is a short hash of the canonical split, the same whatever the rooting.
- **`sptree.tsv`** has one row per species-tree branch with the columns `branch clade label genes_present genes_testable genes_recovered frac_recovered`, which answers questions such as how many genes recover the D+E stem.
- Every TSV starts with `#` lines giving the command and the twig and toytree versions.

### Batches and parallelism (`core/batch.py`)

- **Inputs.** `-i` and `-a` take files, directories or quoted globs, which twig expands itself so thousands of files never hit the shell's argument limit. `-M/--manifest` takes a TSV `locus tree [alignment]` instead.
- **Locus names** come from file names through `--locus-regex` (default `^([^.]+)`), so `OG0001.treefile` pairs with `OG0001.nt.fa`. Duplicate locus names are an error, and genes without an alignment and orphan alignments are reported.
- **Parallel runs** (`-j/--jobs`) use `run_with_pool_iter`. Each job receives the locus, the species-tree text, the config and the mapper arguments. Workers build their `SpeciesIndex` and `NameMapper` once, disable toytree logging after import, and catch per-gene exceptions, so one bad gene never stops the batch. The parent writes reports in input order.

### Options

| group | flag | default | help |
|---|---|---|---|
| Input/Output | `-s, --sptree path` | required | labeled species tree (braces, NHX, FigTree or PAML marks) |
| | `-i, --input path+` | | gene trees (files, directories or quoted globs) |
| | `-a, --alignments path+` | | alignments paired to trees by locus (requires -O) |
| | `-M, --manifest path` | | TSV `<locus> <tree> [alignment]` |
| | `-o, --out path` | stdout | all labeled trees in one multi-newick file |
| | `-O, --outdir path` | | one tree and alignment per locus, plus reports |
| | `-F, --format str+` | from preset, else `hyphy` | writers |
| | `--locus-regex str` | `^([^.]+)` | locus name from file name |
| Taxon Mapping | shared group | | |
| Rooting and Matching | `-R, --root str` | `auto` | auto, outgroup, as-is, mad |
| | `-r, --outgroups str+` | | ordered outgroups (comma-joined for clades) |
| | `--root-fail str` | `unrooted` | unrooted, skip, mad |
| | `--match str` | `split` | split, clade |
| | `--merge str` | `unanimous` | internal merged species branches |
| | `--merge-tips str` | `tipward` | terminal merged species branches |
| | `--intraspecific str` | `inherit` | inherit, unlabel |
| Repeats and Paralogs | `--repeats str` | `error` | error, individuals, copies |
| | `-p, --paralogs str` | `skip` | skip, mask, split, error |
| | `--paralog-tips str` | `unlabel` | unlabel, label |
| Masking | `-m, --mask str` | | CLASS=ACTION, repeatable |
| | `--protect`, `--prune-species str+` | | species never pruned; species to prune |
| | `--drop-empty-columns str` | `none` | none, site, codon |
| | `--mismatch str` | `skip` | skip, intersect, error |
| | `--stops str` | `skip` | skip, mask, error |
| Program | `-P, --preset str` | | relax, busted, busted-ph, absrel, paml-bs, csubst |
| | `--test-label`, `--ref-label str` | `Test`, `Reference` | labels used by presets |
| | `--require`, `--limit`, `--require-length str` | | CLASS=N, repeatable |
| | `--min-species int` | `3` | minimum species per gene |
| | `--on-fail str` | `skip` | skip, warn, error |
| Output Format | shared writer group | | `--root-edges`, `--keep-root`, `--internal-names`, `--sanitize-suffix`, `--no-dists`, `--short-len` |
| Runtime | `-c, --check` | off | run and validate, write reports only |
| | `-j, --jobs int` | `1` | worker processes |
| | `--resume`, `-f, --force` | off | skip completed loci; overwrite outputs |
| | `-l, --log-level level` | `INFO` | |

Options belonging to later milestones are hidden under an Advanced group in help until they exist.

## G.3 `twig sample-reduce`, one specimen per species

`tree-label-transfer` handles several samples per species, but for HyPhy the clean route is to reduce first.

### Strategies (`--strategy`)

| strategy | inputs | what happens | trade-off |
|---|---|---|---|
| `prune` (default) | tree, optional alignment | keep one sample per species, prune the rest (lengths summed along the kept path), subset the alignment rows | fastest; the topology still reflects all samples' data, which may differ from a tree inferred from the reduced data |
| `rows` | alignment, optional tree for tree-based selectors | drop extra rows, optionally drop emptied columns, then emit or run a tree-building command | tree and data match while columns stay identical; one tree inference per gene |
| `realign` | unaligned sequences (`-u`), optional alignment and tree for scoring | write the reduced unaligned FASTA, then emit or run an alignment command and a tree command | most self-consistent, since alignment and trimming depend on which sequences are present; most expensive, and columns can change |
| `collapse` | tree, optional alignment | each species whose samples form a clade is reduced to one sample; `--non-monophyletic keep|select|skip` handles the others | keeps evidence of non-monophyly; remaining multi-sample species are handled by `--intraspecific` at transfer |

- **Commands are templates**, never hard-coded, for example `--tree-cmd 'iqtree2 -s {aln} --prefix {prefix} -m MFP -T {threads} --quiet'` or `--tree-cmd 'raxml-ng --search --msa {aln} --prefix {prefix} --model GTR+G'`, and `--align-cmd 'twig align-pre -i {seqs} -o {prefix}.pre.fa && twig align-cds -i {prefix}.pre.fa -o {prefix}.msa.fa'`. Placeholders are `{aln}`, `{seqs}`, `{prefix}`, `{locus}`, `{outdir}` and `{threads}`. The template is split with `shlex` first and only known `{name}` tokens are replaced, so paths with spaces and literal braces (awk) are safe.
- **Running.** By default twig writes `OUTDIR/commands.sh`, one line per locus, for GNU parallel or a scheduler. Running them directly with `--run` (through `safe_popen`) comes in a later milestone.

### Keep-one selectors (`core/keep.py`)

The `-K/--keep-by` chain is applied in order. Each selector scores the remaining candidates (higher is better), and candidates within `--tolerance` of the best (`score >= best - tol * max(abs(best), eps)`, default 0.02) pass to the next selector. Ordinal selectors use exact ties only. `name` is always appended, so results are deterministic.

| selector | needs | score | trade-off |
|---|---|---|---|
| `priority` | `--priority` file | rank in a list of names or `~regex` lines | encodes expert knowledge (for example reference genomes first); prepended automatically when given |
| `scores` | `--scores` TSV | a user value per sample (coverage, completeness) | flexible, only as good as the table |
| `occupancy` | alignment or sequences | count of A, C, G, T, U | maximizes data per site, usually what matters most for HyPhy power |
| `length` | sequences | non-gap characters including ambiguity codes | like occupancy, but counts Ns |
| `overlap` | alignment | mean sites shared with samples of other species | favors sequences covering what other species cover |
| `medoid` | tree | negative summed distance to conspecific samples | picks the most typical individual and avoids misidentified outliers; needs three or more samples |
| `terminal` | tree | negative terminal branch length | avoids long-branch artifacts but selects on the rate being tested, so opt-in only |
| `first` | none | order of appearance | reproducible baseline |
| `random` | `--seed` | seeded draw | sensitivity checks |
| `name` | none | lexicographic | final tie-break |

- **Default chain** (`auto`) is `[priority] [scores] occupancy name`. `medoid` and `terminal` are opt-in, because a branch-length rule biases toward slowly evolving individuals, and so toward the hypothesis being tested.
- **Scope** (`--scope gene|global`). `global` picks one specimen per species across all loci (most loci present, then total occupancy) and uses it wherever it occurs, which keeps genes comparable. It needs consistent specimen IDs and is the default when they exist, otherwise `gene`.
- **Report.** One row per locus and species with the columns `locus species n_samples kept dropped chain scores reason`.

### Options

| group | flag | default | help |
|---|---|---|---|
| Input/Output | `-i`, `-a`, `-u`, `-M`, `--locus-regex` | | as in tree-label-transfer, plus unaligned sequences |
| | `-O, --outdir path` | required | output directory |
| Strategy | `--strategy str` | `prune` | prune, rows, realign, collapse |
| | `-K, --keep-by str+` | `auto` | selector chain |
| | `--tolerance float` | `0.02` | relative tie tolerance |
| | `--scope str` | auto | gene, global |
| | `--priority path`, `--scores path`, `--seed int` | | selector inputs |
| | `--non-monophyletic str` | `keep` | collapse only (keep, select, skip) |
| | `--rename str` | `keep` | name of kept tips (keep, species) |
| | `--drop-empty-columns str` | `none` | none, site, codon |
| External Tools | `--tree-cmd`, `--align-cmd`, `--tree-out str` | | command templates and expected tree path |
| | `-t, --threads int` | `1` | value for {threads} |
| Taxon Mapping | shared group | | |
| Runtime | `-c`, `-j`, `--resume`, `-f`, `-l` | | |

## G.4 Deliberately excluded

- **Collapsing low-support branches.** A collapsed branch disappears with its label, labeled stems are exactly the short, weakly supported branches it would remove, HyPhy accounts for uncertainty through its own estimates, polytomies break PAML (three children at most), and toytree's `collapse_nodes` also collapses every branch shorter than 1e-6. Anyone who wants it can run `twig tree-filter -s` first.
- **Fuzzy assignment of species names.** Similarity scores only produce suggestions.
- **Copying HyPhy quirks.** Degenerate parsimony, Skip with no-MRCA, and Skip with None are documented as differences and not reproduced.

---

# Part H. Testing, documentation and roadmap

## H.1 Testing

- **Layout.** `tests/labeling/` with `conftest.py` (Newick fixtures, and a helper that builds a command's parser from `twig.labeling.parsers`, parses an argument list and calls its `run_*`), `tests/labeling/data/` (the user's species tree, toy genes, alignments and map files), and one test module per library module and per command.
- **Integration with twig.** One test builds the full `twig` parser and runs `twig -h`, one checks that importing `twig.cli.subcommands` does not load toytree, and one runs each command through `twig.cli.main.main()` with simple unquoted arguments.
- **Running.** `pip install toytree pytest`, then `pytest tests/labeling`. No CI file is added unless the user asks. When added, it would run on toytree 3.0.10, 3.0.11 and master.

### `tree-label` on the user's species tree `(O,((M1,(A,(B,C))),(M2,(D,E))));`

Internal nodes are named BC, ABC, X (M1+ABC), DE and Y (M2+DE), and the selection is A, B, D and E.

| rule and options | labeled branches |
|---|---|
| `tips` | A, B, D, E |
| `all` | A, B, D, E, DE |
| `crown` | A, B, D, E |
| `some` | A, B, D, E, BC, ABC, X, DE, Y, and the ingroup root branch |
| `fitch` | A, B, BC, ABC, D, E, DE (cost 3, gains on ABC and DE, loss on C, unique reconstruction) |
| `fitch --stems unlabel` | A, B, BC, D, E |
| `fitch --stems unlabel --losses Rev --background Reference` | the tree in C.2 |
| `all` (HyPhy's "Parsimony") | A, B, D, E, DE, the regression test showing why real parsimony matters |

### `tree-label-transfer` with the C.2 labeled tree, HyPhy writer, `--internal-names none --no-dists`

| gene | input | expected output |
|---|---|---|
| concordant | `(O,((M1,(A,(B,C))),(M2,(D,E))));` | `(O{Reference},(M1{Reference},(A{Test},(B{Test},C{Rev}){Test})){Reference},(M2{Reference},(D{Test},E{Test})){Reference});` |
| second figure, discordant | `(O,((M1,(B,(A,C))),(D,(M2,E))));` | `(O{Reference},(M1{Reference},(B{Test},(A{Test},C{Rev}))){Reference},(D{Test},(M2{Reference},E{Test})){Reference});` with A+C and M2+E `unmatched`, and the A+B+C stem matched but unlabeled |
| B missing | `(O,((M1,(A,C)),(M2,(D,E))));` | C keeps `C{Rev}` (tipward), and A+C matches the unlabeled stem |
| A missing | `(O,((M1,(B,C)),(M2,(D,E))));` | B+C merges the Test crown branch with the unlabeled stem, so it is unlabeled with status `conflict` (the stem uncertainty made explicit) |
| RELAX | concordant, `-P relax` | as concordant, with Rev reported as Unclassified |
| BUSTED-PH | concordant with alignment, `-P busted-ph -m Rev=prune` | `(O{Reference},(M1{Reference},(A{Test},B{Test})){Reference},(M2{Reference},(D{Test},E{Test})){Reference});` with row C removed. Without C, the B+C crown branch merges into B's tip, and A+B matches the unlabeled stem |
| crown versus stem | species tree labeled with `--stems Stem`, transfer with `-m Stem=Test` | the A+B+C and D+E branches become Test |
| unrooted input | `((M1,(A,(B,C))),(M2,(D,E)),O);` | same as concordant |
| two A samples | `(O_1,((M1_1,((A_1,A_2),(B_1,C_1))),(M2_1,(D_1,E_1))));` with `--repeats individuals` | A_1 and A_2 Test, the A_1+A_2 branch unlabeled under the relax preset, plus a multi-sample warning |
| paralog | a tree with two copies of the whole ingroup | skipped with reason `paralog`; with `-p mask` nothing at or below the duplication is labeled |

Further tests cover root-edge policies, sanitization (`AmaranthusPumilus-Sample1` written as `AmaranthusPumilus_Sample1` in both tree and FASTA), name mapping (both of the user's name styles resolve, an Atriplex homonym makes `A_pumilus` ambiguous, precedence), keep-one selectors, glob pairing, resume, identical results with `-j 2` and `-j 1`, and an optional end-to-end test that runs `hyphy relax` when HyPhy is on PATH.

## H.2 Documentation

All pages go in `docs/labeling/`, following the layout of `docs/tree-filter.md` (introduction, command help, topic sections with examples, parallelization).

| page | contents |
|---|---|
| `design.md` | this document |
| `overview.md` | concepts, statuses, what unlabeled means per program, why no support collapsing |
| `name-map.md`, `tree-label.md`, `tree-label-transfer.md`, `sample-reduce.md` | one page per command |
| `hyphy-workflow.md` | recipes for RELAX, BUSTED, BUSTED-PH and aBSREL, crown versus stem, reversions, root edges, short branches |

mkdocs builds these pages without a nav entry. Adding them to the site menu needs the optional `mkdocs.yml` edit in B.2.

## H.3 Roadmap

| milestone | contents | what the user can do afterwards |
|---|---|---|
| M0 foundations | `twig/labeling` skeleton, the six added lines in `main.py` and `subcommands.py`, integration tests, `core/names.py`, `fasta.py`, `selectors.py`, `batch.py`, `treeio.py`, `newick.py`, `sanitize.py`, `twig name-map` | check and freeze the name mapping across every gene |
| M1 tree-label | `rules.py`, `parsimony.py`, stems, losses, background, `--show`, braces, hyphy, newick and table writers, root-edge planning | build and inspect the labeled species tree |
| M2 transfer core | `SpeciesIndex`, restriction and merges, split and clade matching, rooting, classifier, label masks including prune, alignment pairing and writing, presets relax, busted, busted-ph and absrel, reports, `hyphy_commands.sh`, check, resume | RELAX and BUSTED-PH inputs for all single-copy genes |
| M3 repeats and paralogs | `--repeats`, intraspecific handling, duplication detection, skip, mask and split policies, status masks | safe runs on multi-individual and multi-copy data |
| M4 sample-reduce | keep-one selectors, global scope, prune, rows, realign and collapse strategies, command files | one specimen per species by any strategy |
| M5 more writers and speed | paml, nhx, csubst and lists writers, `-j` parallelism, `--run` for templates, SVG drawing, outgroup grafting for split paralogs | PAML, csubst and table-driven tools |

Each milestone ships with its tests and docs and changes no existing twig behavior. M0 to M2 cover everything the user's current analysis needs.

## H.4 Open questions for the user

1. **Command names.** Are `name-map`, `tree-label`, `tree-label-transfer` and `sample-reduce` acceptable? `name-map` and `sample-reduce` are generic and could clash with a future upstream command, so `label-names` and `label-reduce` are alternatives.
2. **mkdocs nav.** May the optional `mkdocs.yml` nav block be added, or should the pages stay unlisted?
3. **Matching default.** Split matching (identical to the literal clade rule on every non-root branch of a correctly rooted gene, and robust to rooting failures) as the default, with `--match clade` available?
4. **Merged branches.** Tips keep their own label when a sister species is missing (`--merge-tips tipward`), and internal merged branches with disagreeing labels become unlabeled (`--merge unanimous`). Agreed?
5. **Root edges.** For HyPhy output, unlabel the merged root branch when its halves disagree (`unlabel`), or prefer the labeled half (`labeled`)?
6. **Reversions.** Their own label (`Rev`, produced by `--losses Rev`) in the species tree, or unlabeled with the derived `nested-unlabeled` class?
7. **Internal names in HyPhy trees.** None by default, or `twn1`, `twn2` names so HyPhy branch results join to twig's branch table?
8. **HyPhy runner.** Stop at writing inputs and `hyphy_commands.sh`, or later add a command that runs HyPhy?

The answers are recorded in H.6.

## H.5 Execution on approval and verification

1. Create `docs/labeling/design.md` with this document's content. No other file is created or changed in this step.
2. Check with `git status` and `git diff --stat` that the only change is that new file.
3. Commit with a message describing the design document, and push with `git push -u origin tree-labeling`.
4. Start M0 only after the user answers H.4. At each milestone, verify by installing toytree in the session (`pip install toytree pytest`), running `pytest tests/labeling`, running `twig -h` and each new command on `tests/labeling/data/`, comparing outputs with the expected trees in H.1, and confirming with `git diff --stat` that existing files changed only as listed in B.2.

## H.6 Answers to the open questions

Recorded on 2026-10-13, after the design document was committed.

1. **Command names.** Keep `name-map`, `tree-label`, `tree-label-transfer` and `sample-reduce`. None of them, nor the alternatives `label-names` and `label-reduce`, appears in the `DISPATCH` of this fork or of upstream eaton-lab/twig. The upstream check required by B.5 was made on 2026-10-13 against upstream commit `9bc5a7e`, which is identical to the fork's `main`.
2. **mkdocs nav.** Keep the pages unlisted for now. `mkdocs.yml` is not edited, so the optional row in B.2 is dropped unless the user asks for it later.
3. **Matching default.** Split matching, with `--match clade` available.
4. **Merged branches.** `--merge-tips tipward` and `--merge unanimous`, as planned.
5. **Root edges.** `--root-edges unlabel`, as planned.
6. **Reversions and nuisance species.** Neither fixed option was chosen, and the user gave a direction instead. Most people cannot know in advance whether a species is a reversion, so twig must not depend on detecting reversions. A user needs the flexibility to mark any species or clade as nuisance for any reason, and how nuisance species are handled must be a choice made per analysis, because it depends on the program. The consequences are listed below.
   - `Rev` in C.2 is one example of a label the user picked and not a built-in. Documentation uses a generic nuisance label in its examples.
   - A species or clade can be marked as nuisance directly with `-L NAME:tips SELECTOR` (precedence `specific`), whether or not parsimony finds a loss there. `--losses` stays as an optional convenience.
   - How the nuisance label is treated is chosen at transfer time with masks, namely `-m NAME=prune`, `-m NAME=unlabel`, `-m NAME=Test` or `-m NAME=Reference`. RELAX sets aside any branch that is neither Test nor Reference. BUSTED-PH has no such set, so the species is either pruned or accepted as background. The preset warnings in G.2 name that consequence for each program.
   - The exact options (for example a `--nuisance` shorthand for `tree-label`, and whether nuisance tips count as unknown state in parsimony) are settled when M1 and M2 are built, and are confirmed with the user at that point.
7. **Internal names.** `--internal-names none` is the default.
8. **HyPhy runner.** twig stops at writing HyPhy inputs and `hyphy_commands.sh`. A command that runs HyPhy is not planned.

## H.7 Decisions made while building M1

Recorded on 2026-10-13, when `tree-label` was built. Where this section differs from the sections above, this section is what the code does.

1. **Tie rule.** The default of `--ambiguous` is `acctran`, not `deltran` as written in G.1. On the species tree the trait can have been gained on the branch above A, B and C and lost on C, or it can have been gained separately on A and on B, and both histories cost 3. The `fitch` row of H.1 and the figure in C.2 show the first history, and only `acctran` produces it. The value `deltran` stays available and labels A, B, D, E and DE. The statement in H.1 that the reconstruction is unique is therefore not right.
2. **Nuisance species.** `--nuisance SELECTOR...` and `--nuisance-label NAME` (default `Nuisance`) mark species or clades as nuisance, with the rule `all`. A nuisance label beats every other label, labels already in the tree included. It is not changed by `--stems`, `--losses`, `--leaf-nodes skip` or `--background`, and its tips have an unknown state in parsimony. Tips given with `--unknown` that are not nuisance take the state that parsimony infers for them, and are labeled with it.
3. **Extra options.** `tree-label` also has `--paml-codes` (names for PAML label numbers), `--sanitize-suffix` and `--names-table` (the table of written tip names), which the table in G.1 does not list.
4. **Branch table.** The columns of the `tree-label` table are `locus idx node is_tip n_tips tips species label kind source dist support written_name written merged_into note`, where `species` appears only with `--select-on species`. The `tips` column is cut after 200 names. The table describes the tree as it is written, so a dissolved root branch has `written` set to no.
5. **Checks before work.** Label names, rules, costs, the root-edge policy and the output paths are checked before any tree is read. Label names must be valid for HyPhy in every output format, so that a labeled tree can always be used for HyPhy later.
6. **Root labels.** A label on the root of the input is ignored with a warning. Rerooting joins the two old root branches by the `--root-edges` policy before toytree moves the root, and the branch that the new root cuts keeps its label on both halves.
7. **Help and examples.** The help of `tree-label` explains each option in plain words and ends with a list of examples that shows every rule and every option on one small tree. The tests run each command and compare its output with the help, and they check that every option has an example. The rule `parsimony` is accepted as a synonym of `fitch`, and the help says what the name means.

## H.8 Decisions made when the review of M1 was fixed

Recorded on 2026-10-14, when the findings of a review of `tree-label` were fixed. Where this section differs from the sections above, this section is what the code does.

1. **Stems.** A stem is the internal branch on which a trait arose, as in the glossary, and not simply a labeled branch whose parent has another label as G.1 says. With `fitch` it is a branch with a gain (state 1 below state 0), so a child of a root in state 1 is not a stem. With `crown` it is the branch that the rule leaves out at the base of each group, so `--stems NAME` labels it, the background never takes it, and the first crown branch is not a stem. With the other rules G.1 applies, and a labeled child of the root is a stem unless every child of the root carries the same label. The label of `--ambiguous` is not a trait label, so an ambiguous branch is never a stem, and its source is `ambiguous`. A stem never gets the background label, whatever `--stems` says.
2. **Stems of labels already in the tree.** When `--stems` is unlabel or a name and `--existing keep` is in effect, the labels kept from the input tree take part, with the definition for the rules other than `fitch` and `crown`, so a tree labeled elsewhere (a PAML tree, or the output of `--stems keep`) can be turned into its crown or stem variant. twig cannot tell whether such a label marks a trait or whether its stems were handled before, so every stem found among them is named in a warning, and the labels that the same run gives to the background, to losses and to nuisance species are never taken for traits. With `--stems keep` these labels are left alone as before. `--losses` needs a parsimony reconstruction, so it never acts on labels kept from the input tree.
3. **Unrooted trees.** G.1 point 7 first said that only unequal costs make the result of `fitch` depend on the root. A branch is labeled by the state of the node below it, so with equal costs too the labels change with the position of the basal node of an unrooted tree, and the same holds for every rule except `tips`, which do not know which tips lie below a branch without a root. When the root has three or more children, every label whose rule is not `tips` gives a warning that names the problem and suggests `-r` with the outgroup, and so does `--stems` on labels kept from the input tree. With `-r` the root has two children and there is no warning.
4. **Losses and the background.** G.1 says that the background skips losses. With `--losses keep` a loss is a branch where the trait is absent, like any branch that the rule did not label, so it gets the background label. `--losses unlabel` keeps losses out of the background, and `--losses NAME` gives them their own label. Holding losses back under `keep` as well would make `keep` and `unlabel` the same, so the code keeps this reading, and the user docs say so.
5. **HyPhy's Skip mode.** The user docs of `tree-label` now describe how HyPhy's Skip mode differs from `--leaf-nodes skip`, as G.1 and G.4 require, and list the HyPhy to twig option table of G.1.

---

# Appendix. Notes on existing twig code (not changed by this work)

The research found these problems in other commands. They are recorded here for the user to report upstream or fix separately, and this work does not touch them.

- tree-filter stores minmap values as text, so `--minmap` crashes.
- Blank lines in tree files reach `toytree.tree` in tree-filter and tree-rooter.
- `assert args.input.is_file` in tree-filter never checks anything, because the method is not called.
- tree-skeleton's `__main__` block is `pass`.
- macse_align imports a parser that does not exist.
- align-pre ignores its join argument.
- `requires-python` in `pyproject.toml` says 3.8 while `subcommands.py` needs 3.10.
- The README and docs describe `macse-*` commands that are no longer in `DISPATCH`.

A later, separate effort could offer `core/names.py` and `core/fasta.py` to the rest of twig to replace the duplicated name and FASTA parsing, but only as its own change with the user's agreement.
