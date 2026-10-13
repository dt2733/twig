#!/usr/bin/env python

import re
from textwrap import dedent
from pathlib import Path
from tempfile import gettempdir
from argparse import ArgumentParser, RawDescriptionHelpFormatter

TMPDIR = gettempdir()


class SingleMetavarHelpFormatter(RawDescriptionHelpFormatter):
    """Render optional args as '-x, --long VALUE' (single metavar copy)."""

    def _format_action_invocation(self, action):
        if not action.option_strings:
            return super()._format_action_invocation(action)
        if action.nargs == 0:
            return ", ".join(action.option_strings)
        default = self._get_default_metavar_for_optional(action)
        args = self._format_args(action, default)
        return f"{', '.join(action.option_strings)} {args}"


def get_parser_csubst(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for relabel tool.
    """
    KWARGS = dict(
    prog="csubst",
    usage="twig csubst -i MSA -t TREE -g FG [options]",
    help="test for molecular convergence in 'csubst' given CDS & NWK",
    formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),
    description=dedent("""
        -------------------------------------------------------------------
        | csubst: A simple wrapper to run and extract results from csubst
        -------------------------------------------------------------------
        |
        |
        |
        -------------------------------------------------------------------
    """),
    epilog=dedent("""
        Examples
        --------
        $ twig csubst -a MSA -t NWK -g FG -o OUT/PRE   # OUT/PRE.cb...
        $ ...
    """)
    )

    # create parser or connect as subparser to cli parser
    if parser:
        KWARGS['name'] = KWARGS.pop("prog")
        parser = parser.add_parser(**KWARGS)
    else:
        KWARGS.pop("help")
        parser = ArgumentParser(**KWARGS)

    # path args
    parser.add_argument("-a", "--alignment", type=Path, metavar="path", required=True, help="input CDS alignment")
    parser.add_argument("-t", "--tree", type=Path, metavar="path", required=True, help="input rooted tree file")
    parser.add_argument("-g", "--foreground", type=Path, metavar="path", help="foreground file")
    # parser.add_argument("-o", "--out", type=Path, metavar="path", help="out prefix; default is input path [{input}]")
    parser.add_argument("-o", "--outdir", type=Path, metavar="path", default="./csubst", help="Out directory will be created at this path [%(default)s]")
    # parser.add_argument("-p", "--prefix", type=str, metavar="str", help="optional outfile prefix. If None the cds filename is used")

    parser.add_argument("-m", "--max-arity", type=int, metavar="int", default=2, help="max combinatorial number of branches (K) [%(default)s]")
    parser.add_argument("-u", "--exhaustive-until", type=int, metavar="int", default=2, help="perform exhaustive (non-heuristic) search up N branch combs [%(default)s]")
    parser.add_argument("-c", "--cutoff-stat", type=str, metavar="str", default="OCNany2spe,2.0|omegaCany2spe,5.0", help="Cutoff stats for searching higher-order branch combs [%(default)s]")

    parser.add_argument("-b", "--foreground-table", action="store_true", help="foreground file is a table (fg_format=2)")
    parser.add_argument("-w", "--fg-exclude-wg", action="store_true", help="sets 'yes' to exclude within-group comparisons")
    parser.add_argument("-s", "--fg-stem-only", action="store_true", help="sets 'yes' to only stem comparisons")

    # other outputs
    parser.add_argument("--write-cbs", action="store_true", help="sets 'yes' to generate combinatorial-branch-site output")
    parser.add_argument("--write-cs", action="store_true", help="sets 'yes' to generate combinatorial-site output")
    parser.add_argument("--write-s", action="store_true", help="sets 'yes' to generate site output")

    # others
    parser.add_argument("-e", "--env", type=Path, metavar="path", default="csubst", help="conda env name where 'csubst' in installed [csubst]")
    parser.add_argument("-j", "--threads", type=int, metavar="int", default=1, help="number of threads")
    parser.add_argument("-v", "--verbose", action="store_true", help="print macse progress info to stderr")
    parser.add_argument("-f", "--force", action="store_true", help="overwrite existing result files in outdir")
    parser.add_argument("-l", "--log-level", type=str, metavar="level", default="INFO", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    # parser.add_argument("-L", "--log-file", type=Path, metavar="path", help="append stderr log to a file")
    return parser


def get_parser_diamond_bl(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for diamond-bl."""
    KWARGS = dict(
        prog="diamond-bl",
        usage="diamond-bl -q query.fa -t target.fa [options]",
        help="search BLASTP hits of a query proteome against a target or target DB",
        formatter_class=lambda prog: SingleMetavarHelpFormatter(prog, width=220, max_help_position=56),
        description=dedent("""
            ----------------------------------------------------------------
            | diamond-bl: simple DIAMOND BLASTP wrapper for tabular output |
            ----------------------------------------------------------------
            | Searches a query protein FASTA against a target FASTA or a
            | prebuilt DIAMOND database, writes tabular hits, and supports
            | optional filtering and annotation joins for downstream analyses.
            ----------------------------------------------------------------
        """),
        epilog=dedent("""
            Examples
            --------
            # basic query-vs-target search to file
            $ diamond-bl -q query.fa -t target.fa -o hits.tsv

            # reuse a prebuilt target DB
            $ diamond-bl -q sample.fa -t arabidopsis.dmnd --target-is-db -o sample_vs_at.tsv

            # apply coverage filters and append annotations
            $ diamond-bl -q sample.fa -t arabidopsis.dmnd --target-is-db --min-qcov 0.5 --min-scov 0.5 --target-annotations at_go.tsv -o annotated_hits.tsv

            # transfer GO terms directly from a GAF/GPA file
            $ diamond-bl -q sample.fa -t Araport11_pep.fa --best-hit-only --target-annotations tair.gaf.gz -o sample_go.tsv

            # write a consensus GO vote table across all query sequences
            $ diamond-bl -q group.fa -t Araport11_pep.fa --best-hit-only --target-annotations tair.gaf.gz --consensus-go-out group_go_consensus.tsv --consensus-go-id OG_9999 --consensus-go-min-support 0.5
        """)
    )

    # create parser or connect as subparser to cli parser
    if parser:
        KWARGS['name'] = KWARGS.pop("prog")
        parser = parser.add_parser(**KWARGS)
    else:
        KWARGS.pop("help")
        parser = ArgumentParser(**KWARGS)

    io = parser.add_argument_group("Input/Output")
    io.add_argument("-q", "--query", type=Path, metavar="path", required=True, help="query protein FASTA")
    io.add_argument("-t", "--target", type=Path, metavar="path", required=True, help="target protein FASTA or DIAMOND DB")
    io.add_argument("-o", "--out", type=Path, metavar="path", default=None, help="write TSV output to file (default: stdout)")
    io.add_argument("-H", "--header", action="store_true", help="write a header row to output TSV")
    io.add_argument("-c", "--compress", action="store_true", help="gzip-compress output (requires --out)")

    search = parser.add_argument_group("Search")
    search.add_argument("-e", "--evalue", type=float, metavar="float", default=1e-3, help="max e-value threshold [%(default)s]")
    search.add_argument("-b", "--bitscore", type=float, metavar="float", default=None, help="min bit-score threshold (overrides --evalue)")
    search.add_argument("-k", "--max-target-seqs", type=int, metavar="int", default=25, help="max candidate target hits DIAMOND reports per query before post-BLAST ranking filters [%(default)s]")
    search.add_argument("--target-is-db", action="store_true", help="treat --target as an existing DIAMOND database file")
    search.add_argument("-u", "--ultra-sensitive", action="store_true", help="enable DIAMOND ultra-sensitive mode")
    search.add_argument("--outfmt", type=str, metavar="field", nargs="*", default=None, help="DIAMOND outfmt fields: omit flag for extended defaults; pass empty list for standard outfmt 6")

    filtering = parser.add_argument_group("Filtering and Ranking")
    filtering.add_argument("--min-qcov", type=float, metavar="float", default=None, help="minimum query coverage fraction (0-1)")
    filtering.add_argument("--min-scov", type=float, metavar="float", default=None, help="minimum subject coverage fraction (0-1)")
    rank = filtering.add_mutually_exclusive_group()
    rank.add_argument("--best-hit-only", action="store_true", help="post-BLAST filter: from hits retained after -k/coverage filters, keep only the top-scoring hit per query")
    rank.add_argument("--top-hits", type=int, metavar="int", default=None, help="post-BLAST filter: from hits retained after -k/coverage filters, keep up to N top-scoring hits per query")

    ann = parser.add_argument_group("Target Annotation Join")
    ann.add_argument("--target-annotations", type=Path, metavar="path", default=None, help="annotation table to append by target ID (TSV, GAF(.gz), GPA(.gz)); GAF/GPA joins use exact ID then fallback from isoform suffix .<digits> to gene ID")
    ann.add_argument("--annotation-format", choices=["auto", "tsv", "gaf", "gpa"], metavar="format", default="auto", help="annotation file format [%(default)s]")
    ann.add_argument("--annotation-key", type=str, metavar="str", default="sseqid", help="join key column for TSV annotations only [%(default)s]")
    ann.add_argument("--annotation-cols", type=str, metavar="str", nargs="*", default=None, help="annotation columns to append (TSV=any non-key cols; GPA=go_* cols; GAF=go_* plus db_object_id/db_object_symbol/db_object_name)")

    consensus = parser.add_argument_group("Consensus GO Vote")
    consensus.add_argument("--consensus-go-out", type=Path, metavar="path", default=None, help="optional TSV path to write consensus GO-term support across distinct query IDs")
    consensus.add_argument("--consensus-go-id", type=str, metavar="str", default=None, help="identifier label for consensus GO rows (e.g., orthogroup ID); required with --consensus-go-out")
    consensus.add_argument("--consensus-go-min-support", type=float, metavar="float", default=0.0, help="keep GO terms with support fraction >= this threshold [%(default)s]")
    consensus.add_argument("--consensus-go-include-failed", action="store_true", help="include GO terms that fail --consensus-go-min-support in consensus output")

    runtime = parser.add_argument_group("Runtime")
    runtime.add_argument("-j", "--threads", type=int, metavar="int", default=4, help="threads used by DIAMOND [%(default)s]")
    runtime.add_argument("-d", "--tmpdir", type=Path, metavar="path", default=TMPDIR, help=f"directory for temporary DB files [%(default)s]")
    runtime.add_argument("-f", "--force", action="store_true", help="overwrite an existing auto-built DB in tmpdir")
    runtime.add_argument("-s", "--save-db", action="store_true", help="keep an auto-built DB after search finishes")
    runtime.add_argument("-B", "--binary", type=Path, metavar="path", default=None, help="path to DIAMOND executable")
    runtime.add_argument("--kwargs-makedb", type=str, nargs="*", default=None, help="additional passthrough argument tokens for 'diamond makedb'")
    runtime.add_argument("--kwargs-blastp", type=str, nargs="*", default=None, help="additional passthrough argument tokens for 'diamond blastp'")
    runtime.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "EXCEPTION"], metavar="level", default="INFO", help="stderr logging level [DEBUG|INFO|WARNING|EXCEPTION] [%(default)s]")
    return parser


def get_parser_diamond_pw(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for relabel tool.
    """
    kwargs = dict(
        prog="diamond-pw",
        usage="%(prog)s paths [args]",
        help="write .tsv blast hits for all pairs of fasta sequence files",
        formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),
        description=dedent("""
            -------------------------------------------------------------------
            | diamond-pw: write .tsv's of all-by-all pairwise diamond blastp  |
            -------------------------------------------------------------------
            | Perform all-by-all blastp on a set of fasta protein sequence
            | files. There are many more options available from the `diamond
            | blast` tool itself. This is intended only as a simple wrapper.  |
            -------------------------------------------------------------------
        """),
        epilog=dedent("""
            Examples
            --------
            $ diamond-pw -d d1.fa d2.fa d3.fa -o blast/ -e 1e-5
            $ diamond-pw -d d1.fa d2.fa d3.fa -o blast/ -e 1e-5 --no-self
            $ diamond-pw -d d1.fa d2.fa d3.fa -o blast/ -e 1e-5 -k 1 -j 10 -t 4

            # Example result dir (blast/) contains
            # d1_d1.tsv  d2_d2.tsv  d3_d3.tsv
            # d1_d2.tsv  d1_d3.tsv  d2_d3.tsv
            # d2_d1.tsv  d3_d1.tsv  d3_d2.tsv
            # sample_map.tsv
        """)
    )

    # create parser or connect as subparser to cli parser
    if parser:
        kwargs['name'] = kwargs.pop("prog")
        parser = parser.add_parser(**kwargs)
    else:
        kwargs.pop("help")
        parser = ArgumentParser(**kwargs)

    # add arguments
    parser.add_argument("-i", "--input", type=Path, metavar="path", nargs="*", required=True, help="fasta sequence files for all-by-all blastp")
    parser.add_argument("-o", "--outdir", type=Path, metavar="path", default="./blast", help="directory to write results (created if necessary) [%(default)s]")
    parser.add_argument("-e", "--evalue", type=float, metavar="float", default=1e-5, help="max evalue to report an alignment [%(default)s]")
    parser.add_argument("-m", "--min-bitscore", type=float, metavar="float", default=None, help="min bit-score to report an alignment (overrides -e) [%(default)s]")
    parser.add_argument("-k", "--max-target-seqs", type=int, metavar="int", default=25, help="max n target seqs to report hits for [%(default)s]")
    parser.add_argument("-t", "--threads", type=int, default=4, metavar="int", help="number of threads per diamond job [%(default)s]")
    parser.add_argument("-j", "--jobs", type=int, default=1, metavar="int", help="number of diamond jobs to run in parallel [%(default)s]")
    parser.add_argument("-d", "--dbdir", type=Path, metavar="path", default=TMPDIR, help="directory for temp .db files [%(default)s]")
    parser.add_argument("-n", "--no-self", action="store_true", help="do not perform self-self search")
    parser.add_argument("-f", "--force", action="store_true", help="overwrite .db files if they exist in dbdir")
    # parser.add_argument("-r", "--relabel-headers", type=int, default=1, help="...")
    # parser.add_argument("-r", "--relabel-taxa", type=int, default=1, help="...")
    # parser.add_argument("-I", "--imap", type=Path, help="map file to translate sample labels.")
    return parser


def get_parser_consensus_fasta(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for consensus-fasta tool."""
    KWARGS = dict(
        prog="consensus-fasta",
        usage="consensus-fasta -i ALIGN.fa -o CONS.fa [options]",
        help="write a consensus sequence from an aligned AA or CDS fasta",
        formatter_class=lambda prog: SingleMetavarHelpFormatter(prog, width=220, max_help_position=56),
        description=dedent("""
            -------------------------------------------------------------------
            | consensus-fasta: write one consensus sequence from an alignment
            -------------------------------------------------------------------
            | Reads an aligned AA or CDS fasta, computes a per-column consensus,
            | and writes one FASTA record to output.
            -------------------------------------------------------------------
        """),
        epilog=dedent("""
            Examples
            --------
            $ twig consensus-fasta -i MSA.fa -o consensus.fa
            $ twig consensus-fasta -i MSA.fa -o consensus.fa -n geneA_consensus
            $ twig consensus-fasta -i MSA.fa -o consensus.fa -a N -m 0.6
            $ twig consensus-fasta -i MSA.fa -o consensus.fa --ignore-chars '-?.NnXx'
        """),
    )

    if parser:
        KWARGS["name"] = KWARGS.pop("prog")
        parser = parser.add_parser(**KWARGS)
    else:
        KWARGS.pop("help")
        parser = ArgumentParser(**KWARGS)

    io = parser.add_argument_group("Input/Output")
    io.add_argument("-i", "--input", type=Path, metavar="path", required=True, help="input aligned fasta (AA or CDS)")
    io.add_argument("-o", "--output", type=Path, metavar="path", required=True, help="output fasta containing one consensus sequence")

    opts = parser.add_argument_group("Consensus")
    opts.add_argument("-n", "--name", type=str, metavar="str", default="consensus", help="header name for consensus sequence [%(default)s]")
    opts.add_argument("-a", "--ambiguous-char", type=str, metavar="char", default=None, help="character used for ties/no-data/low-support (default: auto N or X)")
    opts.add_argument("-m", "--min-majority", type=float, metavar="float", default=0.5, help="minimum top-character fraction required at a site [%(default)s]")
    opts.add_argument("--ignore-chars", type=str, metavar="str", default="-?.NnXx", help="characters ignored when voting consensus [%(default)s]")

    runtime = parser.add_argument_group("Runtime")
    runtime.add_argument("-f", "--force", action="store_true", help="overwrite existing output file")
    runtime.add_argument("-l", "--log-level", type=str, metavar="level", default="INFO", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    return parser


def get_parser_format_fasta(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for format-fasta tool.
    """
    KWARGS = dict(
        prog="format-fasta",
        usage="twig format-fasta FASTA [options]",
        help="format fasta to sort/filter/subset/relabel/modify sequences",
        formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),
        description=dedent("""
            -------------------------------------------------------------------
            | format-fasta: sort/filter/subset/relabel/modify headers or seqs |
            -------------------------------------------------------------------
            | This tool is used to modify a fasta file by performing one or   |
            | more of the following operations in order: (1) sorting; (2)     |
            | filtering; (3) sub-selecting; (4) relabeling; and (5) modifying |
            | headers and sequences. A modified fasta is written to STDOUT.   |
            | A subset of sequences/scaffolds can be returned in any order    |
            | with --subset-idx or --subset-names but note that --sort-alpha  |
            | or --sort-len override the ordering. This tool can be useful    |
            | for modifying genome fasta files or CDS or PEP sequence files.  |
            -------------------------------------------------------------------
        """),
        epilog=dedent("""
            Examples
            --------
            # operations on genome scaffolds
            $ twig format-fasta GENOME --subset 3 --relabel-list A B C > genome.fa
            $ twig format-fasta GENOME --subset-idx {0..3} --relabel-prefix-idx A_ > subgenomeA.fa
            $ twig format-fasta GENOME --subset-names Chr1 Chr2 --relabel-list chromosome_1 chromosome_2 > genome.fa
            $ twig format-fasta GENOME --sort-len --subset 8 --max-width 80 > chromosomes.fa

            # operations on sequences (e.g., cds, pep)
            $ twig format-fasta FA --relabel-prefix-idx spp --map NAME.map.tsv > NAME.fa
            $ twig format-fasta FA --sort-len --min-len 200 --subset 50 > NAME.fa
            $ twig format-fasta FA --max-width 80 --unmask > NAME.fa

            # complex operations: select sequences from a subset of scaffolds
            $ ids=$(twig format-gff GFF --subset-names Chr1 --feat mRNA --export-ids)
            $ twig format-fasta FA --subset-names $ids > Chr1.pep.fa
        """)
    )

    # create parser or connect as subparser to cli parser
    if parser:
        KWARGS['name'] = KWARGS.pop("prog")
        parser = parser.add_parser(**KWARGS)
    else:
        KWARGS.pop("help")
        parser = ArgumentParser(**KWARGS)

    # add arguments
    parser.add_argument("fasta", type=Path, help="a fasta sequence file (can be .gz)")

    sort_options = parser.add_mutually_exclusive_group()
    sort_options.add_argument("--sort-alpha", action="store_true", help="sort chromosomes by name alphanumerically")
    sort_options.add_argument("--sort-len", action="store_true", help="sort chromosomes by length (longest to shortest)")

    relabel_options = parser.add_mutually_exclusive_group()
    relabel_options.add_argument("--relabel-prefix", type=str, metavar="str", default="", help="append prefix label to start of existing labels [%(default)s]")
    relabel_options.add_argument("--relabel-prefix-idx", type=str, metavar="str", default="", help="overwrite labels with {{prefix}}{{counter}} using incrementing counter [%(default)s]")
    relabel_options.add_argument("--relabel-list", type=str, metavar="str", nargs="+", help="overwrite labels with new labels provided as a list [%(default)s]")
    # relabel_options.add_argument("--relabel-from-map", type=str, metavar="str", nargs="+", help="overwrite labels with new labels provided as a list [%(default)s]")
    parser.add_argument("--map", type=Path, metavar="Path", default=None, help="path to write optional relabel translation map (TSV format) [%(default)s]")

    sub_options = parser.add_mutually_exclusive_group()
    sub_options.add_argument("--subset", type=int, metavar="int", default=None, help="subselect the first N sequences after optional sorting")
    sub_options.add_argument("--subset-idx", type=int, metavar="int", nargs="+", help="subselect one or more sequences by 1-based index after optional sorting")
    sub_options.add_argument("--subset-names", type=str, metavar="str", nargs="+", help="subselect one or more sequences by name before optional relabeling")

    parser.add_argument("--min-len", type=int, metavar="int", default=None, help="sequences shorter than this length are excluded (default=None)")
    parser.add_argument("--max-width", type=int, metavar="int", default=None, help="write sequences with max width textwrapping (e.g., 80)")
    parser.add_argument("--unmask", action="store_true", help="convert any soft-masked sequences to upper case")
    parser.add_argument("--revcomp", action="store_true", help="reverse complement DNA sequences, or reverse AA sequences")

    # compress?
    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "EXCEPTION"], metavar="level", default="INFO", help="stderr logging level (DEBUG, INFO, WARNING, ERROR; default=INFO)")
    return parser


def get_parser_format_gff(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for format-gff tool.
    """
    KWARGS = dict(
        prog="format-gff",
        # usage="%(prog)s fasta [options]",
        usage="twig format-gff GFF [options]",
        help="format gff to filter and relabel genes and scaffolds",
        formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),
        description=dedent("""
            -------------------------------------------------------------------
            | format-gff: filter/relabel genome annotation table              |
            -------------------------------------------------------------------
            | This tool is used to ...
            -------------------------------------------------------------------
        """),
        epilog=dedent(r"""
            Examples
            --------
            # select a subset of scaffolds by name
            $ format-gff GFF --subset-names Chr1 > Chr1.gff

            # sort scaffolds by length and select first 10
            $ format-gff GFF --sort-len --subset 10 > Chr1-10.gff

            # relabel sorted scaffold names using str subst, save map of old names to new
            $ format-gff GFF --sort-len --relabel {idx}_Chr{sidx}.{gidx} --relabel-map scaff-map.tsv > relabeled.gff

            # filter to include only gene features
            $ format-gff GFF --features gene > genes.gff

            # subset to genes on Chr1, keep only ID attrs, and relabel IDs with spp prefix
            $ format-gff GFF --subset-names Chr1 --feat gene --attrs ID --relabel-ids spp_{id} > relabeled.gff

            # export GFF data to bed format and remove overlapping/duplicate features
            $ format-gff GFF --subset-names Chr1 --feat mRNA --to bed --filter-duplicates > Chr1-transcripts.bed

            # get all sequence IDs of mRNAs on Chr1
            $ format-gff GFF --subset-names Chr1 --feat mRNA --attrs ID --to ids > Chr1-ids.txt
        """)
    )

    # create parser or connect as subparser to cli parser
    if parser:
        KWARGS['name'] = KWARGS.pop("prog")
        parser = parser.add_parser(**KWARGS)
    else:
        KWARGS.pop("help")
        parser = ArgumentParser(**KWARGS)

    # add arguments
    parser.add_argument("GFF", type=Path, help="a gff/gtf annotation table (can be .gz)")

    sort_options = parser.add_mutually_exclusive_group()
    sort_options.add_argument("--sort-alpha", action="store_true", help="sort scaffolds by name alphanumerically")
    sort_options.add_argument("--sort-len", action="store_true", help="sort scaffolds by length (longest to shortest)")

    # relabel scaffolds. If you relabel in GFF you need to relabel same way in FA file using format-fasta --relabel-map MAP
    relabel_options = parser.add_mutually_exclusive_group()
    relabel_options.add_argument("--relabel", type=str, metavar="str", default="", help=r"relabel scaffs using str subst w/ {scaff},{sidx}. Applies after sorting/subset")
    relabel_options.add_argument("--relabel-names", type=str, metavar="str", default="", help=r"relabel scaffs with list of new names. Applies after sorting/subset")
    relabel_options.add_argument("--relabel-map", type=str, metavar="str", default="", help=r"optional path to write tsv mapping old scaff names to new")

    # select a subset of scaffolds
    sub_options = parser.add_mutually_exclusive_group()
    sub_options.add_argument("--subset", type=int, metavar="int", default=None, help="subselect the first N scaffolds. Applies after sorting")
    sub_options.add_argument("--subset-idx", type=int, metavar="int", nargs="+", help="subselect one or more scaffolds by 1-based index. Applies after sorting")
    sub_options.add_argument("--subset-names", type=str, metavar="str", nargs="+", help="subselect one or more scaffolds by name. Applies before relabeling/sorting")

    parser.add_argument("--features", type=str, metavar="str", nargs="+", help="subselect one or more features (e.g., gene) to keep (default=all)")
    parser.add_argument("--attrs", type=str, metavar="str", nargs="+", help="subselect one or more attributes (e.g., ID) to keep (default=all)")

    # relabel IDs. If you relabel in GFF you need to relabel same way in PEP, CDS files using format-fasta --relabel-map MAP
    parser.add_argument("--relabel-ids", type=str, metavar="str", default="", help=r"relabel attr IDs using str subst w/ {scaff},{sidx},{gidx},{idx},{type}")
    parser.add_argument("--relabel-ids-map", type=str, metavar="str", default="", help=r"optional path to write tsv mapping old IDS to new")

    parser.add_argument("--exclude-comments", action="store_true", help="exclude comments from the output")
    parser.add_argument("--add-to-comments", action="store_true", help="add this cmd to comments section")
    parser.add_argument("--filter-duplicates", action="store_true", help="keep only the highest score, or first, of duplicates coordinate features")

    parser.add_argument("--to", type=str, metavar="{gff, bed, ids}", default="gff", help="output format (default=gff)")
    # export_options.add_argument("--to-gff2", action="store_true", help="write BED format selected features")
    # export_options.add_argument("--to-gff3", action="store_true", help="write BED format selected features")
    # export_options.add_argument("--to-gtf", action="store_true", help="write BED format selected features")

    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "EXCEPTION"], metavar="level", default="INFO", help="stderr logging level (DEBUG, INFO, WARNING, ERROR; default=INFO)")
    return parser


def get_parser_genome_table(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for relabel tool.
    """
    kwargs = dict(
        prog="genome-table",
        usage="twig genome-table FA GFF [options]",
        help="write table with chrom/scaff [names lengths ngenes]",
        formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),
        description=dedent("""
            -------------------------------------------------------------------
            | genome-table: get table of chrom/scaff [names, lengths, ngenes] |
            -------------------------------------------------------------------
            | Parses the genome fasta (and annotation) files to extract       |
            | scaffold names, lengths, and gene numbers and write as a tsv    |
            | to stdout. This table is useful for downstream analysis such as |
            | dotplots and other synteny analyses. If an annotation file is   |
            | not provided then only names and lengths are extracted. Options |
            | can be used to sort and/or subselect scaffolds for inclusion.   |
            -------------------------------------------------------------------
        """),
        epilog=dedent("""
            Examples
            --------
            $ genome-table FA > genome.info.tsv
            $ genome-table FA GFF > genome.info.tsv
            $ genome-table FA GFF --sort-len --subset 8 > genome.info.tsv
            $ genome-table FA GFF --subset-names chr1 chr2 > genome.info.tsv
        """)
    )

    # create parser or connect as subparser to cli parser
    if parser:
        kwargs['name'] = kwargs.pop("prog")
        parser = parser.add_parser(**kwargs)
    else:
        kwargs.pop("help")
        parser = ArgumentParser(**kwargs)

    # add arguments
    parser.add_argument("genome", type=Path, help="fasta genome sequence (can be .gz)")
    parser.add_argument("annotation", type=Path, nargs="?", default=None, help="gff annotation file (optional; can be .gz)")
    sort_options = parser.add_mutually_exclusive_group()
    sort_options.add_argument("--sort-alpha", action="store_true", help="sort chromosomes by name alphanumerically")
    sort_options.add_argument("--sort-len", action="store_true", help="sort chromosomes by length")
    # parser.add_argument("--subset", type=int, metavar="int", default=None, help="write only the first N scaffolds after sorting")
    sub_options = parser.add_mutually_exclusive_group()
    sub_options.add_argument("--subset", type=int, metavar="int", default=None, help="subselect the first N scaffolds after optional sorting")
    sub_options.add_argument("--subset-idx", type=int, metavar="int", nargs="+", help="subselect one or more scaffolds by 1-based index after optional sorting")
    sub_options.add_argument("--subset-names", type=str, metavar="str", nargs="+", help="subselect one or more scaffolds by name")

    parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "EXCEPTION"], metavar="level", default="INFO", help="stderr logging level (DEBUG, INFO, WARNING, ERROR; default=INFO)")
    return parser


def get_parser_align_pre(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for relabel tool.
    """
    KWARGS = dict(
        prog="align-pre",
        usage="align-pre -i CDS -o OUT.nt.fa [options]",
        help="trim/filter CDS and prune isoforms before alignment",
        formatter_class=lambda prog: SingleMetavarHelpFormatter(prog, width=220, max_help_position=56),
        description=dedent("""
            -------------------------------------------------------------------
            | align-pre: trim/filter CDS and prune isoforms
            -------------------------------------------------------------------
            | Runs MACSE trimNonHomologousFragments (unless skipped), then
            | prunes isoforms by selecting one best sequence per isoform group.
            | Isoform groups can be defined globally or per taxon from names.
            -------------------------------------------------------------------
        """),
        epilog=dedent("""
            Examples
            --------
            $ twig align-pre -i CDS -o OUT.nt.fa
            $ twig align-pre -i CDS -o OUT.nt.fa -x '_' 2
            $ twig align-pre -i CDS -o OUT.nt.fa -X isoform_rules.tsv -d '_' -di 0
            $ twig align-pre -i CDS -o OUT.nt.fa -hf 0.3 -hi 0.8 -hc 10
            $ twig align-pre -i CDS -o OUT.nt.fa -A -e '^outgroup_'
        """)
    )

    # create parser or connect as subparser to cli parser
    if parser:
        KWARGS['name'] = KWARGS.pop("prog")
        parser = parser.add_parser(**KWARGS)
    else:
        KWARGS.pop("help")
        parser = ArgumentParser(**KWARGS)

    io = parser.add_argument_group("Input/Output")
    io.add_argument("-i", "--input", type=Path, metavar="path", required=True, help="input CDS fasta (aligned or unaligned)")
    io.add_argument("-o", "--outpath", type=Path, metavar="path", required=True, help="output path for retained nucleotide fasta")
    io.add_argument("-e", "--exclude", type=str, metavar="str", nargs="*", help="exclude names matching one or more regex patterns")
    io.add_argument("-s", "--subsample", type=str, metavar="str", nargs="*", help="keep only names matching one or more regex patterns")

    trim = parser.add_argument_group("MACSE Trimming")
    trim.add_argument("-A", "--skip-macse-trim", action="store_true", help="skip MACSE trimming; use existing .trim_info if available")
    trim.add_argument("-hf", "--min-homology-full", type=float, metavar="float", default=0.1, help="min full-length homology required to keep a sequence [%(default)s]")
    trim.add_argument("-hi", "--min-homology-internal", type=float, metavar="float", default=0.5, help="min internal homology required to keep a sequence [%(default)s]")
    trim.add_argument("-hc", "--min-homology-coverage", type=int, metavar="int", default=3, help="min sequence count meeting homology thresholds [%(default)s]")
    trim.add_argument("-ti", "--min-trim-length-homology-internal", type=int, metavar="int", default=50, help="trim short internal low-homology fragments below this length [%(default)s]")
    trim.add_argument("-te", "--min-trim-length-homology-external", type=int, metavar="int", default=50, help="trim short edge low-homology fragments below this length [%(default)s]")
    trim.add_argument("-mm", "--mem-length", type=int, metavar="int", default=6, help="amino-acid MEM length used by MACSE homology scoring [%(default)s]")

    pruning = parser.add_argument_group("Isoform Pruning and Scope")
    pruning.add_argument("-x", "--isoform-group-split", metavar=("str", "int"), nargs=2, default=None, help="global grouping rule: split names on <delim>; keep text left of idx-th delimiter")
    pruning.add_argument("-X", "--isoform-group-rules", type=Path, metavar="path", help=r"per-taxon TSV (<taxon> <delim> <idx>); group by text left of that idx-th delimiter")
    pruning.add_argument("-d", "--taxon-delim", type=str, metavar="str", help="optional taxon-scope delimiter before isoform grouping; affects per-taxon grouping only")
    pruning.add_argument("-di", "--taxon-delim-idxs", type=int, metavar="int", nargs="+", default=[0], help="token indices used to build taxon scope keys after --taxon-delim split")
    pruning.add_argument("-dj", "--taxon-delim-join", type=str, metavar="str", default="-", help="join string for taxon scope keys (used only when --taxon-delim is set)")
    pruning.add_argument("-ml", "--min-retained-nt-length", type=int, metavar="int", default=0, help="exclude retained isoform reps shorter than this nt length [%(default)s]")
    pruning.add_argument("-mc", "--min-retained-seqs", type=int, metavar="int", default=0, help="write output only if this many sequences remain after isoform pruning [%(default)s]")

    runtime = parser.add_argument_group("Runtime")
    runtime.add_argument("-v", "--verbose", action="store_true", help="print MACSE progress to stderr")
    runtime.add_argument("-f", "--force", action="store_true", help="overwrite existing output files")
    runtime.add_argument("-k", "--keep", action="store_true", help="keep temporary files for debugging")
    runtime.add_argument("-l", "--log-level", type=str, metavar="level", default="INFO", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    return parser


def get_parser_align_cds(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for relabel tool.
    """
    KWARGS = dict(
        prog="align-cds",
        usage="align-cds -i CDS -o OUT.msa.nt.fa [options]",
        help="align CDS using MACSE codon-aware alignment",
        formatter_class=lambda prog: SingleMetavarHelpFormatter(prog, width=220, max_help_position=56),
        description=dedent("""
            -------------------------------------------------------------------
            | align-cds: run MACSE codon-aware CDS alignment
            -------------------------------------------------------------------
            | Calls `macse -prog alignSequences` and writes aligned CDS to
            | --outpath. The resulting aligned CDS result should next be run in
            | `twig align-post` to trim and filter the alignment.
            -------------------------------------------------------------------
        """),
        epilog=dedent("""
            Examples
            --------
            $ twig align-cds -i CDS -o OUT.msa.nt.fa
            $ twig align-cds -i CDS -o OUT.msa.nt.fa -m 100
            $ twig align-cds -i CDS -o OUT.msa.nt.fa --seq-lr '^pseudo_' '^partial_'
            $ twig align-cds -i CDS -o OUT.msa.nt.fa --seq-lr seq_lr_patterns.txt
            $ twig align-cds -i CDS -o OUT.msa.nt.fa -v
        """)
    )
    # create parser or connect as subparser to cli parser
    if parser:
        KWARGS['name'] = KWARGS.pop("prog")
        parser = parser.add_parser(**KWARGS)
    else:
        KWARGS.pop("help")
        parser = ArgumentParser(**KWARGS)

    io = parser.add_argument_group("Input/Output")
    io.add_argument("-i", "--input", type=Path, metavar="path", required=True, help="input CDS fasta (aligned or unaligned)")
    io.add_argument("-o", "--outpath", type=Path, metavar="path", required=True, help="output path for aligned nucleotide fasta")

    align = parser.add_argument_group("Alignment")
    align.add_argument("-m", "--max-refine-iter", type=int, metavar="int", default=-1, help="maximum MACSE refinement iterations (-1 uses MACSE default) [%(default)s]")
    align.add_argument(
        "-s",
        "--seq-lr",
        type=str,
        metavar="str",
        nargs="+",
        help=(
            "designate set of 'less-reliable' seqs by name/regex or in file path"
            # .selectors passed to MACSE -seq_lr; "
            # "provide regex patterns directly or one selector-file path"
        ),
    )

    runtime = parser.add_argument_group("Runtime")
    runtime.add_argument("-v", "--verbose", action="store_true", help="print MACSE progress to stderr")
    runtime.add_argument("-f", "--force", action="store_true", help="overwrite existing output files")
    runtime.add_argument("-l", "--log-level", type=str, metavar="level", default="INFO", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    return parser


# def get_parser_align_post_trim(parser: ArgumentParser | None = None) -> ArgumentParser:
#     """Return a parser for relabel tool.
#     """
#     KWARGS = dict(
#         prog="align-post-trim",
#         usage="align-post-trim -i CDS [options]",
#         help="trim/filter/convert alignments",
#         formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),
#         description=dedent("""
#             -------------------------------------------------------------------
#             | align-post-trim: trim aligned CDS of poorly aligned sites or seqs
#             -------------------------------------------------------------------
#             | Trims rows and columns of an alignment while retaining coding frame.
#             | Runs 'macse -prog TrimAlignment' to trim terminal CDS edges by min
#             | coverage; runs trimal on the AA translation to trim internal gappy
#             | columns and sequences with too few high quality sites; and runs
#             | terrace filter to iteratiely remove sequences with too little overlap
#             | with other sequences. The filtered alignment is written to both
#             | .nt and .aa with frameshift and stop codons masked.
#             -------------------------------------------------------------------
#         """),
#         epilog=dedent("""
#             Examples
#             --------
#             $ twig align-post-trim -i CDS -o PATH/PRE -ac 0.5
#             $ twig align-post-trim -i CDS -o PATH/PRE -t NWK -R
#             $ twig align-post-trim -i CDS -o PATH/PRE

#             # run parallel jobs on many cds files
#             $ parallel -j 10 "twig align-post-trim -i {} ..." ::: CDS/*.msa.nt.fa

#             # full pipeline
#             $ twig align-post-trim -i CDS -o TRIM    # {TRIM}
#             $ twig align-post-trim -i TRIM -o ALN    # {ALN}
#             $ twig align-post-trim -i ALN -o MSA     # {MSA}.nt.fa, {MSA}.aa.fa
#         """)
#     )

#     # create parser or connect as subparser to cli parser
#     if parser:
#         KWARGS['name'] = KWARGS.pop("prog")
#         parser = parser.add_parser(**KWARGS)
#     else:
#         KWARGS.pop("help")
#         parser = ArgumentParser(**KWARGS)

#     # path args
#     parser.add_argument("-i", "--input", type=Path, metavar="path", required=True, help="input aligned CDS")
#     parser.add_argument("-o", "--outprefix", type=Path, metavar="path", required=True, help="out prefix; creates {prefix}.nt.fa and {prefix}.aa.fa")
#     parser.add_argument("-e", "--exclude", type=str, metavar="str", nargs="*", help="optional names or glob to exclude one or more sequences")
#     parser.add_argument("-s", "--subsample", type=str, metavar="str", nargs="*", help="optional names or glob to include only a subset sequences")
#     parser.add_argument("-t", "--tree", type=Path, metavar="path", help="optional newick file to subsample genes present in tree")
#     # options
#     parser.add_argument("-mo", "--min-overlap", type=int, metavar="int", default=0, help="min sites shared by each pair of samples, else iter prune lowest [%(default)s]")
#     parser.add_argument("-ms", "--min-samples", type=int, metavar="int", default=0, help="min samples after filtering else error [%(default)s]")
#     # parser.add_argument("-ml", "--min-length", type=int, metavar="int", default=0, help="min length of non-missing sequence in a sample [%(default)s]")
#     parser.add_argument("-ac", "--aln-trim-ends-min-coverage", type=float, metavar="float", default=0.4, help="trim alignment edges to where a min percent of samples have data [%(default)s]")
#     parser.add_argument("-as", "--aln-trim-window-size", type=int, metavar="int", default=5, help="trim alignment edges using a sliding 'half_window_size' [%(default)s]")
#     parser.add_argument("-if", "--codon-int-fs", type=str, metavar="str", default="NNN", help="codon to sub for internal frame shift [NNN]")
#     parser.add_argument("-ef", "--codon-ext-fs", type=str, metavar="str", default="NNN", help="codon to sub for external frame shift [NNN]")
#     parser.add_argument("-fs", "--codon-final-stop", type=str, metavar="str", default="NNN", help="codon to sub for final stop [NNN]")
#     parser.add_argument("-is", "--codon-int-stop", type=str, metavar="str", default="NNN", help="codon to sub for internal stop [NNN]")

#     parser.add_argument("-r", "--refine-alignment", action="store_true", help="refine alignment")
#     parser.add_argument("-R", "--refine-alignment-if", action="store_true", help="refine alignment only if >=1 sequences are filtered out")
#     parser.add_argument("-ri", "--max-iter-refine-alignment", type=int, metavar="int", default=-1, help="max iterations in refine alignment [%(default)s]")

#     # others
#     parser.add_argument("-v", "--verbose", action="store_true", help="print macse progress info to stderr")
#     parser.add_argument("-f", "--force", action="store_true", help="overwrite existing result files in outdir")
#     parser.add_argument("-k", "--keep", action="store_true", help="keep tmp files (for debugging)")
#     parser.add_argument("-l", "--log-level", type=str, metavar="level", default="INFO", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
#     # parser.add_argument("-L", "--log-file", type=Path, metavar="path", help="append stderr log to a file")
#     return parser


def get_parser_align_post(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for relabel tool.
    """
    KWARGS = dict(
        prog="align-post",
        usage="align-post -i CDS -o OUT [options]",
        help="trim and filter alignments of gappy edges, columns, and sequences",
        formatter_class=lambda prog: SingleMetavarHelpFormatter(prog, width=220, max_help_position=56),
        description=dedent("""
            -------------------------------------------------------------------
            | align-post: trim gappy alignments and low overlap sequences
            -------------------------------------------------------------------
            | This method takes a CDS alignment from `twig align-cds`, first
            | runs MACSE trimAlignment to trim gappy edges, then runs trimal
            | to trim internal gappy columns (default algorithm: automated1).
            | A pairwise overlap filter can further remove low-overlap samples.
            | Final cleaned alignments are exported to NT and AA and remain in
            | coding frame.
            -------------------------------------------------------------------
        """),
        epilog=dedent("""
            Examples
            --------
            $ twig align-post -i CDS -o PATH/PRE
            $ twig align-post -i CDS -o PATH/PRE -ec 0.5 -ew 7
            $ twig align-post -i CDS -o PATH/PRE --skip-trimal
            $ twig align-post -i CDS -o PATH/PRE -t TREE.nwk -ms 8
        """)
    )

    # create parser or connect as subparser to cli parser
    if parser:
        KWARGS['name'] = KWARGS.pop("prog")
        parser = parser.add_parser(**KWARGS)
    else:
        KWARGS.pop("help")
        parser = ArgumentParser(**KWARGS)

    io = parser.add_argument_group("Input/Output")
    io.add_argument("-i", "--input", type=Path, metavar="path", required=True, help="input aligned CDS fasta")
    io.add_argument("-o", "--out-prefix", type=Path, metavar="path", required=True, help="output prefix; writes {prefix}.nt.fa and {prefix}.aa.fa")

    select = parser.add_argument_group("Selection")
    select.add_argument("-e", "--exclude", type=str, metavar="str", nargs="*", help="exclude names matching one or more regex patterns")
    select.add_argument("-s", "--subsample", type=str, metavar="str", nargs="*", help="keep only names matching one or more regex patterns")
    select.add_argument("-t", "--subsample-tree", type=Path, metavar="path", help="keep only names present in this Newick tree")

    edge = parser.add_argument_group("MACSE Edge Trimming")
    edge.add_argument("-ec", "--edge-min-coverage", type=float, metavar="float", default=0.4, help="trim alignment ends to minimum non-gap coverage [%(default)s]")
    edge.add_argument("-ew", "--edge-window-size", type=int, metavar="int", default=5, help="half-window size used for MACSE edge trimming [%(default)s]")

    trimal = parser.add_argument_group("Trimal Filtering")
    trimal.add_argument("-ta", "--trimal-algorithm", type=str, metavar="str", default="automated1", help="trimal auto mode ['automated1','gappyout','strict','strictplus']")
    trimal.add_argument("-tr", "--trimal-res-overlap", type=float, metavar="float", default=0.25, help="minimum non-gap proportion to keep a column [%(default)s]")
    trimal.add_argument("-ts", "--trimal-seq-overlap", type=float, metavar="float", default=0.25, help="minimum kept-site proportion required per sequence [%(default)s]")
    trimal.add_argument("-tx", "--skip-trimal", action="store_true", help="skip trimal column/sequence filtering")

    retain = parser.add_argument_group("Overlap and Retention")
    retain.add_argument("-mo", "--min-pair-overlap", type=int, metavar="int", default=10, help="minimum non-missing site overlap required for every pair [%(default)s]")
    retain.add_argument("-ms", "--min-retained-seqs", type=int, metavar="int", default=4, help="require at least this many sequences after filtering [%(default)s]")
    retain.add_argument("-ml", "--min-retained-nt-length", type=int, metavar="int", default=25, help="require final alignment length >= this nt length [%(default)s]")
    # parser.add_argument("-if", "--codon-int-fs", type=str, metavar="str", default="NNN", help="codon to sub for internal frame shift [NNN]")
    # parser.add_argument("-ef", "--codon-ext-fs", type=str, metavar="str", default="NNN", help="codon to sub for external frame shift [NNN]")
    # parser.add_argument("-fs", "--codon-final-stop", type=str, metavar="str", default="NNN", help="codon to sub for final stop [NNN]")
    # parser.add_argument("-is", "--codon-int-stop", type=str, metavar="str", default="NNN", help="codon to sub for internal stop [NNN]")

    runtime = parser.add_argument_group("Runtime")
    runtime.add_argument("-v", "--verbose", action="store_true", help="print tool progress to stderr")
    runtime.add_argument("-f", "--force", action="store_true", help="overwrite existing output files")
    runtime.add_argument("-k", "--keep", action="store_true", help="keep temporary files for debugging")
    runtime.add_argument("-l", "--log-level", type=str, metavar="level", default="INFO", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    # parser.add_argument("-L", "--log-file", type=Path, metavar="path", help="append stderr log to a file")
    return parser



def get_parser_tree_filter(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for relabel tool.
    """
    KWARGS = dict(
        prog="tree-filter",
        usage="tree-filter [options]",
        help="filter and relabel a set of trees for downstream analyses",
        formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=140, max_help_position=140),
        description=dedent("""
            -------------------------------------------------------------------
            | tree-filter: Filter/Relabel a set of trees for downstream analyses
            -------------------------------------------------------------------
            | The order of [optional] operations is (1) parse names from delim
            | split tip labels; (2) subsample and/or relabel names by imap;
            | (3) exclude outlier edges; (4) collapse outgroups; (5) exclude
            | trees w/ < min-tips; (6) exclude trees w/ > min-copies in any name.
            | Note that if names are parsed by delim, then the imap should contain
            | the parsed names. The labels on output trees will be the same as in
            | the inputs unless --relabel-delim or --relabel-imap is used. But
            | parsing shorter names can be useful to match names to imap even
            | if keeping full names. Filtering stats are written to stderr, which
            | can be redirected to a log file with `2> log.txt`.
            -------------------------------------------------------------------
        """),
        epilog=dedent("""
            Examples
            --------
            # relabel by split-select-join on tip names
            $ twig tree-filter -i NWK -d '-' -di 0 2 -dj '-rd' > relabeled-trees.nwk

            # get only single-copy gene trees
            $ twig tree-filter -i NWK -c 1 > single-copy-trees.nwk

            # get only trees w/ >20 tips
            $ twig tree-filter -i NWK -m 20 > min20-trees.nwk

            # exclude outlier edges
            $ twig tree-filter -i NWK --exclude -ei 4 > cleaned-trees.nwk

            Examples using IMAP/MINMAP
            --------------------------
            # subsample to names in imap
            $ twig tree-filter -i NWK -I IMAP > subsample-trees.nwk

            # subsample and relabel to pop names in imap
            $ twig tree-filter -i NWK -I IMAP -ri > relabeled-trees.nwk

            # subsample, relabel, and keep only one outgroup (best to root before)
            $ twig tree-filter -i NWK -I IMAP -ri --collapse > final-trees.nwk

            Example on a many trees in different paths
            -------------------------------------------
            $ parallel "twig tree-filter -i {} ... > {}.filtered" ::: TREES/*.nwk

            Example IMAP                       Example MINMAP
            --------------------               ----------------
            sampleA     pop1                   pop1       1
            sampleB     pop1                   pop2       1
            sampleC     pop2                   outgroup   1
            sampleD     pop2
            sampleE     outgroup
            sampleF     outgroup
        """)
    )
    # create parser or connect as subparser to cli parser
    if parser:
        KWARGS['name'] = KWARGS.pop("prog")
        parser = parser.add_parser(**KWARGS)
    else:
        KWARGS.pop("help")
        parser = ArgumentParser(**KWARGS)

    # path i/o args
    parser.add_argument("-i", "--input", type=Path, metavar="path", required=True, help="newick or multi-newick trees file")
    parser.add_argument("-o", "--out", type=Path, metavar="path", help="outfile name else printed to stdout")

    # parsing names
    parser.add_argument("-d", "--delim", type=str, metavar="str", help="delimiter to split tip labels")
    parser.add_argument("-di", "--delim-idxs", type=int, metavar="int", nargs="+", default=[0], help="index of delimited name items to keep")
    parser.add_argument("-dj", "--delim-join", type=str, metavar="str", default="-", help="join character on delimited name items")

    # filter on names
    # parser.add_argument("-i", "--include", type=str, metavar="str", nargs="+", help="subset of names to keep")
    parser.add_argument("-I", "--imap", type=Path, metavar="path", help=r"filepath listing names to keep, or assigning (name population)")
    parser.add_argument("-M", "--minmap", type=Path, metavar="path", help=r"filepath listing (population mincov) for filters")

    # filter on tree data
    parser.add_argument("-m", "--min-tips", type=int, metavar="int", default=4, help="min tips after pruning [4]")
    parser.add_argument("-c", "--max-copies", type=int, metavar="int", help="filter trees with >c gene copies from a taxon [None]")
    parser.add_argument("-s", "--min-support", type=int, metavar="int", default=0, help="collapse edges w/ support <s into polytomies [0]")
    parser.add_argument("-e", "--min-edges", type=int, metavar="int", default=0, help="filter trees with <e non-zero edges [0]")
    parser.add_argument("-eo", "--edge-outlier-outgroup", type=float, metavar="float", default=10, help="exclude 'outgroup' population edges if >eo stdev from mean [10]")
    parser.add_argument("-ei", "--edge-outlier-ingroup", type=float, metavar="float", default=5, help="exclude non 'outgroup' population edges if >ei stdev from mean [5]")

    # actions
    parser.add_argument("-rd", "--relabel-delim", action="store_true", help="relabel tips by their delim parsed names")
    parser.add_argument("-ri", "--relabel-imap", action="store_true", help="relabel tips to their imap mapped names")
    # parser.add_argument("-x", "--nexus", action="store_true", help="export in NEXUS format. Retains path/tree names")
    parser.add_argument("-S", "--subsample", action="store_true", help="subsample to include only tips in imap")
    parser.add_argument("-E", "--exclude-outliers", action="store_true", help="exclude tips with outlier edge lengths (>ei or >eo)")
    parser.add_argument("-R", "--require-outgroups", action="store_true", help="require at least one 'outgroup' sample")
    parser.add_argument("-C", "--collapse-outgroups", action="store_true", help="keep only the most distant 'outgroup' (assumes rooted trees)")
    parser.add_argument("-O", "--outgroups", type=Path, metavar="path", help=r"path listing outgroups for -R (overrides IMAP group 'outgroup')")
    parser.add_argument("-l", "--log-level", type=str, metavar="level", default="INFO", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    # parser.add_argument("-L", "--log-file", type=Path, metavar="path", help="append stderr log to a file")
    return parser


def get_parser_tree_rooter(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for relabel tool.
    """
    KWARGS = dict(
        prog="tree-rooter",
        usage="tree-rooter [options]",
        help="Re-root gene trees by outgroup, species tree, and/or MAD",
        formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=140, max_help_position=140),
        description=dedent("""
            -------------------------------------------------------------------
            | tree-rooter: re-root gene trees by outgroups, species tree or MAD
            -------------------------------------------------------------------
            | Root trees by designating one or more outgroups. If a species tree
            | is provided an iterative process is used to try to root on ordered
            | outgroup subclades. If the ingroup (all non-outgroup samples) is
            | not monophyletic after attempted rooting, the tree is not returned.
            | Stats are printed to stderr. Input can be a multinewick with many
            | trees in it. Use -x to return only the trees that could not be
            | rooted, which can be rooted by alternative methods, like --mad,
            | or manually.
            -------------------------------------------------------------------
        """),
        epilog=dedent("""
            Examples
            --------
            # root trees on a specific outgroup
            $ twig tree-rooter -t NWK -r A > rooted-trees.nwk

            # root trees on an ordered set of options of outgroup clades
            $ twig tree-rooter -t NWK -r A B A,B > rooted-trees.nwk

            # root trees using minimal ancestor deviation (MAD)
            $ twig tree-rooter -t NWK --mad > rooted-trees.nwk

            # root trees on outgroups if present and mad root the rest
            $ twig tree-rooter -t NWK -r A B --mad > rooted-trees.nwk

            # root trees on outgroups in order of a rooted species tree
            $ twig tree-rooter -t NWK -s SPTREE > rooted-trees.nwk

            # root trees on a restricted set of outgroups from a rooted species tree
            $ twig tree-rooter -t NWK -s SPTREE -r Z Y > rooted-trees.nwk

            # root trees on outgroups defined in IMAP given a rooted species tree
            $ twig tree-rooter -t NWK -s SPTREE -I IMAP > rooted-trees.nwk

            # return trees that could not be rooted given the options
            $ twig tree-rooter -t NWK -s SPTREE -x > unrooted-trees.nwk

            # parallel process many files individually
            $ parallel -j 10 "twig tree-rooter -i {} -s SPTREE > ROOTED/{}" ::: TREES/*.nwk
        """)
    )
    # create parser or connect as subparser to cli parser
    if parser:
        KWARGS['name'] = KWARGS.pop("prog")
        parser = parser.add_parser(**KWARGS)
    else:
        KWARGS.pop("help")
        parser = ArgumentParser(**KWARGS)

    # path i/o args
    parser.add_argument("-i", "--input", type=Path, metavar="path", required=True, help="newick or multi-newick trees file")
    parser.add_argument("-o", "--out", type=Path, metavar="path", help="outfile name else printed to stdout")

    # parsing names
    parser.add_argument("-d", "--delim", type=str, metavar="str", help="delimiter to split tip labels")
    parser.add_argument("-di", "--delim-idxs", type=int, metavar="int", nargs="+", default=[0], help="index of delimited name items to keep")
    parser.add_argument("-dj", "--delim-join", type=str, metavar="str", default="-", help="join character on delimited name items")

    # rooting options
    parser.add_argument("-s", "--sptree", type=Path, metavar="path", help="rooted species tree")
    parser.add_argument("-r", "--outgroups", type=str, metavar="str", nargs="+", help="list outgroup tip labels")
    parser.add_argument("-I", "--imap", type=Path, metavar="path", help="get outgroup assignment from tabular imap file w/ rows '{sample}\toutgroup'")

    # return option
    # parser.add_argument("-a", "--allow", action="store_true", help="allow ")
    parser.add_argument("-m", "--mad", action="store_true", help="use minimal ancestor deviation to estimate root of remaining unrooted trees")
    parser.add_argument("-x", "--not-rooted", action="store_true", help="return the trees that could not be rooted given the options")
    parser.add_argument("-rd", "--relabel-delim", action="store_true", help="relabel tips by their delim parsed names")
    # logging
    parser.add_argument("-l", "--log-level", type=str, metavar="level", default="INFO", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    # parser.add_argument("-L", "--log-file", type=Path, metavar="path", help="append stderr log to a file")
    return parser


def get_parser_tree_skeleton(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for relabel tool.
    """
    KWARGS = dict(
        prog="tree-skeleton",
        usage="tree-skeleton [options]",
        help="Force gene tree to match species tree topology",
        formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=140, max_help_position=140),
        description=dedent("""
            -------------------------------------------------------------------
            | tree-skeleton
            -------------------------------------------------------------------
            | ...
            -------------------------------------------------------------------
        """),
        epilog=dedent("""
            Examples
            --------
            $ twig tree-skeleton -i NWK1 -s NWK2 > NWK3
            $ twig tree-skeleton -i NWK1 -s NWK2 -d "|" -di 0 > NWK3
        """)
    )
    # create parser or connect as subparser to cli parser
    if parser:
        KWARGS['name'] = KWARGS.pop("prog")
        parser = parser.add_parser(**KWARGS)
    else:
        KWARGS.pop("help")
        parser = ArgumentParser(**KWARGS)

    # path i/o args
    parser.add_argument("-i", "--input", type=Path, metavar="path", required=True, help="newick or multi-newick gene trees file")
    parser.add_argument("-s", "--sptree", type=Path, metavar="path", required=True, help="newick species tree file")
    parser.add_argument("-o", "--out", type=Path, metavar="path", help="outfile name else printed to stdout")

    # parsing names
    parser.add_argument("-d", "--delim", type=str, metavar="str", help="delimiter to split tip labels")
    parser.add_argument("-di", "--delim-idxs", type=int, metavar="int", nargs="+", default=[0], help="index of delimited name items to keep")
    parser.add_argument("-dj", "--delim-join", type=str, metavar="str", default="-", help="join character on delimited name items")

    # ...
    parser.add_argument("-rd", "--relabel-delim", action="store_true", help="relabel tips by their delim parsed names")
    parser.add_argument("-l", "--log-level", type=str, metavar="level", default="CRITICAL", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    return parser


def get_parser_partition_cds(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for relabel tool.
    """
    kwargs = dict(
        prog="partition-cds",
        usage="partition-cds [options]",
        help="write a partition file for phylogenetic inference",
        formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),
        description=dedent("""
            -------------------------------------------------------------------
            | partition-cds: write a partition file for phylogenetic analysis |
            -------------------------------------------------------------------
            | Write a partition file for a CDS alignment to specify substitution
            | models for phylogenetic inference based on codon position and the
            | length of the CDS. Writes to stdout.
            -------------------------------------------------------------------
        """),
        epilog=dedent("""
            Examples
            --------
            $ twig partition-cds -i ALN.fa -n 2 -m GTR+G GTR+G > ALN.partition
            $ twig partition-cds -i ALN.fa -n 3 -m GTR+G GTR+G GTR+G > ALN.partition
        """)
    )

    # create parser or connect as subparser to cli parser
    if parser:
        kwargs['name'] = kwargs.pop("prog")
        parser = parser.add_parser(**kwargs)
    else:
        kwargs.pop("help")
        parser = ArgumentParser(**kwargs)

    # add arguments
    parser.add_argument("-i", "--input", type=Path, metavar="path", required=True, help="fasta CDS sequence")
    parser.add_argument("-n", "--number", choices=[2, 3], type=int, metavar="int", default=3, help="number of partitions (2 or 3) [3]")
    parser.add_argument("-m", "--model", type=str, metavar="str", nargs="+", default=["GTR+G"], help="subst model for each position (1+2,3) or (1,2,3) for 2 or 3 partitions")
    parser.add_argument("-l", "--log-level", type=str, metavar="level", default="INFO", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    return parser


def get_parser_filter_concat(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for relabel tool.
    """
    KWARGS = dict(
        prog="filter-concat",
        usage="filter-concat -i LOCI.fa [...] -o OUT [options]",
        help="filter loci by IMAP/single-copy and write a concatenated supermatrix",
        formatter_class=lambda prog: SingleMetavarHelpFormatter(prog, width=220, max_help_position=56),
        description=dedent("""
            -------------------------------------------------------------------
            | filter-concat: filter loci and concatenate to .fa/.phy/.nex
            -------------------------------------------------------------------
            | Reads aligned fasta loci, parses taxon names from sequence labels,
            | excludes non-single-copy loci, filters loci by sample coverage,
            | and concatenates retained loci into one supermatrix. If IMAP is
            | provided, filters can be applied over IMAP groups/populations.
            -------------------------------------------------------------------
        """),
        epilog=dedent("""
            Examples
            --------
            $ twig filter-concat -i loci/*.fa -o out.fa
            $ twig filter-concat -i loci/*.fa -I imap.tsv -o out.fa
            $ twig filter-concat -i A.fa B.fa C.fa -I imap.tsv -o out.phy -F phylip
            $ twig filter-concat -i loci/*.fa -I imap.tsv -M minmap.tsv -o out.fa
            $ twig filter-concat -i loci/*.fa -I imap.tsv -o out.nex -F nexus -d "|" -di 0
            $ twig filter-concat -i loci/*.fa -o out.fa --partition-file out.partitions.txt
        """)
    )
    # create parser or connect as subparser to cli parser
    if parser:
        KWARGS['name'] = KWARGS.pop("prog")
        parser = parser.add_parser(**KWARGS)
    else:
        KWARGS.pop("help")
        parser = ArgumentParser(**KWARGS)

    io = parser.add_argument_group("Input/Output")
    io.add_argument("-i", "--input", type=Path, metavar="path", nargs="+", required=True, help="input aligned fasta locus files")
    io.add_argument("-I", "--imap", type=Path, metavar="path", help=r"IMAP TSV with rows: <name>\t<pop>")
    io.add_argument("-o", "--output", type=Path, metavar="path", required=True, help="output concatenated alignment file")
    io.add_argument("-F", "--format", choices=["fasta", "phylip", "nexus"], default="fasta", help="output format [%(default)s]")
    io.add_argument("-P", "--partition-file", type=Path, metavar="path", help="optional output file listing per-locus partition ranges")

    parsing = parser.add_argument_group("Taxon Parsing")
    parsing.add_argument("-d", "--delim", type=str, metavar="str", help="delimiter used to parse taxon names from sequence labels")
    parsing.add_argument("-di", "--delim-idxs", type=int, metavar="int", nargs="+", default=[0], help="token indices retained after split by --delim")
    parsing.add_argument("-dj", "--delim-join", type=str, metavar="str", default="-", help="join string used for parsed taxon names")

    filt = parser.add_argument_group("Locus Filters")
    mx = filt.add_mutually_exclusive_group()
    mx.add_argument("-m", "--min-present", type=int, metavar="int", default=4, help="minimum taxa (or IMAP groups) required for a retained locus [%(default)s]")
    mx.add_argument("-M", "--minmap", type=Path, metavar="path", help=r"TSV per-pop minima: <imap-pop>\t<min-samples> (used instead of -m)")

    runtime = parser.add_argument_group("Runtime")
    runtime.add_argument("-f", "--force", action="store_true", help="overwrite existing output file")
    runtime.add_argument("-l", "--log-level", type=str, metavar="level", default="INFO", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    return parser


# ...


# the tree labeling commands define their parsers in twig.labeling.parsers
from twig.labeling.parsers import get_parser_name_map, get_parser_sample_reduce, get_parser_tree_label, get_parser_tree_label_transfer  # noqa: E402,F401
