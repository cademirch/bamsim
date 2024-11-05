# Command-Line Help for `bamsim`

This document contains the help content for the `bamsim` command-line program.

**Command Overview:**

* [`bamsim`↴](#bamsim)
* [`bamsim generate`↴](#bamsim-generate)

## `bamsim`

Command-line interface for the BAM simulator

**Usage:** `bamsim [COMMAND]`

###### **Subcommands:**

* `generate` — Generate BAM files and truth files (default command)



## `bamsim generate`

Generate BAM files and truth files (default command)

**Usage:** `bamsim generate --num-bams <NUM_BAMS> --num-chrs <NUM_CHRS> --chr-length <CHR_LENGTH> --min-depth <MIN_DEPTH> --max-depth <MAX_DEPTH> --min-mean-depth <MIN_MEAN_DEPTH> --proportion <PROPORTION> <OUTDIR>`

###### **Arguments:**

* `<OUTDIR>` — Output directory for BAM files and the truth file

###### **Options:**

* `--num-bams <NUM_BAMS>` — Number of BAM files to generate
* `--num-chrs <NUM_CHRS>` — Number of chromosomes
* `--chr-length <CHR_LENGTH>` — Length of each chromosome
* `--min-depth <MIN_DEPTH>` — Minimum depth of coverage
* `--max-depth <MAX_DEPTH>` — Maximum depth of coverage
* `--min-mean-depth <MIN_MEAN_DEPTH>` — Minimum mean depth for intervals to be included in the truth file
* `--proportion <PROPORTION>` — Minimum proportion of BAM files covering an interval for it to be included in the truth file



<hr/>

<small><i>
    This document was generated automatically by
    <a href="https://crates.io/crates/clap-markdown"><code>clap-markdown</code></a>.
</i></small>

