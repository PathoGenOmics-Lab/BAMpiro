# Changelog

All notable changes to this project are documented in this file. The format is
based on [Keep a Changelog](https://keepachangelog.com/), and the project follows
[Semantic Versioning](https://semver.org/).

## [Unreleased] - targeting 1.1.0

Headline: a **test suite and CI**, and a pipeline that runs correctly on a machine
that is not the authors' cluster.

`manifest.version` already reads `1.1.0`: that is the number this work is heading
for, not a release that happened. Nothing is tagged and no release is published.
See [RELEASING.md](RELEASING.md) for the order the actual release has to follow.

The **container stays at the 1.0.1 image**, pinned by its digest. Nothing here
changes a bundled tool, and the scripts added under `bin/` are staged by Nextflow
from the pipeline directory rather than baked into the image, so there is nothing
to rebuild. A release only needs a new image when the Dockerfile changes.

### Added

- **Test suite** (`tests/`) - around 1,100 tests in three legs, all runnable with
  `tests/run_tests.sh` and none of them needing a container, a reference genome or
  a network connection:
 - `tests/unit/` - the Python under `bin/`: the consensus decision tree, the QC
    verdict engine, every parser, the k-mer liftover and the read filter.
 - `tests/js/` - the report's hand-written ES5 statistics (Spearman,
    Kruskal-Wallis, chi-square, Benjamini-Hochberg), checked against **SciPy and
    statsmodels** reference values rather than a snapshot of our own output, plus
    the integrity of the asset bundle.
 - `tests/pipeline/` - samplesheet validation and a full `-stub-run` of the DAG.
- **A `stub:` block for every process**, so `-stub-run` walks the whole workflow in
  seconds with no data. A test fails if a new process arrives without one.
- **Fixture cohort** (`tests/data/`, 170 kB) - a 5 kb reference with three genes and
  three samples covering paired-end, single-end and a sample split over two runs,
  all derived from one seed. CI regenerates it and fails on any diff.
- **`test` and `test_full` profiles** - the default feature set, and every optional
  branch on so a single run instantiates every process.
- **Continuous integration** (`.github/workflows/ci.yml`) - lint, unit tests on two
  Python versions, the front-end tests, and a stub run against both the minimum
  supported Nextflow and the current release, on every pull request.
- **`slurm` and `garnatxa` execution profiles**. `-profile garnatxa` replaces the
  hand-written `-c cluster.config` a user previously needed, and works when the
  pipeline is pulled straight from GitHub.
- **`--help`**, listing every parameter with its default.
- **`--max_cpus` / `--max_memory` / `--max_time`** - cap every process to what the
  machine can give, so an over-sized request (Kraken2 asks for 80 GB) fits a laptop
  without editing the modules.
- **Community files** - CONTRIBUTING, a Code of Conduct, a security policy, issue
  forms, a pull-request template, `RELEASING.md` and `.zenodo.json`.

### Changed

- **`-profile standard` no longer submits to SLURM.** It is now the portable default
  and runs on the current host, so a fresh clone works anywhere. Use `-profile slurm`
  for a generic cluster or `-profile garnatxa` for the I2SysBio one. **This changes
  what an existing `-profile standard` command does.**
- The base config no longer carries one site's paths: the SLURM QoS, the
  `module load singularity`, the Singularity `--bind` paths and the Kraken database
  location all moved into `conf/garnatxa.config`.
- **`--tsv` is required** and says what to pass; the default was a leftover
  samplesheet name. `--kraken2_db` has no default either, and without one the
  contamination screen is skipped.
- **The container is pinned by digest** rather than the mutable `1.0.1` tag, so two
  runs a month apart cannot silently use different tool versions.
- **An unrecognised `--parameter` now stops the run** with a suggestion, instead of
  being accepted and ignored. The declared set is read from `nextflow.config`, so
  the two cannot drift.
- The startup banner reports the resolved profile, executor, container and Kraken
  state, because where a run lands is the easiest thing to get wrong.

### Fixed

- **The QC report could not be produced on the default settings.** Its three optional
  inputs all fall back to a placeholder, and all three features are off by default,
  so every run handed `QC_REPORT` the same `NO_FILE` path three times and Nextflow
  refused to stage them (`input file name collision`). The placeholders are now
  distinct files under `assets/`.
- **`main.nf` did not compile under the strict Nextflow parser**, the default from
  25.10, and only ran through the deprecated v1 fallback. The samplesheet parsing is
  now a function, every statement lives inside the workflow, and each dynamic
  `publishDir` is a closure. Verified on 24.04.2 and 26.04.6.
- Removed `process.publishDirMode`, which is not a Nextflow directive and only
  produced a warning on every run.
- A blank line in a bedgraph raised `ValueError` in `build_min_unique_len.py`
  instead of being skipped.
- **A GFF attribute was matched as a substring**, so `locus_tag=` also matched
  inside `old_locus_tag=` and `gene=` inside `pseudogene=`, both routine in RefSeq
  and Prokka output. Genes were labelled with an obsolete tag, and the curated
  H37Rv gene map could be silently overridden.
- `mask_profile` summed raw interval lengths for its headline total, so an
  unmerged BED could report above 100% masked.
- Four parsers in `stats_to_legacy.py` wrapped a whole read loop in one
  `try/except`, so one malformed line silently discarded every line after it; a
  partial genome size then propagated into depth, breadth, evenness and every
  density bin.
- Genuinely-zero metrics were written as `NA`, making a failed sample
  indistinguishable from an unmeasured one.
- Genotype classification assumed diploid spellings, so a haploid no-call was
  counted as heterozygous. This matters for `-profile generic`, which is haploid.
- `safe_tabix` had drifted between its copies: the one in `annotation.nf` still
  used a `zgrep` check that cannot tell a truncated VCF from an empty one, and
  wrote a fake index without saying so.
- A missing `FORMAT/DP` reached the consensus as `ADP=.`, which reads as zero
  depth and turns a called SNP into a gap.

### Refactored

Behaviour-preserving throughout, and each step verified rather than assumed:

- **The science came out of the process scripts.** Around 500 lines of bash and
  awk lived inside `modules/*.nf`, where no test could reach it: a stub run
  replaces the script block and nothing can import it. The genotype
  re-validation, the backbone pileup parser, the hom/het bcftools expressions and
  `safe_tabix` now live in `bin/`, with 62 tests over them, including a run
  against a real bcftools. `variants.nf` drops from 485 lines to 344.
- **`qc_report.py` became a package.** 1697 lines and 53 functions are now a
  319-line CLI over `bin/qcreport/` (parsers, metrics, panels, render). The old
  and new versions produce byte-identical output on the demo cohort.
- **The workflow is one sub-workflow per stage.** `main.nf` goes from 556 lines to
  277, with the stage sequence reduced to eleven calls. Equivalence checked by
  diffing execution traces against the previous version: identical task counts and
  identical published outputs on both test profiles.

## [1.0.1] - 2026-07-19

Headline: a consolidated **interactive QC report**, alignment-free **lineage &
drug-resistance typing**, and a restructured documentation set. Rebranded to
**BAMpiro**.

### Added

- **Consolidated interactive QC report** (`bin/qc_report.py`) - a single
  self-contained HTML dashboard (no internet / CDN) built from the whole cohort,
  plus a machine-readable per-sample `qc_flags.tsv` (PASS/WARN/FAIL). Enabled with
  `--make_qc_report` (on by default). It bundles:
 - **SNP dynamics** - per-variant allele-frequency trajectories over time, gene-
    centric, with per-timepoint depth bars, a zoom control, and a series filter by
    samplesheet metadata.
 - **Epistasis** - pairs of variants with concordant / discordant dynamics,
    scored with a permutation *p*-value and Benjamini-Hochberg FDR, in card, matrix
    and table views.
 - **SNP matrix** - an explorable, downloadable site × sample matrix with the
    depth in each cell, samplesheet metadata as column-header levels, and column
    filtering by metadata (also written as a standalone TSV, `--make_snp_matrix`).
 - **Drug resistance** - a sample × drug matrix (worst WHO grade per drug) and a
    per-mutation table, driven by the Pathotypr WHO-catalogue calls.
 - Live-adjustable QC thresholds with presets, a collapsible left contents
    sidebar, per-section **(i)** info popovers, a sample exclusion basket, and
    panels for lineage summary, distributions, correlations, a metric-correlation
    heatmap, QC-space PCA, genome landscape, functional annotation, gene burden,
    variable genes, taxonomic composition (Kraken2: primary taxon / contaminants /
    unclassified, for a read-level contamination check), and aDNA damage.
 - **Report presentation & UX** - a **dark / light theme** (header toggle, follows
    the OS preference, persisted per viewer); a single homogeneous inline icon set;
    a **mobile-responsive** layout (capped frozen column, scrollable matrices, a
    dismissable contents drawer); a **triage-first** ordering (table sorted worst-QC
    first, the Flagged panel directly under the statistics with each failing margin
    shown inline); per-lineage colours from the canonical **mycolorsTB** palette;
    and the pipeline **version** + a **GitHub** link in the header and footer.
- **Alignment-free lineage & drug-resistance typing with Pathotypr**
  (`--run_pathotypr`) - now run **from the container** (bioconda build), replacing
  the external cluster binary. It types straight from the reads with diagnostic
  k-mers, so calls are **reference-agnostic**. Two `split-fastq` passes per sample
  (nested sub-lineage + WHO drug-resistance). The Zenodo marker panels + RF model
  (v1.0.0) and the MTBC-ancestor reference are bundled in the image under
  `/opt/pathotypr/`. `--pathotypr_min_alt` tunes the alt-allele cut-off for
  heteroresistant / minority DR alleles.
- **Dual amino-acid numbering** - an optional canonical re-annotation pass
  (`--annotate_canonical`, `--canonical_snpeff_db`, `--canonical_label`, default
  H37Rv) so the report shows every protein change in both the used-reference and
  the H37Rv / Mycobrowser numbering, with a link to each gene's Mycobrowser locus.
- **Length-aware read masking** - a per-reference `genmap` mappability track
  (cached and reused across runs) that drops reads too short to be uniquely placed;
  a repeat BED is emitted alongside.
- **Virgin (unmasked) consensus** output in parallel with the masked consensus.
- **CRAM output** for published alignments (`--output_cram`), a selectable
  **publish mode** (`copy` / `link`), and a trimmed default output footprint.
- **Optional region-parallel FreeBayes** for faster calling on deep samples.
- **Documentation** - a [MkDocs Material](https://squidfunk.github.io/mkdocs-material/)
  site under `docs/` (introduction, installation, quick start, configuration,
  troubleshooting, QC report, Pathotypr, outputs, downstream phylogeny, and an
  end-to-end Jupyter tutorial), with a GitHub Pages publishing workflow (activated once
  the repo is public) and the README slimmed to a landing page.
- **Metadata-driven cohort filtering & dose analyses** in the report - samplesheet
  categorical columns (e.g. `treatment`) drive a report-wide cohort filter; a numeric
  `dose` column becomes a first-class metric (scatter axes + correlation matrix, with a
  Spearman *r* / *p* read-out), plus two new panels: **Dose × treatment** (per-treatment
  dose distribution + a Kruskal-Wallis rank test) and **Variant × dose** (per-variant
  allele-frequency-vs-dose Spearman scan with Benjamini-Hochberg FDR).
- **Guided tutorial series** on the docs site - a step-by-step *Tutorials* section (first
  run, samplesheet, configuring a run, reading the QC report, lineage & DR typing, building
  a phylogeny, longitudinal analysis, running a non-TB organism) alongside the runnable
  Jupyter notebook.
- **Live example QC report** - a self-contained demo report published on the docs site so
  the interactive dashboard can be explored before running the pipeline.

### Changed

- **Rebranded to BAMpiro** (larger report UI; panels with no input hide themselves).
- The container is referenced by version **tag** and its tool versions are recorded;
  re-pin to the `@sha256:` digest for byte-for-byte reproducible runs.
- Right-sized per-process CPU / memory, streamed the deduplication pipeline, and
  produced lighter intermediates for faster `-resume`.
- Replaced the report's top navigation with the collapsible contents sidebar.
- Nextflow runtime artifacts (logs, cache, reports) are now git-ignored.

### Fixed

- Consensus generation: multiallelic reference leak, spanning-deletion gate,
  contig guard, and a chromosome-aware masker.
- Genotype revalidation (multiallelic ref-leak, het/hom), robust `mpileup` paste,
  multi-contig VCF headers, and the ploidy-1 reference-allele filter.
- Multi-reference stats wiring and sample-id sanitization.
- Stripped stray headers in the FreeBayes exclude step and masked contig-start
  repeats.
- ES5 compatibility slips in the report and the H37Rv snpEff database pre-fetch in
  the image; hardened the canonical-annotation channel and the drug-resistance
  collector; plus six further issues found in a deep audit (multi-contig handling,
  deep-coverage bins, and output-name collisions).

## [1.0.0]

Initial release. A modular, containerized Nextflow (DSL2) pipeline for bacterial
short-read analysis: automated reference indexing and SnpEff database building,
`nucmer` repeat masking, `FastP` + `Kraken2` read QC, `bwa-mem2` mapping with
duplicate marking and multi-run merging, `FreeBayes` variant calling with all-sites
backbone VCFs, consensus generation, `SnpEff` annotation, basic Pathotypr lineage
classification, and a MultiQC report.
