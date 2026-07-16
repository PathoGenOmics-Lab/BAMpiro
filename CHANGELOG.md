# Changelog

All notable changes to this project are documented in this file. The format is
based on [Keep a Changelog](https://keepachangelog.com/), and the project follows
[Semantic Versioning](https://semver.org/).

## [1.0.1] - 2026-07-16

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
    variable genes, and aDNA damage.
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
- **Documentation** - an [mdBook](https://rust-lang.github.io/mdBook/) under `docs/`
  (introduction, installation, quick start, configuration, QC report, Pathotypr,
  outputs) with the README slimmed to a landing page.

### Changed

- **Rebranded to BAMpiro** (larger report UI; panels with no input hide themselves).
- The container is **pinned to a digest** and its tool versions are recorded, for
  reproducible runs.
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
