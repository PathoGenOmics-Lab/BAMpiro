# Introduction

**BAMpiro** is a modular, containerized bioinformatics pipeline built with
**Nextflow (DSL2)**. It is optimized by default for *Mycobacterium tuberculosis*
(TB), but its architecture is **organism-agnostic** and works on **any bacterial
genome** (e.g. *E. coli*, *Salmonella*, *Staphylococcus*) by adjusting a few
parameters.

It automates the whole path from raw reads to annotated variants, consensus
sequences, and a consolidated interactive quality-control report.

## Key features

- **Universal bacterial support** - any reference genome + GFF annotation.
- **Automated reference prep** - indexes genomes and builds SnpEff databases on the fly.
- **Repeat masking** - `nucmer` (MUMmer4) auto-detects and excludes repetitive
  regions from variant calling, plus a length-aware `genmap` read filter.
- **Robust QC** - `FastP` cleaning and `Kraken2` taxonomic contamination checks.
- **Variant calling** - `FreeBayes` with customizable ploidy (1 or 2) and strict filtering.
- **Backbone generation** - "all-sites" VCFs (WT + variants) for phylogenetic supermatrices.
- **Lineage & drug-resistance typing** - alignment-free (k-mer) MTBC lineage + WHO
  drug-resistance genotyping with [Pathotypr](pathotypr.md), reference-agnostic and
  bundled in the container.
- **Interactive QC report** - a single self-contained HTML dashboard (20 panels)
  plus a machine-readable per-sample `qc_flags.tsv`, alongside the classic MultiQC
  report. See [Interactive QC Report](qc-report.md).
- **Dual amino-acid numbering** - optionally re-annotates variants against a
  canonical reference (H37Rv by default) so protein changes are shown in both
  numberings.

## Workflow summary

1. **Reference** - indexing + repeat masking + SnpEff DB building.
2. **QC** - read validation → Kraken2 (taxonomy) → FastP (trimming).
3. **Pathotypr** *(optional)* - alignment-free MTBC lineage (nested sub-lineage)
   **and** WHO drug-resistance typing from reads.
4. **Mapping** - `bwa-mem2` alignment → merge runs → mark duplicates (`samtools`) →
   length-aware read masking (`genmap`).
5. **Variants** - `FreeBayes` calling + `mpileup` for backbone generation.
6. **Consensus** - FASTA generation, masking low-coverage / low-quality sites.
7. **Annotation** - `SnpEff` annotation of main and legacy VCFs (+ optional
   canonical/H37Rv pass for [dual amino-acid numbering](pathotypr.md#dual-amino-acid-numbering-h37rv--mycobrowser)).
8. **Report** - aggregation into the classic MultiQC report **and** the interactive
   [HTML QC report](qc-report.md) with per-sample PASS/WARN/FAIL flags.

---

Next: [Installation](installation.md) · [Quick Start](quickstart.md)
