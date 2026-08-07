# Outputs

The pipeline organizes results by `sampleId`. Below is a breakdown using a sample
named `MP00091` mapped against reference `LENS`. Run/cohort-level deliverables land
at the top of `outdir` and are prefixed with the sample-sheet name (here
`samples_legio`); per-sample files live under their `sampleId` folder. With the
default `nested_output=true`, that folder is a **nested** path (e.g. `MP00091` →
`MP/00/09/1/`, not the flat `MP00091/` drawn below); set `nested_output=false` for the
flat `sampleId` layout shown here.

```text
results_bampiro/
├── samples_legio_qc_report.html          # 🌟 INTERACTIVE consolidated QC report (self-contained HTML)
├── samples_legio_qc_flags.tsv            # Per-sample PASS/WARN/FAIL verdicts
├── samples_legio_summary.tsv             # Cohort metrics table (feeds the report)
├── samples_legio_snp_matrix.tsv          # Master SNP matrix (site × sample AF & depth)
├── samples_legio_gene_burden.tsv         # Per-gene functional burden (cohort)
├── samples_legio_dr.tsv                  # 💊 Drug-resistance calls (pathotypr; only if --run_pathotypr)
│
├── multiqc/
│   └── samples_legio_multiqc_report.html   # 📊 Aggregate Report (QC, Mapping, Variants summary)
│
├── pipeline_info/                        # Nextflow execution timeline / report / trace + software versions
│
├── references/
│   └── LENS/                               # Processed Reference indices & SnpEff DB
│
└── MP00091/                                # 📁 Per-Sample Results Directory
    │
    ├── MP00091.LENS.filtered.bam           # 🧬 Dedup + length-aware-filtered BAM (analysis BAM, published by default)
    ├── MP00091.LENS.filtered.bam.bai       # BAM Index  (or .cram with --output_cram)
    │                                       # (final.bam is published instead when dynamic_read_filter=false, or in addition with --publish_prefilter_bam)
    │
    ├── MP00091.LENS.ann.vcf.gz             # 🎯 MAIN OUTPUT: Annotated Variants (SNPs/Indels)
    ├── MP00091.LENS.ann.vcf.gz.tbi         # Index for the main VCF
    ├── MP00091.LENS.canonical.ann.vcf.gz   # 🧬 Variants re-annotated vs the canonical ref (only if --annotate_canonical)
    │
    ├── MP00091.LENS.all.pos.vcf.gz         # 🦴 BACKBONE: VCF containing ALL positions (WT + Variants)
    │                                       # (Ideal for phylogenetic supermatrices)
    │
    ├── MP00091.LENS.consensus.fasta        # 📝 Consensus Sequence (Fasta generated from VCF)
    │
    ├── MP00091.LENS.freebayes.raw...       # 🧪 RAW VCF: Unfiltered calls (debug/comparison)
    ├── MP00091.LENS.var.homo.SNPs.ann...   # 📂 Split VCFs: Subset of Homozygous SNPs (Annotated)
    ├── MP00091.LENS.var.het.SNPs.ann...    # 📂 Split VCFs: Subset of Heterozygous SNPs (Annotated)
    ├── MP00091.LENS.var.homo.indel...      # 📂 SPLIT VCF: Homozygous Indels only
    │
    ├── lineage/                            # 🧬 Pathotypr typing (only if --run_pathotypr)
    │   └── MP00091__<runId>.pathotypr.lineage_summary.tsv   # -> Pathotypr lineage / sub-lineage call (__<runId> avoids multi-lane collisions)
    │
    └── stats/                              # 📉 Statistics & Logs Folder
        ├── MP00091.LENS.log                # -> LEGACY summary log (Tab-separated metrics)
        ├── MP00091.LENS.dedup.stats        # -> Samtools stats (reads mapped, coverage, etc.)
        ├── MP00091.LENS.filter_mqc.tsv     # -> Length-aware read-filter drop counts (input/kept/dropped; MultiQC table)
        ├── MP00091.LENS.mask_sites.tsv     # -> Specific positions masked due to low confidence
        ├── MP00091...fastp.html/.json      # -> Trimming quality reports
        ├── MP00091...kraken.report         # -> Taxonomic classification report (only with --kraken2_db)
        ├── MP00091.LENS.snpeff.csv         # -> Variant effect statistics
        ├── MP00091__<runId>.dr_mutations.tsv  # -> Pathotypr per-sample DR mutations (only if --run_pathotypr)
        └── Locus_to_exclude_LENS.txt       # -> List of repetitive regions excluded from calling
```

## QC flag codes

`<samplesheet>_qc_flags.tsv` gives each sample a **verdict** (PASS / WARN / FAIL) and the
`flags` it tripped. A sample **FAILs** if it hits any FAIL-level flag (or has no metrics);
otherwise any flag makes it WARN. All thresholds are `--report_*` params and are
live-adjustable inside the HTML report.

| Code | Condition | Param | Level |
| :--- | :--- | :--- | :---: |
| `LOW_DEPTH` | mean depth below | `--report_depth_min` | :octicons-x-circle-fill-16:{ .red } **FAIL** |
| `LOW_BREADTH` | genome breadth below | `--report_breadth_min` | :octicons-x-circle-fill-16:{ .red } **FAIL** |
| `HIGH_MISSING` | missing % above | `--report_missing_max` | :octicons-x-circle-fill-16:{ .red } **FAIL** |
| `NO_DATA` | no QC metrics for the sample | - | :octicons-x-circle-fill-16:{ .red } **FAIL** |
| `MAPPING_LOW` | mapped % below | `--report_mapping_min` | :octicons-alert-fill-16:{ .amber } WARN |
| `HIGH_DUP` | duplication % above | `--report_dup_max` | :octicons-alert-fill-16:{ .amber } WARN |
| `HIGH_IUPAC` | ambiguous/IUPAC % above | `--report_iupac_max` | :octicons-alert-fill-16:{ .amber } WARN |
| `TITV_LOW` | Ti/Tv below | `--report_titv_min` | :octicons-alert-fill-16:{ .amber } WARN |

**FAIL** = a candidate for the exclusion basket (FAIL samples start pre-selected);
**WARN** = review, usually keep. See [Downstream](downstream.md) for applying the basket.

## Directory layout

```text
BAMpiro/
├── main.nf                  # Entry point: reads the samplesheet, then calls each stage in order
├── subworkflows/            # One per pipeline stage, with an explicit take/emit
│   ├── prepare_references.nf  # 1. Reference indexing + SnpEff DB
│   ├── read_qc.nf             # 2-3. Read validation, Kraken2 screen, FastP
│   ├── lineage_typing.nf      # 4. Pathotypr lineage + drug resistance
│   ├── map_reads.nf           # 5. Mapping, merge/markdup, length-aware read filter
│   ├── call_variants.nf       # 6. FreeBayes, backbone, merge
│   ├── make_consensus.nf      # 7. Masked and virgin consensus
│   ├── annotate.nf            # 8. SnpEff (main / legacy)
│   ├── legacy_stats.nf        # 9. Per-sample metrics
│   ├── multiqc_report.nf      # 10. MultiQC
│   ├── cohort_report.nf       # 11. Cohort summary, SNP matrix, QC report
│   └── gene_conversion.nf     # 12. Gene conversion (opt-in)
├── modules/                 # The processes themselves (Nextflow DSL2)
│   ├── qc.nf                # FastP, Kraken, MultiQC, software versions
│   ├── mapping.nf           # BWA-MEM2, MarkDup, length-aware read filter
│   ├── variants.nf          # FreeBayes, Backbone, Merge
│   ├── annotation.nf        # SnpEff (main / legacy / canonical), legacy stats
│   ├── consensus.nf         # Consensus Fasta
│   ├── report.nf            # Cohort summary, SNP matrix, DR collection, QC report
│   ├── pathotypr.nf         # Lineage + drug-resistance typing
│   ├── reference.nf         # Reference Prep
│   ├── gene_conversion.nf   # Paralog map + conversion tract detection
│   └── utils.nf             # Publish-path routing, parameter validation, --help
├── bin/                     # Everything a process actually runs, kept out of the process scripts
│   ├── qc_report.py               # CLI for the interactive self-contained HTML QC report
│   ├── qcreport/                  # The package behind it
│   │   ├── parsers.py             #   Read every input file the report consumes
│   │   ├── metrics.py             #   Thresholds and the PASS/WARN/FAIL verdict
│   │   ├── panels.py              #   SNP dynamics, epistasis, SNP matrix
│   │   └── render.py              #   CSS/JS assets and the HTML assembly
│   ├── stats_to_legacy.py         # Per-sample metrics -> legacy log (parses the pathotypr lineage call)
│   ├── collect_summary.py         # Aggregate per-sample logs -> cohort summary + gene-burden TSVs
│   ├── collect_dr.py              # Aggregate pathotypr DR calls -> run drug-resistance TSV
│   ├── pathotypr_liftover.py      # Alignment-free k-mer coordinate liftover (mapping ref <-> H37Rv)
│   ├── build_snp_matrix.py        # Build the master SNP matrix (site × sample)
│   ├── WGS_fasta_allpos.py        # Consensus FASTA from the all-positions VCF
│   ├── build_min_unique_len.py    # Per-reference mappability track (genmap)
│   ├── filter_reads_mappability.py  # Length-aware read filter
│   ├── paralog_map.py             # Donor/acceptor map + diagnostic sites from the self-alignment
│   ├── gene_conversion.py         # Gene conversion tracts, with read-level breakpoint evidence
│   ├── vcf_filter_rules.py        # The bcftools expressions defining a hom / het call
│   ├── format_snps_for_backbone.awk # Genotype re-validation before the backbone merge
│   ├── backbone_allpos.awk        # mpileup -> one VCF record per reference position
│   ├── safe_tabix                 # Index a VCF, tolerating only a genuinely empty one
│   └── extract_kraken_reads.py    # Vendored from KrakenTools (decontamination)
├── assets/
│   ├── mycolorsTB_nature.tsv       # Canonical MTBC lineage colour palette (report)
│   ├── H37Rv_blindspots.bed        # H37Rv Illumina blind-spots (Zenodo 3701840; --mask_blindspots)
│   ├── H37Rv_blindspots.README.md  # Provenance / derivation of the blind-spots BED
│   └── NO_FILE*                    # Placeholders for the optional QC report inputs
├── tests/                   # See tests/README.md: unit, front-end and stub-run legs
├── conf/                    # Ready-made configs: sites (garnatxa.config -> -profile garnatxa),
│                            #   organisms (tuberculosis.config, organism.config), tests (test*.config)
├── .github/dockerfile/      # Container recipe (bundles pathotypr + Zenodo panels + H37Rv snpEff DB)
├── docs/                    # This documentation (MkDocs Material; built via mkdocs.yml)
└── nextflow.config          # Global configuration & params
```

!!! note "Why the science lives in `bin/`"

    A process script cannot be tested: `-stub-run` replaces it, and nothing can import it. Anything
    with a decision in it - the hom/het thresholds, the genotype re-validation, the pileup parser -
    is a file under `bin/` that the process calls, so the test suite can exercise it directly.
