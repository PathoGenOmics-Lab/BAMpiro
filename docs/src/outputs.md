# Outputs

The pipeline organizes results by `sampleId`. Below is a breakdown using a sample
named `MP00091` mapped against reference `LENS`. Run/cohort-level deliverables land
at the top of `outdir` and are prefixed with the sample-sheet name (here
`samples_legio`); per-sample files live under their `sampleId` folder.

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
    ├── MP00091.LENS.final.bam              # 🧬 Merged, Coordinate-sorted, Deduplicated BAM
    ├── MP00091.LENS.final.bam.bai          # BAM Index  (or .cram with --output_cram)
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
    │   └── MP00091.pathotypr.lineage_summary.tsv   # -> Pathotypr lineage / sub-lineage call
    │
    └── stats/                              # 📉 Statistics & Logs Folder
        ├── MP00091.log                     # -> LEGACY summary log (Tab-separated metrics)
        ├── MP00091.LENS.dedup.stats        # -> Samtools stats (reads mapped, coverage, etc.)
        ├── MP00091.LENS.mask_sites.tsv     # -> Specific positions masked due to low confidence
        ├── MP00091...fastp.html/.json      # -> Trimming quality reports
        ├── MP00091...kraken.report         # -> Taxonomic classification report
        ├── MP00091.LENS.snpeff.csv         # -> Variant effect statistics
        ├── MP00091.dr_mutations.tsv        # -> Pathotypr per-sample DR mutations (only if --run_pathotypr)
        └── Locus_to_exclude_LENS.txt       # -> List of repetitive regions excluded from calling
```

## Directory layout

```text
BAMpiro/
├── bin/                     # Python helpers
│   ├── stats_to_legacy.py         # Per-sample metrics -> legacy log (parses the pathotypr lineage call)
│   ├── collect_summary.py         # Aggregate per-sample logs -> cohort summary + gene-burden TSVs
│   ├── collect_dr.py              # Aggregate pathotypr DR calls -> run drug-resistance TSV
│   ├── build_snp_matrix.py        # Build the master SNP matrix (site × sample)
│   ├── qc_report.py               # Build the interactive self-contained HTML QC report
│   ├── WGS_fasta_allpos.py        # Consensus FASTA from the all-positions VCF
│   ├── build_min_unique_len.py    # Per-reference mappability track (genmap)
│   └── filter_reads_mappability.py  # Length-aware read filter
├── assets/
│   └── mycolorsTB_nature.tsv     # Canonical MTBC lineage colour palette (report)
├── modules/                 # Nextflow DSL2 Modules
│   ├── qc.nf                # FastP, Kraken, MultiQC, software versions
│   ├── mapping.nf           # BWA-MEM2, MarkDup, length-aware read filter
│   ├── variants.nf          # FreeBayes, Backbone, Merge
│   ├── annotation.nf        # SnpEff (main / legacy / canonical), legacy stats
│   ├── consensus.nf         # Consensus Fasta
│   ├── report.nf            # Cohort summary, SNP matrix, DR collection, QC report
│   ├── pathotypr.nf         # Lineage + drug-resistance typing
│   ├── reference.nf         # Reference Prep
│   └── utils.nf             # Publish-path routing / clean publish dir
├── .github/dockerfile/      # Container recipe (bundles pathotypr + Zenodo panels + H37Rv snpEff DB)
├── docs/                    # This documentation (mdBook source in docs/src/)
├── nextflow.config          # Global configuration & params
└── main.nf                  # Main workflow entry point
```
