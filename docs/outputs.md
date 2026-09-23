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
├── samples_legio_deletions.tsv           # Stretches some samples have no reads for (see below)
├── samples_legio_snp_distances.tsv       # SNPs between every two consensus sequences (square matrix)
├── samples_legio_snp_distances_pairs.tsv # The same, one row per pair, with how much each rests on
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
        ├── MP00091.LENS.gene_conversion.tsv       # -> Conversion tracts (only if --find_gene_conversion)
        ├── MP00091.LENS.gene_conversion_loci.tsv  # -> One row per paralog pair looked at, found or not
        ├── MP00091.LENS.depth_windows.tsv  # -> Depth per 1 kb window: mean, share without reads, share callable
        ├── MP00091.LENS.zero_depth.tsv     # -> Every stretch of 50+ bp no read covers
        ├── MP00091.LENS.gene_depth.tsv     # -> Per gene: share read, share callable, depth against the genome
        └── Locus_to_exclude_LENS.txt       # -> List of repetitive regions excluded from calling
```

## SNP matrix

`<samplesheet>_snp_matrix.tsv` has one row per SNP site called in any sample and two columns per
sample, `<sample>|AF` and `<sample>|DP`. A cell says what that sample's reads showed at the site:

| `AF` | `DP` | Meaning |
| :--- | :--- | :--- |
| a fraction | the depth | The alternate allele was called, carried by that fraction of the reads |
| `0` | the depth | The sample was read there and no alternate allele was called: absent, not unknown |
| `NA` | the depth | A call the files keep no read counts for, so its fraction is not known (see below) |
| empty | `0` | No read shows a base at the site in that sample, so its absence says nothing |
| empty | empty | The site is not in that sample's genome: it was mapped against another reference |

The depth of a `0` cell comes from the sample's all-positions VCF, the pileup the consensus is
built from, and counts the reads that show a base there: a read carrying a deletion over the
site says nothing about the allele. That pileup does not apply the caller's mapping-quality floor
(`--freebayes_min_map_qual`), so in a repeat it can read higher than the depth of a called cell.
Read a `0` together with its depth: at three reads it says little, and at or below
`--consensus_min_dp` the consensus leaves the site uncalled.

An MNP is split into the single-base SNPs it is made of, each with the MNP's allele fraction. A
call the variant VCF writes as part of a complex record is taken from the all-positions VCF,
which holds it as a SNP, with its own read counts; an all-positions VCF written before it kept
read counts gives a homozygous call `1.0000` and a heterozygous one `NA`, rather than the
genotype's dosage passed off as a fraction.

Rows are keyed on contig and position, so two references that share a contig name cannot be told
apart; the matrix step warns when it sees one.

## SNP distances

`<samplesheet>_snp_distances.tsv` is the square matrix of SNPs between every two consensus
sequences, `NA` between samples mapped against different references, which share no coordinates.
Two samples differ at a position only where both called a base (`A`, `C`, `G` or `T`) and the
bases differ: a no-call, a masked position or an ambiguity code at a mixed site is never a
difference, whichever sample it is in. These are the distances a tree built on the same
consensus sequences sees.

`<samplesheet>_snp_distances_pairs.tsv` has one row per pair with `compared`, the number of the
reference's variable positions both samples called. A thinly called sample is close to everyone
because it called little that could differ; the report leaves pairs that compared under half of
the variable positions out of its clusters. Off with `--make_snp_distances false`, and absent
without consensus sequences.

## Deletions

`<samplesheet>_deletions.tsv` lists the stretches of the reference that some samples have no reads
for. A stretch without reads means different things depending on the rest of the cohort, so each
one is compared with the other samples mapped to the same reference:

| `class` | Meaning |
| :--- | :--- |
| `private` | One sample lacks it and the others read it: a deletion in that sample |
| `shared` | Several samples lack it and the others read it: often a deletion that defines a lineage |
| `cohort` | At least 90% of the samples on its reference lack it: a repeat no read can be placed on, or a part of the reference none of these genomes has. Nobody's deletion |
| `alone` | The only sample on its reference read deeply enough, so there is nothing to compare it with |

Two samples' stretches are the same deletion when each covers at least half of the other,
measured against the stretch that opened the region. Plain overlap would chain a run of
neighbouring small deletions into one, and let one sample's long deletion swallow a short gap
every sample shares.

Only stretches of at least `--deletion_min_len` bp (200) count, and only in samples whose median
depth reaches `--report_depth_min`: in a thinly read sample, stretches without reads turn up by
chance. The `samples` column gives each carrier's own stretch and `genes` the genes a carrier
reads less than half of.

The per-sample tables behind it are in each sample's `stats/` folder and are useful on their own
(they are written whether or not a report is made). *Callable* there means read by
`--allpos_min_cov` reads or more, where the consensus can call a base; it is also what the
report's *SNPs / kb* divides by.
`gene_depth.tsv` says, for every gene of the GFF, how much of it the sample read and at what depth
against its genome-wide median, which is also where an amplification shows as a depth well above 1.
A deletion shorter than the minimum, or one the variant caller already reported as an indel, is
not in this file.

## Gene conversion

Off unless `--find_gene_conversion true`. Four files, and which one to open depends on the
question. Full detail in [gene conversion](gene-conversion.md).

| File | Where | What it is |
| :--- | :--- | :--- |
| `<samplesheet>_gene_conversion.tsv` | `outdir/` | **Start here.** Every sample's tracts with the cohort's reading on top: how many samples carry each event, which relative it came from, and whether recurrence says the reference rather than the isolates |
| `<sample>.<ref>.gene_conversion.tsv` | per sample, `stats/` | That sample's tracts on their own, before the cohort was consulted |
| `<sample>.<ref>.gene_conversion_loci.tsv` | per sample, `stats/` | One row per paralog pair analysed, reported or not. Written for the cohort pass, and the place to look when you expected a call somewhere and got none |
| `<refId>.paralog_pairs.tsv`, `<refId>.paralog_sites.tsv` | `references/<refId>/` | The donor/acceptor graph of your reference and every position the copies differ at. Computed once per reference and worth a look on its own |

Every one of them opens with `#` lines naming the tool and each setting that moved a number in
it, so a table can be checked against another run a year later. The cohort file carries the
per-sample lines as well, which is what shows a cohort assembled from samples run differently.

!!! note "These loci are excluded from variant calling by design"

    A tract sits inside a repeat, so it will not appear in the SNP matrix or the consensus. That
    is expected, and it is why the stage reads the deduplicated pre-filter alignment rather than
    the one everything else uses.

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
| `LINEAGE_MISMATCH` | typed as a different lineage from the other samples mapped to its reference (it was probably mapped to the wrong genome) | - | :octicons-alert-fill-16:{ .amber } WARN |
| `GROUP_MISMATCH` | further than the cluster threshold from every other sample of its samplesheet group (patient, line, series), in a group whose members are otherwise within it of each other: swapped, mislabelled, contaminated or reinfected | `--snp_cluster_threshold` | :octicons-alert-fill-16:{ .amber } WARN |

**FAIL** = a candidate for the report's exclusion list (FAIL samples start in it);
**WARN** = review, usually keep. See [Downstream](downstream.md) for applying the list.

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
