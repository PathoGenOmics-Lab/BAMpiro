# Configuration

Customize execution by passing parameters on the command line (e.g.
`--threads 16`) or by editing `nextflow.config`.

## Feature toggles (what to turn on / off)

Each of these is a boolean you flip on the command line (e.g. `--run_pathotypr true`).
This is the quick on/off overview; every flag also appears with its default in the full
table below.

**Opt-in — OFF by default, enable with `--<flag> true`:**

| Flag | Turns on |
| :--- | :--- |
| `--run_pathotypr` | Alignment-free lineage + WHO drug-resistance typing ([pathotypr](pathotypr.md)) |
| `--annotate_canonical` | Dual **amino-acid** numbering (re-annotate vs H37Rv) + Mycobrowser gene links |
| `--variant_liftover` | Dual **coordinate** — the H37Rv position per variant (alignment-free k-mer liftover) |
| `--mask_blindspots` | Mask the H37Rv Illumina blind-spots ([Zenodo 3701840](https://zenodo.org/records/3701840)) |
| `--output_cram` | Publish alignments as CRAM instead of BAM |
| `--report_gate` | Fail the run if any sample is flagged **FAIL** |
| `--publish_prefilter_bam` | Also publish the pre-filter `final.bam` |
| `--publish_virgin_allpos_vcf` | Also publish the virgin (`.raw`) per-position VCF |

**On by default — DISABLE with `--<flag> false`:**

| Flag | Turns off |
| :--- | :--- |
| `--make_qc_report` | The interactive HTML QC report + `qc_flags.tsv` |
| `--make_snp_matrix` | The cohort master SNP-matrix TSV |
| `--make_consensus` | The per-sample consensus FASTA |
| `--exclude_repeats` | `nucmer` self-alignment repeat masking |
| `--dynamic_read_filter` | The length-aware (mappability) read filter |
| `--annotate_main_vcf` · `--annotate_legacy_vcfs` | snpEff annotation of the main / split VCFs |
| `--blindspot_liftover` | Lifting the blind-spots onto a non-H37Rv reference (only relevant with `--mask_blindspots`) |
| `--nested_output` | Nested per-sample output folders (e.g. `MP/00/09/1/`) |
| `--publish_allpos_vcf` | Publishing the all-positions (backbone) VCF |
| `--kraken_memory_mapping` | Reading the Kraken2 DB via mmap (RAM-sharing) |

## All parameters

| Category | Parameter | Default | Description |
| :--- | :--- | :--- | :--- |
| **Input/Output** | `--tsv` | `samples_legio.tsv` | Path to the input sample sheet (TSV). |
| | `--outdir` | `results_bampiro` | Directory where results will be saved. |
| | `--threads` | `8` | Max CPUs per process (where applicable). |
| | `--container` | `docker://paururo/bampiro:1.0.1` | Container image. The default is the mutable `1.0.1` Docker Hub **tag** (not a digest); re-pin to `docker://paururo/bampiro@sha256:<digest>` for byte-exact reproducibility, or override. |
| | `--nested_output` | `true` | Nest per-sample folders (e.g. `MP00091` → `MP/00/09/1`). |
| | `--publish_mode` | `copy` | `copy` duplicates outputs into `outdir`; `link` hardlinks them to the work dir. |
| | `--output_cram` | `false` | Publish the alignment as CRAM (~40-50% smaller) instead of BAM. |
| | `--publish_prefilter_bam` | `false` | Also publish the pre-filter `final.bam`. With the length-aware filter on (`--dynamic_read_filter`, the default), `filtered.bam` is the analysis BAM and `final.bam` is a ~redundant second full BAM, so it is dropped by default. |
| | `--publish_allpos_vcf` | `true` | Publish the masked per-position (all-positions) VCF the consensus is built from. Set `false` to drop it and save disk. |
| | `--publish_virgin_allpos_vcf` | `false` | Also publish the virgin (`.raw`) per-position VCF (needs `--publish_allpos_vcf`); redundant with the virgin consensus, so off by default. |
| **Lineage & DR (Pathotypr)** | `--run_pathotypr` | `false` | Enable alignment-free MTBC lineage + WHO drug-resistance typing. |
| | `--pathotypr_bin` | `pathotypr` | Executable name (on `PATH` inside the container). |
| | `--pathotypr_ref` | `/opt/pathotypr/reference.fasta` | MTBC-ancestor FASTA the markers are defined on (bundled). |
| | `--pathotypr_markers` | `/opt/pathotypr/lineage_markers.tsv` | Zenodo lineage markers (bundled). |
| | `--pathotypr_dr_markers` | `/opt/pathotypr/dr_markers.tsv` | Zenodo WHO drug-resistance markers (bundled). |
| | `--pathotypr_min_alt` | `95` | Min alt-allele % for a DR call. Lower it (10-25) to catch heteroresistant / minority alleles. |
| **Dual AA numbering** | `--annotate_canonical` | `false` | Re-annotate variants against a canonical snpEff DB for dual (H37Rv) amino-acid numbering. |
| | `--canonical_snpeff_db` | `Mycobacterium_tuberculosis_h37rv` | Canonical snpEff genome for the second annotation. |
| | `--canonical_label` | `H37Rv` | How that numbering is labelled in the report. |
| **Canonical coords & masking** | `--variant_liftover` | `false` | Fill the H37Rv **coordinate** per variant via the alignment-free k-mer liftover (report SNP tables); works for any reference. |
| | `--canonical_ref` | `/opt/pathotypr/reference.fasta` | FASTA the canonical (H37Rv-colinear) coordinates are defined on — the liftover / blind-spot **target**. |
| | `--mask_blindspots` | `false` | Add the H37Rv Illumina blind-spots ([Zenodo 3701840](https://zenodo.org/records/3701840)) to the exclusion mask. |
| | `--blindspot_bed` | `assets/H37Rv_blindspots.bed` | The blind-spots BED (H37Rv / NC_000962.3 coordinates). |
| | `--blindspot_liftover` | `true` | Lift the blind-spots onto the run reference (k-mer liftover) so the mask is correct on **any** reference. |
| **QC Report** | `--make_qc_report` | `true` | Build the interactive HTML QC report + `qc_flags.tsv`. |
| | `--make_snp_matrix` | `true` | Also emit the master SNP-matrix TSV. |
| | `--report_gate` | `false` | Fail the run if any sample is flagged **FAIL**. |
| | `--report_depth_min` | `10` | Gate: min mean depth (all `report_*` cut-offs are editable live in the report). |
| | `--report_breadth_min` | `90` | Gate: min breadth %. |
| | `--report_missing_max` | `10` | Gate: max missing %. |
| | `--report_dup_max` | `40` | Gate: max duplication %. |
| | `--report_iupac_max` | `2` | Gate: max IUPAC (ambiguous) %. |
| | `--report_mapping_min` | `80` | Gate: min mapped %. |
| **QC & Filter** | `--kraken2_db` | *(path)* | Path to the Kraken2 database directory. |
| | `--kraken_memory_mapping` | `true` | Read the Kraken2 DB via mmap so parallel tasks share it in RAM (large peak-memory saving; slightly slower per task on a cold page cache). |
| | `--fastp_min_length` | `35` | Discard reads shorter than this length. |
| **Mapping/Backbone**| `--allpos_min_cov` | `30` | Minimum coverage to call a site "WT" (otherwise "NC"). |
| | `--allpos_max_depth` | `10000` | Max depth for mpileup to avoid memory issues. |
| | `--allpos_min_bq` | `20` | Minimum base quality for backbone calling. |
| **Variant Calling** | `--freebayes_ploidy` | `2` | Ploidy (1 for haploid, 2 for mixed/diploid). |
| | `--freebayes_min_map_qual`| `30` | Min mapping quality to use a read. |
| | `--freebayes_min_base_qual`| `20` | Min base quality to use a base. |
| | `--freebayes_min_alt_fraction`| `0.05` | Min fraction of alt reads to propose a variant. |
| | `--freebayes_min_alt_count`| `2` | Min alt-supporting reads to propose a variant. |
| | `--bcftools_sort_mem` | `2G` | Max memory (`-m`) for `bcftools sort` when building the all-positions (`all.pos`) VCF. |
| **VAF & Filters** | `--hom_threshold` | `0.90` | Frequency ≥ 0.90 is called **Homozygous**. |
| | `--het_min_frac` | `0.10` | Frequency between 0.10 and 0.90 is **Heterozygous**. |
| | `--filter_min_dp` | `30` | Minimum depth required to call a variant. |
| | `--min_alt_fwd/rev` | `2` | Min variant supporting reads in FWD and REV strands. |
| **Consensus** | `--make_consensus` | `true` | Generate a consensus FASTA for each sample. |
| | `--consensus_min_dp` | `7` | Depth threshold below which a base becomes "No Call". |
| | `--consensus_mask_char` | `X` | Character for masked/low-quality sites. |
| | `--consensus_nocall_char`| `-` | Character for no-coverage sites (gaps). |
| | `--consensus_wrap` | `80` | FASTA line width (bp per line) for the output consensus (`0` = single line). |
| **Flags** | `--exclude_repeats` | `true` | Mask self-aligned repetitive regions from reference. |
| | `--annotate_main_vcf` | `true` | Run SnpEff on the main VCF. |
| | `--annotate_legacy_vcfs`| `true` | Run SnpEff on split VCFs (homo/het/indel). |

> **More knobs.** `nextflow.config` also exposes the length-aware read filter
> (`genmap_*`, `dynamic_read_filter`, `filter_keep_unmapped`, `filter_strict_contigs`,
> `filter_max_drop_pct`, `mappability_*`), the reference-bias consensus gate
> (`consensus_ref_*`, `consensus_max_ref_*`), the virgin/unmasked consensus (`keep_virgin_consensus`,
> `virgin_full_freebayes`), the region-parallel FreeBayes (`freebayes_parallel*`),
> and the backbone/masking options (`allpos_*`, `mask_baq_dropouts`). See the
> commented `params { … }` block for the full, annotated list.
