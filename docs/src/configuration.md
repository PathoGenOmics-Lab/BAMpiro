# Configuration

Customize execution by passing parameters on the command line (e.g.
`--threads 16`) or by editing `nextflow.config`.

| Category | Parameter | Default | Description |
| :--- | :--- | :--- | :--- |
| **Input/Output** | `--tsv` | `samples_legio.tsv` | Path to the input sample sheet (TSV). |
| | `--outdir` | `results_bampiro` | Directory where results will be saved. |
| | `--threads` | `8` | Max CPUs per process (where applicable). |
| | `--container` | *(pinned digest)* | Container image. Defaults to a pinned `paururo/bampiro` digest for reproducibility. |
| | `--nested_output` | `true` | Nest per-sample folders (e.g. `MP001` → `MP/00/1`). |
| | `--publish_mode` | `copy` | `copy` duplicates outputs into `outdir`; `link` hardlinks them to the work dir. |
| | `--output_cram` | `false` | Publish the alignment as CRAM (~40-50% smaller) instead of BAM. |
| **Lineage & DR (Pathotypr)** | `--run_pathotypr` | `false` | Enable alignment-free MTBC lineage + WHO drug-resistance typing. |
| | `--pathotypr_bin` | `pathotypr` | Executable name (on `PATH` inside the container). |
| | `--pathotypr_ref` | `/opt/pathotypr/reference.fasta` | MTBC-ancestor FASTA the markers are defined on (bundled). |
| | `--pathotypr_markers` | `/opt/pathotypr/lineage_markers.tsv` | Zenodo lineage markers (bundled). |
| | `--pathotypr_dr_markers` | `/opt/pathotypr/dr_markers.tsv` | Zenodo WHO drug-resistance markers (bundled). |
| | `--pathotypr_rf_model` | `/opt/pathotypr/rf_model.pathotypr` | Zenodo pre-trained RF lineage model (bundled). |
| | `--pathotypr_min_alt` | `95` | Min alt-allele % for a DR call. Lower it (10-25) to catch heteroresistant / minority alleles. |
| **Dual AA numbering** | `--annotate_canonical` | `false` | Re-annotate variants against a canonical snpEff DB for dual (H37Rv) numbering. |
| | `--canonical_snpeff_db` | `Mycobacterium_tuberculosis_h37rv` | Canonical snpEff genome for the second annotation. |
| | `--canonical_label` | `H37Rv` | How that numbering is labelled in the report. |
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
| | `--fastp_min_length` | `35` | Discard reads shorter than this length. |
| **Mapping/Backbone**| `--allpos_min_cov` | `30` | Minimum coverage to call a site "WT" (otherwise "NC"). |
| | `--allpos_max_depth` | `10000` | Max depth for mpileup to avoid memory issues. |
| | `--allpos_min_bq` | `20` | Minimum base quality for backbone calling. |
| **Variant Calling** | `--freebayes_ploidy` | `2` | Ploidy (1 for haploid, 2 for mixed/diploid). |
| | `--freebayes_min_map_qual`| `30` | Min mapping quality to use a read. |
| | `--freebayes_min_base_qual`| `20` | Min base quality to use a base. |
| | `--freebayes_min_alt_fraction`| `0.05` | Min fraction of alt reads to propose a variant. |
| | `--freebayes_min_alt_count`| `2` | Min alt-supporting reads to propose a variant. |
| **VAF & Filters** | `--hom_threshold` | `0.90` | Frequency ≥ 0.90 is called **Homozygous**. |
| | `--het_min_frac` | `0.10` | Frequency between 0.10 and 0.90 is **Heterozygous**. |
| | `--filter_min_dp` | `30` | Minimum depth required to call a variant. |
| | `--min_alt_fwd/rev` | `2` | Min variant supporting reads in FWD and REV strands. |
| **Consensus** | `--make_consensus` | `true` | Generate a consensus FASTA for each sample. |
| | `--consensus_min_dp` | `7` | Depth threshold below which a base becomes "No Call". |
| | `--consensus_mask_char` | `X` | Character for masked/low-quality sites. |
| | `--consensus_nocall_char`| `-` | Character for no-coverage sites (gaps). |
| **Flags** | `--exclude_repeats` | `true` | Mask self-aligned repetitive regions from reference. |
| | `--annotate_main_vcf` | `true` | Run SnpEff on the main VCF. |
| | `--annotate_legacy_vcfs`| `true` | Run SnpEff on split VCFs (homo/het/indel). |

> **More knobs.** `nextflow.config` also exposes the length-aware read filter
> (`genmap_*`, `dynamic_read_filter`, `mappability_*`), the reference-bias consensus
> gate (`consensus_ref_*`), the virgin/unmasked consensus (`keep_virgin_consensus`,
> `virgin_full_freebayes`), the region-parallel FreeBayes (`freebayes_parallel*`),
> and the backbone/masking options (`allpos_*`, `mask_baq_dropouts`). See the
> commented `params { … }` block for the full, annotated list.
