# Configuration

Customize execution by passing parameters on the command line (e.g.
`--threads 16`) or by editing `nextflow.config`.

`nextflow run main.nf --help` prints every parameter with its default and exits, so the
table below is never the only source of truth. A parameter BAMpiro does not know is
**rejected** before any work starts:

```console
$ nextflow run main.nf --tsv samples.tsv --freebayes_ploydi 1
Unrecognised parameter(s): --freebayes_ploydi
  --freebayes_ploydi: did you mean --freebayes_ploidy?
Run with --help to list every parameter.
```

Previously a misspelled flag was accepted, ignored, and the run finished with the default.

## Execution profiles

`-profile` decides **where** and **how** tasks execute. Nothing site-specific is baked into
the defaults, so a fresh clone (or `nextflow run PathoGenOmics-Lab/BAMpiro`) runs correctly
anywhere. Pick one executor profile and add as many option profiles as you like,
comma-separated and with no spaces.

| Profile | Sets | Notes |
| :--- | :--- | :--- |
| `standard` | `executor = local` | The default, applied when you pass no `-profile`. Runs on the machine you launched from. |
| `local` | `executor = local` | An explicit alias for `standard`, for a command that says what it does. |
| `slurm` | `executor = slurm` | Any SLURM cluster. No queue, account or module system assumed; add yours with `-c`. |
| `garnatxa` | `conf/garnatxa.config` | The I2SysBio / University of Valencia cluster: QoS per step, `module load singularity`, bind paths, the shared Kraken2 DB, a mappability cache and a work directory on scratch. |
| `docker` | Docker on, Singularity off | Combine it with an executor profile: `-profile local,docker`. |
| `generic` | ploidy 1, MTBC features off | A non-tuberculosis organism (see [below](#example-configs)). |
| `test` · `test_full` | `conf/test.config` (+ `test_full`) | The bundled fixture cohort, sized for `nextflow run main.nf -profile test -stub-run`. |

```bash
-profile local,docker      # a laptop, with Docker
-profile slurm             # any SLURM cluster
-profile garnatxa          # the I2SysBio cluster, with its paths already set
-profile slurm,generic     # a SLURM cluster, non-tuberculosis organism
```

!!! warning "On a cluster, choose the executor profile deliberately"

    Because the default executor is `local`, launching from a cluster login node with no
    `-profile` runs the **whole pipeline on the login node** instead of submitting it. The
    startup banner prints the profile, executor, container and Kraken state it resolved:
    check the `Executor` line before walking away.

Adapting the pipeline to another cluster is covered in
[Configuring a run](tutorials/configuring-a-run.md#on-an-hpc-cluster-make-your-files-visible-to-the-container),
with `conf/garnatxa.config` as the worked template.

## Feature toggles (what to turn on / off)

Each of these is a boolean you flip on the command line (e.g. `--run_pathotypr true`).
This is the quick on/off overview; every flag also appears with its default in the full
table below.

**Opt-in - OFF by default, enable with `--<flag> true`:**

| Flag | Default | Turns on |
| :--- | :---: | :--- |
| `--run_pathotypr` | :octicons-x-circle-fill-16:{ .red } off | Alignment-free lineage + WHO drug-resistance typing ([pathotypr](pathotypr.md)) |
| `--annotate_canonical` | :octicons-x-circle-fill-16:{ .red } off | Dual **amino-acid** numbering (each SNP lifted to H37Rv, then annotated there) + Mycobrowser gene links |
| `--variant_liftover` | :octicons-x-circle-fill-16:{ .red } off | Dual **coordinate** - the H37Rv position per variant (alignment-free k-mer liftover) |
| `--mask_blindspots` | :octicons-x-circle-fill-16:{ .red } off | Mask the H37Rv Illumina blind-spots ([Zenodo 3701840](https://zenodo.org/records/3701840)) |
| `--output_cram` | :octicons-x-circle-fill-16:{ .red } off | Publish alignments as CRAM instead of BAM |
| `--report_gate` | :octicons-x-circle-fill-16:{ .red } off | Fail the run if any sample is flagged **FAIL** |
| `--publish_prefilter_bam` | :octicons-x-circle-fill-16:{ .red } off | Also publish the pre-filter `final.bam` |
| `--publish_virgin_allpos_vcf` | :octicons-x-circle-fill-16:{ .red } off | Also publish the virgin (`.raw`) per-position VCF |

**On by default - DISABLE with `--<flag> false`:**

| Flag | Default | Turns off |
| :--- | :---: | :--- |
| `--make_qc_report` | :octicons-check-circle-fill-16:{ .green } on | The interactive HTML QC report + `qc_flags.tsv` |
| `--make_snp_matrix` | :octicons-check-circle-fill-16:{ .green } on | The cohort master SNP-matrix TSV |
| `--make_indel_matrix` | :octicons-check-circle-fill-16:{ .green } on | The cohort master indel-matrix TSV |
| `--run_mnv` | :octicons-check-circle-fill-16:{ .green } on | Codon-level changes with get_MNV (per-sample `mnv/`, `<samplesheet>_mnv.tsv`, the SNP matrix's `codon_change`) |
| `--make_consensus` | :octicons-check-circle-fill-16:{ .green } on | The per-sample consensus FASTA |
| `--exclude_repeats` | :octicons-check-circle-fill-16:{ .green } on | `nucmer` self-alignment repeat masking |
| `--dynamic_read_filter` | :octicons-check-circle-fill-16:{ .green } on | The length-aware (mappability) read filter |
| `--annotate_main_vcf` · `--annotate_legacy_vcfs` | :octicons-check-circle-fill-16:{ .green } on | snpEff annotation of the main / split VCFs |
| `--blindspot_liftover` | :octicons-check-circle-fill-16:{ .green } on | Lifting the blind-spots onto a non-H37Rv reference (only relevant with `--mask_blindspots`) |
| `--nested_output` | :octicons-check-circle-fill-16:{ .green } on | Nested per-sample output folders (e.g. `MP/00/09/1/`) |
| `--publish_allpos_vcf` | :octicons-check-circle-fill-16:{ .green } on | Publishing the all-positions (backbone) VCF |
| `--kraken_memory_mapping` | :octicons-check-circle-fill-16:{ .green } on | Reading the Kraken2 DB via mmap (RAM-sharing) |

## Example configs

Two documented, ready-to-use **organism** configs live in
[`conf/`](https://github.com/PathoGenOmics-Lab/BAMpiro/tree/indel-mask/conf) - apply either with `-c`.
(The same directory also holds `garnatxa.config`, the **site** config behind
`-profile garnatxa`, and the `test*` configs behind `-profile test`.)

=== "M. tuberculosis (default)"

    BAMpiro's built-in defaults are already TB-tuned, so a plain run is a valid TB run.
    [`conf/tuberculosis.config`](https://github.com/PathoGenOmics-Lab/BAMpiro/blob/indel-mask/conf/tuberculosis.config)
    additionally turns on the **full recommended MTBC feature set** - lineage +
    drug-resistance typing, dual H37Rv amino-acid numbering, the H37Rv coordinate
    liftover, and blind-spot masking (all bundled in the container):

    ```bash
    nextflow run main.nf --tsv samples.tsv --outdir results -profile local,docker -c conf/tuberculosis.config
    ```

    *(Running with no `-c` is still a valid TB run - it just leaves those optional features off.)*

=== "Another organism"

    BAMpiro is organism-agnostic, but its defaults are TB-tuned (diploid calling for
    mixed infections + the MTBC-only features). Add the built-in **`generic`** profile,
    which sets **ploidy 1** and forces the MTBC-only features off, and copy the commented
    template [`conf/organism.config`](https://github.com/PathoGenOmics-Lab/BAMpiro/blob/indel-mask/conf/organism.config)
    for the rest (Kraken DB, a canonical numbering, gate thresholds):

    ```bash
    nextflow run main.nf --tsv samples.tsv --outdir results -profile local,docker,generic -c my_organism.config
    ```

    What `generic` changes versus the TB defaults:

    | Setting | TB default | `generic` | Why |
    | :--- | :--- | :--- | :--- |
    | `freebayes_ploidy` | `2` | **`1`** | TB calls diploid to detect mixed infections; a clonal bacterium is haploid |
    | `run_pathotypr` · `variant_liftover` · `mask_blindspots` · `annotate_canonical` | *(off)* | **off** | all *M. tuberculosis*-specific (H37Rv lineage/DR, blind-spots, coordinate liftover) |

    Everything else is already generic: point the samplesheet `refFasta` / `refGff` at
    your genome, set `--kraken2_db` and the `taxId` column for your species, and the
    report's TB-only panels (lineage, drug resistance) simply **self-hide** when there is
    no such data. To keep a canonical amino-acid numbering for your organism, set
    `--annotate_canonical true` with `--canonical_snpeff_db` / `--canonical_label` - but
    note **only the H37Rv SnpEff DB is bundled** in the image, so for another organism
    either add its genome to the image and rebuild, or make sure a bare
    `canonical_snpeff_db` id can be downloaded by SnpEff at run time (needs internet).

## All parameters

| Category | Parameter | Default | Description |
| :--- | :--- | :--- | :--- |
| **Input/Output** | `--tsv` | *(required, no default)* | Path to the input sample sheet (TSV). Omitting it stops the run at once with a message naming the expected columns. |
| | `--outdir` | `results_bampiro` | Directory where results will be saved. |
| | `--threads` | `8` | CPUs for the **BWA-MEM2 mapping** steps only; every other process has its own fixed `cpus` (see [Resource requirements](installation.md#resource-requirements)). |
| | `--container` | `docker://paururo/bampiro@sha256:bef4375…` | Container image, **pinned by digest** (the 1.1.0-rc2 image) so two runs months apart cannot pick up a rebuilt tag. Override it to use your own copy; the resolved value is recorded in the report's provenance footer. See [the container](installation.md#the-container). |
| | `--nested_output` | `true` | Nest per-sample folders (e.g. `MP00091` → `MP/00/09/1`). |
| | `--publish_mode` | `copy` | `copy` duplicates outputs into `outdir`; `link` hardlinks them to the work dir. |
| | `--output_cram` | `false` | Publish the alignment as CRAM (~40-50% smaller) instead of BAM. |
| | `--publish_prefilter_bam` | `false` | Also publish the pre-filter `final.bam`. With the length-aware filter on (`--dynamic_read_filter`, the default), `filtered.bam` is the analysis BAM and `final.bam` is a ~redundant second full BAM, so it is dropped by default. |
| | `--publish_allpos_vcf` | `true` | Publish the masked per-position (all-positions) VCF the consensus is built from. Set `false` to drop it and save disk. |
| | `--publish_virgin_allpos_vcf` | `false` | Also publish the virgin (`.raw`) per-position VCF (needs `--publish_allpos_vcf`); redundant with the virgin consensus, so off by default. |
| **Resource ceilings** | `--max_cpus` | *(none)* | Cap **every** process request to what the machine can actually give, e.g. `--max_cpus 4 --max_memory 8.GB`. Without a cap the requests stand as written, which is what you want on a cluster but not on a laptop (Kraken2 alone asks for 80 GB). |
| | `--max_memory` | *(none)* | As above, for memory. Nextflow units: `8.GB`, `16.GB`. |
| | `--max_time` | *(none)* | As above, for wall time. Nextflow units: `30.m`, `2.h`. |
| **Lineage & DR (Pathotypr)** | `--run_pathotypr` | `false` | Enable alignment-free MTBC lineage + WHO drug-resistance typing. |
| | `--pathotypr_bin` | `pathotypr` | Executable name (on `PATH` inside the container). |
| | `--pathotypr_ref` | `/opt/pathotypr/reference.fasta` | MTBC-ancestor FASTA the markers are defined on (bundled). |
| | `--pathotypr_markers` | `/opt/pathotypr/lineage_markers.tsv` | Zenodo lineage markers (bundled). |
| | `--pathotypr_dr_markers` | `/opt/pathotypr/dr_markers.tsv` | Zenodo WHO drug-resistance markers (bundled). |
| | `--pathotypr_min_alt` | `95` | Min alt-allele % for a DR call. Lower it (10-25) to catch heteroresistant / minority alleles. |
| **Dual AA numbering** | `--annotate_canonical` | `false` | Annotate each SNP against a canonical snpEff DB at its lifted canonical position, for dual (H37Rv) amino-acid numbering and gene names on any reference. |
| | `--canonical_snpeff_db` | `Mycobacterium_tuberculosis_h37rv` | Canonical snpEff genome for the second annotation. |
| | `--canonical_chrom` | `Chromosome` | The name that database gives its chromosome, written for a single-contig `--canonical_ref`; empty keeps the FASTA's names. |
| | `--canonical_label` | `H37Rv` | How that numbering is labelled in the report. |
| **Canonical coords & masking** | `--variant_liftover` | `false` | Fill the H37Rv **coordinate** per variant via the alignment-free k-mer liftover (report SNP tables); every reference of the run from its own sequence, every contig of it. |
| | `--canonical_ref` | `/opt/pathotypr/reference.fasta` | FASTA the canonical (H37Rv-colinear) coordinates are defined on - the liftover / blind-spot **target**. |
| | `--mask_blindspots` | `false` | Add the H37Rv Illumina blind-spots ([Zenodo 3701840](https://zenodo.org/records/3701840)) to the exclusion mask. |
| | `--blindspot_bed` | `assets/H37Rv_blindspots.bed` | The blind-spots BED (H37Rv / NC_000962.3 coordinates). |
| | `--blindspot_liftover` | `true` | Lift the blind-spots onto the run reference, every contig of it (k-mer liftover), so the mask is correct on **any** reference. |
| **QC Report** | `--make_qc_report` | `true` | Build the interactive HTML QC report + `qc_flags.tsv`. |
| | `--make_snp_matrix` | `true` | Also emit the master SNP-matrix TSV. |
| | `--depth_window` | `1000` | Window (bp) of each sample's depth profile, behind the deletions and SNPs per callable kb. |
| | `--deletion_min_len` | `200` | Shortest stretch without reads reported as a deletion. |
| | `--make_snp_distances` | `true` | Pairwise SNP distances between the consensus sequences, and the report's relatedness page. |
| | `--snp_cluster_threshold` | `12` | SNPs within which the report starts drawing clusters (movable there), and beyond which a sample is far from its group. 12 is the usual *M. tuberculosis* transmission cut; set your organism's own. |
| | `--report_gate` | `false` | Fail the run if any sample is flagged **FAIL**. |
| | `--report_depth_min` | `10` | Gate: min mean depth (all `report_*` cut-offs are editable live in the report). |
| | `--report_breadth_min` | `90` | Gate: min breadth %. |
| | `--report_missing_max` | `10` | Gate: max missing %. |
| | `--report_dup_max` | `40` | Gate: max duplication %. |
| | `--report_iupac_max` | `2` | Gate: max IUPAC (ambiguous) %. |
| | `--report_mapping_min` | `80` | Gate: min mapped %. |
| **Cluster** | `--extra_binds` | *(empty)* | Extra host paths a container must be able to see, comma separated and no spaces. Only the site profiles bind anything, so only they use it. An unbound path does not fail at submission: it fails inside the task, as a file that is plainly there on the login node and absent to the tool that opens it. |
| **QC & Filter** | `--kraken2_db` | *(none)* | Path to a Kraken2 database directory. **No default**: without one the contamination screen is skipped and the report's taxonomy panel hides itself. `-profile garnatxa` sets the shared database on that cluster. |
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
| **Indels & codons** | `--make_indel_matrix` | `true` | Emit `<samplesheet>_indel_matrix.tsv`: every indel a sample calls with a PASS (the rules above), each sample's AF, depth and filter. |
| | `--run_mnv` | `true` | Read the SNPs that share a codon whole, on the reads that carry them ([get_MNV](https://github.com/PathoGenOmics-Lab/get_MNV)): the amino-acid change a per-base annotation misses. |
| | `--mnv_min_reads` | `2` | Reads that must carry every change of a codon for it to count as one change in `<samplesheet>_mnv.tsv`. |
| | `--mnv_gff_features` | `CDS` | GFF feature types get_MNV reads codons from. |
| | `--mnv_translation_table` | `11` | NCBI genetic code get_MNV translates with (11: bacteria and archaea). |
| | `--mnv_container` | `quay.io/biocontainers/get_mnv` 1.1.5, by digest | The image get_MNV runs in until the pipeline image carries it; empty runs it in the pipeline image. |
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
