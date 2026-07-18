# Quick Start

By default BAMpiro assumes *M. tuberculosis* settings (ploidy = 2 to detect mixed
infections). Lineage/DR typing is **off** by default - enable it with
`--run_pathotypr true`. **Working with a different organism?** Add `-profile
standard,generic` (ploidy 1 + the MTBC-only features off) — see the
[example configs](configuration.md#example-configs) (a ready-made TB config and a
non-TB template).

```bash
nextflow run main.nf \
    --tsv samples.tsv \
    --outdir results_bampiro \
    -profile standard
```

`--tsv` is **required** — there is no usable default samplesheet, so omitting it fails
with `Samplesheet not found: samples_legio.tsv` (a leftover placeholder). `-profile
standard` targets a **SLURM cluster**; on a laptop / VM / non-SLURM host use
`-profile local` (or `-profile local,docker` for Docker) — see
[Running without SLURM](installation.md#running-without-slurm).

Enable alignment-free lineage / drug-resistance typing and dual amino-acid numbering:

```bash
nextflow run main.nf \
    --tsv samples.tsv --outdir results_bampiro -profile standard \
    --run_pathotypr true --annotate_canonical true
```

When it finishes, open `results_bampiro/<samplesheet>_qc_report.html` - the
consolidated [interactive QC report](qc-report.md). See [Outputs](outputs.md) for
the full result layout and [Configuration](configuration.md) for every parameter.

> **Prefer a hands-on walkthrough?** The [Jupyter tutorial](tutorial/bampiro_tutorial.ipynb) runs the same
> journey end-to-end and lets you explore an example cohort's outputs with
> `pandas` / `matplotlib` — no pipeline run needed.

## Samplesheet

Create a Tab-Separated Value (TSV) file (e.g. `samples.tsv`) with the following columns:

| Column | Description |
| :--- | :--- |
| `sampleId` | Unique identifier for the sample (e.g. `Sample_A`). |
| `runId` | (Optional) Sequencing run ID. |
| `r1` | Path to Read 1 (FastQ). |
| `r2` | Path to Read 2 (FastQ). Leave empty for Single-End. |
| `refId` | Identifier for the reference genome (e.g. `H37Rv`, `Ecoli_K12`). |
| `refFasta` | Path to the reference FASTA file. |
| `refGff` | Path to the reference GFF file (for SnpEff annotation). |
| `taxId` | (Optional) NCBI TaxID for Kraken2 filtering (e.g. `1773` for TB, `562` for E. coli). |

**Example `samples.tsv`:**

```tsv
sampleId	runId	r1	r2	refId	refFasta	refGff	taxId
T00001	RUN1	/data/T1_R1.fq.gz	/data/T1_R2.fq.gz	H37Rv	/refs/tb.fa	/refs/tb.gff	1773
ECOLI_1	RUN2	/data/EC_R1.fq.gz	/data/EC_R2.fq.gz	K12	/refs/ecoli.fa	/refs/ecoli.gff	562
```

Reads, `refFasta` and `refGff` may be **plain or gzipped** (`.gz` is auto-detected).

### Reference requirements

- **`refGff` is required** (a mandatory column), even if you don't care about annotation —
  it builds the per-sample SnpEff database. If you truly have no annotation, give a minimal
  GFF3 stub: a `##gff-version 3` line plus one `region`/`gene` feature whose seqid matches
  the FASTA header.
- The GFF must be **GFF3**, and its **first column (seqid) must be identical to the FASTA
  header id** (up to the first whitespace), or the SnpEff build fails / annotates nothing.
  Check with:

  ```bash
  grep '^>' ref.fa | head          # FASTA contig ids
  cut -f1 ref.gff | grep -v '^#' | sort -u   # GFF seqids  (must match)
  ```

### One sample against several references

Give the **same `sampleId`** different `refId` / `refFasta` / `refGff` rows to map one
sample against several references; outputs are produced per `(sample, reference)`. (Merging,
below, groups on `sampleId` **and** `refId`, so multi-reference rows don't merge together.)

## Handling multiple runs per sample (merging)

BAMpiro automatically handles multiple sequencing runs (e.g. different lanes or
re-sequencing) for the same biological sample.

- **How to trigger merging:** assign the **same `sampleId`** to multiple rows in
  your TSV file.
- **The logic:** the pipeline processes QC and mapping for each run independently
  (in parallel) and then **merges** all BAM files for that `sampleId` before the
  deduplication and variant-calling steps.

**Example of merging 2 runs into 1 sample:**

```tsv
sampleId    runId   r1                 r2                 refId ...
Sample_A    Run_L1  /data/A_L1_R1.fq   /data/A_L1_R2.fq   H37Rv ...
Sample_A    Run_L2  /data/A_L2_R1.fq   /data/A_L2_R2.fq   H37Rv ...
```
