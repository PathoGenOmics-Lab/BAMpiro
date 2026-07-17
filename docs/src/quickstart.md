# Quick Start

By default BAMpiro assumes *M. tuberculosis* settings (ploidy = 2 to detect mixed
infections). Lineage/DR typing is **off** by default - enable it with
`--run_pathotypr true`.

```bash
nextflow run main.nf \
    --tsv samples.tsv \
    --outdir results_bampiro \
    -profile standard
```

`--tsv` is **required** — there is no usable default samplesheet, so omitting it fails
with `Samplesheet not found: samples_legio.tsv` (a leftover placeholder). The only
profile defined is `standard` (which runs on SLURM, the default executor).

Enable alignment-free lineage / drug-resistance typing and dual amino-acid numbering:

```bash
nextflow run main.nf \
    --tsv samples.tsv --outdir results_bampiro -profile standard \
    --run_pathotypr true --annotate_canonical true
```

When it finishes, open `results_bampiro/<samplesheet>_qc_report.html` - the
consolidated [interactive QC report](qc-report.md). See [Outputs](outputs.md) for
the full result layout and [Configuration](configuration.md) for every parameter.

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
