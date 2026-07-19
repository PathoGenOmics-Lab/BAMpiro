# The samplesheet

The samplesheet is BAMpiro's one required input: a tab-separated table with **one row per `(sample, run, reference)`**. This page walks you through the columns, a worked example, how to merge runs, and the optional metadata columns that unlock extra panels in the QC report. Est. time: ~10 minutes.

If you have not run the pipeline at all yet, start with [Your first run](first-run.md).

## Required columns

Create a plain TSV (e.g. `samples.tsv`) — real tabs between fields, one header row — with these columns:

| Column | Required | Description |
| :--- | :---: | :--- |
| `sampleId` | ✔ | Unique sample identifier (e.g. `TB-P1-d0`). |
| `runId` | | Sequencing run / lane ID. |
| `r1` | ✔ | Path to Read 1 (FastQ). |
| `r2` | | Path to Read 2 (FastQ). Leave empty for single-end. |
| `refId` | ✔ | Reference genome identifier (e.g. `H37Rv`). |
| `refFasta` | ✔ | Path to the reference FASTA. |
| `refGff` | ✔ | Path to the reference GFF3 (builds the per-sample SnpEff database). |
| `taxId` | | NCBI TaxID for Kraken2 filtering (e.g. `1773` for TB). |

Reads, `refFasta` and `refGff` may be **plain or gzipped** (`.gz` is auto-detected).

!!! warning "`refGff` is mandatory"

    Even if you don't care about annotation, `refGff` is a required column — and its
    first column (seqid) must match the FASTA header id. If you have no annotation, pass a
    minimal GFF3 stub. See [Reference requirements](../quickstart.md#reference-requirements)
    for the stub and the id-matching check.

## A worked example

=== "Minimal"

    ```tsv
    sampleId	runId	r1	r2	refId	refFasta	refGff	taxId
    TB-P1-d0	RUN1	/data/TB-P1-d0_R1.fq.gz	/data/TB-P1-d0_R2.fq.gz	H37Rv	/refs/H37Rv.fa.gz	/refs/H37Rv.gff.gz	1773
    ```

    Single-end? Leave the `r2` field empty (keep the tab).

=== "With metadata"

    ```tsv
    sampleId	r1	r2	refId	refFasta	refGff	taxId	patient	timepoint	treatment	dose	site
    TB-P1-d0	/data/TB-P1-d0_R1.fq.gz	/data/TB-P1-d0_R2.fq.gz	H37Rv	/refs/H37Rv.fa.gz	/refs/H37Rv.gff.gz	1773	TB-P1	0	HRZE	360	Madrid
    TB-P1-d60	/data/TB-P1-d60_R1.fq.gz	/data/TB-P1-d60_R2.fq.gz	H37Rv	/refs/H37Rv.fa.gz	/refs/H37Rv.gff.gz	1773	TB-P1	60	HRZE	360	Madrid
    ```

    The trailing columns are the optional metadata described [below](#optional-metadata) —
    they are the flavour of the bundled 17-sample MTBC demo cohort.

## Merging multiple runs per sample

To combine several sequencing runs (extra lanes, re-sequencing) into one biological sample, give them the **same `sampleId`** on separate rows. BAMpiro processes each run's QC and mapping in parallel, then **merges** the BAMs before deduplication and variant calling.

```tsv
sampleId	runId	r1	r2	refId	refFasta	refGff	taxId
TB-P1-d0	Run_L1	/data/d0_L1_R1.fq.gz	/data/d0_L1_R2.fq.gz	H37Rv	/refs/H37Rv.fa.gz	/refs/H37Rv.gff.gz	1773
TB-P1-d0	Run_L2	/data/d0_L2_R1.fq.gz	/data/d0_L2_R2.fq.gz	H37Rv	/refs/H37Rv.fa.gz	/refs/H37Rv.gff.gz	1773
```

## One sample against several references

Give the **same `sampleId`** different `refId` / `refFasta` / `refGff` rows to map one sample against multiple references — outputs are produced per `(sample, reference)`. Merging groups on `sampleId` **and** `refId`, so these multi-reference rows are kept apart rather than merged.

## Optional metadata

Here is the payoff. The QC report reads richer context **straight from this samplesheet** — columns are matched **by name (case-insensitive)**, so there is nothing to configure: add a column and the matching panel lights up (and hides itself when absent). The pipeline's own columns (`r1`, `r2`, `refId`, `refFasta`, `refGff`, `taxId`, `runId`) are ignored, so this is purely additive to a normal samplesheet.

| Add a column like… | And you unlock |
| :--- | :--- |
| **time** — `timepoint`, `day`, `date`, `week`, `passage`… | The x-axis of the **SNP dynamics** allele-frequency trajectories. |
| **group / series** — `patient`, `series`, `host`, `subject`, `samples`… | Connects samples into one longitudinal series for **SNP dynamics** and **epistasis**. |
| any other categorical — `site`, `ward`, `batch`… | A **cohort-filter dropdown** in the toolbar (restrict the whole view to one value) and a SNP-matrix header level. |
| **`treatment`** (categorical) | The cohort filter **plus** the **Dose × treatment** test (per-group dose distribution + Kruskal–Wallis). |
| **`dose`** (numeric) | A first-class **metric** (selectable on scatter axes and the correlation matrix, with a Spearman *r*) that drives the **Variant × dose** association scan. |

!!! note "When the dynamics panels appear"

    SNP dynamics and epistasis need **both** a time and a group column **and** per-sample
    VCFs. The metadata keys on your first column (`sampleId`). For the complete column table
    and every panel a column feeds, see [Interactive QC Report](../qc-report.md).

## What's next

You have a valid samplesheet. Next, tune the run itself — profiles, references, pathotypr, and QC thresholds — in [Configuring a run](configuring-a-run.md).
