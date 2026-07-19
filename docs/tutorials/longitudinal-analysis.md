# Longitudinal & within-host analysis

When you sample the same host (or passage line) at several timepoints, BAMpiro can track how each
variant's **allele frequency** moves over time and flag pairs of variants that rise and fall together.
This tutorial shows you how to set the run up for that and how to read the two panels it unlocks —
**SNP dynamics** and **Epistasis**. Allow about **15 minutes**; the bundled 17-sample demo cohort
(patients sampled over several timepoints) is exactly this shape, so you can click along with the
[interactive demo report](../examples/qc_report_demo.html){ target="_blank" }.

!!! note "What you need"

    Per-sample **VCFs** (the pipeline produces these) **and** a samplesheet with a **time** column and a
    **group** column. With both present the panels appear automatically; without them they self-hide.

## 1. Design the samplesheet

The report connects samples into trajectories purely from the metadata columns — no extra flags. Add:

- a **time** column — any of `timepoint`, `day`, `date`, `week`, `month`, `passage`, `generation`, `visit`… — the **x-axis** of every trajectory, and
- a **group** column — `patient`, `samples`, `series`, `host`, `subject`, `cluster`… — that ties the samples of one host/line into a single series.

```tsv
sampleId    r1              r2              refId   refFasta        refGff          taxId   patient   timepoint
TB-P1-d0    …/P1-d0_R1.fq   …/P1-d0_R2.fq   H37Rv   …/H37Rv.fa.gz   …/H37Rv.gff.gz  1773    TB-P1     0
TB-P1-d60   …/P1-d60_R1.fq  …/P1-d60_R2.fq  H37Rv   …/H37Rv.fa.gz   …/H37Rv.gff.gz  1773    TB-P1     60
TB-P1-d180  …/P1-d180_R1.fq …/P1-d180_R2.fq H37Rv   …/H37Rv.fa.gz   …/H37Rv.gff.gz  1773    TB-P1     180
```

Here `patient` groups the three samples into one series and `timepoint` orders them. See
[The samplesheet](samplesheet.md) for how every metadata column is matched by name.

## 2. Read the SNP dynamics

Open the report and go to **SNP dynamics** (under *Variants over time*). For every connected series it
draws each variant's **allele-frequency trajectory over time**, with the **per-timepoint read depth**
beneath so you can weigh a move against its support. Events are marked as the frequency crosses key
thresholds:

- **emergence** — a variant appears / climbs,
- **fixation** — it reaches near-100 %,
- **loss** — it falls back to the reference,
- **non-synonymous** — the change alters an amino acid (outlined), the ones most likely under selection.

Use the **series filter** to focus on one host, and **zoom** to inspect a busy region. A variant that
rises and fixes across a treatment course is a textbook candidate for positive selection — cross-check it
against the [Drug resistance](lineage-and-dr-typing.md) calls.

## 3. Read the Epistasis panel

**Epistasis** looks for **pairs of variants whose trajectories co-vary** within a series:

- **concordant** — they rise and fall together (candidate linkage or co-selection), or
- **discordant** — one rises as the other falls (competing lineages / clonal interference).

Each pair is scored by the **Pearson correlation** of the two trajectories, a **permutation *p*-value**
(how often shuffling the timepoints reaches that \|r\| by chance), and a **Benjamini–Hochberg FDR *q***
across every reported pair. Views: cards, matrix or table.

!!! warning "Recurrence is what makes a pair 'strong'"

    With only a few timepoints a **single** series can't beat chance on its own — so a pattern that
    **recurs across independent series** is what drives a pair to **strong** (q ≤ 0.05). Read every pair
    as a **candidate** linked / co-selected / competing set of SNPs, **not** proof of a functional
    interaction.

## 4. Complementary views

- **SNP matrix** — the full site × sample allele-frequency matrix; filter by gene / position or by a
  metadata column, and export it as a TSV for your own analysis.
- **Variant × dose** — if your samplesheet also carries a numeric `dose` column, the report correlates
  each variant's allele frequency with dose across the cohort (Spearman + FDR). See
  [Reading the QC report](reading-the-qc-report.md#7-the-dose-analyses).

The exhaustive panel reference — every control and caveat — is in the
[Interactive QC Report](../qc-report.md).

## Where to go next

- Turn a longitudinal cohort into a tree with [Building a phylogeny](building-a-phylogeny.md).
- Run the numbers yourself: the [Jupyter notebook](../tutorial/bampiro_tutorial.ipynb) opens the SNP
  matrix and metadata in `pandas`.
- Back to the [Tutorials overview](index.md).
