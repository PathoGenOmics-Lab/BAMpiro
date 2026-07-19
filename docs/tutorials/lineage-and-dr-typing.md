# Lineage & drug-resistance typing

Turn your mapped reads into MTBC lineage assignments and a WHO-catalogue drug-resistance screen, then read the results without over-reading them. Enabling typing is two flags; this page is really about interpreting what comes back. **Time: ~10 minutes** (plus the run itself).

!!! warning "A genomic screen, not a clinical DST"

    Drug-resistance calls here are a **genomic screen against the WHO catalogue — not a clinical drug-susceptibility test.** Treat every call as a flag for review, not a diagnosis. Allele frequency and depth matter, and the **absence** of a marker is *not* proof of susceptibility. Never make a treatment decision from this report alone.

## What pathotypr does

[Pathotypr](../pathotypr.md) types **directly from the reads** using diagnostic k-mers, so it does two jobs in one alignment-free pass:

- **Lineage** — nested MTBC sub-lineage classification, which also colours the per-lineage panels in the [QC report](../qc-report.md).
- **Drug resistance** — WHO-catalogue marker calling, aggregated into a sample × drug matrix.

Because it reads k-mers straight from the FASTQs, it is **reference-agnostic**: lineage and DR calls stay correct even when your samples were mapped to a non-H37Rv reference. The marker panels and MTBC-ancestor reference are bundled in the container, so there is nothing extra to download.

## 1. Enable typing

Add `--run_pathotypr true` to your run:

```bash
nextflow run main.nf \
  --tsv samples.tsv \
  --outdir results_bampiro \
  --run_pathotypr true \
  -profile local,docker
```

That is enough to get lineage and DR calls for the whole cohort (the bundled 17-sample MTBC demo included).

!!! tip "Catching minority alleles"

    DR calling is gated by `--pathotypr_min_alt` (default `95`, i.e. **near-fixed variants only**). Lower it — roughly `10`–`25` — to also surface **heteroresistant / minority** alleles:

    ```bash
    --run_pathotypr true --pathotypr_min_alt 15
    ```

    Lower thresholds are more sensitive but noisier; weigh each new call by its allele frequency and depth. See [Configuration](../configuration.md) for every flag, including `--pathotypr_ref` / `--pathotypr_markers` / `--pathotypr_dr_markers` if you supply your own panels.

## 2. Read the Drug resistance panel

Open `<samplesheet>_qc_report.html` and go to the **Drug resistance** panel (under *Variants over time*). It has two views:

- a **sample × drug matrix**, coloured by the **worst WHO grade** seen for each drug, and
- a **per-mutation table** underneath.

Every call carries its **WHO catalogue grade**. Read them like this:

| Grade | Meaning | Read as |
| :--- | :--- | :--- |
| **1** Assoc w R | Associated with resistance | **Resistance-associated** |
| **2** Assoc w R – Interim | Associated (interim) | **Resistance-associated** |
| **3** Uncertain | Uncertain significance | Uncertain — review |
| **4** Not assoc w R – Interim | Not associated (interim) | Not associated |
| **5** Not assoc w R | Not associated | Not associated |

Only grades **1–2** count as "resistant" in the report and in any resistant-sample tally. Grade **3** is genuinely uncertain — flag it for a human, don't score it either way.

!!! note "The underlying table"

    The same calls are written to `<samplesheet>_dr.tsv` with columns `sample, drug, gene, mutation, grade, marker_name, af, dp`. Filter to `grade` 1–2 for a conservative "resistant" set, and always keep `af` (allele frequency) and `dp` (depth) in view when you triage.

## 3. Dual amino-acid numbering

For MTBC work you usually want mutations reported in **canonical H37Rv** coordinates as well as your mapping reference's. Add:

```bash
--run_pathotypr true --annotate_canonical true
```

The report then shows each amino-acid change in **both** numberings — your used reference *and* the canonical label (`--canonical_label`, default `H37Rv`) — and links every gene to its **Mycobrowser** locus tag (`Rv…`). This is what lets you cross-check a call against the literature even when you mapped to a non-H37Rv reference.

!!! note "Non-MTB organisms"

    Turn `--annotate_canonical` **off** for non-MTB organisms, or point `--canonical_snpeff_db` at another snpEff genome. The H37Rv numbering is meaningful only for MTBC. See [Lineage & Drug-Resistance Typing](../pathotypr.md) for the k-mer liftover and coordinate details.

## What's next

You have lineage and DR calls, and you know how far to trust them. Next, use your QC-passing samples to [**build a phylogeny**](building-a-phylogeny.md) from the cohort's consensus sequences.
