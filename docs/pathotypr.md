# Lineage & Drug-Resistance Typing (Pathotypr)

`--run_pathotypr true` runs
[**Pathotypr**](https://github.com/PathoGenOmics-Lab/pathotypr) - alignment-free
(k-mer) MTBC lineage + WHO drug-resistance genotyping - **directly from the
container** (a bioconda build; the old hard-coded cluster binary is gone). Because
it types straight from the reads using diagnostic k-mers, it is
**reference-agnostic**: lineage and DR calls are correct even when your samples were
mapped to a non-H37Rv reference.

Each sample gets two `split-fastq` passes:

1. **Lineage** - nested sub-lineage classification against the bundled lineage
   markers → `<sample>__<runId>.pathotypr.lineage_summary.tsv` (the `__<runId>` infix
   keeps a multi-lane / multi-reference sample from colliding; feeds the per-lineage
   panel and the lineage colours).
2. **Drug resistance** - WHO-catalogue markers, gated by `--pathotypr_min_alt`
   (default `95` = near-fixed only; lower to ~10-25 to also catch heteroresistant /
   minority alleles) → per-sample DR mutations, aggregated into
   `<samplesheet>_dr.tsv` and shown in the **Drug resistance** panel of the
   [QC report](qc-report.md) (a sample × drug matrix with the worst WHO grade per
   drug, plus a per-mutation table).

Each call carries its **WHO catalogue grade**: **1) Assoc w R** and **2) Assoc w R –
Interim** are resistance-associated (the report and any "resistant" count use grades
1–2); **3) Uncertain**; **4) Not assoc w R – Interim** and **5) Not assoc w R** are not
associated. `<samplesheet>_dr.tsv` columns: `sample, drug, gene, mutation, grade,
marker_name, af, dp`.

!!! danger "Responsible use — a genomic screen, not a clinical DST"

    This is a **genomic screen against the WHO catalogue, not a clinical DST result.**
    Treat it as a flag for review, not a diagnosis; a call's allele frequency (`af`) and
    depth (`dp`) matter, and absence of a marker is not proof of susceptibility. When
    filtering, restrict to grades 1–2 for "resistant".

The marker panels and pre-trained model (Zenodo v1.0.0, DOI
[10.5281/zenodo.19210044](https://doi.org/10.5281/zenodo.19210044)) and the
MTBC-ancestor reference are **bundled in the image** under `/opt/pathotypr/`;
override the reference or marker panels with `--pathotypr_ref` / `--pathotypr_markers`
/ `--pathotypr_dr_markers` if you supply your own. (BAMpiro types lineage by k-mer
**markers** — `split-fastq --nested-classification` — not by the RF `predict` model, so
the bundled `rf_model.pathotypr` is currently unused and no flag exposes it.) See
[Configuration](configuration.md) for every flag.

## Dual amino-acid numbering (H37Rv / Mycobrowser)

With `--annotate_canonical true`, each sample's variants are re-annotated with a
canonical snpEff database (`--canonical_snpeff_db`, default H37Rv) *in addition to*
your mapping reference. The report then shows every amino-acid change in **both**
numberings (the used reference + `--canonical_label`, default `H37Rv`) and links each
gene to its Mycobrowser locus (`Rv…`). The **amino-acid** number is exact when your
mapping reference shares (or is lifted to) H37Rv coordinates; the canonical
**coordinate** shown beside it is provided independently and alignment-free by the
k-mer liftover below, so the SNP tables carry the H37Rv position even for a non-H37Rv
reference. Turn `--annotate_canonical` off for non-MTB organisms, or point
`--canonical_snpeff_db` at another snpEff genome.

## Reference-agnostic coordinates (k-mer liftover)

`--variant_liftover true` translates every variant position into canonical (H37Rv)
coordinates for the report's SNP tables **without a whole-genome alignment and without
assuming shared coordinates**, so it works for *any* MTBC reference. It is done by
`bin/pathotypr_liftover.py` (pure Python, `lift --global-chain`): every k-mer unique in
both the mapping reference and `--canonical_ref` becomes an anchor, the anchors are
chained into collinear (and inverted) blocks, and each position is placed by
interpolating between its flanking anchors — then the placement is **verified against
the actual sequence**. The method is **"correct or absent"**: SNP sites, indels and
inversions are followed, and any position it cannot pin (a repeat / low-complexity
desert, or a rearrangement below the k-mer resolution) is *dropped* rather than
mis-mapped. It replaces the Picard/bcftools liftover for the report coordinate (the
amino-acid numbering still comes from the snpEff re-annotation above). See
`bin/pathotypr_liftover.README.md` for the algorithm, options and limits.

## Blind-spot masking (H37Rv problematic sites)

`--mask_blindspots true` adds the H37Rv Illumina **blind spots** — the repetitive /
low-mappability positions that are unreliable to call (Zenodo record
[3701840](https://zenodo.org/records/3701840), shipped as
`assets/H37Rv_blindspots.bed`, NC_000962.3 coordinates) — to each reference's
exclusion mask, so those sites are dropped from variant calling, consensus and the
report. Because the BED is in H37Rv coordinates, `--blindspot_liftover true` (default)
lifts it onto the run reference with the **same k-mer liftover**, treating
`--canonical_ref` as the FASTA the blind-spots are defined on (default: the
H37Rv-colinear MTBC ancestor bundled at `/opt/pathotypr/reference.fasta`; point it at a
true H37Rv FASTA for an exact lift), so the mask is correct even when the reference is
**not** H37Rv. Set `--blindspot_liftover false` — which appends the BED directly and
ignores `--canonical_ref` — only when the reference already shares H37Rv coordinates.
The repetitive fraction is already covered by the pipeline's own repeat /
mappability masking — the blind-spots add the non-repetitive problematic sites.
