# Lineage & Drug-Resistance Typing (Pathotypr)

`--run_pathotypr true` runs
[**Pathotypr**](https://github.com/PathoGenOmics-Lab/pathotypr) — alignment-free
(k-mer) MTBC lineage + WHO drug-resistance genotyping — **directly from the
container** (a bioconda build; the old hard-coded cluster binary is gone). Because
it types straight from the reads using diagnostic k-mers, it is
**reference-agnostic**: lineage and DR calls are correct even when your samples were
mapped to a non-H37Rv reference.

Each sample gets two `split-fastq` passes:

1. **Lineage** — nested sub-lineage classification against the bundled lineage
   markers → `<sample>.pathotypr.lineage_summary.tsv` (feeds the per-lineage panel
   and the lineage colours).
2. **Drug resistance** — WHO-catalogue markers, gated by `--pathotypr_min_alt`
   (default `95` = near-fixed only; lower to ~10–25 to also catch heteroresistant /
   minority alleles) → per-sample DR mutations, aggregated into
   `<samplesheet>_dr.tsv` and shown in the **Drug resistance** panel of the
   [QC report](qc-report.md) (a sample × drug matrix with the worst WHO grade per
   drug, plus a per-mutation table).

The marker panels and pre-trained model (Zenodo v1.0.0, DOI
[10.5281/zenodo.19210044](https://doi.org/10.5281/zenodo.19210044)) and the
MTBC-ancestor reference are **bundled in the image** under `/opt/pathotypr/`;
override any of them with `--pathotypr_ref` / `--pathotypr_markers` /
`--pathotypr_dr_markers` / `--pathotypr_rf_model` if you supply your own. See
[Configuration](configuration.md) for every flag.

## Dual amino-acid numbering (H37Rv / Mycobrowser)

With `--annotate_canonical true`, each sample's variants are re-annotated with a
canonical snpEff database (`--canonical_snpeff_db`, default H37Rv) *in addition to*
your mapping reference. The report then shows every amino-acid change in **both**
numberings (the used reference + `--canonical_label`, default `H37Rv`) and links each
gene to its Mycobrowser locus (`Rv…`). Turn it off for non-MTB organisms, or point
`--canonical_snpeff_db` at another snpEff genome. This is only meaningful when the
mapping reference shares the canonical reference's coordinates.
