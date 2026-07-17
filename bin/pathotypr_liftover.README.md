# pathotypr_liftover.py — k-mer coordinate liftover between MTBC references

Translates a set of genome positions from one MTBC reference to another **without a whole-genome
alignment and without assuming the references share coordinates**, using `pathotypr classify`:
for a position it builds the flanking-context k-mer and locates it in the target genome, which reports
the target coordinate (`snp_position`) next to the source one (`ref_position`).

Two halves, with `pathotypr classify` (in the container) in between:

```bash
# A -> B liftover of the positions in POS (BED or 1-based list), context built on reference A:
python3 pathotypr_liftover.py markers POS  A.fasta  -o markers.tsv
printf 'genome\tpath\nB\tB.fasta\n' > genomes.tsv          # --tsv_genomes: --fasta-genomes alone trips a
pathotypr classify --tsv_pos markers.tsv --ref_fasta A.fasta \
        --tsv_genomes genomes.tsv --output classify_out --kmer-size 21   # required-args bug in pathotypr 0.1.0
python3 pathotypr_liftover.py apply classify_out --out-map map.tsv --out-bed lifted.bed --contig B \
        --offset 1 --source-fasta A.fasta --target-fasta B.fasta --kmer-size 21   # uniqueness guard
```

**Always pass `--source-fasta` and `--target-fasta`.** A k-mer that recurs in either genome makes the lift
ambiguous (pathotypr keeps the last occurrence -> a WRONG coordinate; seen on real MTBC at a duplicated
segment around H37Rv 2,300,000 and 832,200). The guard drops those, so `apply` only ever emits a coordinate
that is *correct*. Validated on the real H37Rv vs MTBC-ancestor pair: with the guard, 0 wrong coordinates
(≈96 % of positions lift; the rest are dropped, not mis-mapped), and the DR sites (rpoB 761155, katG
2155168, gyrA 7570, rrs 1473246) all map correctly.

Flag style is mixed in pathotypr 0.1.0: `--tsv_pos` / `--ref_fasta` / `--tsv_genomes` keep underscores,
but `--kmer-size` / `--fasta-genomes` use dashes. The classify main output is written to the `--output`
name exactly (no extension); its columns are `genome  k-mer  k-merPOS  SNPgenome  SNPreference  lineage`.

`map.tsv` = `src_pos <TAB> tgt_pos`; `lifted.bed` = the positions collapsed into intervals in B's coords.

## The two BAMpiro use cases

- **(A) SNP tables — variant position in H37Rv.** Lift the run's variant positions (mapping-reference
  coords) **to H37Rv**: `A = mapping reference`, `B = H37Rv`. `snp_position` is the H37Rv coordinate; feed
  it as the report's `pos_h37rv` (the same field the `--vcfs-h37rv` path fills). Alignment-free, replaces
  the Picard/bcftools liftover.
- **(B) Blind-spots mask on any reference.** Lift `assets/H37Rv_blindspots.bed` **from H37Rv to the run's
  reference**: `A = H37Rv`, `B = mapping reference`. The resulting BED (target coords) goes into the
  per-reference exclusion — so the mask works even when the reference is not H37Rv.

## Calibration and limits (be honest)

- **Coordinate convention:** `snp_position = kmer_start + k/2` is 0-based-ish; run an **identity liftover**
  (`A == B`) once and set `--offset` so `tgt_pos == src_pos`. For k=21 it is typically `--offset 1`.
- **A position does NOT lift if** its flanking k-mer is not unique in B (repeats — i.e. exactly the
  repetitive blind-spots) or is broken by a nearby A/B difference. `apply` reports the lift/ambiguous
  counts. For the blind-spots the repetitive fraction is already covered by the pipeline's own
  nucmer + genmap masking on the target reference; this liftover adds the non-repetitive problematic sites.
- **k-mer size** must match between `markers`/`classify` (default 21).
