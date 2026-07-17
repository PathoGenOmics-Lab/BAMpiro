# pathotypr_liftover.py — k-mer coordinate liftover between MTBC references

Translates a set of genome positions from one MTBC reference to another **without a whole-genome
alignment and without assuming the references share coordinates**, using `pathotypr classify`:
for a position it builds the flanking-context k-mer and locates it in the target genome, which reports
the target coordinate (`snp_position`) next to the source one (`ref_position`).

Two halves, with `pathotypr classify` (in the container) in between:

```bash
# A -> B liftover of the positions in POS (BED or 1-based list), context built on reference A:
python3 pathotypr_liftover.py markers POS  A.fasta  -o markers.tsv
pathotypr classify --tsv_pos markers.tsv --ref_fasta A.fasta --fasta_genomes B.fasta -o classify_out
python3 pathotypr_liftover.py apply classify_out.tsv --out-map map.tsv \
        --out-bed lifted.bed --contig <B_contig> --offset 1
```

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
