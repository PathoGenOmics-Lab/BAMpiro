# pathotypr_liftover.py — k-mer coordinate liftover between MTBC references

Translates a set of genome positions from one MTBC reference to another **without a whole-genome
alignment and without assuming the references share coordinates**: for a position it takes the
flanking-context k-mer and finds it in the target genome; the target position is the equivalent coordinate.

## `lift` — the recommended method (self-contained, synteny-anchored)

```bash
# lift the positions in POS (BED or 1-based list) from reference A onto reference B:
python3 pathotypr_liftover.py lift POS A.fasta B.fasta --out-map map.tsv --out-bed lifted.bed --contig B --kmer-size 21
```

Builds **anchors** from k-mers unique in both genomes, then for a position whose k-mer **recurs** (a repeat)
picks the occurrence **consistent with the surrounding anchors (synteny)** instead of dropping it — so it
never mis-maps and recovers most repeat positions. Pure Python, no pathotypr call. Validated on the real
H37Rv / MTBC-ancestor pair: **≈99 % of positions lift with zero wrong coordinates**, it tracks indels
(before an insertion → identity, after → shifted), and the DR sites (rpoB 761155, katG 2155168, gyrA 7570,
rrs 1473246) all map correctly. Tune repeat resolution with `--max-shift` / `--min-margin`.

## `markers` + `apply` — the `pathotypr classify` alternative (Rust, for very large sets)

```bash
python3 pathotypr_liftover.py markers POS A.fasta -o markers.tsv
printf 'genome\tpath\nB\tB.fasta\n' > genomes.tsv            # --tsv_genomes: --fasta-genomes alone trips a
pathotypr classify --tsv_pos markers.tsv --ref_fasta A.fasta \
        --tsv_genomes genomes.tsv --output classify_out --kmer-size 21   # required-args bug in pathotypr 0.1.0
python3 pathotypr_liftover.py apply classify_out --out-map map.tsv --contig B \
        --offset 1 --source-fasta A.fasta --target-fasta B.fasta --kmer-size 21   # ALWAYS pass both fastas
```

classify reports only one occurrence per k-mer, so `apply` **drops** repeats (doesn't synteny-resolve them);
`--source-fasta`/`--target-fasta` are required so a non-unique k-mer is dropped rather than mis-mapped.

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
