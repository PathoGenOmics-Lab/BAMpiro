# pathotypr_liftover.py — k-mer coordinate liftover between MTBC references

Translates a set of genome positions from one MTBC reference to another **without a whole-genome
alignment and without assuming the references share coordinates**: for a position it takes the
flanking-context k-mer and finds it in the target genome; the target position is the equivalent coordinate.

## `lift` — the recommended method (self-contained, synteny-anchored)

```bash
# lift the positions in POS (BED or 1-based list) from reference A onto reference B:
python3 pathotypr_liftover.py lift POS A.fasta B.fasta --out-map map.tsv --out-bed lifted.bed --contig B \
        --kmer-size 21 --global-chain
```

**`--global-chain` (recommended, what BAMpiro uses)** builds a whole-genome coordinate map: every k-mer
unique in BOTH genomes is an anchor (~97 % of the MTBC genome — anchors roughly every base), chained into a
collinear order (LIS), and any position is placed by **interpolating between its flanking anchors**. So a
position is mapped by its NEIGHBOURS, not by its own k-mer — SNP sites and other differing positions lift
correctly, and indels are followed. Benchmarked on the real H37Rv / MTBC-ancestor pair vs the per-position
method:

| | per-position (default) | `--global-chain` |
|---|---|---|
| colinear, 500 positions | 99.0 % lift | **100 % lift** |
| SNP AT the query site (e.g. DR rpoB 761155) | **dropped** | **placed correctly** |
| +25 bp indel | tracked | tracked |
| wrong coordinates | 0 | 0 |

Cost: builds the map by scanning each genome once — ≈6 s and ≈2.1 GB for a 4.4 Mb genome (a per-reference,
one-time step). Without `--global-chain`, `lift` uses the lighter per-position method: anchors from
unique-in-both k-mers, and repeats resolved by the nearest anchors (`--max-shift` / `--min-margin`), but a
SNP at the queried site drops the position.

### Anchor-chain design — "correct or absent"

A position is placed **only** when it is bracketed by two consecutive anchors of a chain across a **colinear**
gap (source span == target span, so no indel lies between them). It takes the chain whose bracketing gap is
**tightest**; a forward gap that itself contains reverse anchors is *not* interpolated (it would be crossing an
inverted block). So SNP sites and colinear regions are exact, and anything ambiguous **drops** rather than
receiving a smeared or extrapolated coordinate. Validated on the real H37Rv / ancestor pair (DR sites, 100 %
lift, 0 wrong) and on synthetic constructions of every failure mode below (18/18 scenarios, 0 wrong coords).

- **`--sample N` (FracMinHash, memory).** Keep only ~1/N of the k-mers as anchors — deterministically, hashing
  the **canonical** code (`hash(min(kmer, revcomp(kmer))) % N == 0`), so a k-mer *and its reverse complement*
  share one keep/drop decision and forward **and** reverse (inversion) anchors both co-survive at rate 1/N. On
  the 4.4 Mb pair `--sample 10` drops the anchor map from **2.1 GB to ~300 MB** (7×), still **0 wrong** and
  100 % lift. BAMpiro uses `--sample 10` for the **blind-spots mask** (memory) and **no sampling** for the
  **variant-coordinate** lift (max anchor density near indels).
- **Inversions / rearrangements (automatic).** A REVERSE chain of reverse-complement anchors (longest
  *decreasing* run) places a position inside an inverted block by the reflected mapping `s+e−P+1`. It works
  under sampling (canonical hashing above) and for small inversions (the forward chain refuses to cross a gap
  containing reverse anchors, so boundary positions drop instead of mis-mapping). A parity correction makes it
  exact for even k too. Verified on 100 kb and 400 bp inversions, with `--sample 10`/`30`, and with k=20.
- **Indels — deletions AND insertions (automatic + optional RD BED).** A gap whose source and target spans
  differ holds an indel; its interior (the ±k breakpoint shadow, which has no clean equivalent) is **dropped**,
  never smeared across the indel — in either direction (target lost *or* gained sequence). A forward gap where
  the source advanced ≥ `--rd-min` (default 50) bp more than the target is additionally reported as an **RD
  deletion** to `--rd-out FILE` (source coords). Verified on 300 bp insertion / deletion and a 5 kb RD block.

The one blind spot is inherent to k-mers: a rearrangement **shorter than k** (e.g. a sub-21 bp inversion) has
no k-mer inside it and preserves the flanking gap length, so it is invisible and its few interior positions
take the colinear (identity) coordinate. Keep k below any structural feature you must resolve.

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
- **At an A/B indel breakpoint** the k-mers spanning it don't anchor, so the ±k window around the breakpoint
  is not colinear and those positions **drop** (correct-or-absent) rather than getting a smeared coordinate;
  `--sample` widens this dropped window slightly. Away from indels the coordinate is exact. The only case that
  mis-maps is a rearrangement shorter than k (invisible to k-mers) — see the design note above.
- **k-mer size** must match between `markers`/`classify` (default 21).
