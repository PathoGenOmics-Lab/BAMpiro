# H37Rv Illumina "blind spots" mask

`H37Rv_blindspots.bed` — genome positions in the *M. tuberculosis* H37Rv reference
(`NC_000962.3`, 4,411,532 bp) where short-read Illumina calling is unreliable and that
should be excluded from SNP analysis.

- **Source:** Zenodo record [3701840](https://zenodo.org/records/3701840) — *"Supplemental
  Table 7: Blind Spots and their attributes"* (`S7-seq-attributes.csv`, one row per H37Rv
  position with `pal,homopolymer,GCrich,repetitive,blindspot` flags).
- **Derivation:** all positions with `blindspot == 1`, collapsed into intervals.
- **Coordinates:** `NC_000962.3`, **0-based half-open** BED (`contig  start  end  name`).
- **Extent:** 159,659 positions = **3.62 %** of the genome, in **2,292** intervals.

Regenerate:

```bash
curl -sL "https://zenodo.org/records/3701840/files/S7-seq-attributes.csv?download=1" -o S7.csv
awk -F, 'NR>1 && $6==1{print $1}' S7.csv | \
  awk 'BEGIN{OFS="\t"} NR==1{s=p=$1;next} $1==p+1{p=$1;next}
       {print "NC_000962.3",s-1,p,"blindspot"; s=p=$1}
       END{print "NC_000962.3",s-1,p,"blindspot"}' > H37Rv_blindspots.bed
```

## Use in BAMpiro

Enable it with `--mask_blindspots true` (opt-in): the intervals are added to the per-reference
exclusion, so they are dropped from variant calling, the consensus, and shown in the report's
*mask regions* panel.

By default (`blindspot_liftover = true`) the BED is **lifted onto the run's reference with pathotypr**
(`bin/pathotypr_liftover.py` -> `pathotypr classify`, alignment-free k-mer liftover), so it is correct
even when the reference does **not** share H37Rv coordinates. Set `blindspot_liftover = false` to append
the BED directly (only correct when the reference already shares H37Rv coordinates, e.g. H37Rv or the
MTBC ancestor). `canonical_ref` is the FASTA the blind-spot coordinates are defined on (default: the
container's MTBC-ancestor reference, which is H37Rv-colinear).
