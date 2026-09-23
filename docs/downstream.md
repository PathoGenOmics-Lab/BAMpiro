# From outputs to a phylogeny

BAMpiro stops at per-sample variants + a cohort SNP matrix; here is how to turn those
into an alignment and a tree. Pick the route that matches what you published.

## Route A — consensus FASTAs → tree

Each sample's `*.consensus.fasta` is a whole-genome consensus (masked at low-confidence
sites). Concatenate the ones you want to keep and build a tree:

```bash
# one multi-FASTA of all kept samples (same reference!)
cat results/*/*.consensus.fasta > cohort.aln.fasta
iqtree2 -s cohort.aln.fasta -m GTR+G -B 1000 -T AUTO
```

## Route B — all-positions VCFs → SNP alignment → tree

The `*.all.pos.vcf.gz` "backbone" VCFs carry every position (WT + variant), so merging
them gives a gap-free SNP alignment (better for recombination-aware / ML tools):

```bash
bcftools merge results/*/*.all.pos.vcf.gz -Oz -o merged.vcf.gz && bcftools index merged.vcf.gz
snp-sites -c -o cohort.snps.fasta <(bcftools view merged.vcf.gz)   # variable sites only
iqtree2 -s cohort.snps.fasta -m GTR+G -B 1000 -T AUTO
```

**Mask problematic sites** before tree-building: drop the columns listed in each
sample's `stats/*.mask_sites.tsv` and the reference `Locus_to_exclude_*.txt` (repeats /
mappability), and — for MTBC — the H37Rv blind-spots (`assets/H37Rv_blindspots.bed`).

## Applying the report's exclusion list

The interactive report exports two files from the exclusion list (on its *Sample QC* page):

- **`keep_list.txt`** — the sampleIds that survived (one per line).
- **`exclusion.tsv`** — the dropped samples and the reason (which QC flag).

Filter your alignment to the survivors before building the tree, e.g.:

```bash
# keep only the survivors' consensus FASTAs
while read s; do cat results/*/"$s".*.consensus.fasta; done < keep_list.txt > kept.aln.fasta

# ...or subset the master SNP matrix to the kept samples' AF/DP columns
python3 - <<'PY'
import pandas as pd
keep = set(open('keep_list.txt').read().split())
m = pd.read_csv('results/<samplesheet>_snp_matrix.tsv', sep='\t')
meta = [c for c in m.columns if '|' not in c]
cols = meta + [c for c in m.columns if c.split('|')[0] in keep]
m[cols].to_csv('snp_matrix.kept.tsv', sep='\t', index=False)
PY
```

The [Jupyter tutorial](tutorial/bampiro_tutorial.ipynb) shows how to load these TSVs and
explore them with pandas.
