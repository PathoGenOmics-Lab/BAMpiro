# Building a phylogeny

BAMpiro stops at per-sample variants and a cohort SNP matrix — it does not build trees itself. This page takes you from those outputs to a maximum-likelihood tree, dropping your failing samples first. Budget about 15 minutes of setup, plus tree-inference runtime that scales with your cohort.

!!! note "Tools you'll need"
    The steps below use `bcftools`, `snp-sites` and `iqtree2` — standard phylogenetics tools that live **outside** the pipeline. Install them yourself (e.g. via conda/mamba or your cluster modules). For the complete recipe with every variation, follow [From outputs to a phylogeny](../downstream.md); this tutorial is the guided short path.

## 1. Start from the survivors, not the whole cohort

A tree built on low-quality genomes invents branches. So before anything else, drop the samples that failed QC.

If you worked through [Reading the QC report](reading-the-qc-report.md), you exported the exclusion basket to two files:

- **`keep_list.txt`** — the sampleIds that survived, one per line.
- **`exclusion.tsv`** — the dropped samples and the QC flag that sank them.

Everything downstream filters against `keep_list.txt`, so keep it beside your results.

!!! tip "No report export yet?"
    You can regenerate `keep_list.txt` from `<samplesheet>_qc_flags.tsv` — keep the `PASS` (and reviewed `WARN`) rows. FAIL samples start pre-selected in the basket for exactly this reason.

## 2. Pick a route and assemble the alignment

Two pipeline outputs can seed a tree. Use whichever you published — for the demo MTBC cohort either works.

=== "Consensus FASTAs (Route A)"

    Each sample's `*.consensus.fasta` is a whole-genome consensus, masked at low-confidence sites. Concatenate only the survivors:

    ```bash
    # one multi-FASTA of the kept samples (all mapped to the SAME reference!)
    while read s; do cat results/*/"$s".*.consensus.fasta; done < keep_list.txt > kept.aln.fasta
    ```

=== "All-positions VCFs (Route B)"

    The `*.all.pos.vcf.gz` "backbone" VCFs carry every position (WT + variant), so merging them yields a gap-free SNP alignment — the better choice for recombination-aware / ML tools:

    ```bash
    bcftools merge results/*/*.all.pos.vcf.gz -Oz -o merged.vcf.gz && bcftools index merged.vcf.gz
    snp-sites -c -o cohort.snps.fasta <(bcftools view merged.vcf.gz)   # variable sites only
    ```

!!! warning "Same reference, always"
    Every sample in the alignment must have been mapped against the **same** reference genome. Mixing references produces meaningless columns.

## 3. Mask the problematic sites

Repeats, poor-mappability regions and low-confidence positions create false SNPs. Before tree-building, drop:

- the columns in each sample's `stats/*.mask_sites.tsv`,
- the reference repeat regions in `Locus_to_exclude_*.txt`,
- and — for MTBC only — the H37Rv Illumina blind-spots in `assets/H37Rv_blindspots.bed`.

The exact filtering commands live in [From outputs to a phylogeny](../downstream.md) — apply them to your alignment before the next step rather than duplicating them here.

## 4. Infer the tree

Feed the masked alignment to IQ-TREE with a GTR+G model and 1000 ultrafast bootstraps:

```bash
# Route A: iqtree2 -s kept.aln.fasta   -m GTR+G -B 1000 -T AUTO
# Route B: iqtree2 -s cohort.snps.fasta -m GTR+G -B 1000 -T AUTO
iqtree2 -s kept.aln.fasta -m GTR+G -B 1000 -T AUTO
```

The resulting `*.treefile` is a Newick tree you can open in FigTree, iTOL or `ggtree`. Colour the tips with the pipeline's lineage calls (see [Lineage & drug-resistance typing](lineage-and-dr-typing.md)) to read structure at a glance.

!!! tip "Sanity-check the topology"
    For the 17-sample demo cohort, samples of the same MTBC lineage should cluster together. If they don't, revisit your masking and confirm every input mapped to the same reference.

## Where to go next

You've completed the tutorial path — from a first run to a typed, filtered, phylogeny-ready cohort. To go deeper:

- [From outputs to a phylogeny](../downstream.md) — the full recipe, including every masking and matrix-subsetting variation.
- [Outputs](../outputs.md) — a complete map of every file the pipeline writes.
- The bundled **Jupyter notebook** shows how to load the SNP matrix and QC TSVs and explore them with pandas.
