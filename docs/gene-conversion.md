# Gene conversion

Gene conversion copies a stretch of one paralog onto another without a reciprocal exchange. The
acceptor locus stops carrying its own alleles and carries the donor's over a bounded **tract**.
Against a reference that never saw the event, that reads as a run of variants which are not
random: they are exactly the donor's bases, in order, between two breakpoints.

This is off by default. Turn it on with:

```bash
nextflow run main.nf --tsv samples.tsv --outdir results -profile local,docker \
  --find_gene_conversion true
```

!!! warning "It needs `--exclude_repeats true`, which is the default"

    The paralog map is read from the reference self-alignment that repeat masking already
    computes. With `--exclude_repeats false` that alignment is never produced and the stage
    finds nothing to work with.

## Why this needs to bypass the rest of the pipeline

Conversion happens in paralogous sequence, which is what most of BAMpiro spends effort removing:

| Default | What it does | Why it is a problem here |
| :--- | :--- | :--- |
| `exclude_repeats = true` | Drops repeats from variant calling | The tract sits inside a repeat |
| `dynamic_read_filter = true` | Drops reads that cannot be placed uniquely | Those are the reads carrying the evidence |

Both are right for variant calling. So this stage does not use the filtered BAM or the masked
VCF: it reads the **deduplicated, pre-filter** alignment directly, and it does **not** filter on
mapping quality, because a paralogous region is exactly where an aligner assigns MAPQ 0.

That also means the tract it reports is, by construction, in a region the rest of the run has
excluded. A tract will not appear in the SNP matrix or the consensus, and that is expected.

## How it works

**1. The paralog map.** `nucmer --maxmatch` already aligns the reference against itself for repeat
masking, and the masking step collapses the result to a flat list of intervals, discarding which
copy aligns to which. That pairing is recovered here, and `show-snps` on the same alignment gives
the positions where the two copies differ.

**2. Diagnostic sites.** Only those differing positions can carry evidence. Everywhere else the two
copies are identical and a read is uninformative by construction, so the analysis works over that
list rather than over the whole locus.

**3. Tracts.** A candidate is a run of *consecutive diagnostic sites* where the reads carry the
donor allele. Consecutive means adjacent in the site list, not in base coordinates: the bases
between two diagnostic sites are the same in both copies and cannot testify either way.

**4. Not believing it.** This is the part that decides whether the output is useful. Reads that
mismap from the donor produce the same per-site picture as a real conversion. Three things
separate them, and all three are reported rather than collapsed into one score:

| Column | Real conversion | Mismapping |
| :--- | :--- | :--- |
| `donor_af_outside` | Near 0: the tract is bounded | Elevated: donor alleles everywhere in the locus |
| `donor_af_in` | Near 1 in a clonal sample | Intermediate, and much the same at every site |
| `breakpoint_reads` | Above 0 | Always 0 |

A site with too little depth to genotype is **undetermined**, and undetermined testifies in
neither direction. It neither joins a tract nor ends one, and it stays out of the
`donor_af_outside` average. Treating it as "the acceptor allele is here" would mean one ordinary
coverage dip inside a clean tract splits it into fragments that each fall below
`--gconv_min_sites`, and the tract disappears. `n_undetermined` reports how many sites inside a
tract were skipped that way, so a tract does not look more solid than the data supports.

`breakpoint_reads` is the one a mismapping cannot fake. It counts single molecules carrying donor
alleles on one side of a breakpoint and acceptor alleles on the other, in cis. A mismapped read is
a donor read: it carries donor alleles everywhere it reaches and never crosses back.

## Verdicts

| Verdict | Meaning |
| :--- | :--- |
| `gene_conversion` | A read crosses a breakpoint in cis, or the tract is bounded with the donor allele near-fixed and supported by reads spanning several sites |
| `mismapping` | Donor alleles are present outside the tract as well, so it is not bounded |
| `ambiguous` | Reported, but not called. The `reason` column says which test it failed |

The most common `ambiguous` case is a tract covering **every** diagnostic site of the locus. With
no site outside it, there is nothing the donor alleles are bounded by, and a whole locus replaced
by its paralog is indistinguishable from a locus whose reads all arrived from its paralog. Absence
of contrary evidence is not evidence, so it is not called.

## Output

Per sample, next to the other per-sample files:

```
<sample>.<ref>.gene_conversion.tsv
```

and one cohort-level `<samplesheet>_gene_conversion.tsv`. Columns: `sample`, `pair_id`, `contig`,
`donor`, `verdict`, `reason`, `start`, `end`, `span_bp`, `n_sites`, `n_sites_outside`,
`donor_af_in`, `donor_af_outside`, `n_undetermined`, `min_depth`, `cis_reads`,
`breakpoint_reads`, `donor_only_reads`.

The paralog map itself is published under `references/<refId>/` as
`<refId>.paralog_pairs.tsv` and `<refId>.paralog_sites.tsv`, and is worth a look on its own: it is
the donor/acceptor graph of your reference.

## Tuning

| Parameter | Default | What it controls |
| :--- | :--- | :--- |
| `--gconv_min_identity` | `90.0` | Ignore paralog pairs below this % identity |
| `--gconv_min_paralog_length` | `200` | Ignore paralog pairs shorter than this |
| `--gconv_min_af` | `0.7` | Donor-allele fraction for a site to join a tract |
| `--gconv_min_sites` | `3` | Diagnostic sites needed before a run is called a tract |
| `--gconv_min_depth` | `5` | Informative depth needed at a site |
| `--gconv_min_bq` | `13` | Base-quality floor when reading an allele off a read |

Lower `--gconv_min_af` to catch a conversion present in only part of the population, at the cost of
more `ambiguous` calls. Raise `--gconv_min_sites` on a reference with many close paralogs.

## What it does not do

Breakpoints are located **to diagnostic-site resolution**, not to the base. A tract is bounded
between the last site carrying the acceptor allele and the first carrying the donor's; the actual
crossover is somewhere in between. Pinning it down needs local reassembly of the region, which
this does not do.

It also reports each paralog pair independently. In a family where three or more copies are
mutually similar, the same tract can appear against more than one donor, and choosing between them
needs evidence this analysis does not gather.
