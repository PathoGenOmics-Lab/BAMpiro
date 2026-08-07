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

**3. The model.** There are three ways an acceptor site can show the donor's base, and only one
of them is a conversion:

- it was converted, along with its neighbours, as one tract;
- it mutated to that base on its own, which for a single site is entirely ordinary;
- the read carrying it came from the donor, which in a paralogous region is what aligners do.

All three are in the comparison. A tract is an interval of diagnostic sites carrying donor bases;
the fraction of reads that arrived from the donor is a property of the locus, so it is fitted and
marginalised out rather than thresholded on. Every possible interval is evaluated exactly, on
every value of that fraction, and what comes out is a **log10 Bayes factor** for a conversion over
the best of the other two explanations, plus a posterior distribution over the breakpoints.

Three things follow from writing it this way, none of which the earlier allele-fraction version
could do:

- **Depth is evidence.** 7 donor reads out of 10 and 700 out of 1000 both give a fraction of 0.7,
  and they are not remotely the same evidence.
- **Base quality is evidence.** A Q30 base counts for more than a Q14 one instead of both merely
  clearing a floor, and a base under the floor is worth a little rather than nothing.
- **The molecule is the unit, not the site.** At 100 bp reads and diagnostic sites 30 bp apart it
  is the same fragments voting at both, so counting sites counts the same evidence twice. The
  likelihood is a product over reads, which is also what makes a read spanning a breakpoint decide
  the call on its own: it carries donor bases on one side and acceptor bases on the other, and a
  read that came from the donor carries donor bases everywhere it reaches.

**4. What sets the size a tract has to reach.** For the same set of sites, "converted" and
"mutated independently" have identical likelihoods. What separates them is their priors: one
conversion event against *k* independent substitutions of probability `--gconv_mut_rate` each. At
one site the substitution wins, at two it is close, and by three the conversion wins by orders of
magnitude. The old hand-picked `min_sites = 3` is where that lands, except now it is derived, it
moves with the spacing of the diagnostic sites, and it can be argued with by changing a rate.

That rate is **measured at the locus, not taken from the genome**. PE/PPE genes are
hypervariable, and they are hypervariable in exactly the positions where the copies already
differ, so a gene carrying a scattering of donor-matching bases on its own is a gene where a run
of two or three of them is unremarkable. Every diagnostic site outside the candidate tract is a
direct observation of that rate, and `--gconv_mut_rate` becomes a floor rather than the answer.
Only isolated substituted sites count: a run of them is what a conversion looks like, so counting
runs would let a second genuine tract talk the first one down. The rate used is reported per row
as `mut_rate`.

**5. How much of the read pool carries it.** A clonal conversion in a two-copy family shows up in
every read. Two ordinary things pull that down, and from allele data they are the same thing: a
mixed infection, and a gene family with more than two members whose unconverted copies contribute
reads to the acceptor. The second is not an edge case, and this is worth stating plainly:

!!! warning "21% of H37Rv's paralog pairs have a relative closer than their own donor"

    Reads from that third copy land on the acceptor carrying unconverted bases, so a perfectly
    clonal conversion can show an allele fraction near 0.5 and nothing is wrong. The fraction is
    fitted and reported as `tract_af` rather than counted against the tract, with a floor at 0.25:
    below about a fifth of the reads there is nothing left to separate a minority conversion from
    index hopping or a contaminating sample.

## Verdicts

| Verdict | Meaning |
| :--- | :--- |
| `gene_conversion` | The Bayes factor clears `--gconv_min_bf` against both alternatives |
| `mismapping` | The locus is explained by a fitted fraction of reads arriving from the donor, with nothing left for a tract to account for |
| `ambiguous` | Reported, but not called. The `reason` column says what came closest |
| `coverage_shift` | The acceptor lost its reads to the donor over a run of sites. Consistent with a conversion longer than the library insert, and equally with a deletion. See below |

A tract covering **every** diagnostic site of the locus is a special case that needs no special
handling: "the whole locus was converted" and "every read here came from the donor" predict
exactly the same bases, so their likelihoods are equal and the Bayes factor collapses to the ratio
of their priors on its own. The model reports the ambiguity rather than having to be told about it.

## Reading the numbers

| Column | What it answers |
| :--- | :--- |
| `log10_bf` | Is this a conversion, rather than a coincidence of substitutions or a pile of donor reads? Above 3 is decisive |
| `log10_bf_vs_null` | Are these sites really carrying the donor's bases? This is where depth, base quality and read linkage show up |
| `post_conv` | The same as `log10_bf`, read through `--gconv_prior` |
| `tract_af` | What fraction of the reads carry the tract. Below 1 means a mixed infection or a third copy of the family; nothing in short reads tells those apart |
| `mut_rate` | The per-site substitution rate this locus turned out to have, which is what the tract had to beat |
| `mismap_frac` | What fraction of the reads here the model had to assume came from the donor |
| `start_ci`, `end_ci` | Where the breakpoints are, to diagnostic-site resolution, as a 95% credible interval |

`log10_bf` is the weaker of two things, so it is the one to threshold on: thin data leave the
sites themselves in doubt, and deep data settle the sites but cannot settle whether a run of donor
bases is a conversion or a coincidence. A tract with a huge `log10_bf_vs_null` and a small
`log10_bf` is one the reads are certain about and the biology is not.

The descriptive columns are still there and still worth reading, because they are what you can go
and check in the BAM: `donor_af_in` (how fixed the donor allele is inside the tract),
`donor_af_outside` (whether it also turns up outside), `breakpoint_reads` (molecules carrying both
in cis) and `n_undetermined` (sites inside the tract with too little depth to genotype). None of
them decides anything any more.

## Output

Per sample, next to the other per-sample files:

```
<sample>.<ref>.gene_conversion.tsv
```

and one cohort-level `<samplesheet>_gene_conversion.tsv`. Columns: `sample`, `pair_id`, `contig`,
`donor`, `verdict`, `reason`, `start`, `end`, `span_bp`, `post_conv`, `log10_bf`,
`log10_bf_vs_null`, `mismap_frac`, `start_ci`, `end_ci`, `n_sites`, `n_sites_outside`,
`n_undetermined`, `donor_af_in`, `donor_af_outside`, `min_depth`, `cis_reads`,
`breakpoint_reads`, `donor_only_reads`.

The paralog map itself is published under `references/<refId>/` as
`<refId>.paralog_pairs.tsv` and `<refId>.paralog_sites.tsv`, and is worth a look on its own: it is
the donor/acceptor graph of your reference.

## Tuning

| Parameter | Default | What it controls |
| :--- | :--- | :--- |
| `--gconv_min_identity` | `90.0` | Ignore paralog pairs below this % identity |
| `--gconv_min_paralog_length` | `200` | Ignore paralog pairs shorter than this |
| `--gconv_min_bf` | `3.0` | log10 Bayes factor before a tract is called a conversion |
| `--gconv_report_bf` | `1.0` | log10 Bayes factor below which a tract is not written out |
| `--gconv_prior` | `0.01` | Prior that a given paralog pair carries a tract |
| `--gconv_mut_rate` | `0.0003` | FLOOR on the chance a site carries the donor base by plain substitution. The locus raises it when it turns out to be hypervariable |
| `--gconv_min_tract_af` | `0.25` | Read fraction a tract needs before it is called |
| `--gconv_mean_tract_bp` | `1000` | Mean of the exponential prior on tract length |
| `--gconv_max_tract_bp` | `10000` | Longest tract considered |
| `--gconv_max_tracts` | `2` | Conversion tracts looked for per paralog pair |
| `--gconv_min_mismap` | `0.2` | Fitted donor-read fraction reported as mismapping |
| `--gconv_min_sites` | `3` | Diagnostic sites a pair needs before it is analysed at all |
| `--gconv_min_depth` | `5` | Depth below which a site is undetermined in the summaries |
| `--gconv_min_bq` | `13` | Base-quality floor when reading an allele off a read |

`--gconv_min_bf` is the one to reach for. 3 is "decisive" on the usual scale; drop it to 2 to see
more candidates, and read `log10_bf_vs_null` alongside to see which kind of doubt is behind each.

`--gconv_mut_rate` is a floor on the per-site chance of an independent substitution to exactly the
donor's base, roughly the genome's SNP rate against its reference divided by three. Each locus
raises it if its own sites say so, so on a divergent isolate the adjustment happens on its own;
raise the floor if you want every locus treated as hypervariable. `--gconv_mean_tract_bp` only breaks ties between overlapping intervals of
similar likelihood; it is not a filter, and a tract far longer than it will still be called if the
reads say so.

!!! note "Performance"

    The search is exhaustive and costs O(sites² x reads) per paralog pair. Above 150 diagnostic
    sites the breakpoints are placed on a coarser grid, which costs breakpoint resolution and
    nothing else; the `reason` column says so when it happens. Short intervals are always kept at
    full resolution, since those are the ones a coarse grid would mangle.

## What it was measured on

Not a toy. The model was tuned and checked against the real H37Rv genome and its **411 real
paralog pairs**, with isolates simulated at 12x and 30x carrying position-dependent quality decay,
errors the Phred score does not admit to, indels, PCR duplicates removed by `samtools markdup`,
uneven coverage with dropouts, 1200 background substitutions and 3% reads from a near neighbour,
then mapped and deduplicated exactly as the pipeline does.

| | Result |
| :--- | :--- |
| Conversion tracts implanted across 23 real paralog pairs, at 2 to 20 diagnostic sites and 25% to 100% frequency | **21 of 23 found** at the default threshold, at 30x and at 12x alike |
| Breakpoints on the clonal tracts | exact, with the credible interval covering the truth |
| False positives over the 388 pairs with nothing implanted | **0**, at every threshold down to `log10_bf` 1 |
| A negative control isolate: 411 real pairs, no conversion anywhere | **no conversion call at all** |
| 22 isolated substitutions that happen to match the donor's base | **0 called** |
| Hypervariable loci, 10% of diagnostic sites substituted | **3% called**, against 39% with a fixed genome-average substitution rate, with the same 10 of 10 real tracts found either way |

The two tracts that were missed sat at 25% and 40% frequency, which is what the floor under
`tract_af` is for.

## The blind spot: a tract longer than the insert

A conversion makes the acceptor identical to the donor over the tract. Once that tract is **longer
than the library insert**, a read pair falling entirely inside it has no unique anchor left, and
the aligner assigns it to one copy arbitrarily. The acceptor's tract loses coverage and the
donor's matching region gains the same reads.

Measured on a simulated 99.7% paralog pair with a 1.4 kb tract and a 320 bp insert:

| Region | Converted sample | Unconverted control |
| :--- | ---: | ---: |
| Acceptor, inside the tract | **36.7x** | 83.8x |
| Donor, matching region | **127.4x** | 83.3x |
| Acceptor, outside the tract | 80.1x | 81.1x |

The reads did not disappear, they moved: the two loci together carry 164x against 167x in the
control. Every allele-based signal this tool depends on evaporates exactly when the conversion is
most complete, which is why the same detector finds a 350 bp tract at 98% identity without
difficulty and misses a 1.4 kb one at 99.7%.

That failure would otherwise be silent, so a run of diagnostic sites where the acceptor is
depleted while the donor carries the reads is reported with the verdict `coverage_shift`. **It is
not a conversion call.** A deletion of the acceptor produces exactly the same picture, and telling
them apart needs evidence short reads do not carry: longer reads spanning the whole tract, or a
read-depth analysis over the paralog pair as a whole.

## What it does not do

Breakpoints are located **to diagnostic-site resolution**, not to the base. A tract is bounded
between the last site carrying the acceptor allele and the first carrying the donor's; the actual
crossover is somewhere in between. Pinning it down needs local reassembly of the region, which
this does not do.

It also reports each paralog pair independently. In a family where three or more copies are
mutually similar, the same tract can appear against more than one donor, and choosing between them
needs evidence this analysis does not gather.
