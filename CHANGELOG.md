# Changelog

All notable changes to this project are documented in this file. The format is
based on [Keep a Changelog](https://keepachangelog.com/), and the project follows
[Semantic Versioning](https://semver.org/).

## [Unreleased] - targeting 1.1.0

### Added

- **Two guards against a sample mapped to a genome it does not belong to.** In a 185-sample
  cohort, 22 samples annotated as one lineage were really another, were routed to that lineage's
  reference, carried about 2,000 SNPs where their line mates carried 5, and produced 507 of the
  run's 522 gene-conversion tracts. Tracts per sample against genome-wide SNP count correlated at
  r = 0.98. Nothing failed: the reference was a valid genome, the reads mapped to it, and every
  metric was computed correctly about the wrong comparison.

  `LINEAGE_MISMATCH` flags a sample whose k-mer lineage call disagrees with the rest of the
  samples sharing its reference. K-mer typing does not use the reference, so the two are
  independent evidence, and the majority stands in for what should have been there without
  needing to know the reference's own lineage.

  `divergent_sample` is a cohort verdict that refuses to read a sample's tracts as conversion when
  it calls them in too many loci at once, since conversion is local and a diverged genotype
  carries the donor base at every paralogous locus by inheritance. The threshold is
  `--max-locus-frac`, default 0.1: on the cohort that surfaced this the mislabelled samples called
  tracts in 16 to 29 percent of the stretches their reference reports anything in, and every
  correctly mapped sample in at most 1 percent. It is counted in distinct stretches of the
  acceptor rather than in rows or pairs, because one tract is emitted once per candidate donor
  and a gene family reports one stretch through a pair per relative, and it takes precedence over
  the corroboration rule, which diverged samples would otherwise satisfy by sharing their
  artefacts with each other at identical coordinates.

### Fixed

- **Gene conversion was called from reads that could not carry it, and the cohort multiplied
  it.** On a 185-sample cohort, 1,520 of the 1,579 calls were rows the samples themselves had
  left `ambiguous`, promoted because another sample had called the same stretch outright. Those
  rows carried the donor's bases in 2 to 5% of the reads; the model's tract fraction stops at
  20%, so it fitted 20% to all of them and reported Bayes factors for a fraction the reads did not
  have. The outright calls vouching for them came from samples at 0.2x to 2x depth, contaminated
  cultures the QC fails outright, whose tracts were read by one to four molecules: one read with
  the donor's bases at four sites is log10 BF 7 on base quality alone. Reads from a third copy of
  a gene family reproduce the same stretch in every library mapped to the same reference, which is
  exactly what the corroboration rule took for independent confirmation.

  A tract most of whose sites are under `--gconv_min_depth`, or whose reads carry the donor's
  bases below half the model's thinnest tract (10%), is now `ambiguous` with a reason that quotes
  the reads rather than the floor, and neither kind takes part in corroboration. The cohort pass
  applies the same check, so re-running only that step corrects a run. On that cohort the calls
  drop from 38 events to 4, two of them in the mixed cultures the QC also fails.

- **The cohort's recurrence rule was measured against the wrong cohort.** It divided by every
  sample in the run, so on a cohort mapped against two references an artefact of one could never
  reach the 90% that marks a `reference_artifact`: 113 of 185 samples is 61%. It now divides by
  the samples mapped to the same reference, leaving out the samples that are diverged from it.
  The pooled locus statistics were keyed on the pair number alone, which pooled pair 61 of one
  reference with pair 61 of the other; they are keyed on the contig as well. The reason on a
  corroborated row quoted the Bayes factor as short of the threshold on rows at log10 BF 5.6 and a
  threshold of 3; it now says what the row fell short of.

- **Kraken counted runs and read species, so a clean cohort looked contaminated.** One report
  per run became one row per run, so 185 samples gave 225 rows, and a sample's "primary taxon"
  was its top species: Kraken2 leaves most *M. tuberculosis* reads on the complex, so the
  species held about 5% and every sample fell under the 90% line. The report summarised 225 of
  185 samples as possibly contaminated, while the cultures that were really another organism sat
  unnoticed among them. A sample's reports are now merged by their reads, each sample is placed
  in the clade its reads actually sit in, and contamination is measured as the share of the
  classified reads in the clade most samples are dominated by: below 90% a sample is mixed,
  below 50% it is another organism.

- **The lineage read-out called every typed sample mixed.** It counted the samples carrying a
  lineage breakdown, which every typed sample does, and reported 176 mixed-lineage samples
  above a table that correctly counted none. It now counts the `MIXED` flag.

- **The report dropped `LINEAGE_MISMATCH` as soon as it loaded.** Loading the page, and every
  threshold change, rebuilt the flags from the thresholds alone, so a flag decided from evidence
  the page does not hold vanished: the report never showed a wrong-reference sample and its
  verdicts could disagree with `qc_flags.tsv`. Flags the page cannot recompute are now kept, and
  the report says what was compared (*types as A4; the rest of E1ASM0035 types as L7*).

- **A sample's collection date was the day the pipeline ran.** The summary's `date` is
  `DAT_OUT`, stamped at processing, and the temporal panel read it as a collection date: a whole
  cohort "sampled in 2026, a span of 0 years". The date now comes only from a samplesheet column
  that names one (`collection_date`, `sampling_date`, `date`, `year`...).

- **The reference-bias screen accused samples mapped to their own ancestor.** Its cut is
  relative to the cohort median, and where that median is a handful of SNPs, 2 SNPs against 5
  made 18 passing samples "reference-bias suspects". Below 10 SNPs per callable Mb (50 SNPs when
  there is no density) the screen is off, and the panel says so and points at the samples far
  above the rest instead.

- **The resistance read-out counted lineage markers as resistance.** pncA H57D is in every
  *M. bovis*, so every A4 sample of a cohort was "resistant to PZA" and the mutations that
  arose during the experiment were a minority of what the panel reported. A grade 1-2 mutation
  carried by at least 90% of a lineage's samples is now listed as that lineage's marker, apart
  from the rest.

- **Smaller read-out and drawing faults.** The correlation read-outs headlined pairs that are one
  quantity measured twice (*Mapped %* against *Unmapped %* at r = -1); the gene-conversion panel
  had no verdict for `divergent_sample`, the most common one on the cohort that introduced it, so
  its 2,430 rows read "verdict not recognised"; the metric scatter drew its top values above the
  axis; the genome landscape carried a fixed purple key over a red heatmap; and a printed report
  came out under the grey backdrop of the mobile sidebar.

- **A GFF attribute spelled `nan` was read as a gene name.** A GFF built from a table writes
  pandas' NaN as four ordinary characters, and the attribute cascade in `gconv_annotate` only
  asked whether a key was PRESENT. On a real MTBC reference 3,001 CDS lines carry `Name=nan`
  beside a perfectly good `ID=Rv0001_1-1524`, and none of them carry `locus_tag` at all, so
  `Name` won: every tract in a 185-sample cohort came out labelled `nan`, with amino-acid changes
  reported as `nan:P1443A`. Nothing downstream could tell it from a gene really called that.

  `_attr` now treats a placeholder value as an absent one, so the cascade falls through to the
  ID. `nan`, `NaN`, `NA`, `None`, `null`, `n/a`, `.`, `-` and the empty string are all compared
  case-insensitively, and a real symbol such as `dnaN` still wins over the ID.

### Changed

- **The QC report reads as seven pages, starting from what the run found.** On a 185-sample
  cohort it was one scroll of 27 panels about 150,000 pixels long: the SNP dynamics alone drew
  800 trajectory cards, the resistance panel listed all 4,986 calls, the table of all samples
  had 48 columns with a bar in every cell, and the flags were codes. It now opens on a
  *Summary* that says what the run shows in sentences with their numbers (the samples to
  exclude and why, those that are another organism or were mapped to the wrong reference, the
  resistance mutations beyond each lineage's own, the alleles that swept, the gene conversion
  called), each linking to the page that holds the evidence and saying what it does not prove.
  The rest is split into *Sample QC*, *Genome & genes*, *Variants over time*, *Resistance*,
  *Gene conversion* and *Diagnostics*, with a badge in the sidebar where a page needs attention.

  A flag reads as what went wrong and the value against its rule (*Low depth 2.2x, needs 10x*),
  and the codes stay in every export. The table of all samples opens on the 14 metrics a
  verdict is made of and colours only the cells that tripped a flag. The resistance panel groups
  the calls by mutation and opens on grades 1-2; the contamination panel is one row per sample,
  flagged ones first. The big panels draw a first screenful and keep the rest one click away:
  24 trajectory cards at a time, 40 gene chips, 12 parallel genes. Printing lays out every page.

- **Garnatxa submits four times faster.** `submitRateLimit` goes from 50 to 200 jobs a minute,
  because it binds whenever the tasks are shorter than the interval between submissions. A
  185-sample cohort instantiates about 4,300 tasks, 740 of them `ANNOTATE_LEGACY_VCF` at two or
  three seconds each, so at the old rate those alone spent a quarter of an hour doing nothing but
  calling sbatch. It is a courtesy limit toward the scheduler rather than a correctness one, and
  the comment says to lower it again if submissions start being refused.

### Fixed

- **The dedup pipe asked for 160 percent of its own memory, and the resulting OOM was invisible.**
  `samtools sort -m` is per thread and `MERGE_AND_MARKDUP` keeps two sorts alive at once, a name
  sort feeding fixmate and a coordinate sort feeding markdup. Each was sized at 80 percent of the
  task's allocation, so together they asked for 13.1 GB of an 8 GB request. A single sequencing
  run produces a BAM too small to fill those buffers, which is why it stayed hidden until the
  first merged sample, whose four input BAMs gave a 2.8 GB stream.

  The second half was worse. The kernel kills one stage of a pipe, its stdout closes mid-BGZF
  block, and the next stage reports a short read and exits 1, so an out-of-memory is
  indistinguishable from a corrupt file and `errorStrategy` retried nothing. The pipe's
  PIPESTATUS is now inspected and a signal death re-raised as 137, which the existing policy
  already retries with double the memory. A genuine non-zero exit still fails fast.

### Changed

- **Kraken's memory is a parameter and scales with the attempt.** It was the one fixed `memory`
  directive left in the pipeline, at 56 GB, and on a 185-sample cohort that put 61 of these jobs
  behind SLURM's `QOSMaxMemoryPerUser` while the cluster itself had room. The right number is a
  property of the database, not of the pipeline: with `--memory-mapping` the resident set is the
  part of the index actually probed, measured at 34.7-39.7 GB against a 133 GB standard database
  and well under 20 GB against the capped 16 GB one. `--kraken_memory` now sets it, defaulting to
  24 GB, and it doubles on each retry like every other memory directive here, so an underestimate
  costs one retry rather than the run.


- **Every task on Garnatxa goes to the same QoS by default, and the queue is a parameter.** Two
  of the three steps that asked for a longer QoS were measured at about twenty seconds each, so
  they never needed it; the third only did because Kraken was reading a 133 GB database. `--qos`
  and `--qos_heavy` now set the queue for the ordinary and the long steps, both defaulting to
  `short`, and the long steps' walltime comes down from two hours to one. Declared in the root
  config rather than in the profile, because `main.nf` validates command-line parameters against
  that block and one declared only in a profile is rejected.


- **Kraken filtering compresses with threads and asks for what it uses.** The step ended with two
  sequential `gzip` calls, which are single-threaded: on a 1.5 million pair sample that is 72
  seconds spent compressing 480 MB per mate while eleven of the twelve reserved cores sit idle.
  `bgzip -@` does the same work in 2.8 seconds. BGZF is valid gzip, the decompressed bytes are
  identical and the file is slightly smaller, and fastp reads it whole, which is the check that
  matters because BGZF is a multi-member gzip and a reader that stops at the first member would
  truncate silently. All of it measured against the pinned container.

  The reservation is sized from the same measurements: three tasks averaged 102 percent CPU over
  ten minutes on a 12-CPU request and peaked at 39.7 GB against an 80 GB one, so they now ask for
  8 CPUs and 56 GB. A smaller reservation is scheduled sooner, which is most of what a run waits
  for on a shared queue.

- **Tasks on Garnatxa declare a walltime.** No process set `time`, so every job inherited the QoS
  default of six hours, and SLURM's backfill scheduler can only place a job in a gap longer than
  the walltime it asked for. A 200 ms collection step therefore queued behind everything rather
  than filling the next gap. The default is now 30 minutes and the long steps get two hours, both
  scaled by `task.attempt` so a timeout is retried with more rather than failing the run.

### Fixed

- **Lineage and drug-resistance typing ran on a fraction of each library.** FastP is run with
  `--merge`, so a pair whose mates overlap leaves as a single merged fragment in the orphan
  stream rather than as r1/r2. `MAP_READS` has taken all three read streams from the start;
  `LINEAGE_TYPING` took two, so pathotypr typed only whatever happened to stay unmerged. How
  much that is depends on the insert size, not on the pipeline: 8-19 percent in the cohort that
  surfaced it, and on a purpose-built fully overlapping library r1 and r2 reached the typer
  holding zero reads each while the merged stream held 3,000.

  Nothing failed and no count read zero anywhere: a sample typed on nothing reports
  `Unclassified`, which is also what a genuinely unclassifiable sample reports. Any lineage
  call or resistance panel from an earlier run was computed on a subset of the evidence and is
  worth regenerating.

  The three streams are now concatenated and typed as one sample. They cannot be passed as
  three `-i` arguments, because `--paired` requires an even number of files and without it the
  names collapse to one duplicate sample. Both risks of concatenating were measured against the
  pinned binary rather than assumed: pathotypr reads multi-member gzip (r1 alone falls below
  `--min-depth` and yields no call, while the concatenation yields the same counts as
  `--paired`, which it could not if it read only the first member), and `--paired` only groups
  files, so the same pair passed as one concatenated file gives identical ref and alt counts.

### Changed

- **Drug-resistance catalogue upgraded to Zenodo v1.0.2, which corrects v1.0.0.** The WHO
  catalogue grades every variant-drug pair separately; v1.0.0 instead gave each variant a
  single drug inherited from its gene and took the grade from whichever catalogue row it
  found first. Amikacin therefore did not occur anywhere in v1.0.0 and could never be
  reported, the *rrs* aminoglycoside determinants were attributed to streptomycin, and the
  *inhA* promoter variants to isoniazid alone. Beyond the 89 relabelled rows, 15,969 rows
  carried the wrong grade: the upgrade moves 9,439 into grades 1-2 and none out of them,
  growing the reportable marker set by 59 percent.

  Detection is unchanged, and so are MDR and pre-XDR assignment. Results reporting
  amikacin, kanamycin, capreomycin, ethionamide, linezolid, streptomycin or delamanid from
  an image built before this change should be regenerated. An image ships v1.0.0 if
  `pipeline_info/software_versions.txt` records the `dr_markers` SHA-256
  `9769864774a11ed325d05b0543932782f6cd8c6ee444998c616e44c5763576c2`.

  v1.0.2 publishes two coordinate frames; the bundled one is the MTBC-ancestor frame, which
  is the frame of the reference bundled beside it. All 102,216 of its REF alleles match that
  reference and the H37Rv file mismatches 609 of them at identical coordinates, so the wrong
  choice would have dropped those markers with no error. The lineage markers and the
  pre-trained model are byte-identical to v1.0.0.

### Added

- **Every run records which resistance catalogue produced its calls.** `DUMP_VERSIONS` now
  emits the pathotypr binary version, the catalogue version and the SHA-256 of the marker
  file actually read. Which catalogue a set of calls came from was not recoverable from the
  calls themselves, and two catalogues that disagree on 15,969 rows can produce reports that
  look identical.

- **A composite drug label now counts towards every drug it names.** Catalogue v1.0.2 writes
  a variant graded for several drugs as `AMI_KAN_CAP` or `INH_ETH`. The resistance matrix
  used to give such a label a column of its own, leaving the columns of the drugs it covers
  empty, so a sample whose only isoniazid evidence was an *inhA* promoter variant read as
  clean under INH. The calls table and the TSV still carry the label as the catalogue wrote it.

### Fixed

- **A detected variant with no WHO grade no longer looks like a clean drug.** Both an
  ungraded call and a genuine no-call rendered as an empty cell in the resistance matrix.
  v1.0.0 contained 6,056 ungraded rows, ten of which are grade 1-2 under v1.0.2. Ungraded
  calls now carry a mark and a legend entry of their own; only a no-call is blank.

- **The container build verifies its Zenodo assets by checksum.** The comment claimed it
  failed loudly on a truncated download, but the only guard was `test -s`, which a partial
  9 MB file passes. `sha256sum --strict` is used rather than plain `-c`, because without
  `--strict` a malformed digest line is skipped with a warning and the build still succeeds,
  which was confirmed by building with one.

- **The pathotypr binary is pinned alongside its data.** It was installed unpinned while its
  marker files were pinned to a Zenodo record, so a rebuild paired whichever binary bioconda
  served that day with a fixed catalogue. The image published before this change carries
  binary 1.0.2 against catalogue v1.0.0, which is exactly the skew v1.0.2 renumbered the data
  files to prevent.

### Added

- **Whether the donor kept its own bases.** Gene conversion is non-reciprocal by
  definition, and until now that was only asserted in the definition. The donor's
  reads are read as alleles over the same sites, and `donor_swap_af` says how much
  of the donor carries the ACCEPTOR's bases. Past `--gconv_reciprocal_af` both
  copies changed, which is an unequal crossover: the verdict is
  `reciprocal_exchange`. The swap has to be BOUNDED, high over the tract and low
  outside it: reads from the unconverted acceptor that the aligner placed at the
  donor carry the acceptor's bases everywhere, and reading the tract's sites
  alone called a donor that had changed nothing an exchange.
- **How many copies are contributing reads at a locus.** `locus_cn` and
  `expected_af` say how many times the genome's depth a locus runs at and what a
  clonal conversion of one copy would reach there, so a fraction that reads as a
  minority event can be recognised as one whole copy of several. It explains a
  diluted fraction and does not discriminate one, so no verdict moves on it.
- **Which copy the conversion is on** (`--gconv_outgroup`, off by default). A sample
  whose acceptor carries the donor's base has either changed or not, and the site
  alone cannot tell which: where the REFERENCE's acceptor carries a derived allele,
  a sample carrying the donor's base is holding the ancestral state and has
  converted nothing. An outgroup aligned against the reference labels every
  diagnostic site, and a tract whose sites mostly say ancestral is reported as
  `reference_derived` instead of as a conversion.
- **What a tract does to the genes it lands on.** With the reference's GFF, each
  tract reports the genes it covers, how many of its copied bases are synonymous and
  how many are not, and the amino-acid changes as `Rv0001:K2Q`. Deletion markers are
  left out, being a frameshift question rather than a codon one.
- **Gene conversion detection** (`--find_gene_conversion`, off by default). Finds
  tracts where one paralog has been copied onto another, which on a reference that
  never saw the event reads as a run of variants that are exactly the donor's
  bases between two breakpoints.
 - The donor/acceptor map is recovered from the reference self-alignment the
    repeat-masking step already computes and then discards, so no alignment is
    redone; `show-snps` on the same file gives the positions where the two copies
    differ, which are the only ones that can carry evidence.
 - The stage deliberately bypasses the pipeline's own masking. Repeat exclusion
    and the length-aware read filter remove paralogous sequence and multi-mapping
    reads, which is right for variant calling and removes exactly this signal, so
    the analysis reads the deduplicated **pre-filter** alignment and does not
    filter on mapping quality.
 - The hard part is not finding tracts but deciding whether to believe them. There
    are three ways an acceptor site can show the donor's base and only one of them
    is a conversion: it was converted with its neighbours, it mutated to that base
    on its own, or the read carrying it came from the donor. All three are weighed
    against each other by an explicit model, which evaluates every possible tract
    exactly, weights each base by its quality, treats the molecule rather than the
    site as the unit of evidence, and fits the fraction of reads that arrived from
    the donor instead of testing an average against a threshold. It reports a log10
    Bayes factor and a posterior over the breakpoints.
 - Two components exist because a run over the REAL H37Rv genome, with its 411 real
    paralog pairs and an isolate simulated with position-dependent quality decay,
    indels, duplicates, uneven coverage and 3% contamination, showed the model
    getting the wrong answer for a reason no toy dataset would have shown. The
    fraction of reads that carry a tract is fitted rather than counted against it,
    because 21% of H37Rv's paralog pairs have a relative closer than their own donor
    and that relative's unconverted reads land on the acceptor. And the substitution
    rate a run of donor bases has to beat is measured at the locus rather than taken
    from the genome, because a hypervariable gene is one where a short run is
    unremarkable: 39% of such loci were called conversions with a fixed rate
    against 3% with the locus's own, and the same real tracts were found either way.
 - Measured on that benchmark: 21 of 23 implanted tracts found at 30x and at 12x
    alike, breakpoints exact on the clonal ones, and no false positive over the 388
    pairs with nothing implanted. A negative-control isolate produced no conversion
    call at all.
 - Two things that used to need special cases now fall out of the arithmetic. A
    tract covering the whole locus predicts the same bases as every read having
    come from the donor, so its Bayes factor collapses to the ratio of the priors
    on its own. And the number of sites a tract needs is set by how implausible
    that many independent substitutions would be, so it is derived from a rate
    rather than chosen, and it moves with the spacing of the diagnostic sites.
 - The tracts whose reads have moved to the donor are reported as `coverage_shift`.
    Once a tract is longer than the library insert a read pair falling inside it has
    no unique anchor left, so the acceptor loses its coverage to the donor and every
    allele-based signal evaporates exactly when the conversion is most complete. It
    is not a conversion call: a deletion of the acceptor looks the same.
 - The paralog map attaches each diagnostic site to EVERY pair that spans it. A gene
    family produces overlapping and nested alignments on purpose, so stopping at the
    first one cost 140 sites on H37Rv and left 12 pairs with none at all, among them
    a 5.8 kb paralog at 98% identity that the caller then skipped for having nothing
    to work with. Where two nested alignments compete for one acceptor position, the
    offset between the copies says which pair it belongs to.
 - `coverage_shift` requires the reads to have MOVED, not merely to be scarce: the
    acceptor below its own level elsewhere in the locus AND the donor above its own.
    The absolute depth floor it used before fired on any low-coverage paralog, which
    is how a clean negative control reported two shifts over a locus running at 4x
    throughout.
 - **Deletions between the copies count as markers.** A stretch the donor does not
    have still has a coordinate in the acceptor, so a read either carries a base
    there or spans it with a deletion. A run of missing bases is one marker rather
    than one per base, because it happened once. That is 423 more markers on H37Rv,
    11% on top of the substitutions, and each is worth more than a substitution: the
    headline Bayes factor is capped by how implausible it is that the markers arose
    independently, and two copies losing the same bases is far longer odds than two
    copies mutating to the same base. Measured against conversions that copy the
    donor's deletions as a real one would, 15 of 15 implanted tracts were found
    against 14, and the one that changed sides is a 99.4% pair with seven diagnostic
    sites in the whole locus.
 - **A gene family's several views of one event collapse into one, with a source
    named where the reads can name it.** Those rows are two different things wearing
    the same shape, and the donor coordinates tell them apart: naming the same
    stretch of donor is redundancy, naming different places is a real choice of
    source. It is settled by evidence per marker, which on the real genome separates
    a true source from a bystander by 19.1 against 1.4, and honestly fails to
    separate four relatives sitting within 0.17 of each other. Across the benchmark:
    14 sources named correctly, 1 reported ambiguous, none named wrongly, and 27
    rows collapsing to 19 events.
 - **The cohort is read together, not one sample at a time.** Two questions cannot
    be answered from a single sample however good the model is. A tract in one
    sample of fifty is a finding; the same tract at the same coordinates in all
    fifty is the reference being wrong there, or the aligner doing it to everybody,
    and those are reported as `reference_artifact`. And a signal too weak for one
    sample to commit to is a different proposition when the same tract is called
    outright in another, because contamination and index hopping do not reproduce a
    specific tract across independent libraries. Validated on the real genome: six
    isolates sharing one tract had all twelve of its rows demoted while their
    isolate-specific rows were left alone, and on an eight-isolate cohort three
    sub-threshold tracts were corroborated, all three of them real.
 - **The checks that found the defects are in the suite**, not run once and thrown
    away. Every defect in this stage came from a randomised run: a prior of 0 taking
    log(0), positions out of order answered confidently and wrongly, a tract reported
    outside its own credible interval, a CIGAR offset. Each became a fixed-case
    regression test while the generator that found it was discarded, so the next
    defect of the same shape waited for someone to go looking. The generators now run
    on every change, and each is verified to still catch the defect it was built
    from. None of them needs real data: real data would test whether the model
    describes biology, these test whether the code does what the model says.
 - **Every output records the settings that produced it**, as `#` header lines. A
    results table whose verdicts depend on seventeen settings and does not say what
    they were cannot be checked against another run or reproduced a year later. The
    cohort file carries the per-sample lines too, which is what shows a cohort
    assembled from samples run differently.
 - A **Gene conversion** panel in the QC report leads with the Bayes factor and puts
    the observable evidence beside it, one row per event by default, with how many
    samples of the cohort carry each and the settings that produced the verdicts.
    See [gene conversion](docs/gene-conversion.md), and the output reference for
    which of the four files answers which question.

Headline: a **test suite and CI**, and a pipeline that runs correctly on a machine
that is not the authors' cluster.

`manifest.version` already reads `1.1.0`: that is the number this work is heading
for, not a release that happened. Nothing is tagged and no release is published.
See [RELEASING.md](RELEASING.md) for the order the actual release has to follow.

The **container stays at the 1.0.1 image**, pinned by its digest. Nothing here
changes a bundled tool, and the scripts added under `bin/` are staged by Nextflow
from the pipeline directory rather than baked into the image, so there is nothing
to rebuild. A release only needs a new image when the Dockerfile changes.

### Added

- **Test suite** (`tests/`) - around 1,100 tests in three legs, all runnable with
  `tests/run_tests.sh` and none of them needing a container, a reference genome or
  a network connection:
 - `tests/unit/` - the Python under `bin/`: the consensus decision tree, the QC
    verdict engine, every parser, the k-mer liftover and the read filter.
 - `tests/js/` - the report's hand-written ES5 statistics (Spearman,
    Kruskal-Wallis, chi-square, Benjamini-Hochberg), checked against **SciPy and
    statsmodels** reference values rather than a snapshot of our own output, plus
    the integrity of the asset bundle.
 - `tests/pipeline/` - samplesheet validation and a full `-stub-run` of the DAG.
- **A `stub:` block for every process**, so `-stub-run` walks the whole workflow in
  seconds with no data. A test fails if a new process arrives without one.
- **Fixture cohort** (`tests/data/`, 170 kB) - a 5 kb reference with three genes and
  three samples covering paired-end, single-end and a sample split over two runs,
  all derived from one seed. CI regenerates it and fails on any diff.
- **`test` and `test_full` profiles** - the default feature set, and every optional
  branch on so a single run instantiates every process.
- **Continuous integration** (`.github/workflows/ci.yml`) - lint, unit tests on two
  Python versions, the front-end tests, and a stub run against both the minimum
  supported Nextflow and the current release, on every pull request.
- **`slurm` and `garnatxa` execution profiles**. `-profile garnatxa` replaces the
  hand-written `-c cluster.config` a user previously needed, and works when the
  pipeline is pulled straight from GitHub.
- **`--help`**, listing every parameter with its default.
- **`--max_cpus` / `--max_memory` / `--max_time`** - cap every process to what the
  machine can give, so an over-sized request (Kraken2 asks for 80 GB) fits a laptop
  without editing the modules.
- **Community files** - CONTRIBUTING, a Code of Conduct, a security policy, issue
  forms, a pull-request template, `RELEASING.md` and `.zenodo.json`.

### Fixed

- **`--flag false` now means false on every supported Nextflow.** Up to Nextflow 24
  a command-line parameter was coerced to the type of its config default; from
  Nextflow 26 it arrives as the string `"false"`, which is truthy in Groovy. All 15
  boolean flags were affected, so a run asked to skip the QC report, the consensus
  or the SNP matrix produced them anyway, and the change arrived with a Nextflow
  upgrade rather than with anything in this repository. Every Groovy-side test now
  goes through `asBool()`, and the stub-run tests check both directions.
- **`--find_gene_conversion` without `--exclude_repeats` stops the run.** The paralog
  map comes from the self-alignment repeat masking computes, so with masking off the
  stage produced a header-only file: an empty result that reads like "no conversion
  anywhere" and means "this never ran".

### Changed

- **`-profile standard` no longer submits to SLURM.** It is now the portable default
  and runs on the current host, so a fresh clone works anywhere. Use `-profile slurm`
  for a generic cluster or `-profile garnatxa` for the I2SysBio one. **This changes
  what an existing `-profile standard` command does.**
- The base config no longer carries one site's paths: the SLURM QoS, the
  `module load singularity`, the Singularity `--bind` paths and the Kraken database
  location all moved into `conf/garnatxa.config`.
- **`--tsv` is required** and says what to pass; the default was a leftover
  samplesheet name. `--kraken2_db` has no default either, and without one the
  contamination screen is skipped.
- **The container is pinned by digest** rather than the mutable `1.0.1` tag, so two
  runs a month apart cannot silently use different tool versions.
- **An unrecognised `--parameter` now stops the run** with a suggestion, instead of
  being accepted and ignored. The declared set is read from `nextflow.config`, so
  the two cannot drift.
- The startup banner reports the resolved profile, executor, container and Kraken
  state, because where a run lands is the easiest thing to get wrong.

### Fixed

- **The QC report could not be produced on the default settings.** Its three optional
  inputs all fall back to a placeholder, and all three features are off by default,
  so every run handed `QC_REPORT` the same `NO_FILE` path three times and Nextflow
  refused to stage them (`input file name collision`). The placeholders are now
  distinct files under `assets/`.
- **`main.nf` did not compile under the strict Nextflow parser**, the default from
  25.10, and only ran through the deprecated v1 fallback. The samplesheet parsing is
  now a function, every statement lives inside the workflow, and each dynamic
  `publishDir` is a closure. Verified on 24.04.2 and 26.04.6.
- Removed `process.publishDirMode`, which is not a Nextflow directive and only
  produced a warning on every run.
- A blank line in a bedgraph raised `ValueError` in `build_min_unique_len.py`
  instead of being skipped.
- **A GFF attribute was matched as a substring**, so `locus_tag=` also matched
  inside `old_locus_tag=` and `gene=` inside `pseudogene=`, both routine in RefSeq
  and Prokka output. Genes were labelled with an obsolete tag, and the curated
  H37Rv gene map could be silently overridden.
- `mask_profile` summed raw interval lengths for its headline total, so an
  unmerged BED could report above 100% masked.
- Four parsers in `stats_to_legacy.py` wrapped a whole read loop in one
  `try/except`, so one malformed line silently discarded every line after it; a
  partial genome size then propagated into depth, breadth, evenness and every
  density bin.
- Genuinely-zero metrics were written as `NA`, making a failed sample
  indistinguishable from an unmeasured one.
- Genotype classification assumed diploid spellings, so a haploid no-call was
  counted as heterozygous. This matters for `-profile generic`, which is haploid.
- `safe_tabix` had drifted between its copies: the one in `annotation.nf` still
  used a `zgrep` check that cannot tell a truncated VCF from an empty one, and
  wrote a fake index without saying so.
- A missing `FORMAT/DP` reached the consensus as `ADP=.`, which reads as zero
  depth and turns a called SNP into a gap.

### Refactored

Behaviour-preserving throughout, and each step verified rather than assumed:

- **The science came out of the process scripts.** Around 500 lines of bash and
  awk lived inside `modules/*.nf`, where no test could reach it: a stub run
  replaces the script block and nothing can import it. The genotype
  re-validation, the backbone pileup parser, the hom/het bcftools expressions and
  `safe_tabix` now live in `bin/`, with 62 tests over them, including a run
  against a real bcftools. `variants.nf` drops from 485 lines to 344.
- **`qc_report.py` became a package.** 1697 lines and 53 functions are now a
  319-line CLI over `bin/qcreport/` (parsers, metrics, panels, render). The old
  and new versions produce byte-identical output on the demo cohort.
- **The workflow is one sub-workflow per stage.** `main.nf` goes from 556 lines to
  277, with the stage sequence reduced to eleven calls. Equivalence checked by
  diffing execution traces against the previous version: identical task counts and
  identical published outputs on both test profiles.

## [1.0.1] - 2026-07-19

Headline: a consolidated **interactive QC report**, alignment-free **lineage &
drug-resistance typing**, and a restructured documentation set. Rebranded to
**BAMpiro**.

### Added

- **Consolidated interactive QC report** (`bin/qc_report.py`) - a single
  self-contained HTML dashboard (no internet / CDN) built from the whole cohort,
  plus a machine-readable per-sample `qc_flags.tsv` (PASS/WARN/FAIL). Enabled with
  `--make_qc_report` (on by default). It bundles:
 - **SNP dynamics** - per-variant allele-frequency trajectories over time, gene-
    centric, with per-timepoint depth bars, a zoom control, and a series filter by
    samplesheet metadata.
 - **Epistasis** - pairs of variants with concordant / discordant dynamics,
    scored with a permutation *p*-value and Benjamini-Hochberg FDR, in card, matrix
    and table views.
 - **SNP matrix** - an explorable, downloadable site × sample matrix with the
    depth in each cell, samplesheet metadata as column-header levels, and column
    filtering by metadata (also written as a standalone TSV, `--make_snp_matrix`).
 - **Drug resistance** - a sample × drug matrix (worst WHO grade per drug) and a
    per-mutation table, driven by the Pathotypr WHO-catalogue calls.
 - Live-adjustable QC thresholds with presets, a collapsible left contents
    sidebar, per-section **(i)** info popovers, a sample exclusion basket, and
    panels for lineage summary, distributions, correlations, a metric-correlation
    heatmap, QC-space PCA, genome landscape, functional annotation, gene burden,
    variable genes, taxonomic composition (Kraken2: primary taxon / contaminants /
    unclassified, for a read-level contamination check), and aDNA damage.
 - **Report presentation & UX** - a **dark / light theme** (header toggle, follows
    the OS preference, persisted per viewer); a single homogeneous inline icon set;
    a **mobile-responsive** layout (capped frozen column, scrollable matrices, a
    dismissable contents drawer); a **triage-first** ordering (table sorted worst-QC
    first, the Flagged panel directly under the statistics with each failing margin
    shown inline); per-lineage colours from the canonical **mycolorsTB** palette;
    and the pipeline **version** + a **GitHub** link in the header and footer.
- **Alignment-free lineage & drug-resistance typing with Pathotypr**
  (`--run_pathotypr`) - now run **from the container** (bioconda build), replacing
  the external cluster binary. It types straight from the reads with diagnostic
  k-mers, so calls are **reference-agnostic**. Two `split-fastq` passes per sample
  (nested sub-lineage + WHO drug-resistance). The Zenodo marker panels + RF model
  (v1.0.0) and the MTBC-ancestor reference are bundled in the image under
  `/opt/pathotypr/`. `--pathotypr_min_alt` tunes the alt-allele cut-off for
  heteroresistant / minority DR alleles.
- **Dual amino-acid numbering** - an optional canonical re-annotation pass
  (`--annotate_canonical`, `--canonical_snpeff_db`, `--canonical_label`, default
  H37Rv) so the report shows every protein change in both the used-reference and
  the H37Rv / Mycobrowser numbering, with a link to each gene's Mycobrowser locus.
- **Length-aware read masking** - a per-reference `genmap` mappability track
  (cached and reused across runs) that drops reads too short to be uniquely placed;
  a repeat BED is emitted alongside.
- **Virgin (unmasked) consensus** output in parallel with the masked consensus.
- **CRAM output** for published alignments (`--output_cram`), a selectable
  **publish mode** (`copy` / `link`), and a trimmed default output footprint.
- **Optional region-parallel FreeBayes** for faster calling on deep samples.
- **Documentation** - a [MkDocs Material](https://squidfunk.github.io/mkdocs-material/)
  site under `docs/` (introduction, installation, quick start, configuration,
  troubleshooting, QC report, Pathotypr, outputs, downstream phylogeny, and an
  end-to-end Jupyter tutorial), with a GitHub Pages publishing workflow (activated once
  the repo is public) and the README slimmed to a landing page.
- **Metadata-driven cohort filtering & dose analyses** in the report - samplesheet
  categorical columns (e.g. `treatment`) drive a report-wide cohort filter; a numeric
  `dose` column becomes a first-class metric (scatter axes + correlation matrix, with a
  Spearman *r* / *p* read-out), plus two new panels: **Dose × treatment** (per-treatment
  dose distribution + a Kruskal-Wallis rank test) and **Variant × dose** (per-variant
  allele-frequency-vs-dose Spearman scan with Benjamini-Hochberg FDR).
- **Guided tutorial series** on the docs site - a step-by-step *Tutorials* section (first
  run, samplesheet, configuring a run, reading the QC report, lineage & DR typing, building
  a phylogeny, longitudinal analysis, running a non-TB organism) alongside the runnable
  Jupyter notebook.
- **Live example QC report** - a self-contained demo report published on the docs site so
  the interactive dashboard can be explored before running the pipeline.

### Changed

- **Rebranded to BAMpiro** (larger report UI; panels with no input hide themselves).
- The container is referenced by version **tag** and its tool versions are recorded;
  re-pin to the `@sha256:` digest for byte-for-byte reproducible runs.
- Right-sized per-process CPU / memory, streamed the deduplication pipeline, and
  produced lighter intermediates for faster `-resume`.
- Replaced the report's top navigation with the collapsible contents sidebar.
- Nextflow runtime artifacts (logs, cache, reports) are now git-ignored.

### Fixed

- Consensus generation: multiallelic reference leak, spanning-deletion gate,
  contig guard, and a chromosome-aware masker.
- Genotype revalidation (multiallelic ref-leak, het/hom), robust `mpileup` paste,
  multi-contig VCF headers, and the ploidy-1 reference-allele filter.
- Multi-reference stats wiring and sample-id sanitization.
- Stripped stray headers in the FreeBayes exclude step and masked contig-start
  repeats.
- ES5 compatibility slips in the report and the H37Rv snpEff database pre-fetch in
  the image; hardened the canonical-annotation channel and the drug-resistance
  collector; plus six further issues found in a deep audit (multi-contig handling,
  deep-coverage bins, and output-name collisions).

## [1.0.0]

Initial release. A modular, containerized Nextflow (DSL2) pipeline for bacterial
short-read analysis: automated reference indexing and SnpEff database building,
`nucmer` repeat masking, `FastP` + `Kraken2` read QC, `bwa-mem2` mapping with
duplicate marking and multi-run merging, `FreeBayes` variant calling with all-sites
backbone VCFs, consensus generation, `SnpEff` annotation, basic Pathotypr lineage
classification, and a MultiQC report.
