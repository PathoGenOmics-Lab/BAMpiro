# Lineage & Drug-Resistance Typing (Pathotypr)

`--run_pathotypr true` runs
[**Pathotypr**](https://github.com/PathoGenOmics-Lab/pathotypr) - alignment-free
(k-mer) MTBC lineage + WHO drug-resistance genotyping - **directly from the
container** (a bioconda build; the old hard-coded cluster binary is gone). Because
it types straight from the reads using diagnostic k-mers, it is
**reference-agnostic**: lineage and DR calls are correct even when your samples were
mapped to a non-H37Rv reference.

Each sample gets two `split-fastq` passes:

1. **Lineage** - nested sub-lineage classification against the bundled lineage
   markers → `<sample>__<runId>.pathotypr.lineage_summary.tsv` (the `__<runId>` infix
   keeps a multi-lane / multi-reference sample from colliding; feeds the per-lineage
   panel and the lineage colours).
2. **Drug resistance** - WHO-catalogue markers, gated by `--pathotypr_min_alt`
   (default `95` = near-fixed only; lower to ~10-25 to also catch heteroresistant /
   minority alleles) → per-sample DR mutations, aggregated into
   `<samplesheet>_dr.tsv` and shown in the **Drug resistance** panel of the
   [QC report](qc-report.md) (a sample × drug matrix with the worst WHO grade per
   drug, plus a per-mutation table).

Each call carries its **WHO catalogue grade**: **1) Assoc w R** and **2) Assoc w R –
Interim** are resistance-associated (the report and any "resistant" count use grades
1–2); **3) Uncertain**; **4) Not assoc w R – Interim** and **5) Not assoc w R** are not
associated. `<samplesheet>_dr.tsv` columns: `sample, drug, gene, mutation, grade,
marker_name, af, dp`.

!!! danger "Responsible use — a genomic screen, not a clinical DST"

    This is a **genomic screen against the WHO catalogue, not a clinical DST result.**
    Treat it as a flag for review, not a diagnosis; a call's allele frequency (`af`) and
    depth (`dp`) matter, and absence of a marker is not proof of susceptibility. When
    filtering, restrict to grades 1–2 for "resistant".

!!! warning "Catalogue version matters, and the results do not record it"
    The resistance catalogue is versioned separately from BAMpiro. Images built before
    catalogue **v1.0.2** bundled v1.0.0, which assigned each variant a single drug inherited
    from its gene instead of grading each variant-drug pair the way the WHO catalogue does.
    Amikacin did not occur anywhere in v1.0.0, so it could never be reported; the *rrs*
    aminoglycoside determinants were attributed to streptomycin, the *inhA* promoter variants
    to isoniazid alone, and 15,969 rows carried a grade belonging to a different drug. The
    upgrade moves 9,439 rows into grades 1-2 and none out of them.

    Detection is unchanged: the same variants are found either way, and MDR and pre-XDR
    assignment is identical. Results that report **amikacin, kanamycin, capreomycin,
    ethionamide, linezolid, streptomycin or delamanid** from an older image should be
    regenerated. Which catalogue a run used is printed in
    `pipeline_info/software_versions.txt` as `pathotypr markers`, together with the
    SHA-256 of the file actually read.

The marker panels and pre-trained model (Zenodo v1.0.2, DOI
[10.5281/zenodo.21915539](https://doi.org/10.5281/zenodo.21915539)) and the
MTBC-ancestor reference are **bundled in the image** under `/opt/pathotypr/`;
override the reference or marker panels with `--pathotypr_ref` / `--pathotypr_markers`
/ `--pathotypr_dr_markers` if you supply your own. (BAMpiro types lineage by k-mer
**markers** — `split-fastq --nested-classification` — not by the RF `predict` model, so
the bundled `rf_model.pathotypr` is currently unused and no flag exposes it.) See
[Configuration](configuration.md) for every flag.

## Dual amino-acid numbering (H37Rv / Mycobrowser)

With `--annotate_canonical true`, each sample's SNPs are also annotated against a
canonical snpEff database (`--canonical_snpeff_db`, default H37Rv), *in addition to*
your mapping reference. Each SNP is first moved to its H37Rv position with the k-mer
liftover below (`bin/lift_vcf.py`): its REF becomes H37Rv's base, its allele is
complemented where your reference lies reversed against H37Rv, and a SNP whose allele is
H37Rv's own base is left out, since in H37Rv numbering nothing changed. So the gene and
the amino-acid change are H37Rv's on **any** MTBC reference, whatever it calls and
numbers its genes. The report then shows every amino-acid change in **both** numberings
(the used reference + `--canonical_label`, default `H37Rv`), links each gene to its
Mycobrowser locus (`Rv…`), and matches trajectories to the resistance catalogue in
H37Rv's gene names and numbering.

The database names H37Rv's chromosome `Chromosome` (`--canonical_chrom`), which is the
name the lifted records are written under. For another organism, point
`--canonical_snpeff_db` at its snpEff genome, `--canonical_ref` at that genome's FASTA,
and `--canonical_chrom` at the name the database gives its chromosome (empty keeps the
FASTA's own names, for a multi-chromosome genome whose names already match). Turn
`--annotate_canonical` off for non-MTB organisms that have no canonical genome.

## Reference-agnostic coordinates (k-mer liftover)

`--variant_liftover true` translates every variant position into canonical (H37Rv)
coordinates for the report's SNP tables **without a whole-genome alignment and without
assuming shared coordinates**, so it works for *any* MTBC reference. It is done by
`bin/pathotypr_liftover.py` (pure Python, `lift --global-chain`): every k-mer unique in
both the mapping reference and `--canonical_ref` becomes an anchor, the anchors are
chained into collinear (and inverted) blocks, and each position is placed by
interpolating between its flanking anchors — then the placement is **verified against
the actual sequence**. The method is **"correct or absent"**: SNP sites, indels and
inversions are followed, and any position it cannot pin (a repeat / low-complexity
desert, or a rearrangement below the k-mer resolution) is *dropped* rather than
mis-mapped. `--annotate_canonical` runs it too. See `bin/pathotypr_liftover.README.md`
for the algorithm, options and limits.

Every reference of the run is lifted from its own sequence (`LIFT_VARIANTS`, one task
per reference), and every contig of it: each contig of a draft assembly is chained on
its own, so it is placed wherever and in whichever orientation it lies in H37Rv. The
maps name the contig, and the report looks each variant up by its contig and position.

Where a position cannot be interpolated — an indel between its two anchors, or no
shared unique k-mer for longer than ~2k (a SNP-dense or repeated stretch) — the pipeline
aligns the stretch between the two anchors (`--align-gaps`). The anchors pin both ends,
so the alignment cannot drift to another copy of a repeat. A position is placed when a
gap-free, ≥90%-identical run of the alignment joins it to one of the anchors and its own
context agrees; one inside a stretch H37Rv does not have (an IS copy, a region H37Rv
lost) is reported as **not in H37Rv**, which the report shows beside the variant; and
one that could sit on either side of an indel in a repeat is left without a coordinate.

On the two references of a 185-sample cohort, checked against minimap2's alignment of
each to H37Rv: every position of 27 resistance genes lifts to minimap2's coordinate,
99.8% of 20,000 random positions and of the report's variant sites do, and not one of
those that differ fits H37Rv worse than minimap2's coordinate — they are exact ties in
identical repeat units, bases at the junction of a structural difference, or places
where minimap2 aligned another copy. About 0.2% stay without a coordinate, in PE_PGRS-type
repeats and rearranged stretches that minimap2 mostly cannot place either. Cut into ten
contigs with six of them reversed, the same reference lifts just as well.

## Blind-spot masking (H37Rv problematic sites)

`--mask_blindspots true` adds the H37Rv Illumina **blind spots** — the repetitive /
low-mappability positions that are unreliable to call (Zenodo record
[3701840](https://zenodo.org/records/3701840), shipped as
`assets/H37Rv_blindspots.bed`, NC_000962.3 coordinates) — to each reference's
exclusion mask, so those sites are dropped from variant calling, consensus and the
report. Because the BED is in H37Rv coordinates, `--blindspot_liftover true` (default)
lifts it onto the run reference, every contig of it, with the **same k-mer liftover**, treating
`--canonical_ref` as the FASTA the blind-spots are defined on (default: the
H37Rv-colinear MTBC ancestor bundled at `/opt/pathotypr/reference.fasta`; point it at a
true H37Rv FASTA for an exact lift), so the mask is correct even when the reference is
**not** H37Rv. Set `--blindspot_liftover false` — which appends the BED directly and
ignores `--canonical_ref` — only when the reference already shares H37Rv coordinates.
The repetitive fraction is already covered by the pipeline's own repeat /
mappability masking — the blind-spots add the non-repetitive problematic sites.
