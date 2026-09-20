/*
 * Gene conversion detection (opt-in: --find_gene_conversion).
 *
 * A conversion copies a stretch of one paralog onto another, so the acceptor stops carrying its own
 * alleles and carries the donor's over a bounded tract. That happens in paralogous sequence, which
 * is precisely what the rest of the pipeline spends effort removing:
 *
 *   - `exclude_repeats` drops repeats from variant calling
 *   - `dynamic_read_filter` drops reads that cannot be placed uniquely
 *
 * Both are right for variant calling and fatal here, so this stage bypasses them: it works from the
 * paralog map (recovered from the self-alignment the masking step already computes) and from the
 * DEDUPLICATED, PRE-FILTER BAM.
 */

include { ALIGN_OUTGROUP; PARALOG_MAP; FIND_GENE_CONVERSION; COLLECT_GENE_CONVERSION } from '../modules/gene_conversion'
include { asBool } from '../modules/utils'

workflow GENE_CONVERSION {

    take:
    ref_delta      // [refId, self_aln.delta] from PREPARE_REFERENCE
    dedup_bam      // [sampleId, refId, bam, bai, ref_fa, exclude_txt] BEFORE the mappability filter
    ref_gff        // [refId, gff] so a tract can name the genes it landed on
    tsv_name

    main:
    def cohort = Channel.empty()
    def tracts = Channel.empty()

    if (asBool(params.find_gene_conversion)) {
        // The outgroup is one sequence for the whole run: a genome ancestral to, or outside, the
        // clade being studied. Absent, the placeholder rides along and every polarity column
        // comes out empty rather than guessed at.
        def no_outgroup = file("${projectDir}/assets/NO_FILE_OUTGROUP")
        def anc_delta = params.gconv_outgroup
            ? ALIGN_OUTGROUP(dedup_bam
                  .map { sId, rId, bam, bai, fa, excl -> tuple(rId, fa) }
                  .unique { it[0] }
                  .map { rId, fa -> tuple(rId, fa, file(params.gconv_outgroup)) }).delta
            : ref_delta.map { rId, d -> tuple(rId, no_outgroup) }

        // The reference FASTA travels with the delta, because the delta only names where it
        // used to be. See the note on PARALOG_MAP's input.
        def ref_fa_by_id = dedup_bam
            .map { sId, rId, bam, bai, fa, excl -> tuple(rId, fa) }
            .unique { it[0] }
        def maps = PARALOG_MAP(ref_delta.join(anc_delta).join(ref_fa_by_id))

        // Fan the per-reference site list out to every sample mapped against that reference.
        def no_gff = file("${projectDir}/assets/NO_FILE_GCONV_GFF")
        def gff_by_ref = ref_gff.map { rId, gff -> tuple(rId, gff ? file(gff) : no_gff) }
        def gconv_in = dedup_bam
            .map { sId, rId, bam, bai, fa, excl -> tuple(rId, sId, bam, bai, fa) }
            .combine(maps.sites, by: 0)
            .combine(gff_by_ref, by: 0)
            .map { rId, sId, bam, bai, fa, sites, gff -> tuple(sId, rId, bam, bai, sites, fa, gff) }

        def found = FIND_GENE_CONVERSION(gconv_in)
        tracts = found.tracts
        cohort = COLLECT_GENE_CONVERSION(
            tracts.map { sId, rId, tsv -> tsv }.collect().ifEmpty([]),
            found.loci.map { sId, rId, tsv -> tsv }.collect().ifEmpty([]),
            tsv_name).cohort
    }

    emit:
    per_sample = tracts
    cohort_tsv = cohort
}
