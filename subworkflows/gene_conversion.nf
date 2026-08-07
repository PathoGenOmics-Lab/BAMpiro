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

include { PARALOG_MAP; FIND_GENE_CONVERSION; COLLECT_GENE_CONVERSION } from '../modules/gene_conversion'
include { asBool } from '../modules/utils'

workflow GENE_CONVERSION {

    take:
    ref_delta      // [refId, self_aln.delta] from PREPARE_REFERENCE
    dedup_bam      // [sampleId, refId, bam, bai, ref_fa, exclude_txt] BEFORE the mappability filter
    tsv_name

    main:
    def cohort = Channel.empty()
    def tracts = Channel.empty()

    if (asBool(params.find_gene_conversion)) {
        def maps = PARALOG_MAP(ref_delta)

        // Fan the per-reference site list out to every sample mapped against that reference.
        def gconv_in = dedup_bam
            .map { sId, rId, bam, bai, fa, excl -> tuple(rId, sId, bam, bai) }
            .combine(maps.sites, by: 0)
            .map { rId, sId, bam, bai, sites -> tuple(sId, rId, bam, bai, sites) }

        tracts = FIND_GENE_CONVERSION(gconv_in).tracts
        cohort = COLLECT_GENE_CONVERSION(
            tracts.map { sId, rId, tsv -> tsv }.collect().ifEmpty([]), tsv_name).cohort
    }

    emit:
    per_sample = tracts
    cohort_tsv = cohort
}
