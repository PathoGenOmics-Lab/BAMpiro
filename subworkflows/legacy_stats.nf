/*
===============================================================================
 Stage 9: Legacy Stats Generation
===============================================================================
 One per-sample log combining the VCF, the BAM stats, the FastP JSONs and the
 optional pathotypr lineage summary.
*/

include { GENERATE_LEGACY_STATS } from '../modules/annotation'

workflow LEGACY_STATS {

    take:
    ref_bundle      // (refId, ref_fa, indices, exclude_txt)
    pe_json         // (sId, fastp json) from FASTP_PE
    se_json         // (sId, fastp json) from FASTP_SE
    // BAM stats already carry (sampleId, refId) as vals -- do NOT re-parse the filename with
    // tokenize('.'), which truncates any dotted refId (e.g. NC_000962.3) and breaks the join.
    bam_stats_ch    // [sId, rId, stats]
    vcf_for_stats   // (sId, rId, vcf)
    patho_results   // (sId, lineage summary); empty when pathotypr is off

    main:

    // 9. Legacy Stats Generation
    // Prepare Reference Index channel
    def ref_fai_ch = ref_bundle.map { rId, fa, indices, excl ->
        def fai_file = indices.find { it.name.endsWith('.fai') }
        tuple(rId, fai_file)
    }

    // Prepare FastP JSONs (grouped by sample). A multi-lane sample has one JSON per lane; pass them ALL
    // so GENERATE_LEGACY_STATS aggregates read counts (SUM) and quality rates (read/base-weighted) instead
    // of silently keeping a single lane. Grouped on the real sampleId val (do NOT re-derive it from the
    // filename: a sampleId containing '__' would be truncated by split('__')). Sorted for a stable -resume.
    def json_ch = pe_json.mix(se_json)                // [sId, json]
        .groupTuple()                                 // -> [sId, [json_lane1, json_lane2, ...]]
        .map { sId, jsons -> tuple(sId, jsons.toSorted { it.name }) }

    // Join logic for Stats
    def vcf_bam_joined = vcf_for_stats.join(bam_stats_ch, by: [0,1]) // [sId, rId, vcf, stats]
    // combine (not join) on sampleId: one json per sample must fan out to EVERY (sample,ref)
    // row, else a sample mapped to >1 reference silently loses all but one reference's stats.
    def vcf_bam_fastp = vcf_bam_joined.combine(json_ch, by: 0) // [sId, rId, vcf, stats, json]

    // Join with Reference Index (Using COMBINE to reuse reference)
    def ready_no_patho = vcf_bam_fastp
        .map { sId, rId, vcf, stats, json -> tuple(rId, sId, vcf, stats, json) }
        .combine(ref_fai_ch)
        .filter { it[0] == it[5] } // Filter where sample RefID == FAI RefID
        .map { rId, sId, vcf, stats, json, rId2, fai -> tuple(sId, rId, json, stats, vcf, fai) }

    // Join with Pathotypr (Using remainder to handle missing files safely)
    def final_stats_input = ready_no_patho
        .join(patho_results, by: 0, remainder: true)
        .map { list ->
            // list size < 7 means missing left-side data (e.g. mapping failed)
            if (list.size() < 7) return null

            def sId   = list[0]
            def rId   = list[1]
            def json  = list[2]
            def stats = list[3]
            def vcf   = list[4]
            def fai   = list[5]
            def patho = list[6] // Can be null

            // Fallback placeholder file (see the assets/NO_FILE* note in section 11a)
            def real_patho = patho ? patho : file("${projectDir}/assets/NO_FILE")
            tuple(sId, rId, json, stats, vcf, fai, real_patho)
        }
        .filter { it != null }


    def stats_out = GENERATE_LEGACY_STATS(final_stats_input)

    emit:
    legacy_log = stats_out.legacy_log   // one log per (sample, reference)
}
