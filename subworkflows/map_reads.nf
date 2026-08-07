/*
===============================================================================
 Stage 5: Mapping, per-sample merge and the optional length-aware read filter
===============================================================================
 Ends on the BAM that variant calling and the consensus are built from.
*/

include { MAPPING_PE; MAPPING_SE; MERGE_AND_MARKDUP; FILTER_READS } from '../modules/mapping'
include { BUILD_MAPPABILITY } from '../modules/reference'
include { asBool } from '../modules/utils'

workflow MAP_READS {

    take:
    pe_reads      // (sId, runId, r1, r2, refId, taxId) from FASTP_PE
    pe_orphans    // (sId, runId, rSE, refId, taxId): PE mates orphaned by FastP
    se_reads      // (sId, runId, rSE, refId, taxId) from FASTP_SE
    ref_bundle    // (refId, ref_fa, indices, exclude_txt)
    expectedMap   // "sampleId||refId" -> number of BAMs to expect, for the groupTuple groupKey

    main:

    // 5. Mapping (BWA)
    def map_pe  = pe_reads.map { sId, runId, r1, r2, refId, taxId -> tuple(refId, "PE", sId, runId, r1, r2, taxId) }
    def map_se1 = pe_orphans.map { sId, runId, rSE, refId, taxId -> tuple(refId, "SE", sId, runId, rSE, null, taxId) }
    def map_se2 = se_reads.map { sId, runId, rSE, refId, taxId -> tuple(refId, "SE", sId, runId, rSE, null, taxId) }
    def map_all = map_pe.mix(map_se1).mix(map_se2)

    // Join reads with references (matching RefID)
    def joined = map_all.combine(ref_bundle).filter { it[0] == it[7] }.map { refId, mode, sId, runId, r1, r2, taxId, refId2, ref_fa, indices, exclude_txt -> tuple(refId, mode, sId, runId, r1, r2, taxId, ref_fa, indices, exclude_txt) }
    def branched = joined.branch { pe: it[1] == "PE"; se: it[1] == "SE" }

    def map_pe_in = branched.pe.map { refId, mode, sId, runId, r1, r2, taxId, ref_fa, indices, exclude_txt -> tuple(refId, sId, runId, r1, r2, taxId, ref_fa, indices, exclude_txt) }
    def map_se_in = branched.se.map { refId, mode, sId, runId, r1, r2, taxId, ref_fa, indices, exclude_txt -> tuple(refId, sId, runId, r1, taxId, ref_fa, indices, exclude_txt) }

    def mapped_pe = MAPPING_PE(map_pe_in)
    def mapped_se = MAPPING_SE(map_se_in)
    def all_bams = mapped_pe.bam.mix(mapped_se.bam)

    // Group BAMs by SampleID + RefID for merging
    def bams_grouped = all_bams.map { sId, rId, bam, ref_fa, exclude_txt ->
            def key = "${sId}||${rId}"
            def nExp = expectedMap[key]
            tuple(groupKey(key, nExp as int), sId, rId, bam, ref_fa, exclude_txt)
        }.groupTuple(by: 0).map { gk, sIds, rIds, bams, ref_fas, excls -> tuple(sIds[0], rIds[0], bams, ref_fas[0], excls[0]) }

    def final_bams = MERGE_AND_MARKDUP(bams_grouped)

    // Optional length-aware read filter: drop reads too short to map uniquely at their locus
    // (short-read false positives in near-repeats), then call/consensus on the filtered BAM.
    // Feature off -> call on the dedup BAM unchanged.
    def vbase
    def filter_stats = Channel.empty()
    if (asBool(params.dynamic_read_filter)) {
        def mapp = BUILD_MAPPABILITY( ref_bundle.map { rId, fa, idx, excl -> tuple(rId, fa) } )
        def track_bed = mapp.track.join(mapp.repeat_bed, by: 0)                    // (rId, npz, repeat_bed)
        def filt_in = final_bams.final_bam
            .map { sId, rId, bam, bai, fa, excl -> tuple(rId, sId, bam, bai, fa, excl) }
            .combine(track_bed, by: 0)                                             // fan the per-reference track+bed to each sample
            .map { rId, sId, bam, bai, fa, excl, npz, bed -> tuple(sId, rId, bam, bai, fa, excl, npz, bed) }
        def filt = FILTER_READS(filt_in)
        vbase = filt.filtered_bam                                                  // exclude_txt now = nucmer + genmap repeats
        filter_stats = filt.stats
    } else {
        vbase = final_bams.final_bam
    }

    emit:
    variant_base = vbase                    // (sId, rId, bam, bai, ref_fa, exclude_txt): the BAM to call on
    dedup_bam    = final_bams.final_bam     // the ORIGINAL dedup BAM, before any read filtering
    stats        = final_bams.stats         // (sId, rId, dedup.stats)
    filter_mqc   = filter_stats             // empty when dynamic_read_filter is off
}
