/*
===============================================================================
 Stage 4: Pathotypr typing
===============================================================================
*/

include { RUN_PATHOTYPR_PE; RUN_PATHOTYPR_SE } from '../modules/pathotypr'

workflow LINEAGE_TYPING {

    take:
    pe_reads   // (sId, runId, r1, r2, refId, taxId) from FASTP_PE
    se_reads   // (sId, runId, rSE, refId, taxId) from FASTP_SE

    main:

    // 4. Pathotypr (Optional lineage + drug-resistance typing; k-mer, reference-agnostic)
    def patho_results    = Channel.empty()   // lineage summary per sample
    def patho_dr_results = Channel.empty()   // drug-resistance summary per sample

    if (params.run_pathotypr) {
        // Reshape channels for Pathotypr
        def patho_pe_ch = pe_reads.map { sId, runId, r1, r2, refId, taxId -> tuple(refId, sId, runId, r1, r2) }
        def patho_se_ch = se_reads.map { sId, runId, rSE, refId, taxId -> tuple(refId, sId, runId, rSE) }

        def run_pe = RUN_PATHOTYPR_PE(
            patho_pe_ch,
            file(params.pathotypr_ref),
            file(params.pathotypr_markers),
            file(params.pathotypr_dr_markers),
            params.pathotypr_bin
        )

        def run_se = RUN_PATHOTYPR_SE(
            patho_se_ch,
            file(params.pathotypr_ref),
            file(params.pathotypr_markers),
            file(params.pathotypr_dr_markers),
            params.pathotypr_bin
        )

        // Combine summaries for later statistics / the QC report
        patho_results    = run_pe.summary.mix(run_se.summary)
        patho_dr_results = run_pe.dr_mutations.mix(run_se.dr_mutations)
    }

    emit:
    summary      = patho_results      // (sId, lineage summary); empty when pathotypr is off
    dr_mutations = patho_dr_results   // (sId, dr mutations); empty when pathotypr is off
}
