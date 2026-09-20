/*
===============================================================================
 Stage 4: Pathotypr typing
===============================================================================
*/

include { RUN_PATHOTYPR_PE; RUN_PATHOTYPR_SE } from '../modules/pathotypr'
include { asBool } from '../modules/utils'

workflow LINEAGE_TYPING {

    take:
    pe_reads   // (sId, runId, r1, r2, refId, taxId) from FASTP_PE
    pe_orphans // (sId, runId, rSE, refId, taxId) from FASTP_PE: merged pairs + unpaired mates
    se_reads   // (sId, runId, rSE, refId, taxId) from FASTP_SE

    main:

    // 4. Pathotypr (Optional lineage + drug-resistance typing; k-mer, reference-agnostic)
    def patho_results    = Channel.empty()   // lineage summary per sample
    def patho_dr_results = Channel.empty()   // drug-resistance summary per sample

    if (asBool(params.run_pathotypr)) {
        // Reshape channels for Pathotypr
        // Pair each PE run with its own orphan stream before typing. FastP --merge puts the
        // overlapping fragments there, and they are most of a typical library, so r1+r2 alone is a
        // minority of the reads. MAP_READS has taken all three streams from the start; this is the
        // same read set, reaching the typer.
        //
        // failOnMismatch, because the failure this replaces was silent. FASTP_PE always emits an
        // se_combined file, so a run with no partner here means the assumption has changed, and a
        // plain join would answer that by dropping the sample from lineage typing without a word.
        def orphan_by_run = pe_orphans.map { sId, runId, rSE, refId, taxId -> tuple(tuple(sId, runId), rSE) }
        def patho_pe_ch = pe_reads
            .map { sId, runId, r1, r2, refId, taxId -> tuple(tuple(sId, runId), refId, sId, runId, r1, r2) }
            .join(orphan_by_run, failOnMismatch: true, failOnDuplicate: true)
            .map { key, refId, sId, runId, r1, r2, rSE -> tuple(refId, sId, runId, r1, r2, rSE) }
        def patho_se_ch = se_reads.map { sId, runId, rSE, refId, taxId -> tuple(refId, sId, runId, rSE) }

        def run_pe = RUN_PATHOTYPR_PE(
            patho_pe_ch,
            params.pathotypr_ref,
            params.pathotypr_markers,
            params.pathotypr_dr_markers,
            params.pathotypr_bin
        )

        def run_se = RUN_PATHOTYPR_SE(
            patho_se_ch,
            params.pathotypr_ref,
            params.pathotypr_markers,
            params.pathotypr_dr_markers,
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
