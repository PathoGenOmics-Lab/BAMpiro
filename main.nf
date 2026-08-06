#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
===============================================================================
 BAMpiro Pipeline - Short Read Mapping, Variant Calling, Lineage & Reporting
===============================================================================
*/

// --- MODULE IMPORTS ---
include { PREPARE_REFERENCE; SNPEFF_BUILD_DB; BUILD_MAPPABILITY } from './modules/reference'
include { VALIDATE_RAW_READS_PE; VALIDATE_RAW_READS_SE; KRAKEN_FILTER_PE; KRAKEN_FILTER_SE; FASTP_PE; FASTP_SE; MULTIQC; DUMP_VERSIONS } from './modules/qc'
include { RUN_PATHOTYPR_PE; RUN_PATHOTYPR_SE } from './modules/pathotypr'
include { MAPPING_PE; MAPPING_SE; MERGE_AND_MARKDUP; FILTER_READS } from './modules/mapping'
include { CALL_FREEBAYES; CALL_BACKBONE; MERGE_VCFS; CALL_FREEBAYES_RAW } from './modules/variants'
include { CONSENSUS_FASTA } from './modules/consensus'
// Aliases for the parallel virgin (unmasked) consensus path (a DSL2 process runs once per name)
include { CALL_BACKBONE as CALL_BACKBONE_RAW; MERGE_VCFS as MERGE_VCFS_RAW } from './modules/variants'
include { CONSENSUS_FASTA as CONSENSUS_FASTA_RAW } from './modules/consensus'
include { ANNOTATE_LEGACY_VCF; ANNOTATE_MAIN_VCF; ANNOTATE_CANONICAL; GENERATE_LEGACY_STATS } from './modules/annotation'
include { COLLECT_SUMMARY; COLLECT_DR; LIFT_VARIANTS; QC_REPORT; SNP_MATRIX } from './modules/report'
include { cleanStr; nullish; sanitizeId } from './modules/utils'

/* ----------------------------- Helpers ----------------------------- */
// cleanStr / nullish / sanitizeId are shared with the modules; imported from modules/utils.nf above.
//
// Everything below is a FUNCTION, and every statement lives inside `workflow`. That is what the
// strict parser (the default from Nextflow 25.10) requires: a script may declare includes,
// processes, workflows and functions, but may not run statements at the top level.

def inferRunId(String r1) {
    def name = new File(r1).getName()
    name = name.replaceAll(/(\.fastq|\.fq)(\.gz|\.bz2)?$/,'')
    name = name.replaceAll(/(\.gz|\.bz2)$/,'')
    return sanitizeId(name)
}

/* ----------------------------- TSV Parsing ----------------------------- */

// Read and validate the samplesheet, returning the parsed cohort as
//   [refMap, refGffMap, expectedMap, peList, seList, hasTax]
// Throws with every problem listed at once, so the user fixes them in a single pass.
def parseSamplesheet(String tsvPath) {

def tsvFile = new File(tsvPath)
if (!tsvFile.exists()) throw new RuntimeException("Samplesheet not found: ${tsvPath}")

def lines = tsvFile.readLines()
if (lines.size() < 2) throw new RuntimeException("Samplesheet has a header but no data rows: ${tsvPath}")

def header = lines[0].replace('\r','').split('\t')*.trim()
def col = [:]
header.eachWithIndex { h,i -> col[h]=i }

// Validate the WHOLE samplesheet before running anything, collecting every problem so the
// user can fix them all in one pass instead of rerunning after each first error.
def errors = []

['sampleId','r1','refId','refFasta','refGff'].each { req ->
    if (!col.containsKey(req)) errors << "header is missing the required column '${req}'"
}
if (errors) {
    log.error "Samplesheet ${tsvPath} is unusable:\n" + errors.collect{ "  - ${it}" }.join('\n')
    throw new RuntimeException("Samplesheet header invalid (${errors.size()} problem(s)); see the list above.")
}

def refMap    = [:]
def refGffMap = [:]
def expectedMap = [:].withDefault{0}
def seenRuns = [:].withDefault{0}
def peList = []
def seList = []
def hasTax = false

lines.drop(1).eachWithIndex { raw, idx ->
    def rowNum = idx + 2   // header is line 1, so the first data row is line 2
    def line = raw.replace('\r','')
    if (!line.trim()) return
    if (line.startsWith('#')) return

    // split(-1) only keeps trailing empties when the tabs are physically present; pad short rows so a
    // row that omits trailing columns becomes empty fields (reported below) instead of crashing on an
    // out-of-bounds Java-array read.
    def cells = line.split('\t', -1).toList()
    def p = cells.size() < header.size() ? cells + ([''] * (header.size() - cells.size())) : cells

    def sampleId = sanitizeId(p[col.sampleId])
    def r1Str    = cleanStr(p[col.r1])
    def r2Str    = col.containsKey('r2') ? cleanStr(p[col.r2]) : null
    def refId    = sanitizeId(p[col.refId])
    def refFasta = cleanStr(p[col.refFasta])
    def refGff   = cleanStr(p[col.refGff])
    def taxId    = col.containsKey('taxId') ? cleanStr(p[col.taxId]) : null
    def runId    = col.containsKey('runId') ? cleanStr(p[col.runId]) : null

    // Required fields: record every missing one for this row, then skip the row.
    def missingCols = []
    if (!sampleId) missingCols << 'sampleId'
    if (!r1Str)    missingCols << 'r1'
    if (!refId)    missingCols << 'refId'
    if (!refFasta) missingCols << 'refFasta'
    if (!refGff)   missingCols << 'refGff'
    if (missingCols) { errors << "row ${rowNum}: empty required field(s): ${missingCols.join(', ')}"; return }

    if (nullish(r2Str)) r2Str = null
    def mode = (r2Str ? "PE" : "SE")

    if (nullish(taxId)) taxId = null   // treat NA / N/A / null / . / blank as "no taxId" (as r2 already does)
    if (taxId != null) hasTax = true

    if (!runId) runId = inferRunId(r1Str)
    runId = sanitizeId(runId)
    if (!runId) runId = "run"

    // Guarantee a unique run token per (sampleId, refId): lane-split inputs with no runId column
    // infer the same basename for every lane -> identical BAM names collide in the merge group.
    // Kept as two statements: the strict (v2) parser used by Nextflow >= 25.10 rejects an
    // assignment used as an expression, so `def n = (map[k] = map[k] + 1)` fails to compile.
    def runKey = "${sampleId}||${refId}||${runId}"
    seenRuns[runKey] = seenRuns[runKey] + 1
    def runN = seenRuns[runKey]
    if (runN > 1) runId = "${runId}_${runN}"

    // Input files: record every missing path for this row, then skip the row.
    def missingFiles = []
    if (!new File(r1Str).exists()) missingFiles << "r1=${r1Str}"
    if (mode == "PE" && !new File(r2Str).exists()) missingFiles << "r2=${r2Str}"
    if (!new File(refFasta).exists()) missingFiles << "refFasta=${refFasta}"
    if (!new File(refGff).exists())   missingFiles << "refGff=${refGff}"
    if (missingFiles) { errors << "row ${rowNum} (sample ${sampleId}): file(s) not found: ${missingFiles.join('; ')}"; return }

    if (!refMap.containsKey(refId)) {
        refMap[refId] = refFasta
        refGffMap[refId] = refGff
    } else if (refMap[refId] != refFasta) {
        errors << "row ${rowNum}: refId '${refId}' points to a different refFasta than an earlier row (${refMap[refId]} vs ${refFasta})"
        return
    }

    def key = "${sampleId}||${refId}"
    expectedMap[key] = expectedMap[key] + (mode == "PE" ? 2 : 1)

    if (mode == "PE") {
        peList << [sampleId, runId, file(r1Str), file(r2Str), refId, taxId]
    } else {
        seList << [sampleId, runId, file(r1Str), refId, taxId]
    }
}

if (errors) {
    log.error "Samplesheet ${tsvPath} has ${errors.size()} problem(s):\n" + errors.collect{ "  - ${it}" }.join('\n')
    throw new RuntimeException("Samplesheet validation failed with ${errors.size()} problem(s); fix the rows listed above and rerun.")
}
if (peList.isEmpty() && seList.isEmpty()) {
    throw new RuntimeException("Samplesheet ${tsvPath} produced no usable samples (every data row was blank or commented out).")
}

return [ refMap: refMap, refGffMap: refGffMap, expectedMap: expectedMap,
         peList: peList, seList: seList, hasTax: hasTax ]
}

/* ----------------------------- Main Workflow ----------------------------- */

workflow {

    /* ----------------------------- Configuration ----------------------------- */

    // Checked first: everything below dereferences params.tsv, so without this the user gets
    // "Argument of `file()` function cannot be null" instead of being told what to pass.
    if (!params.tsv) {
        throw new RuntimeException(
            "--tsv is required: a Tab-Separated samplesheet with one row per (sample, run, reference).\n" +
            "  Columns: sampleId, runId, r1, r2, refId, refFasta, refGff, taxId\n" +
            "  Example: nextflow run main.nf --tsv samples.tsv --outdir results -profile local,docker\n" +
            "  See docs/tutorials/samplesheet.md")
    }

    def final_outdir = file(params.outdir).toAbsolutePath().toString()
    // Base name of the TSV ("samples.tsv" -> "samples"), used to name the cohort-level outputs.
    def tsv_name = file(params.tsv).simpleName
    def multiqc_report_filename = "${tsv_name}_multiqc_report"

    def sheet       = parseSamplesheet(params.tsv as String)
    def refMap      = sheet.refMap
    def refGffMap   = sheet.refGffMap
    def expectedMap = sheet.expectedMap

    // Check if Kraken DB exists
    def KRAKEN_ENABLED = false
    if (params.kraken2_db && sheet.hasTax) {
        def db = new File(params.kraken2_db as String)
        if (db.exists()) KRAKEN_ENABLED = true
        else log.warn "Kraken DB not found at: ${params.kraken2_db} -> Kraken disabled"
    }

    // Executor and profile are shown because WHERE the run lands is the easiest thing to get wrong:
    // with no -profile everything runs on the current host, which on a cluster is the login node.
    log.info """
================================================================
 BAMpiro Pipeline 🧛‍♂️🧬  v${workflow.manifest.version}
================================================================
TSV              : ${params.tsv}
Output Absolute  : ${final_outdir}
Report Name      : ${multiqc_report_filename}.html
Profile(s)       : ${workflow.profile}
Executor         : ${workflow.session.config.navigate('process.executor') ?: 'local'}
Container        : ${params.container}
Kraken2          : ${KRAKEN_ENABLED ? params.kraken2_db : 'disabled'}
Annotate Legacy  : ${params.annotate_legacy_vcfs}
Run Pathotypr    : ${params.run_pathotypr}
================================================================
"""

    /* ----------------------------- Channels ----------------------------- */

    def ref_in_ch = Channel.fromList(refMap.collect { k, v -> [k, file(v)] })
                           .map { refId, fasta -> tuple(refId, fasta) }

    def snpeff_db_in = Channel.fromList(refMap.keySet().collect { rid ->
        [rid, file(refMap[rid]), file(refGffMap[rid])]
    }).map { rid, fa, gff -> tuple(rid, fa, gff) }

    def pe_reads_ch = Channel.fromList(sheet.peList).map { sId, runId, r1, r2, refId, taxId -> tuple(sId, runId, r1, r2, refId, taxId) }
    def se_reads_ch = Channel.fromList(sheet.seList).map { sId, runId, r1, refId, taxId -> tuple(sId, runId, r1, refId, taxId) }

    // 1. Reference Preparation
    def ref_bundle = PREPARE_REFERENCE(ref_in_ch)
    def snpeff_db  = SNPEFF_BUILD_DB(snpeff_db_in)

    // 2. Read Validation
    def checked_pe = VALIDATE_RAW_READS_PE(pe_reads_ch)
    def checked_se = VALIDATE_RAW_READS_SE(se_reads_ch)
    def pe_final = checked_pe.reads
    def se_final = checked_se.reads

    // Channel to collect Kraken reports if enabled
    def ch_kraken_reports = Channel.empty()

    if (KRAKEN_ENABLED) {
        def pe_with_tax = checked_pe.reads.filter { it[5] != null }
        def pe_no_tax   = checked_pe.reads.filter { it[5] == null }
        def se_with_tax = checked_se.reads.filter { it[4] != null }
        def se_no_tax   = checked_se.reads.filter { it[4] == null }
        
        def krk_pe = KRAKEN_FILTER_PE(pe_with_tax, params.kraken2_db)
        def krk_se = KRAKEN_FILTER_SE(se_with_tax, params.kraken2_db)
        
        pe_final = krk_pe.reads.mix(pe_no_tax)
        se_final = krk_se.reads.mix(se_no_tax)

        // Accumulate reports for MultiQC
        ch_kraken_reports = krk_pe.report.mix(krk_se.report)
    }

    // 3. FastP QC
    def fastp_pe = FASTP_PE(pe_final)
    def fastp_se = FASTP_SE(se_final)
    
    // 4. Pathotypr (Optional lineage + drug-resistance typing; k-mer, reference-agnostic)
    def patho_results    = Channel.empty()   // lineage summary per sample
    def patho_dr_results = Channel.empty()   // drug-resistance summary per sample

    if (params.run_pathotypr) {
        // Reshape channels for Pathotypr
        def patho_pe_ch = fastp_pe.pe_reads.map { sId, runId, r1, r2, refId, taxId -> tuple(refId, sId, runId, r1, r2) }
        def patho_se_ch = fastp_se.se_reads.map { sId, runId, rSE, refId, taxId -> tuple(refId, sId, runId, rSE) }

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

    // 5. Mapping (BWA)
    def map_pe  = fastp_pe.pe_reads.map { sId, runId, r1, r2, refId, taxId -> tuple(refId, "PE", sId, runId, r1, r2, taxId) }
    def map_se1 = fastp_pe.se_reads.map { sId, runId, rSE, refId, taxId -> tuple(refId, "SE", sId, runId, rSE, null, taxId) }
    def map_se2 = fastp_se.se_reads.map { sId, runId, rSE, refId, taxId -> tuple(refId, "SE", sId, runId, rSE, null, taxId) }
    def map_all = map_pe.mix(map_se1).mix(map_se2)

    // Join reads with references (matching RefID)
    def joined = map_all.combine(ref_bundle.bundle).filter { it[0] == it[7] }.map { refId, mode, sId, runId, r1, r2, taxId, refId2, ref_fa, indices, exclude_txt -> tuple(refId, mode, sId, runId, r1, r2, taxId, ref_fa, indices, exclude_txt) }
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

    // 6. Variant Calling
    // Optional length-aware read filter: drop reads too short to map uniquely at their locus
    // (short-read false positives in near-repeats), then call/consensus on the filtered BAM.
    // Feature off -> call on the dedup BAM unchanged.
    def vbase
    def filter_stats = Channel.empty()
    if (params.dynamic_read_filter) {
        def mapp = BUILD_MAPPABILITY( ref_bundle.bundle.map { rId, fa, idx, excl -> tuple(rId, fa) } )
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

    // Run FreeBayes and Backbone parallel
    def fb_out = CALL_FREEBAYES(vbase)
    def bb_in  = vbase.map { sId, rId, bam, bai, ref_fa, exclude_txt -> tuple(sId, rId, bam, bai, ref_fa) }
    def bb_out = CALL_BACKBONE(bb_in)
    
    // Merge specific variants with backbone for All-Positions VCF (outLabel "" = masked path)
    def merge_in = fb_out.snps.join(bb_out.backbone, by: [0,1]).join(bb_out.header, by: [0,1])
        .map { sId, rId, sv, st, bv, bt, hdr -> tuple(sId, rId, sv, st, bv, bt, hdr, "") }
    def vcf_ch = MERGE_VCFS(merge_in)

    // 7. Consensus Generation (Optional)
    def masked_consensus = Channel.empty()   // captured for the cohort QC report (section 11)
    if (params.make_consensus) {
        // Prepare inputs: VCF + Reference + Mask sites
        // Note: Script is called from bin/ directly in the module
        // refmeta from vbase so the masked consensus gets the combined (nucmer + genmap) exclude
        def refmeta = vbase.map { sId, rId, bam, bai, ref_fa, exclude_txt -> tuple(sId, rId, ref_fa, exclude_txt) }
        def allpos_mask = vcf_ch.allpos.join(fb_out.mask_sites, by: [0,1])
        def cons_in = allpos_mask.join(refmeta, by: [0,1]).map { sId, rId, vcf_gz, tbi, mask, ref_fa, exclude_txt -> tuple(sId, rId, vcf_gz, tbi, mask, ref_fa, exclude_txt, "") }

        masked_consensus = CONSENSUS_FASTA(cons_in).fasta

        // 7b. Virgin (unmasked) consensus from the ORIGINAL dedup BAM, in parallel with the masked one.
        if (params.keep_virgin_consensus && params.dynamic_read_filter) {
            def raw_bb_in = final_bams.final_bam.map { sId, rId, bam, bai, ref_fa, excl -> tuple(sId, rId, bam, bai, ref_fa) }
            def bb_raw = CALL_BACKBONE_RAW(raw_bb_in)
            // Virgin SNPs: re-call FreeBayes on the raw bam (shows pre-filter variants) or reuse masked SNPs.
            def raw_snps = params.virgin_full_freebayes ? CALL_FREEBAYES_RAW(final_bams.final_bam).snps : fb_out.snps
            def merge_raw = raw_snps.join(bb_raw.backbone, by: [0,1]).join(bb_raw.header, by: [0,1])
                .map { sId, rId, sv, st, bv, bt, hdr -> tuple(sId, rId, sv, st, bv, bt, hdr, ".raw") }
            def vcf_raw = MERGE_VCFS_RAW(merge_raw)
            def refmeta_raw = final_bams.final_bam.map { sId, rId, bam, bai, ref_fa, excl -> tuple(sId, rId, ref_fa, excl) }
            def cons_raw = vcf_raw.allpos.join(fb_out.mask_sites, by: [0,1]).join(refmeta_raw, by: [0,1])
                .map { sId, rId, vcf_gz, tbi, mask, ref_fa, excl -> tuple(sId, rId, vcf_gz, tbi, mask, ref_fa, excl, ".raw") }
            CONSENSUS_FASTA_RAW(cons_raw)
        }
    }

    // 8. Annotation
    // Annotate Legacy split VCFs (Homo/Het/Indel/Raw)
    def freebayes_ann = Channel.empty()   // annotated freebayes.raw (AD + snpEff ANN) -> SNP dynamics
    if (params.annotate_legacy_vcfs) {
        def leg_inputs = Channel.empty()
            .mix(fb_out.homo_snp.map   { sId, rId, vcf -> tuple(rId, sId, "var.homo.SNPs", vcf) })
            .mix(fb_out.het_snp.map    { sId, rId, vcf -> tuple(rId, sId, "var.het.SNPs", vcf) })
            .mix(fb_out.homo_indel.map { sId, rId, vcf -> tuple(rId, sId, "var.homo.indel", vcf) })
            .mix(fb_out.raw_fb.map     { sId, rId, vcf, tbi -> tuple(rId, sId, "freebayes.raw", vcf) })

        def ann_leg_in = leg_inputs.combine(snpeff_db.db, by: 0)
            .map { rId, sId, lbl, vcf, cfg, dat -> tuple(sId, rId, lbl, vcf, cfg, dat) }

        def ann_leg = ANNOTATE_LEGACY_VCF(ann_leg_in)
        // Keep the annotated freebayes VCF: it carries per-alt AD (true AF) AND gene/effect (ANN),
        // so the dynamics panel gets continuous allele frequencies without losing gene annotation.
        freebayes_ann = ann_leg.out
            .filter { sId, rId, vcf, tbi -> vcf.name.contains('freebayes.raw') }
            .map { sId, rId, vcf, tbi -> vcf }
    }

    // Annotate Main VCF
    def vcf_for_stats = Channel.empty()
    def ch_snpeff_stats = Channel.empty()
    
    if (params.annotate_main_vcf) {
        def main_pre = vcf_ch.main_vcf.map { sId, rId, vcf, tbi -> tuple(rId, sId, vcf, tbi) }
        def ann_main_in = main_pre.combine(snpeff_db.db, by: 0)
            .map { rId, sId, vcf, tbi, cfg, dat -> tuple(sId, rId, vcf, tbi, cfg, dat) }
        
        def ann_out = ANNOTATE_MAIN_VCF(ann_main_in)
        
        vcf_for_stats = ann_out.vcf_ann.map { sId, rId, vcf, tbi -> tuple(sId, rId, vcf) }
        ch_snpeff_stats = ann_out.stats // Collect SnpEff stats CSV for MultiQC
    } else {
        vcf_for_stats = vcf_ch.allpos.map { sId, rId, vcf, tbi -> tuple(sId, rId, vcf) }
    }

    // 9. Legacy Stats Generation
    // Prepare Reference Index channel
    def ref_fai_ch = ref_bundle.bundle.map { rId, fa, indices, excl -> 
        def fai_file = indices.find { it.name.endsWith('.fai') }
        tuple(rId, fai_file)
    }

    // Prepare FastP JSONs (grouped by sample). A multi-lane sample has one JSON per lane; pass them ALL
    // so GENERATE_LEGACY_STATS aggregates read counts (SUM) and quality rates (read/base-weighted) instead
    // of silently keeping a single lane. Grouped on the real sampleId val (do NOT re-derive it from the
    // filename: a sampleId containing '__' would be truncated by split('__')). Sorted for a stable -resume.
    def json_ch = fastp_pe.json.mix(fastp_se.json)   // [sId, json]
        .groupTuple()                                 // -> [sId, [json_lane1, json_lane2, ...]]
        .map { sId, jsons -> tuple(sId, jsons.toSorted { it.name }) }

    // BAM stats already carry (sampleId, refId) as vals -- do NOT re-parse the filename with
    // tokenize('.'), which truncates any dotted refId (e.g. NC_000962.3) and breaks the join.
    def bam_stats_ch = final_bams.stats // [sId, rId, stats]

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

    
    def legacy_stats = GENERATE_LEGACY_STATS(final_stats_input)

    // 10. MultiQC Report
    // Collect all relevant metrics from previous processes.
    // FastP JSONs and Kraken reports are reference-independent (QC of the raw reads), so a sample
    // mapped to >1 reference emits identically-named copies from parallel tasks. Deduplicate by
    // filename before collecting, otherwise MultiQC hits a fatal input-name collision.
    def fastp_json_mqc = fastp_pe.json.mix(fastp_se.json).map { sId, json -> json }.unique { it.name }
    def kraken_mqc     = ch_kraken_reports.unique { it.name }
    // Provenance: dump the container's tool versions once and surface them in the report.
    def versions_mqc   = DUMP_VERSIONS().mqc
    def qc_collection = Channel.empty()
        .mix(fastp_json_mqc)                          // FastP
        .mix(final_bams.stats.map { s, r, st -> st }) // Samtools (drop the (sId,rId) key -> bare path)
        .mix(kraken_mqc)                              // Kraken
        .mix(ch_snpeff_stats)                         // SnpEff
        .mix(filter_stats)                            // Length-aware read filter (drop rate)
        .mix(versions_mqc)                            // Software versions
        .collect()

    MULTIQC(qc_collection, multiqc_report_filename)

    // 11. Cohort-level outputs (QC report + master SNP matrix).
    // Per-sample annotated VCFs, shared by both: prefer the annotated freebayes VCF (per-alt AD ->
    // continuous AF + gene/effect); fall back to the annotated main VCF (GT-based AF) if legacy
    // annotation is off. reference name(s) the samples were mapped against.
    def report_vcfs = (params.annotate_legacy_vcfs ? freebayes_ann
                                                   : vcf_for_stats.map { sId, rId, vcf -> vcf })
                      .collect().ifEmpty([])
    def ref_name = refMap.keySet().join(',')

    // 11a. Consolidated QC report: aggregate every per-sample legacy log into one summary TSV, then
    // render the self-contained interactive HTML + PASS/WARN/FAIL flags.
    if (params.make_qc_report) {
        def all_logs = legacy_stats.legacy_log.collect()
        def summ = COLLECT_SUMMARY(all_logs, tsv_name)
        def cons_files = masked_consensus.map { sId, rId, fa -> fa }.collect().ifEmpty([])
        // Reference-level extras (cohort report -> take the reference bundle; single-ref is the norm):
        // GFF enables the per-gene SNP-density panel, the nucmer/repeat BED the masked-regions panel.
        def report_ref  = refGffMap.keySet().toList().first()   // deterministic: the first reference
        def report_gff  = file(refGffMap[report_ref])
        def report_mask = ref_bundle.bundle.filter { rId, fa, idx, excl -> rId == report_ref }
                                           .map { rId, fa, idx, excl -> excl }.first()
        // Provenance footer: pinned container digest + reference(s).
        def provenance  = (["container=${params.container}"] + refMap.keySet().collect { "reference=${it}" })
                          .collect { "\"${it}\"" }.join(' ')
        // SNP dynamics: the samplesheet is the metadata source (auto-detects time/group columns; the
        // panel self-hides if absent).
        def report_meta = file(params.tsv)
        // Optional QC_REPORT inputs use a per-input placeholder from assets/. They must have DISTINCT
        // names: Nextflow refuses to stage two inputs of the same task under one filename, and all three
        // of these are off by default, so a single shared "NO_FILE" fails the run with an input file name
        // collision. They are real empty files (a dangling symlink breaks under stageInMode 'copy') and
        // are referenced through projectDir so a stray NO_FILE in the launch directory is never picked up.
        // QC_REPORT recognises them by the NO_FILE* prefix.
        // Drug-resistance calls (pathotypr DR run -> one run TSV); placeholder when pathotypr is off.
        def dr_report = params.run_pathotypr
            ? COLLECT_DR(patho_dr_results.map { sId, f -> f }.collect().ifEmpty([]), tsv_name).dr
            : file("${projectDir}/assets/NO_FILE_DR")
        // Canonical-reference-annotated VCFs for the dual amino-acid numbering (off unless annotate_canonical).
        // vcf_for_stats is the per-sample (sId,rId,vcf) channel; its positions match report_vcfs, so the
        // report's sample+position merge finds each variant's canonical amino-acid change.
        def report_vcfs_h37rv = params.annotate_canonical
            ? ANNOTATE_CANONICAL(vcf_for_stats, params.canonical_snpeff_db).out
                             .map { sId, vcf -> vcf }.collect().ifEmpty([])
            : file("${projectDir}/assets/NO_FILE_H37RV")
        // Kraken2 per-sample reports (deduped) -> Taxonomic composition panel; empty when Kraken is off.
        def report_kraken = ch_kraken_reports.unique { it.name }.collect().ifEmpty([])
        // Alignment-free canonical COORDINATE per variant via pathotypr (alternative to the --vcfs-h37rv path):
        // lift the run's variant positions (mapping-reference coords) onto H37Rv and hand the map to the report.
        def report_ref_fa = ref_bundle.bundle.filter { rId, fa, idx, excl -> rId == report_ref }
                                             .map { rId, fa, idx, excl -> fa }.first()
        def pos_liftover  = params.variant_liftover
            ? LIFT_VARIANTS(report_vcfs, report_ref_fa, report_ref, tsv_name).map
            : file("${projectDir}/assets/NO_FILE_LIFTOVER")
        QC_REPORT(summ.summary, summ.gene_burden, cons_files, report_gff, report_mask,
                  report_meta, report_vcfs, report_vcfs_h37rv, pos_liftover, dr_report, report_kraken, provenance, tsv_name)
    }

    // 11b. Master SNP matrix: rows = SNP sites, columns = reference/annotation + per-sample AF & depth.
    if (params.make_snp_matrix) {
        SNP_MATRIX(report_vcfs, ref_name, tsv_name)
    }
}
