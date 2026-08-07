#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
===============================================================================
 BAMpiro Pipeline - Short Read Mapping, Variant Calling, Lineage & Reporting
===============================================================================
*/

// --- SUB-WORKFLOW IMPORTS ---
// One sub-workflow per numbered stage; each owns the processes it runs and the channel plumbing
// that feeds them, so a stage can be read and changed on its own. The processes themselves still
// live in modules/.
include { PREPARE_REFERENCES } from './subworkflows/prepare_references'   // 1. Reference preparation
include { READ_QC }            from './subworkflows/read_qc'              // 2-3. Validation, Kraken, FastP
include { LINEAGE_TYPING }     from './subworkflows/lineage_typing'       // 4. Pathotypr
include { MAP_READS }          from './subworkflows/map_reads'            // 5. Mapping, merge, read filter
include { CALL_VARIANTS }      from './subworkflows/call_variants'        // 6. Variant calling
include { MAKE_CONSENSUS }     from './subworkflows/make_consensus'       // 7. Consensus
include { ANNOTATE }           from './subworkflows/annotate'             // 8. Annotation
include { LEGACY_STATS }       from './subworkflows/legacy_stats'         // 9. Legacy stats
include { MULTIQC_REPORT }     from './subworkflows/multiqc_report'       // 10. MultiQC
include { COHORT_REPORT }      from './subworkflows/cohort_report'        // 11. Cohort outputs
include { GENE_CONVERSION }    from './subworkflows/gene_conversion'      // 12. Gene conversion (opt-in)

// --- MODULE IMPORTS ---
include { cleanStr; nullish; sanitizeId; validateParams; paramsHelp } from './modules/utils'

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

    if (params.help) {
        log.info paramsHelp("${projectDir}/nextflow.config", workflow.manifest.version)
        return
    }

    // A misspelled --flag is otherwise accepted, ignored, and the run finishes with the default.
    validateParams(params, "${projectDir}/nextflow.config")

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

    /* ----------------------------- Stages ----------------------------- */
    // Each stage lives in subworkflows/ and keeps its own conditionals, so a branch that is off
    // still emits the Channel.empty() placeholder the next stage mixes in.

    // 1. Reference Preparation
    def refs = PREPARE_REFERENCES(ref_in_ch, snpeff_db_in)

    // 2. Read Validation + 3. FastP QC (with the optional Kraken filter in between)
    def reads = READ_QC(pe_reads_ch, se_reads_ch, KRAKEN_ENABLED)

    // 4. Pathotypr (Optional lineage + drug-resistance typing; k-mer, reference-agnostic)
    def typing = LINEAGE_TYPING(reads.pe_reads, reads.se_reads)

    // 5. Mapping (BWA), per-sample merge, and the optional length-aware read filter
    def bams = MAP_READS(reads.pe_reads, reads.pe_orphans, reads.se_reads, refs.bundle, expectedMap)

    // 6. Variant Calling
    def variants = CALL_VARIANTS(bams.variant_base)

    // 12. Gene conversion (opt-in). Reads the PRE-FILTER bam and the reference self-alignment,
    // because the masking the rest of the pipeline applies removes exactly this signal.
    GENE_CONVERSION(refs.delta, bams.dedup_bam, tsv_name)

    // 7. Consensus Generation (Optional)
    def consensus = MAKE_CONSENSUS(bams.variant_base, bams.dedup_bam,
                                   variants.allpos, variants.snps, variants.mask_sites)

    // 8. Annotation
    def annotated = ANNOTATE(variants.homo_snp, variants.het_snp, variants.homo_indel,
                             variants.raw_fb, variants.main_vcf, variants.allpos, refs.snpeff)

    // 9. Legacy Stats Generation
    def legacy = LEGACY_STATS(refs.bundle, reads.pe_json, reads.se_json,
                              bams.stats, annotated.stats_vcf, typing.summary)

    // 10. MultiQC Report
    MULTIQC_REPORT(reads.pe_json, reads.se_json, bams.stats, reads.kraken_reports,
                   annotated.snpeff_stats, bams.filter_mqc, multiqc_report_filename)

    // 11. Cohort-level outputs (QC report + master SNP matrix).
    COHORT_REPORT(legacy.legacy_log, consensus.fasta, refs.bundle, reads.kraken_reports,
                  typing.dr_mutations, annotated.stats_vcf, annotated.freebayes_vcf,
                  refMap, refGffMap, tsv_name)
}
