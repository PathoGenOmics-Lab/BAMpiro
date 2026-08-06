nextflow.enable.dsl=2

// Import centralized functions for path generation and file classification
include { getSavePath; getSampleDir } from './utils'

/* ====================================================================
    ANNOTATION MODULES
    Contains: SnpEff annotation processes and Final Legacy Stats generation
==================================================================== */

process ANNOTATE_CANONICAL {
    // Re-annotate a sample's variants against a CANONICAL reference snpEff DB (H37Rv by default; set
    // params.canonical_snpeff_db for another organism) so the report can show the amino-acid change in
    // both the used-reference and the canonical numbering. Existing ANN is stripped first so only the
    // canonical annotation remains (qc_report.py reads the first ANN). Positions must line up with the
    // canonical reference (true when the mapping reference shares its coordinates).
    tag "AnnCanonical: ${sampleId}"
    publishDir "${params.outdir}/${getSampleDir(sampleId, params)}", mode: params.publish_mode, saveAs: { filename -> getSavePath(filename, params) }
    cpus 2
    memory { 6.GB * task.attempt }

    input:
    tuple val(sampleId), val(refId), path(vcf_in)
    val canonical_db

    output:
    tuple val(sampleId), path("${sampleId}.${refId}.canonical.ann.vcf.gz"), emit: out

    script:
    """
    set -euo pipefail
    # Drop any existing ANN, then annotate against the canonical genome DB. refId is in the name so a
    # sample mapped to more than one reference does not produce a file-name collision downstream.
    ( bcftools annotate -x INFO/ANN "${vcf_in}" -Ov 2>/dev/null || zcat -f "${vcf_in}" ) > stripped.vcf
    snpEff ann -v ${canonical_db} stripped.vcf 2> ${sampleId}.${refId}.canonical.snpeff.log \\
      | awk 'BEGIN{FS="\\t"; OFS="\\t"} /^#/{print;next} /^\\[/{next} NF>=8{print}' \\
      | bgzip -c > ${sampleId}.${refId}.canonical.ann.vcf.gz
    tabix -f -p vcf ${sampleId}.${refId}.canonical.ann.vcf.gz 2>/dev/null || : > ${sampleId}.${refId}.canonical.ann.vcf.gz.tbi
    """

    stub:
    """
    touch ${sampleId}.${refId}.canonical.ann.vcf.gz
    """
}

process ANNOTATE_LEGACY_VCF {
    tag "AnnLegacy: ${sampleId}"
    // Use getSampleDir for nested output support
    publishDir "${params.outdir}/${getSampleDir(sampleId, params)}", mode: params.publish_mode, saveAs: { filename -> getSavePath(filename, params) }
    cpus 2
    memory { 6.GB * task.attempt }

    input:
    tuple val(sampleId), val(refId), val(label), path(vcf_in), path(snpeff_config), path(data_dir)

    output:
    tuple val(sampleId), val(refId), path("${sampleId}.${refId}.${label}.ann.vcf.gz"), path("${sampleId}.${refId}.${label}.ann.vcf.gz.tbi"), emit: out
    path("${sampleId}.${refId}.${label}.snpeff.stderr.log"), emit: log

    shell:
    """
    set -euo pipefail
    
    # Function to safely index VCFs, handling potential empty files
    safe_tabix () {
      local gz="\$1"; local idx="\${gz}.tbi"
      set +e; tabix -f -p vcf "\$gz"; st=\$?; set -e
      if [ \$st -ne 0 ]; then 
        # Check if file has variants or is just header
        if zgrep -vq '^#' "\$gz"; then exit \$st; else : > "\$idx"; fi
      fi
    }

    # SnpEff requires the data directory structure to be present locally
    if [ ! -e data ]; then ln -s !{data_dir} data; fi

    # Run SnpEff
    # -c: Config file
    # -v: Verbose (captured in log)
    # The output is piped to awk to remove SnpEff summary lines and keep valid VCF lines
    snpEff ann -c !{snpeff_config} -v !{refId} !{vcf_in} 2> !{sampleId}.!{refId}.!{label}.snpeff.stderr.log \\
      | awk 'BEGIN{FS="\\t"; OFS="\\t"} /^#/ {print; next} /^\\[/ {next} NF>=8 {print}' \\
      | bgzip -c > !{sampleId}.!{refId}.!{label}.ann.vcf.gz

    safe_tabix !{sampleId}.!{refId}.!{label}.ann.vcf.gz
    """

    // label must stay in the name: main.nf picks the freebayes.raw VCF by filename
    stub:
    """
    touch ${sampleId}.${refId}.${label}.ann.vcf.gz
    touch ${sampleId}.${refId}.${label}.ann.vcf.gz.tbi
    touch ${sampleId}.${refId}.${label}.snpeff.stderr.log
    """
}

process ANNOTATE_MAIN_VCF {
    tag "AnnMain: ${sampleId}"
    // Use getSampleDir for nested output support
    publishDir "${params.outdir}/${getSampleDir(sampleId, params)}", mode: params.publish_mode, saveAs: { filename -> getSavePath(filename, params) }
    cpus 2
    memory { 8.GB * task.attempt }

    input:
    tuple val(sampleId), val(refId), path(vcf_gz), path(vcf_tbi), path(snpeff_config), path(data_dir)

    output:
    tuple val(sampleId), val(refId), path("${sampleId}.${refId}.ann.vcf.gz"), path("${sampleId}.${refId}.ann.vcf.gz.tbi"), emit: vcf_ann
    path("${sampleId}.${refId}.snpeff.csv"), emit: stats
    path("${sampleId}.${refId}.snpeff.stderr.log"), emit: log

    shell:
    """
    set -euo pipefail
    
    safe_tabix () {
      local gz="\$1"; local idx="\${gz}.tbi"
      set +e; tabix -f -p vcf "\$gz"; st=\$?; set -e
      if [ \$st -ne 0 ]; then if zgrep -vq '^#' "\$gz"; then exit \$st; else : > "\$idx"; fi; fi
    }

    # Link SnpEff database directory
    if [ ! -e data ]; then ln -s !{data_dir} data; fi

    # Run SnpEff with CSV stats generation (used by MultiQC)
    snpEff ann -c !{snpeff_config} -v !{refId} !{vcf_gz} -csvStats !{sampleId}.!{refId}.snpeff.csv 2> !{sampleId}.!{refId}.snpeff.stderr.log \\
      | awk 'BEGIN{FS="\\t"; OFS="\\t"} /^#/ {print; next} /^\\[/ {next} NF>=8 {print}' \\
      | bgzip -c > !{sampleId}.!{refId}.ann.vcf.gz

    safe_tabix !{sampleId}.!{refId}.ann.vcf.gz
    """

    stub:
    """
    touch ${sampleId}.${refId}.ann.vcf.gz
    touch ${sampleId}.${refId}.ann.vcf.gz.tbi
    touch ${sampleId}.${refId}.snpeff.csv
    touch ${sampleId}.${refId}.snpeff.stderr.log
    """
}

process GENERATE_LEGACY_STATS {
    tag "Stats: ${sampleId}"
    // Use getSampleDir. getSavePath automatically places .log files into the 'stats/' subfolder.
    publishDir "${params.outdir}/${getSampleDir(sampleId, params)}", mode: params.publish_mode, saveAs: { filename -> getSavePath(filename, params) }
    cpus 1
    
    input:
    tuple val(sampleId), val(refId), path(fastp_json), path(bam_stats), path(vcf), path(ref_fai), path(pathotypr_report)

    output:
    // refId in the name so a sample mapped to >1 reference does not emit two identically-named logs
    // that collide when COLLECT_SUMMARY .collect()s them (the sample_id lives INSIDE the log content,
    // which collect_summary.py parses, so renaming the file is safe).
    path("${sampleId}.${refId}.log"), emit: legacy_log

    script:
    // Logic to handle optional Pathotypr report
    // "NO_FILE" is a placeholder passed by main.nf if pathotypr didn't run
    def patho_arg = (pathotypr_report.name != "NO_FILE") ? "--pathotypr-report ${pathotypr_report}" : ""

    """
    # Run the Python script (located in the bin/ directory)
    python3 ${projectDir}/bin/stats_to_legacy.py \\
        --sample ${sampleId} \\
        --fastp-json ${fastp_json} \\
        --bam-stats ${bam_stats} \\
        --vcf ${vcf} \\
        --ref-fai ${ref_fai} \\
        ${patho_arg}
    mv "${sampleId}.log" "${sampleId}.${refId}.log"
    """

    stub:
    """
    touch ${sampleId}.${refId}.log
    """
}
