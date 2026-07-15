nextflow.enable.dsl=2

// Import centralized functions for path generation and file classification
include { getSavePath; getSampleDir } from './utils'

/* ====================================================================
    CONSENSUS MODULES
    Contains: Generation of Consensus Fasta from VCF and Backbone
==================================================================== */

process CONSENSUS_FASTA {
    tag "Consensus: ${sampleId}"
    
    // Use getSampleDir for nested output support
    publishDir "${params.outdir}/${getSampleDir(sampleId, params)}", mode: 'copy', saveAs: { filename -> getSavePath(filename, params) }
    
    cpus 1
    memory '4 GB'

    input:
    // Input tuple: SampleID, RefID, VCF (all positions), TBI, Mask Sites, Ref Fasta, Excluded Regions
    // Note: The python script 'WGS_fasta_allpos.py' is assumed to be in the 'bin/' folder of the project
    // outLabel is "" for the normal (masked) consensus and ".raw" for the virgin one.
    tuple val(sampleId), val(refId), path(allpos_vcf_gz), path(allpos_tbi), path(mask_sites), path(ref_fa), path(exclude_txt), val(outLabel)

    output:
    tuple val(sampleId), val(refId), path("${sampleId}.${refId}${outLabel}.consensus.fasta"), emit: fasta
    path("${sampleId}.${refId}${outLabel}.consensus.log"), emit: log

    script:
    """
    set -euo pipefail

    # Reference-bias gate: emit the reference base at a monomorphic site only with enough
    # base-supporting depth (RO+AO) and low alt fraction, else N. Disable via consensus_ref_bias_gate=false.
    REF_MIN_DP=${params.consensus_ref_min_dp}
    MAX_REF_ALTFRAC=${params.consensus_max_ref_altfrac}
    MAX_REF_DEL_FRAC=${params.consensus_max_ref_del_frac}
    if [[ "${params.consensus_ref_bias_gate}" != "true" ]]; then
        REF_MIN_DP=0
        MAX_REF_ALTFRAC=1.0
        MAX_REF_DEL_FRAC=1.0
    fi

    # Virgin (unmasked) pass: drop the region-exclude and mask-sites so nothing is masked (X).
    EXCL="${exclude_txt}"
    MASK="${mask_sites}"
    if [ "${outLabel}" = ".raw" ]; then : > empty_mask.tsv; EXCL=empty_mask.tsv; MASK=empty_mask.tsv; fi

    # Run the consensus generation script
    # This script merges the reference, the backbone (all positions), and the variants
    # while masking low-confidence areas.
    python3 ${projectDir}/bin/WGS_fasta_allpos.py \\
      --vcf ${allpos_vcf_gz} \\
      --reference ${ref_fa} \\
      --exclude \$EXCL \\
      --mask-sites \$MASK \\
      --min-dp ${params.consensus_min_dp} \\
      --ref-min-dp \$REF_MIN_DP \\
      --max-ref-altfrac \$MAX_REF_ALTFRAC \\
      --max-ref-min-alt ${params.consensus_max_ref_min_alt} \\
      --max-ref-del-frac \$MAX_REF_DEL_FRAC \\
      --output ${sampleId}.${refId}${outLabel}.consensus.fasta \\
      --wrap ${params.consensus_wrap} \\
      --mask-char ${params.consensus_mask_char} \\
      --nocall-char ${params.consensus_nocall_char} \\
      > ${sampleId}.${refId}${outLabel}.consensus.log 2>&1

    # Rename the Fasta header to match the Sample ID (instead of the reference ID)
    sed -i "s/^>.*/>${sampleId}${outLabel}/" ${sampleId}.${refId}${outLabel}.consensus.fasta
    """
}
