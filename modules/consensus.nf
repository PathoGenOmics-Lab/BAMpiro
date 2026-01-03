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
    tuple val(sampleId), val(refId), path(allpos_vcf_gz), path(allpos_tbi), path(mask_sites), path(ref_fa), path(exclude_txt)

    output:
    tuple val(sampleId), val(refId), path("${sampleId}.${refId}.consensus.fasta"), emit: fasta
    path("${sampleId}.${refId}.consensus.log"), emit: log

    script:
    """
    set -euo pipefail

    # Run the consensus generation script
    # This script merges the reference, the backbone (all positions), and the variants
    # while masking low-confidence areas.
    python3 ${projectDir}/bin/WGS_fasta_allpos.py \\
      --vcf ${allpos_vcf_gz} \\
      --reference ${ref_fa} \\
      --exclude ${exclude_txt} \\
      --mask-sites ${mask_sites} \\
      --min-dp ${params.consensus_min_dp} \\
      --output ${sampleId}.${refId}.consensus.fasta \\
      --wrap ${params.consensus_wrap} \\
      --mask-char ${params.consensus_mask_char} \\
      --nocall-char ${params.consensus_nocall_char} \\
      > ${sampleId}.${refId}.consensus.log 2>&1

    # Rename the Fasta header to match the Sample ID (instead of the reference ID)
    sed -i "s/^>.*/>${sampleId}/" ${sampleId}.${refId}.consensus.fasta
    """
}
