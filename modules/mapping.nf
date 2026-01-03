nextflow.enable.dsl=2

// Import centralized functions for path generation and file classification
include { getSavePath; getSampleDir } from './utils'

/* ====================================================================
    MAPPING MODULES
    Contains: BWA Alignment and BAM Processing (Merge/Sort/MarkDup)
==================================================================== */

process MAPPING_PE {
    tag "Map PE: ${sampleId}"
    // We typically don't publish these intermediate BAMs to save space.
    // If you wanted to publish them, you would use:
    // publishDir "${params.outdir}/${getSampleDir(sampleId, params)}", mode: 'copy', saveAs: { filename -> getSavePath(filename, params) }
    
    cpus { params.threads as int }
    memory '32 GB'

    input:
    tuple val(refId), val(sampleId), val(runId), path(r1), path(r2), val(taxId), path(ref_fa), path(indices), path(exclude_txt)

    output:
    // We pass ref_fa and exclude_txt to keep the channel context for the next step
    tuple val(sampleId), val(refId), path("${sampleId}__${runId}.pe.sorted.bam"), path(ref_fa), path(exclude_txt), emit: bam

    shell:
    '''
    set -euo pipefail
    
    # Align PE reads using BWA-MEM2
    # -R: Adds Read Group header (Critical for GATK/FreeBayes)
    bwa-mem2 mem -t !{task.cpus} \
      -R "@RG\\tID:!{sampleId}.!{runId}.PE\\tSM:!{sampleId}\\tPL:ILLUMINA" \
      !{ref_fa} !{r1} !{r2} | \
      samtools sort -@ !{task.cpus} -o !{sampleId}__!{runId}.pe.sorted.bam -
    '''
}

process MAPPING_SE {
    tag "Map SE: ${sampleId}"
    // Intermediate BAMs are usually not published
    // publishDir "${params.outdir}/${getSampleDir(sampleId, params)}", mode: 'copy', saveAs: { filename -> getSavePath(filename, params) }

    cpus { params.threads as int }
    memory '32 GB'

    input:
    tuple val(refId), val(sampleId), val(runId), path(reads_se), val(taxId), path(ref_fa), path(indices), path(exclude_txt)

    output:
    tuple val(sampleId), val(refId), path("${sampleId}__${runId}.se.sorted.bam"), path(ref_fa), path(exclude_txt), emit: bam

    shell:
    '''
    set -euo pipefail
    
    # Align SE reads using BWA-MEM2
    bwa-mem2 mem -t !{task.cpus} \
      -R "@RG\\tID:!{sampleId}.!{runId}.SE\\tSM:!{sampleId}\\tPL:ILLUMINA" \
      !{ref_fa} !{reads_se} | \
      samtools sort -@ !{task.cpus} -o !{sampleId}__!{runId}.se.sorted.bam -
    '''
}

process MERGE_AND_MARKDUP {
    tag "Dedup: ${sampleId}"
    
    // Use getSampleDir for nested output support
    // This handles putting .stats in the stats/ folder and keeping the final .bam in the sample root
    publishDir "${params.outdir}/${getSampleDir(sampleId, params)}", mode: 'copy', saveAs: { filename -> getSavePath(filename, params) }
    
    // --- OOM (Out of Memory) Protection Strategy ---
    // If the process fails with exit code 137 (OOM), it retries with more memory
    cpus 4
    memory { 32.GB * task.attempt }
    time { 6.h * task.attempt }
    
    errorStrategy { task.exitStatus in [137, 140, 143, 134, 139] ? 'retry' : 'finish' }
    maxRetries 3
    // -----------------------------------------------

    input:
    tuple val(sampleId), val(refId), path(bams), path(ref_fa), path(exclude_txt)

    output:
    tuple val(sampleId), val(refId), path("${sampleId}.${refId}.final.bam"), path("${sampleId}.${refId}.final.bam.bai"), path(ref_fa), path(exclude_txt), emit: final_bam
    path("${sampleId}.${refId}.dedup.stats"), emit: stats
    
    shell:
    '''
    set -euo pipefail
    
    # Capture the list of input BAM files
    set -- !{bams}
    
    # 1. Merge Logic
    # If there is only one BAM file, just rename it. If more, merge them.
    if [[ $# -eq 1 ]]; then
      cp "$1" merged.bam
    else
      samtools merge -@ !{task.cpus} -f merged.bam "$@"
    fi

    # 2. Sorting & Deduplication Pipeline
    # Pipeline: Name Sort -> Fixmate -> Coord Sort -> Markdup
    
    # Calculate memory limit for samtools sort to prevent OOM
    # Uses 80% of allocated memory divided by CPUs
    mem_per_thread=$(python3 -c "import math; print(int(!{task.memory.toMega()} * 0.8 / !{task.cpus}))")M

    # A. Sort by Name (Required for Fixmate)
    samtools sort -n -@ !{task.cpus} -m $mem_per_thread -o merged.name.bam merged.bam
    
    # B. Fixmate (Fills in mate coordinates/flags)
    samtools fixmate -m merged.name.bam merged.fixmate.bam
    
    # C. Sort by Coordinate (Required for Markdup)
    samtools sort -@ !{task.cpus} -m $mem_per_thread -o merged.coord.bam merged.fixmate.bam
    
    # D. Mark Duplicates (Final step)
    samtools markdup -r -@ !{task.cpus} merged.coord.bam !{sampleId}.!{refId}.final.bam
    
    # 3. Final Indexing and Stats
    samtools index !{sampleId}.!{refId}.final.bam
    samtools stats !{sampleId}.!{refId}.final.bam > !{sampleId}.!{refId}.dedup.stats
    
    # 4. Cleanup intermediate files
    rm -f merged.bam merged.name.bam merged.fixmate.bam merged.coord.bam
    '''
}
