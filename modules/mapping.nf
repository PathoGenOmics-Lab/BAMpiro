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
    // publishDir "${params.outdir}/${getSampleDir(sampleId, params)}", mode: params.publish_mode, saveAs: { filename -> getSavePath(filename, params) }
    
    cpus { params.threads as int }
    // A 4.4 Mbp bacterial index is tiny; 6 GB is generous. Escalates on the rare OOM retry.
    memory { 6.GB * task.attempt }

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
      samtools sort -@ !{task.cpus} -l 1 -o !{sampleId}__!{runId}.pe.sorted.bam -
    '''
}

process MAPPING_SE {
    tag "Map SE: ${sampleId}"
    // Intermediate BAMs are usually not published
    // publishDir "${params.outdir}/${getSampleDir(sampleId, params)}", mode: params.publish_mode, saveAs: { filename -> getSavePath(filename, params) }

    cpus { params.threads as int }
    // A 4.4 Mbp bacterial index is tiny; 6 GB is generous. Escalates on the rare OOM retry.
    memory { 6.GB * task.attempt }

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
      samtools sort -@ !{task.cpus} -l 1 -o !{sampleId}__!{runId}.se.sorted.bam -
    '''
}

process MERGE_AND_MARKDUP {
    tag "Dedup: ${sampleId}"
    
    // Use getSampleDir for nested output support
    // This handles putting .stats in the stats/ folder and keeping the final .bam in the sample root
    publishDir "${params.outdir}/${getSampleDir(sampleId, params)}", mode: params.publish_mode, saveAs: { filename -> getSavePath(filename, params) }
    
    // --- OOM (Out of Memory) Protection Strategy ---
    // If the process fails with exit code 137 (OOM), it retries with more memory
    cpus 4
    // Sort/markdup over a few-hundred-MB bacterial BAM; 8 GB is plenty and doubles on OOM retry.
    memory { 8.GB * task.attempt }
    time { 6.h * task.attempt }
    
    errorStrategy { task.exitStatus in [137, 140, 143, 134, 139] ? 'retry' : 'finish' }
    maxRetries 3
    // -----------------------------------------------

    input:
    tuple val(sampleId), val(refId), path(bams), path(ref_fa), path(exclude_txt)

    output:
    tuple val(sampleId), val(refId), path("${sampleId}.${refId}.final.bam"), path("${sampleId}.${refId}.final.bam.bai"), path(ref_fa), path(exclude_txt), emit: final_bam
    tuple val(sampleId), val(refId), path("${sampleId}.${refId}.dedup.stats"), emit: stats
    
    shell:
    '''
    set -euo pipefail
    
    # Capture the list of input BAM files
    set -- !{bams}

    # 1. Merge Logic
    # One BAM -> use it in place (no full-file copy). More -> merge them.
    if [[ $# -eq 1 ]]; then
      MERGED="$1"
    else
      samtools merge -@ !{task.cpus} -f merged.bam "$@"
      MERGED=merged.bam
    fi

    # 2. Sorting & Deduplication Pipeline (streamed -- no intermediate BAMs on disk)
    # Name Sort -> Fixmate -> Coord Sort -> Markdup, piped end to end. The sort stages emit
    # uncompressed BAM (-l 0) so we don't (de)compress throwaway intermediates or round-trip
    # them through disk. fixmate -m keeps the MC/ms tags the length-aware read filter needs.

    # samtools sort per-thread memory: 80% of the allocation split across threads
    mem_per_thread=$(python3 -c "print(int(!{task.memory.toMega()} * 0.8 / !{task.cpus}))")M

    samtools sort -n -@ !{task.cpus} -m $mem_per_thread -l 0 "$MERGED" \
      | samtools fixmate -m - - \
      | samtools sort -@ !{task.cpus} -m $mem_per_thread -l 0 - \
      | samtools markdup -r -@ !{task.cpus} - !{sampleId}.!{refId}.final.bam

    # 3. Final Indexing and Stats
    samtools index !{sampleId}.!{refId}.final.bam
    samtools stats !{sampleId}.!{refId}.final.bam > !{sampleId}.!{refId}.dedup.stats

    # 4. Cleanup (only the merge intermediate exists, and only in the multi-BAM case)
    rm -f merged.bam
    '''
}

process FILTER_READS {
    // Length-aware per-read-PAIR mappability filter. Drops reads too short to be uniquely
    // placed at their locus (short-read false positives in near-repeats); a fragment is kept
    // if EITHER mate anchors uniquely (concordant-pair rescue). One streaming pass, no re-sort.
    tag "Filter: ${sampleId}"

    publishDir "${params.outdir}/${getSampleDir(sampleId, params)}", mode: params.publish_mode, saveAs: { filename -> getSavePath(filename, params) }

    cpus 4
    memory { 4.GB * task.attempt }

    input:
    // final_bam tuple + the per-reference mappability track and repeat BED (joined by refId in main.nf)
    tuple val(sampleId), val(refId), path(bam), path(bai), path(ref_fa), path(exclude_txt), path(mul_npz), path(repeat_bed)

    output:
    // Same 6-field shape as MERGE_AND_MARKDUP.final_bam, but exclude_txt is now nucmer + genmap
    // repeat regions combined -> the always-repetitive interior is masked (X) downstream.
    tuple val(sampleId), val(refId),
          path("${sampleId}.${refId}.filtered.bam"), path("${sampleId}.${refId}.filtered.bam.bai"),
          path(ref_fa), path("${sampleId}.${refId}.exclude.txt"), emit: filtered_bam
    path("${sampleId}.${refId}.filter_mqc.tsv"), emit: stats

    shell:
    '''
    set -euo pipefail

    # Combine the nucmer repeat exclude with the genmap always-repetitive interior (-> X in consensus)
    if [[ "!{params.mappability_mask_consensus}" == "true" ]]; then
        cat !{exclude_txt} !{repeat_bed} > !{sampleId}.!{refId}.exclude.txt
    else
        cp !{exclude_txt} !{sampleId}.!{refId}.exclude.txt
    fi

    KU=""
    if [[ "!{params.filter_keep_unmapped}" == "true" ]]; then KU="--keep-unmapped"; fi
    SC=""
    if [[ "!{params.filter_strict_contigs}" == "true" ]]; then SC="--strict-contigs"; fi

    # Drops-only over the coordinate-sorted BAM -> order preserved -> just re-index (no sort)
    samtools view -h -@ !{task.cpus} !{bam} \
      | python3 !{projectDir}/bin/filter_reads_mappability.py \
          --mul !{mul_npz} --kmin !{params.genmap_min_k} --sentinel !{params.genmap_infinity} $KU $SC \
      | samtools view -b -@ !{task.cpus} -o !{sampleId}.!{refId}.filtered.bam -
    samtools index -@ !{task.cpus} !{sampleId}.!{refId}.filtered.bam

    a=$(samtools view -c !{bam})
    b=$(samtools view -c !{sampleId}.!{refId}.filtered.bam)
    d=$((a-b))
    pct=$(python3 -c "print(f'{100*$d/$a:.2f}') if $a else print('0.00')")

    # QC: MultiQC custom-content table (one row per sample; rows merge under one section)
    {
      printf '# id: read_filter\n'
      printf '# section_name: Length-aware read filter\n'
      printf '# description: Reads dropped as unable to map uniquely at their locus (near-repeat short reads).\n'
      printf '# plot_type: table\n'
      printf 'Sample\tinput_reads\tkept_reads\tdropped_reads\tdropped_pct\n'
      printf '%s\t%s\t%s\t%s\t%s\n' "!{sampleId}.!{refId}" "$a" "$b" "$d" "$pct"
    } > !{sampleId}.!{refId}.filter_mqc.tsv

    # Flag samples that drop an abnormal fraction of reads (contamination / short reads / wrong ref)
    if [ "$(python3 -c "print(1 if $pct > !{params.filter_max_drop_pct} else 0)")" = "1" ]; then
        echo "[FILTER_READS] WARNING: !{sampleId}.!{refId} dropped ${pct}% of reads (> !{params.filter_max_drop_pct}%) -- check contamination / read length / reference" >&2
    fi
    '''
}
