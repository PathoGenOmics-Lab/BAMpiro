nextflow.enable.dsl=2

// Import centralized functions for path generation and file classification
include { getSavePath; getSampleDir } from './utils'

/* ====================================================================
    PATHOTYPR MODULES
    Alignment-free (k-mer) MTBC lineage + WHO drug-resistance genotyping.
    Reference-agnostic: the diagnostic k-mers are built from the MTBC-ancestor
    reference the markers are defined on, so calls do NOT depend on the reference
    the sample was mapped against. Two split-fastq passes per sample:
      - lineage markers (nested sub-lineage) -> emit: summary   (${prefix}_summary.tsv:
        genome / lineage:count / major_lineage; parsed by stats_to_legacy.py)
      - WHO drug-resistance markers          -> emit: dr_mutations (${prefix}_<s>_mutations.tsv:
        pos/ref/alt/ref_count/alt_count/alt_fraction/lineage_path[=drug;resistance;marker;grade;gene;mutation];
        aggregated by collect_dr.py for the report's Drug-resistance panel)
    pathotypr v1.0.2 CLI verified against the built image: split-fastq -i --reference --markers
    --output-prefix --paired --nested-classification --min-alt-percent --threads.
==================================================================== */

process RUN_PATHOTYPR_PE {
    tag "PathotyprPE: ${sampleId}"
    publishDir path: { "${params.outdir}/${getSampleDir(sampleId, params)}" }, mode: params.publish_mode, saveAs: { filename -> getSavePath(filename, params) }
    cpus 4
    memory '8 GB'

    input:
    tuple val(refId), val(sampleId), val(runId), path(r1), path(r2), path(orphans)
    // val, not path, and that is the whole point. These three live INSIDE the image, bundled
    // from Zenodo at build time. Declaring them as `path` tells Nextflow to stage them FROM THE
    // HOST, so it adds a bind for their parent and Singularity refuses:
    //
    //   mount source /opt/pathotypr doesn't exist
    //
    // As values they are interpolated into the script unchanged and resolved inside the
    // container, where they do exist. An override pointing at your own markers has to live under
    // a path the profile binds, since nothing auto-mounts a value.
    val ref_fasta_pathotypr
    val lineage_markers
    val dr_markers
    val pathotypr_bin

    output:
    // runId in the file names so a multi-lane (or multi-reference) sample does not emit identically-named
    // files that collide when COLLECT_DR .collect()s them; the tuple key stays sampleId for the joins.
    tuple val(sampleId), path("${sampleId}__${runId}.pathotypr.lineage_summary.tsv"), emit: summary
    tuple val(sampleId), path("${sampleId}__${runId}.dr_mutations.tsv"),              emit: dr_mutations
    path("${sampleId}__${runId}.pathotypr.*"), emit: results

    script:
    """
    set -euo pipefail
    # Everything the sample has, in one file. FastP is run with --merge, so on a library whose
    # pairs overlap MOST of the reads leave as merged fragments in the orphan stream rather than
    # as r1/r2. Typing only r1 and r2 therefore typed a small and variable slice of each sample:
    # 8-19% of the reads in the cohort that found this. Mapping had all three streams all along.
    #
    # Concatenated rather than passed as three -i arguments, because pathotypr refuses those:
    # --paired wants an even number of files, and without it the three names collapse to one
    # duplicate sample. Concatenating is sound in both directions that could bite here, and both
    # were measured against this binary rather than assumed:
    #   - gzip members: cat produces a multi-member file, and a reader using GzDecoder rather than
    #     MultiGzDecoder would take the first member and report no error. Feeding r1 alone is below
    #     --min-depth and yields no call at all, while the concatenation yields the same counts as
    #     --paired, which it could not do if it were reading r1 only.
    #   - pairing: --paired only groups files into a sample. The same pair passed as one
    #     concatenated single-end file gives byte-identical ref_count, alt_count and alt_fraction.
    cat ${r1} ${r2} ${orphans} > ${sampleId}__${runId}.all.fq.gz
    # --no-auto-paired so a sampleId that happens to contain _R1 or _1 is not read as half of a pair.
    # Lineage (nested sub-lineage classification) -> \${prefix}_summary.tsv
    ${pathotypr_bin} split-fastq -i ${sampleId}__${runId}.all.fq.gz --no-auto-paired \\
        --reference ${ref_fasta_pathotypr} --markers ${lineage_markers} \\
        --nested-classification --output-prefix ${sampleId}__${runId}.pathotypr.lineage --threads ${task.cpus}
    # Drug resistance (WHO catalogue markers) -> \${prefix}_<sample>_mutations.tsv
    ${pathotypr_bin} split-fastq -i ${sampleId}__${runId}.all.fq.gz --no-auto-paired \\
        --reference ${ref_fasta_pathotypr} --markers ${dr_markers} \\
        --output-prefix ${sampleId}__${runId}.pathotypr.dr --min-alt-percent ${params.pathotypr_min_alt} --threads ${task.cpus}
    rm -f ${sampleId}__${runId}.all.fq.gz
    # Fix the DR detail file name to carry our sampleId+runId (pathotypr names the sample after the FASTQ)
    cp "\$(ls ${sampleId}__${runId}.pathotypr.dr_*_mutations.tsv | head -1)" ${sampleId}__${runId}.dr_mutations.tsv
    """

    // The lineage summary also satisfies the pathotypr.* glob (emit: results), so two files cover all three outputs.
    stub:
    """
    touch ${sampleId}__${runId}.pathotypr.lineage_summary.tsv
    touch ${sampleId}__${runId}.dr_mutations.tsv
    """
}

process RUN_PATHOTYPR_SE {
    tag "PathotyprSE: ${sampleId}"
    publishDir path: { "${params.outdir}/${getSampleDir(sampleId, params)}" }, mode: params.publish_mode, saveAs: { filename -> getSavePath(filename, params) }
    cpus 4
    memory '8 GB'

    input:
    tuple val(refId), val(sampleId), val(runId), path(r1)
    // val, not path, and that is the whole point. These three live INSIDE the image, bundled
    // from Zenodo at build time. Declaring them as `path` tells Nextflow to stage them FROM THE
    // HOST, so it adds a bind for their parent and Singularity refuses:
    //
    //   mount source /opt/pathotypr doesn't exist
    //
    // As values they are interpolated into the script unchanged and resolved inside the
    // container, where they do exist. An override pointing at your own markers has to live under
    // a path the profile binds, since nothing auto-mounts a value.
    val ref_fasta_pathotypr
    val lineage_markers
    val dr_markers
    val pathotypr_bin

    output:
    // runId in the file names (see RUN_PATHOTYPR_PE) so multi-lane / multi-ref samples don't collide in COLLECT_DR
    tuple val(sampleId), path("${sampleId}__${runId}.pathotypr.lineage_summary.tsv"), emit: summary
    tuple val(sampleId), path("${sampleId}__${runId}.dr_mutations.tsv"),              emit: dr_mutations
    path("${sampleId}__${runId}.pathotypr.*"), emit: results

    script:
    """
    set -euo pipefail
    ${pathotypr_bin} split-fastq -i ${r1} \\
        --reference ${ref_fasta_pathotypr} --markers ${lineage_markers} \\
        --nested-classification --output-prefix ${sampleId}__${runId}.pathotypr.lineage --threads ${task.cpus}
    ${pathotypr_bin} split-fastq -i ${r1} \\
        --reference ${ref_fasta_pathotypr} --markers ${dr_markers} \\
        --output-prefix ${sampleId}__${runId}.pathotypr.dr --min-alt-percent ${params.pathotypr_min_alt} --threads ${task.cpus}
    cp "\$(ls ${sampleId}__${runId}.pathotypr.dr_*_mutations.tsv | head -1)" ${sampleId}__${runId}.dr_mutations.tsv
    """

    stub:
    """
    touch ${sampleId}__${runId}.pathotypr.lineage_summary.tsv
    touch ${sampleId}__${runId}.dr_mutations.tsv
    """
}
