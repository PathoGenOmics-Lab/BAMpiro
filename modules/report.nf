include { asBool; getSampleDir; getSavePath; refPairs } from './utils'
nextflow.enable.dsl=2

/* ====================================================================
    REPORT MODULES
    Cohort-level consolidated QC report (ported from sBAMpiro):
      COLLECT_SUMMARY -> one summary TSV from all per-sample legacy logs
      QC_REPORT       -> self-contained interactive HTML + qc_flags.tsv
==================================================================== */

process COLLECT_SUMMARY {
    tag "Summary"
    publishDir "${params.outdir}", mode: params.publish_mode
    cpus 1
    memory '2 GB'

    input:
    path(logs)              // every ${sampleId}.log (collected)
    val(basename)

    output:
    path("${basename}_summary.tsv"),     emit: summary
    path("${basename}_gene_burden.tsv"), emit: gene_burden

    script:
    """
    set -euo pipefail
    python3 ${projectDir}/bin/collect_summary.py \\
        -o ${basename}_summary.tsv \\
        --gene-burden-out ${basename}_gene_burden.tsv \\
        ${logs}
    # collect_summary only writes the gene-burden file when snpEff ANN is present in the logs;
    # guarantee it exists (header only) so the channel/QC_REPORT always have a file.
    [ -f ${basename}_gene_burden.tsv ] || printf 'gene\\thigh\\tmoderate\\tlow\\tmodifier\\tdominant_effect\\tn_samples\\ttotal_impactful\\tstart\\tend\\n' > ${basename}_gene_burden.tsv
    """

    stub:
    """
    touch ${basename}_summary.tsv
    touch ${basename}_gene_burden.tsv
    """
}

process COLLECT_DR {
    tag "DR calls"
    publishDir "${params.outdir}", mode: params.publish_mode
    cpus 1
    memory '2 GB'

    input:
    path(dr_mutations)      // every sample's ${sampleId}__${runId}.dr_mutations.tsv (pathotypr DR run)
    val(basename)

    output:
    path("${basename}_dr.tsv"), emit: dr

    script:
    """
    set -euo pipefail
    python3 ${projectDir}/bin/collect_dr.py -o ${basename}_dr.tsv ${dr_mutations}
    """

    stub:
    """
    touch ${basename}_dr.tsv
    """
}

process LIFT_VARIANTS {
    // Whole-genome anchor-chain liftover of ONE reference's variant positions onto the canonical (H37Rv-
    // colinear) reference -> contig-aware map (src_contig src_pos tgt_contig tgt_pos strand) for the report's
    // SNP tables (--pos-liftover) and for ANNOTATE_CANONICAL. One task per reference: every reference is lifted
    // from its own sequence, and only the positions on its own contigs. Lifting a second reference's positions
    // with the first one's sequence put 96% of them at a wrong coordinate. Alignment-free, no shared-coordinate
    // assumption; every position (incl. SNP sites) is placed by interpolation between flanking unique anchors,
    // on each of the reference's contigs -> correct or absent. Where it cannot interpolate (an indel between two
    // anchors, or no anchor for longer than 2k) the stretch between the anchors is aligned (--align-gaps), and a
    // position the canonical genome lacks is written as such. Pure Python (bin/pathotypr_liftover.py lift
    // --global-chain), about 25 s and 2.5 GB for a 4.4 Mb genome.
    tag "LiftVariants: ${refId}"
    cpus 1
    memory { 6.GB * task.attempt }

    input:
    tuple val(refId), path(ref_fa)   // the mapping reference FASTA (k-mer context source)
    path(vcfs)                       // the report's per-sample VCFs, of every reference (may be .gz)

    output:
    path("${refId}.pos_liftover.tsv"), emit: map

    script:
    """
    set -euo pipefail
    # This reference's contigs, and the variant positions on them. A VCF of another reference names none.
    grep '^>' ${ref_fa} | sed 's/^>//; s/[[:space:]].*//' | sort -u > contigs.txt
    zcat -f ${vcfs} | awk 'BEGIN{FS=OFS="\\t"} NR==FNR {c[\$1]=1; next} !/^#/ && (\$1 in c) && \$2 ~ /^[0-9]+\$/ {print \$1, \$2}' contigs.txt - \\
        | sort -u > positions.txt
    if [ -s positions.txt ]; then
        python3 ${projectDir}/bin/pathotypr_liftover.py lift positions.txt ${ref_fa} ${params.canonical_ref} \\
            --out-map ${refId}.pos_liftover.tsv --kmer-size 21 --global-chain --align-gaps
    else
        printf 'src_contig\\tsrc_pos\\ttgt_contig\\ttgt_pos\\tstrand\\n' > ${refId}.pos_liftover.tsv
    fi
    """

    stub:
    """
    printf 'src_contig\\tsrc_pos\\ttgt_contig\\ttgt_pos\\tstrand\\n' > ${refId}.pos_liftover.tsv
    """
}

process QC_REPORT {
    tag "QC report"
    publishDir "${params.outdir}", mode: params.publish_mode
    cpus 1
    // The report now holds every sample's depth summary and the distances as well; a retry after
    // running out gets more rather than failing again at the same size.
    memory { 4.GB * task.attempt }

    input:
    path(summary)
    path(gene_burden)
    path(consensus)         // all masked consensus FASTAs (may be empty)
    path(gffs, stageAs: 'ref_gff??/*')    // every reference's GFF3, in the order of ref_ids -> its genes under its samples
    path(masks, stageAs: 'ref_mask??/*')  // every reference's repeat/exclude list, same order -> masked-regions panel
    path(fais, stageAs: 'ref_fai??/*')    // every reference's FASTA index, same order -> the axis genes and masks are placed on
    val(ref_ids)                          // the reference ids, in that order
    path(metadata)          // samplesheet/metadata TSV -> SNP dynamics panel (auto-detects time+group)
    path(vcfs)              // per-sample annotated VCFs -> per-SNP allele frequencies for dynamics
    path(vcfs_h37rv)        // per-sample VCFs annotated vs H37Rv -> dual amino-acid numbering (may be NO_FILE)
    path(pos_liftover)      // mapping_pos<TAB>h37rv_pos map (pathotypr liftover) -> canonical COORDINATE per variant (may be NO_FILE)
    path(dr_report)         // run drug-resistance calls TSV (collect_dr) -> Drug resistance panel (may be NO_FILE)
    path(gene_conversion)   // cohort gene-conversion tracts TSV (COLLECT_GENE_CONVERSION) -> Gene conversion panel (may be empty)
    path(kraken_reports)    // per-sample Kraken2 .report files -> Taxonomic composition panel (may be empty)
    path(depth_profiles, stageAs: 'depth/*')   // per-sample DEPTH_PROFILE tables -> deletions + SNPs per callable kb (may be empty)
    path(snp_distances)     // SNP_DISTANCES pairs TSV -> relatedness page + GROUP_MISMATCH flag (may be NO_FILE)
    path(snp_matrix)        // SNP_MATRIX TSV -> what each series gained since its first time point (may be NO_FILE)
    path(indel_matrix)      // INDEL_MATRIX TSV -> the indel view of the matrix panel (may be NO_FILE)
    path(mnv)               // MNV_TABLE TSV -> the SNPs that are one codon-level change (may be NO_FILE)
    val(provenance)         // pre-quoted provenance tokens (container=..., reference=...)
    val(basename)

    output:
    path("${basename}_qc_report.html"), emit: html
    path("${basename}_qc_flags.tsv"),   emit: flags
    path("${basename}_deletions.tsv"),  emit: deletions, optional: true

    script:
    def gate_arg = asBool(params.report_gate) ? "--gate" : ""
    def cons_arg = consensus ? "--consensus ${consensus}" : ""
    def palette  = "${projectDir}/assets/mycolorsTB_nature.tsv"
    """
    set -euo pipefail
    # Optional inputs self-hide their panel when absent/empty (the GFF, BED and index readers are tolerant).
    LC_ARG=""; [ -f "${palette}" ] && LC_ARG="--lineage-colors ${palette}"
    MD_ARG=""; [ -s "${metadata}" ] && MD_ARG="--metadata ${metadata}"
    VCF_ARG=""; [ -n "${vcfs}" ] && VCF_ARG="--vcfs ${vcfs}"
    # NO_FILE* are the per-input placeholders main.nf stages when a feature is off (assets/NO_FILE_*).
    VH_ARG="";  case "${vcfs_h37rv}" in ""|NO_FILE*) ;; *) VH_ARG="--vcfs-h37rv ${vcfs_h37rv}";; esac
    PL_ARG="";  case "${pos_liftover}" in ""|NO_FILE*) ;; *) [ -s "${pos_liftover}" ] && PL_ARG="--pos-liftover ${pos_liftover}";; esac
    DR_ARG="";  case "${dr_report}" in ""|NO_FILE*) ;; *) [ -s "${dr_report}" ] && DR_ARG="--dr-report ${dr_report}";; esac
    # A header-only TSV (the stage ran but found nothing) parses to no tracts, so the panel self-hides.
    GC_ARG="";  case "${gene_conversion}" in ""|NO_FILE*) ;; *) [ -s "${gene_conversion}" ] && GC_ARG="--gene-conversion ${gene_conversion}";; esac
    KRK_ARG=""; [ -n "${kraken_reports}" ] && KRK_ARG="--kraken ${kraken_reports}"
    DEL_ARG=""; [ -n "${depth_profiles}" ] && DEL_ARG="--depth-profiles ${depth_profiles} --out-deletions ${basename}_deletions.tsv"
    DIST_ARG=""; case "${snp_distances}" in ""|NO_FILE*) ;; *) [ -s "${snp_distances}" ] && DIST_ARG="--snp-distances ${snp_distances}";; esac
    MX_ARG="";  case "${snp_matrix}" in ""|NO_FILE*) ;; *) [ -s "${snp_matrix}" ] && MX_ARG="--snp-matrix ${snp_matrix}";; esac
    IX_ARG="";  case "${indel_matrix}" in ""|NO_FILE*) ;; *) [ -s "${indel_matrix}" ] && IX_ARG="--indel-matrix ${indel_matrix}";; esac
    MNV_ARG=""; case "${mnv}" in ""|NO_FILE*) ;; *) [ -s "${mnv}" ] && MNV_ARG="--mnv ${mnv}";; esac
    python3 ${projectDir}/bin/qc_report.py \\
        --summary ${summary} \\
        ${cons_arg} \\
        --gene-burden ${gene_burden} \\
        --ref-gff ${refPairs(ref_ids, gffs)} \\
        --ref-mask ${refPairs(ref_ids, masks)} \\
        --ref-fai ${refPairs(ref_ids, fais)} \\
        \$LC_ARG \$MD_ARG \$VCF_ARG \$VH_ARG \$PL_ARG \$DR_ARG \$GC_ARG \$KRK_ARG \$DEL_ARG \$DIST_ARG \$MX_ARG \$IX_ARG \$MNV_ARG \\
        --cluster-snps ${params.snp_cluster_threshold} \\
        --min-dp ${params.consensus_min_dp} \\
        --deletion-min-len ${params.deletion_min_len} \\
        --deletion-min-depth ${params.report_depth_min} \\
        --aa2-label "${params.canonical_label}" \\
        --provenance ${provenance} \\
        --version "${workflow.manifest.version}" \\
        --out-html ${basename}_qc_report.html \\
        --out-flags ${basename}_qc_flags.tsv \\
        --title "BAMpiro QC report" \\
        --depth-min ${params.report_depth_min} \\
        --breadth-min ${params.report_breadth_min} \\
        --missing-max ${params.report_missing_max} \\
        --dup-max ${params.report_dup_max} \\
        --iupac-max ${params.report_iupac_max} \\
        --mapping-min ${params.report_mapping_min} \\
        ${gate_arg}
    """

    stub:
    """
    touch ${basename}_qc_report.html
    touch ${basename}_qc_flags.tsv
    """
}

process DEPTH_PROFILE {
    // What one sample's reads cover, reduced from its all-positions VCF to windows, the stretches
    // no read covers, and a per-gene table. The report reads them against the rest of the cohort,
    // which is what tells a deletion from a part of the reference none of these genomes has.
    tag "Depth profile: ${sampleId}"
    publishDir path: { "${params.outdir}/${getSampleDir(sampleId, params)}" }, mode: params.publish_mode, saveAs: { filename -> getSavePath(filename, params) }
    cpus 1
    memory { 2.GB * task.attempt }

    input:
    tuple val(sampleId), val(refId), path(allpos_vcf), path(gff)

    output:
    tuple val(sampleId), val(refId), path("${sampleId}.${refId}.depth_windows.tsv"), path("${sampleId}.${refId}.zero_depth.tsv"), emit: profile
    // Absent when the GFF has neither gene nor CDS features.
    tuple val(sampleId), val(refId), path("${sampleId}.${refId}.gene_depth.tsv"), optional: true, emit: genes

    script:
    """
    set -euo pipefail
    python3 ${projectDir}/bin/depth_profile.py \\
        --vcf ${allpos_vcf} \\
        --gff ${gff} \\
        --reference ${refId} \\
        --window ${params.depth_window} \\
        --callable-dp ${params.allpos_min_cov} \\
        --out-prefix ${sampleId}.${refId}
    """

    stub:
    """
    printf '# sample=${sampleId}\\ncontig\\tstart\\tend\\tmean_dp\\tzero_frac\\tcallable_frac\\n' > ${sampleId}.${refId}.depth_windows.tsv
    printf '# sample=${sampleId}\\ncontig\\tstart\\tend\\tlength\\n' > ${sampleId}.${refId}.zero_depth.tsv
    """
}

process SNP_DISTANCES {
    // Pairwise SNP distances between the consensus sequences: two samples differ where both called
    // a base and the bases differ, so a gap, a masked position or a mixed site never counts. Taken
    // from the consensus because it already carries every masking decision a tree built on it sees.
    tag "SNP distances"
    publishDir "${params.outdir}", mode: params.publish_mode
    cpus 1
    memory { 4.GB * task.attempt }

    input:
    path(manifest)                      // sample<TAB>reference<TAB>consensus file name
    path(fastas, stageAs: 'cons/*')     // the masked consensus FASTAs
    val(basename)

    output:
    path("${basename}_snp_distances.tsv"),       emit: square
    path("${basename}_snp_distances_pairs.tsv"), emit: pairs

    script:
    """
    set -euo pipefail
    python3 ${projectDir}/bin/snp_distances.py \\
        --manifest ${manifest} \\
        --dir cons \\
        -o ${basename}_snp_distances.tsv \\
        --pairs ${basename}_snp_distances_pairs.tsv
    """

    stub:
    """
    printf 'sample\\n' > ${basename}_snp_distances.tsv
    printf 'sample_a\\tsample_b\\treference\\tsnps\\tcompared\\tvariable_sites\\n' > ${basename}_snp_distances_pairs.tsv
    """
}

process SNP_MATRIX {
    tag "SNP matrix"
    publishDir "${params.outdir}", mode: params.publish_mode
    // Reading every sample's all-positions VCF is the cost: about a second and a half per
    // 4.4 Mb genome, read in parallel.
    cpus 4
    memory { 4.GB * task.attempt }

    input:
    path(vcfs)              // all per-sample annotated VCFs
    // Every sample's all-positions VCF (may be empty). It has a record for each reference position,
    // so a sample with no call at a site reads AF 0 at the depth it was read, where a blank also
    // meant "never read". Staged in a folder of its own: with annotate_legacy_vcfs and
    // annotate_main_vcf both off the VCFs above ARE these files, and one task cannot stage two
    // inputs under one name.
    path(depth_vcfs, stageAs: 'depth/*')
    path(mnv)               // the run's codon-level table (MNV_TABLE) -> codon_change columns (may be NO_FILE)
    val(reference)          // reference id the samples were mapped against
    val(basename)

    output:
    path("${basename}_snp_matrix.tsv"), emit: matrix

    script:
    """
    set -euo pipefail
    DEPTH_ARG=""; [ -n "${depth_vcfs}" ] && DEPTH_ARG="--depth-vcfs ${depth_vcfs}"
    MNV_ARG="";  case "${mnv}" in ""|NO_FILE*) ;; *) [ -s "${mnv}" ] && MNV_ARG="--mnv ${mnv}";; esac
    python3 ${projectDir}/bin/build_snp_matrix.py \\
        --vcfs ${vcfs} \\
        --reference "${reference}" \\
        \$DEPTH_ARG \$MNV_ARG \\
        --threads ${task.cpus} \\
        -o ${basename}_snp_matrix.tsv
    """

    stub:
    """
    touch ${basename}_snp_matrix.tsv
    """
}

process INDEL_MATRIX {
    tag "Indel matrix"
    publishDir "${params.outdir}", mode: params.publish_mode
    // Like SNP_MATRIX, the cost is reading every sample's all-positions VCF, in parallel.
    cpus 4
    memory { 4.GB * task.attempt }

    input:
    path(vcfs)                              // every sample's indel VCF (annotated when the run annotates)
    path(depth_vcfs, stageAs: 'depth/*')    // every sample's all-positions VCF: the depth where a sample calls none
    val(reference)
    val(basename)

    output:
    path("${basename}_indel_matrix.tsv"), emit: matrix

    script:
    """
    set -euo pipefail
    DEPTH_ARG=""; [ -n "${depth_vcfs}" ] && DEPTH_ARG="--depth-vcfs ${depth_vcfs}"
    python3 ${projectDir}/bin/build_indel_matrix.py \\
        --vcfs ${vcfs} \\
        --reference "${reference}" \\
        \$DEPTH_ARG \\
        --threads ${task.cpus} \\
        -o ${basename}_indel_matrix.tsv
    """

    stub:
    """
    touch ${basename}_indel_matrix.tsv
    """
}
