nextflow.enable.dsl=2

/* ====================================================================
    REFERENCE MODULES
    Contains: Reference indexing/preparation and SnpEff DB building
==================================================================== */

process PREPARE_REFERENCE {

    tag "Ref: ${refId}"

    // We keep custom logic here because we NEED to publish the index files (.fai, .bwt, etc.),
    // which the global 'getSavePath' function would filter out.
    publishDir path: { "${params.outdir}/references/${refId}" }, mode: params.publish_mode, saveAs: { filename ->
        // Hide the raw copy of reference.fa to save space (the original input already exists),
        // BUT keep it when publishing CRAM so the outputs are self-decodable.
        if (filename == "reference.fa" && !params.output_cram) return null
        return filename
    }

    cpus 4
    memory '8 GB'

    input:
    tuple val(refId), path(fasta)

    output:
    tuple val(refId), path("reference.fa"), path("reference.fa.*"), path("Locus_to_exclude_${refId}.txt"), emit: bundle

    shell:
    '''
    set -euo pipefail

    # 1) Handle Input Fasta (Gunzip if necessary)
    if [[ "!{fasta}" == *.gz ]]; then
        gzip -cd "!{fasta}" > reference.fa
    else
        cp "!{fasta}" reference.fa
    fi

    # 2) Generate Indexes
    bwa-mem2 index reference.fa
    samtools faidx reference.fa

    # 3) Generate Exclusion File (Repeat Masking via self-alignment)
    out="Locus_to_exclude_!{refId}.txt"
    echo -e "Chrom\tStart\tEnd\tTag\tComment" > "$out"

    if [[ "!{params.exclude_repeats}" == "true" ]]; then

        nucmer --maxmatch --nosimplify reference.fa reference.fa -p self_aln
        show-coords -r -c -l -T self_aln.delta > coords.tsv

        # Parse nucmer coordinates to identify repetitive regions
        awk -v RID="!{refId}" 'BEGIN{OFS="\t"}
        {
            # Skip header/blank lines and enforce numeric coords in columns 1..4
            if ($1=="" || $2=="" || $3=="" || $4=="") next
            if ($1 ~ /[^0-9]/ || $2 ~ /[^0-9]/ || $3 ~ /[^0-9]/ || $4 ~ /[^0-9]/) next

            # Normalize coordinates so Start <= End
            s1=$1+0; e1=$2+0; if (s1>e1){tmp=s1; s1=e1; e1=tmp}
            s2=$3+0; e2=$4+0; if (s2>e2){tmp=s2; s2=e2; e2=tmp}

            # TAGS: usually last 2 columns (TAGS1 TAGS2). If only one tag column exists, reuse it.
            t1=""; t2=""
            if (NF>=13) { t1=$(NF-1); t2=$NF }
            else if (NF==12) { t1=$12; t2=$12 }

            if (t1=="" || t1==".") t1=RID
            if (t2=="" || t2==".") t2=RID

            # Drop exact self-diagonal hits: same contig AND same interval
            if (t1==t2 && s1==s2 && e1==e2) next

            # Emit BOTH sides of the alignment
            print t1, s1, e1, "", ""
            print t2, s2, e2, "", ""
        }' coords.tsv \
        | sort -k1,1 -k2,2n \
        | awk 'BEGIN{OFS="\t"}
               NR==1 { c=$1; cs=$2; ce=$3; next }
               {
                   if ($1==c && $2<=ce+1) {
                       if ($3>ce) ce=$3
                   } else {
                       print c, cs, ce, "", ""
                       c=$1; cs=$2; ce=$3
                   }
               }
               END { if (NR>0) print c, cs, ce, "", "" }' \
        >> "$out"
    fi

    # H37Rv Illumina "blind spots" (Zenodo 3701840): add them to the exclusion. The BED is 0-based
    # half-open in NC_000962.3 coords; rewrite to the exclusion's 1-based-inclusive coords (start+1, end)
    # and to THIS reference's contig name. Opt-in and only correct for an H37Rv-coordinate reference
    # (H37Rv / the MTBC ancestor share coordinates). Flows to variant calling, consensus and the report.
    if [[ "!{params.mask_blindspots}" == "true" && -s "!{params.blindspot_bed}" ]]; then
        CONTIG=$(head -1 reference.fa | sed 's/^>//; s/[[:space:]].*//')
        if [[ "!{params.blindspot_liftover}" == "true" && -s "!{params.canonical_ref}" ]]; then
            # Whole-genome anchor-chain liftover of the H37Rv blind-spots onto THIS reference, so the mask is
            # correct without assuming shared coordinates: every position is placed by interpolation between
            # flanking unique anchors (handles SNP sites + indels; anchor-desert positions drop, never mis-map).
            python3 !{projectDir}/bin/pathotypr_liftover.py lift "!{params.blindspot_bed}" "!{params.canonical_ref}" reference.fa \
                --out-bed bs_lifted.bed --contig "$CONTIG" --kmer-size 21 --global-chain --sample 10
            awk 'BEGIN{OFS="\t"} $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ {print $1, $2+1, $3, "blindspot", ""}' bs_lifted.bed >> "$out"
        else
            # reference already shares H37Rv coordinates: append the blind-spots directly (contig rewrite + 1-based)
            awk -v C="$CONTIG" 'BEGIN{OFS="\t"} $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ {print C, $2+1, $3, "blindspot", ""}' "!{params.blindspot_bed}" >> "$out"
        fi
    fi
    '''

    stub:
    """
    touch reference.fa
    # main.nf picks the .fai out of this index list, so it has to be part of the glob
    touch reference.fa.amb reference.fa.ann reference.fa.bwt.2bit.64 reference.fa.pac reference.fa.0123 reference.fa.fai
    touch Locus_to_exclude_${refId}.txt
    """
}


process BUILD_MAPPABILITY {
    // Length-aware mappability track: for each position, the smallest read length at
    // which the k-mer starting there is genome-unique (min_unique_len). Feeds FILTER_READS.
    tag "Mappability: ${refId}"

    // storeDir persists the track and SKIPS this (expensive) step if it already exists -> computed
    // once per reference across runs. Keyed by refId + the genmap params + the reference size, so
    // swapping the FASTA for a given refId (different genome / added-removed contigs) recomputes it
    // instead of silently reusing a stale track. NOTE: a same-byte-size edit (e.g. a single-base
    // substitution) won't change the key -- clear the dir manually in that rare case.
    storeDir "${params.mappability_dir}/${refId}_k${params.genmap_min_k}-${params.genmap_max_k}_s${params.genmap_step}_E${params.genmap_errors}r${params.genmap_error_rate}_tp${params.genmap_tail_policy}_inf${params.genmap_infinity}_sz${ref_fa.size()}"

    cpus 8
    memory '16 GB'

    input:
    tuple val(refId), path(ref_fa)

    output:
    tuple val(refId), path("${refId}.min_unique_len.npz"), emit: track
    tuple val(refId), path("Locus_to_exclude_mappability_${refId}.txt"), emit: repeat_bed

    shell:
    '''
    set -euo pipefail
    samtools faidx !{ref_fa}

    # genmap refuses to write into an existing index dir (matters on -resume / retry)
    rm -rf gmidx
    genmap index -F !{ref_fa} -I gmidx

    # (k,E)-mappability sweep across candidate read lengths; value 1 == unique (both strands).
    # E is fixed (genmap_errors) unless genmap_error_rate>0, in which case it scales with k.
    BGS=""
    for K in $(seq !{params.genmap_min_k} !{params.genmap_step} !{params.genmap_max_k}); do
        if [ "!{params.genmap_error_rate}" != "0" ]; then
            E=$(python3 -c "print(max(1, round(${K} * !{params.genmap_error_rate})))")
        else
            E=!{params.genmap_errors}
        fi
        genmap map -K ${K} -E ${E} -T !{task.cpus} -I gmidx -O map_K${K} -bg
        BGS="${BGS} ${K}:map_K${K}.bedgraph"
    done

    # Collapse the sweep into one uint16 per-contig min_unique_len array (65535 = never unique),
    # and emit the always-repetitive interior as a 1-based BED for consensus masking.
    python3 !{projectDir}/bin/build_min_unique_len.py \
        --fai !{ref_fa}.fai \
        --bedgraphs ${BGS} \
        --sentinel !{params.genmap_infinity} \
        --tail-policy !{params.genmap_tail_policy} \
        --mask-bed Locus_to_exclude_mappability_!{refId}.txt \
        --mask-window !{params.genmap_max_k} \
        --out-prefix !{refId}
    '''

    stub:
    """
    touch ${refId}.min_unique_len.npz
    touch Locus_to_exclude_mappability_${refId}.txt
    """
}


process SNPEFF_BUILD_DB {
    tag "SnpEff DB: ${refId}"
    publishDir path: { "${params.outdir}/references/${refId}/snpeff" }, mode: params.publish_mode
    cpus 1
    memory '8 GB'
    
    input:
    tuple val(refId), path(fasta), path(gff)
    
    output:
    tuple val(refId), path("snpEff.config"), path("data"), emit: db
    
    shell:
    '''
    set -euo pipefail
    
    # Create SnpEff directory structure
    mkdir -p data/!{refId}
    
    # Handle Fasta (Gunzip or Copy)
    if [[ "!{fasta}" == *.gz ]]; then 
        gzip -cd "!{fasta}" > data/!{refId}/sequences.fa
    else 
        cp "!{fasta}" data/!{refId}/sequences.fa
    fi

    # Handle GFF (Gunzip or Copy)
    if [[ "!{gff}"   == *.gz ]]; then 
        gzip -cd "!{gff}"   > data/!{refId}/genes.gff
    else 
        cp "!{gff}"   data/!{refId}/genes.gff
    fi
    
    # Generate Config File
    cat > snpEff.config <<EOF
data.dir = ./data
!{refId}.genome : !{refId}
EOF
    
    # Build Database
    snpEff build -c snpEff.config -gff3 -noCheckCds -noCheckProtein -v !{refId}
    '''

    stub:
    """
    touch snpEff.config
    mkdir -p data/${refId}
    """
}
