nextflow.enable.dsl=2

/* ====================================================================
    REFERENCE MODULES
    Contains: Reference indexing/preparation and SnpEff DB building
==================================================================== */

process PREPARE_REFERENCE {

    tag "Ref: ${refId}"

    publishDir "${params.outdir}/references/${refId}", mode: 'copy', saveAs: { filename ->
        // We do not publish the 'reference.fa' copy to save space, only indices and exclusion files
        if (filename == "reference.fa") return null
        return filename
    }

    cpus 4
    memory '16 GB'

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
    '''
}


process SNPEFF_BUILD_DB {
    tag "SnpEff DB: ${refId}"
    publishDir "${params.outdir}/references/${refId}/snpeff", mode: 'copy'
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
}
