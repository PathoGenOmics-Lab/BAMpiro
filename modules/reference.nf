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
    
    # 1. Handle Input Fasta (Gunzip if necessary)
    if [[ "!{fasta}" == *.gz ]]; then 
        gzip -cd "!{fasta}" > reference.fa
    else 
        cp "!{fasta}" reference.fa
    fi

    # 2. Generate Indexes
    bwa-mem2 index reference.fa
    samtools faidx reference.fa
    
    # 3. Generate Exclusion File (Repeat Masking)
    out="Locus_to_exclude_!{refId}.txt"
    echo -e "Chrom\tStart\tEnd\tTag\tComment" > "$out"

    if [[ "!{params.exclude_repeats}" == "true" ]]; then
      # Run Nucmer to find self-alignments
      nucmer --maxmatch --nosimplify reference.fa reference.fa -p self_aln
      show-coords -r -c -l -T self_aln.delta > coords.txt
      
      # Parse coordinates to find repeats (excluding self-identity on diagonals)
      awk -v RID="!{refId}" 'BEGIN{OFS="\t"} NR>4 {
        if ($1!="" && $2!="" && $3!="" && $4!="" && $1 !~ /[^0-9]/ && $2 !~ /[^0-9]/ && $3 !~ /[^0-9]/ && $4 !~ /[^0-9]/) {
          if ($1 != $3) {
            chrom=$12; if(chrom=="" || chrom==".") chrom=RID;
            s=$1; e=$2; if(s>e){t=s;s=e;e=t};
            print chrom, s, e, "", ""
          }
        }
      }' coords.txt | sort -k1,1 -k2,2n | awk 'BEGIN{OFS="\t"}
             NR==1{c=$1; cs=$2; ce=$3; next}
             { if($1==c && $2<=ce+1){ if($3>ce) ce=$3 } else { print c, cs, ce, "", ""; c=$1; cs=$2; ce=$3 } }
             END{ if(NR>0) print c, cs, ce, "", "" }' >> "$out"
    fi
    '''
}
