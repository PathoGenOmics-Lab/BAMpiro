#!/usr/bin/env python3
import sys
import json
import argparse
import datetime
import gzip
import os

# List of keys expected by the legacy system
# Added LIN_COUNTS at the end for the specific counts
ORDERED_KEYS = [
    'DAT_OUT', 'GID_IN', 'FQ_IN', 'PHE_OUT', 'FV_IN', 'FV_OUT', 
    'FQC_IN', 'TRD_OUT', 'RDL_OUT', 'LNG_OUT', 'TRM_IN', 'TRM_OUT',
    'ASV_OUT', 'FSV_OUT', 'RSV_OUT', 'ADP_OUT', 'SQP_IN', 'SIN_OUT', 
    'PM_OUT', 'PA_OUT', 'PD_OUT', 'MG_OUT', 'BWA_IN', 'SM1_IN', 
    'MD_IN', 'SM2_IN', 'IDX1_IN', 'RTC_IN', 'IDR_IN', 'BAMF_IN', 
    'IDX2_IN', 'IDX3_IN', 'SFS_IN', 'RDC_OUT', 'DUP_OUT', 'PP_OUT', 
    'MAP_OUT', 'BQC_IN', 'COV_OUT', 'STD_OUT', 'NC_OUT', 'N7_OUT', 
    'VC_IN', 'BPF_OUT', 'VFL_OUT', 'VAR_OUT', 'SNP_OUT', 'HOM_OUT', 
    'HET_OUT', 'HOMO_IND_OUT', 'ANN_IN', 'LIN_OUT', 'LIN_COUNTS'
]

def parse_args():
    parser = argparse.ArgumentParser(description="Generate Legacy TB Pipeline Log")
    parser.add_argument('--sample', required=True, help="Sample ID")
    parser.add_argument('--fastp-json', required=True, help="Path to fastp.json")
    parser.add_argument('--bam-stats', required=True, help="Path to samtools stats output")
    parser.add_argument('--vcf', required=True, help="Path to final VCF")
    parser.add_argument('--ref-fai', required=True, help="Path to reference .fai index")
    parser.add_argument('--pathotypr-report', required=False, help="Path to Pathotypr summary TSV output")
    return parser.parse_args()

def get_genome_size(fai_file):
    total = 0
    try:
        with open(fai_file, 'r') as f:
            for line in f:
                parts = line.split('\t')
                if len(parts) >= 2:
                    total += int(parts[1])
    except:
        pass
    return total if total > 0 else 4411532 

def get_bam_stats(stats_file):
    mapped = 0
    total = 0
    try:
        with open(stats_file, 'r') as f:
            for line in f:
                if line.startswith('SN\treads mapped:'):
                    mapped = int(line.split('\t')[2])
                elif line.startswith('SN\traw total sequences:'):
                    total = int(line.split('\t')[2])
    except:
        pass
    return total, mapped

def analyze_vcf(vcf_file):
    stats = {'snps': 0, 'indels': 0, 'homo_total': 0, 'het_total': 0, 'homo_indels': 0}
    opener = gzip.open if vcf_file.endswith('.gz') else open
    try:
        with opener(vcf_file, 'rt') as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                if len(parts) < 10: continue
                ref = parts[3]
                alt = parts[4]
                if alt in ['.', '*']: continue
                fmt = parts[8].split(':')
                if 'GT' not in fmt: continue
                gt_idx = fmt.index('GT')
                sample = parts[9].split(':')
                if len(sample) <= gt_idx: continue
                gt = sample[gt_idx]
                if gt in ['./.', '0/0', '0|0']: continue
                
                is_snp = (len(ref) == 1 and len(alt) == 1)
                is_homo = False
                alleles = gt.replace('|', '/').split('/')
                if len(alleles) == 2 and alleles[0] == alleles[1]:
                    is_homo = True
                
                if is_snp: stats['snps'] += 1
                else:
                    stats['indels'] += 1
                    if is_homo: stats['homo_indels'] += 1
                if is_homo: stats['homo_total'] += 1
                else: stats['het_total'] += 1
    except: pass
    return stats

def get_pathotypr_data(report_file):
    """
    Parses the Pathotypr summary TSV.
    Header usually: genome <tab> lineage:count <tab> major_lineage
    Returns a dict with 'major' and 'counts'.
    """
    results = {'major': 'NA', 'counts': 'NA'}
    
    if not report_file or not os.path.exists(report_file) or "NO_LINEAGE_FILE" in report_file:
        return results
        
    try:
        with open(report_file, 'r') as f:
            # Read header
            header_line = f.readline().strip()
            if not header_line: return results
            headers = header_line.split('\t')
            
            # Identify columns
            try:
                idx_major = headers.index('major_lineage')
            except ValueError:
                idx_major = -1
            
            try:
                idx_count = headers.index('lineage:count')
            except ValueError:
                idx_count = -1

            # Read data line
            data_line = f.readline().strip()
            if data_line:
                cols = data_line.split('\t')
                
                # Extract Major Lineage
                if idx_major != -1 and len(cols) > idx_major:
                    results['major'] = cols[idx_major]
                
                # Extract Lineage Counts
                if idx_count != -1 and len(cols) > idx_count:
                    results['counts'] = cols[idx_count]
                    
    except Exception as e:
        sys.stderr.write(f"Warning: Error parsing Pathotypr file: {e}\n")
        
    return results

def main():
    args = parse_args()
    data = {k: "NA" for k in ORDERED_KEYS} 
    
    # 1. Metadata & Fastp
    data['DAT_OUT'] = str(datetime.date.today())
    data['GID_IN'] = args.sample
    data['FQ_IN'] = f"{args.sample}_R1/R2" 
    
    g_size = get_genome_size(args.ref_fai)
    r_len = 0.0

    try:
        with open(args.fastp_json, 'r') as f:
            j = json.load(f)
            before = j.get('summary', {}).get('before_filtering', {})
            after = j.get('summary', {}).get('after_filtering', {})
            dup = j.get('duplication', {})
            
            data['TRD_OUT'] = str(before.get('total_reads', 'NA'))
            data['TRM_OUT'] = str(after.get('total_reads', 'NA'))
            r_len = float(after.get('read1_mean_length', 0))
            data['RDL_OUT'] = f"{r_len:.1f}"
            data['DUP_OUT'] = f"{float(dup.get('rate', 0))*100:.2f}"
            data['PHE_OUT'] = '33'
            data['FV_OUT'] = 'FASTQ_SUCCESS'
            
            if g_size > 0:
                data['LNG_OUT'] = f"{float(after.get('total_bases', 0))/g_size:.2f}"
    except: pass

    # 2. Mapping
    tot, mapped = get_bam_stats(args.bam_stats)
    data['MAP_OUT'] = str(mapped)
    if g_size > 0 and r_len > 0:
        data['COV_OUT'] = f"{(mapped * r_len) / g_size:.2f}"

    # 3. Variants
    v = analyze_vcf(args.vcf)
    data['SNP_OUT'] = str(v['snps'])
    data['HOM_OUT'] = str(v['homo_total'])
    data['HET_OUT'] = str(v['het_total'])
    data['VAR_OUT'] = str(v['snps'] + v['indels'])
    data['HOMO_IND_OUT'] = str(v['homo_indels'])
    
    # 4. Pathotypr Data
    if args.pathotypr_report:
        p_data = get_pathotypr_data(args.pathotypr_report)
        data['LIN_OUT'] = p_data['major']
        data['LIN_COUNTS'] = p_data['counts']

    # Write output
    with open(f"{args.sample}.log", 'w') as out:
        for key in ORDERED_KEYS:
            out.write(f"{key}\t{data[key]}\n")

if __name__ == "__main__":
    main()
