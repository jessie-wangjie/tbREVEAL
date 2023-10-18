#!/usr/bin/env python

import argparse
import pysam
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import pandas as pd 
from utils.base import *
import subprocess
import sys
import psycopg2
from pysam import VariantFile
import re

def get_cryptic_dinucleotide_info():
    dinucleotide_info_dict = {}
    query = '''
    SELECT file_registry_id$, chr, central_nucleotides_start FROM integrase_cryptic_attachment_site
    WHERE integrase = 'Bxb1'
    '''
    cur.execute(query)
    result = cur.fetchall()
    window_size = 50
    for row in result:
        cas_id, chr, dinucleotide_pos = row
        dinucleotide_info_dict[cas_id] = (chr, dinucleotide_pos, dinucleotide_pos-window_size, dinucleotide_pos+window_size)
    return(dinucleotide_info_dict)

# function to check cryptic database to see if any of the translocations occur near any cryptic dinucleotide pos
def find_matching_cryptic_dinucleotide(chromosome, position, data):
    for key, value in data.items():

        if value[0] == ('chr'+str(chromosome)) and int(value[2]) <= int(position) <= int(value[3]):
            return(key, value[0], int(position), int(value[2]), int(value[1]), int(value[3]))
    return(None,None,None,None,None,None)

def parse_vcf_file(vcf_file,dinucleotide_info):
    vcf_in = VariantFile(vcf_file)
    translocation_info_list = set()
    for rec in vcf_in.fetch():
        ref_chrom = rec.chrom
        ref_pos = rec.pos
        alt_allele = rec.alleles[1]
        
        if '[' in alt_allele or ']' in alt_allele:
            match = re.search(r'[\[\]]chr(\d+):(\d+)', alt_allele)
            if match:
                num_reads_supporting = rec.info["PAIR_COUNT"]
                alt_chrom = match.group(1)
                alt_pos = match.group(2)
                cas_id, alt_chrom, alt_breakpoint_pos, window_lower, dinucleotide_pos, window_upper = find_matching_cryptic_dinucleotide(alt_chrom, alt_pos, dinucleotide_info)
                if cas_id is not None:
                    distance_from_dinucleotide = alt_breakpoint_pos - dinucleotide_pos
                    if ref_chrom == alt_chrom:
                        translocation_type = 'intra-chromosomal'
                    else:
                        translocation_type = 'inter-chromosomal'
                    translocation_info_list.add((cas_id, ref_chrom, ref_pos, alt_chrom, alt_breakpoint_pos, window_lower, dinucleotide_pos, window_upper, distance_from_dinucleotide, num_reads_supporting, translocation_type))
    return(translocation_info_list)

def create_final_df(translocation_list):
    column_names = [
    "cas_id", 
    "ref_chrom", 
    "ref_pos", 
    "alt_chrom", 
    "alt_pos", 
    "window_lower", 
    "dinucleotide_pos", 
    "window_upper", 
    "distance_from_dinucleotide", 
    "read_pairs_supporting",
    "translocation_type"
    ]
    df = pd.DataFrame(translocation_list, columns=column_names)
    return(df)


# def retrieve_evidence_reads(bam_evidence_file, translocation_list):
#     bamfile = pysam.AlignmentFile(bam_evidence_file, "rb")
#     for read in bamfile.fetch('chr1', 100, 120):
#         print(read)

#     bamfile.close()
                

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process VCF file from Manta and output CSV with translocation information.')
    parser.add_argument('--vcf_file', type=str, help='Path to the VCF file')
    # parser.add_argument('--bam_evidence_file', type=str, help='Path to the BAM evidence file')
    parser.add_argument('--output_csv', type=str, help='Output CSV file name')
    
    args = parser.parse_args()
    
    cryptic_site_dinucleotide_info = get_cryptic_dinucleotide_info()
    translocation_list = parse_vcf_file(args.vcf_file, cryptic_site_dinucleotide_info)
    translocation_df = create_final_df(translocation_list)
    translocation_df.to_csv(args.output_csv, index=False)


        