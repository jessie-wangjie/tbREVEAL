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
        if value[0] == (str(chromosome)) and int(value[2]) <= int(position) <= int(value[3]):
            return(key, value[0], int(position), int(value[2]), int(value[1]), int(value[3]))
    return(None,None,None,None,None,None)

def parse_vcf(vcf_file_path,dinucleotide_info):

    translocation_info_list = set()

    with open(vcf_file_path, 'r') as vcf_file:
        for line in vcf_file:
            # Skip header lines
            if line.startswith('#'):
                continue
        
            fields = line.strip().split('\t')
            ref_chrom = fields[0]
            ref_pos = fields[1]
            info = fields[7]
            
            

            if "<TRA>" in fields[4]: # Check for translocation in ALT field
                # Extract CHR2 and CHR2_POS
                for field in info.split(';'):
                    if field.startswith("CHR2="):
                        alt_chrom = field.split('=')[1]
                    elif field.startswith("CHR2_POS="):
                        alt_pos = field.split('=')[1]
                
                print(f"Translocation between {ref_chrom}:{ref_pos} and {alt_chrom}:{alt_pos}")

                cas_id, alt_chrom, alt_breakpoint_pos, window_lower, dinucleotide_pos, window_upper = find_matching_cryptic_dinucleotide(alt_chrom, alt_pos, dinucleotide_info)
                if cas_id is not None:
                    print(cas_id)
                    distance_from_dinucleotide = alt_breakpoint_pos - dinucleotide_pos
                    if ref_chrom == alt_chrom:
                        translocation_type = 'intra-chromosomal'
                    else:
                        translocation_type = 'inter-chromosomal'

                    translocation_info_list.add((cas_id, ref_chrom, ref_pos, alt_chrom, alt_breakpoint_pos, window_lower, dinucleotide_pos, window_upper, distance_from_dinucleotide, translocation_type))
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
    "translocation_type"
    ]
    df = pd.DataFrame(translocation_list, columns=column_names)
    return(df)
                

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse VCF file to get translocations.")
    parser.add_argument('--vcf', type=str, required=True, help="Path to the VCF file")
    parser.add_argument('--output_csv', type=str, required=True, help="Path to the VCF file")

    cryptic_dinucleotide_info = get_cryptic_dinucleotide_info()

    args = parser.parse_args()
    translocation_list = parse_vcf(args.vcf, cryptic_dinucleotide_info)
    df = create_final_df(translocation_list)
    df.to_csv(args.output_csv)