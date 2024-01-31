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
import re


import os
import pandas as pd
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def compute_integration_percentage(target_info, alignment_dir, sample_name):
    os.makedirs(f"{sample_name}_attL_extracted_reads", exist_ok=True)
    os.makedirs(f"{sample_name}_attR_extracted_reads", exist_ok=True)
    os.makedirs(f"{sample_name}_beacon_extracted_reads", exist_ok=True)
    
    target_info_df = pd.read_csv(target_info)
    integration_dict = {}

    p_reg_sequence = 'GTGGTTTGTCTGGTCAACCACCGCG'
    p_prime_sequence = 'CTCAGTGGTGTACGGTACAAACCCA'

    def get_umi(read):
        return read.query_name
    
    def calculate_soft_clipped_indices(read):
        start, end = 0, len(read.seq)
        for i, (operation, length) in enumerate(read.cigartuples):
            if i == 0 and operation == 4:
                start += length
            elif operation == 4:
                end -= length
        return start, end
    
    for index, row in target_info_df.iterrows():
        alignment_file = f"{alignment_dir}/{row['id']}_alignment.sam"
        
        attL_records = {}
        attR_records = {}
        beacon_records = {}

        full_attb_sequence = row['beacon'][20:len(row['beacon'])-21]

        categories = {
            'ambiguous_attL':set(), 'ambiguous_attR':set(),
            'partial_attL': set(), 'complete_attL': set(), 'cargo_attL': set(),
            'partial_attR': set(), 'complete_attR': set(), 'cargo_attR': set(),
            'ambiguous_beacon': set(), 'partial_beacon': set(), 'complete_beacon': set(), 'wt': set()
        }

        with pysam.AlignmentFile(alignment_file, "r") as alignment:
            for read in alignment:
                if read.flag in [2048, 2064, 4]:  # Supplementary, Secondary, or Unaligned
                    continue

                umi = get_umi(read)
                start, end = calculate_soft_clipped_indices(read)
                seq = read.seq[start:end]
                qual = read.query_qualities[start:end]
                
                seq_record = SeqRecord(Seq(seq), id=read.qname, description="", letter_annotations={"phred_quality": qual})
                alignment_start, alignment_end = read.reference_start + 1, read.reference_end + 1

                
                rules = [
                    (read.reference_name == 'attL_amplicon' and 'CAS' in row['id'] and alignment_start < 32 and alignment_end > 54, 
                     ['ambiguous_attL']),
                    (read.reference_name == 'attL_amplicon' and 'CAS' in row['id'] and alignment_start < 32 and alignment_end > 69 and p_prime_sequence in seq, 
                     ['complete_attL']),
                    (read.reference_name == 'attL_amplicon' and 'CAS' in row['id'] and alignment_start < 32 and alignment_end > 79 and p_prime_sequence in seq, 
                     ['complete_attL', 'cargo_attL']),

                    (read.reference_name == 'attR_amplicon' and 'CAS' in row['id'] and alignment_start < 30 and alignment_end > 57, 
                     ['ambiguous_attR']),
                    (read.reference_name == 'attR_amplicon' and 'CAS' in row['id'] and alignment_start < 20 and alignment_end > 57 and p_reg_sequence in seq, 
                     ['complete_attR']),
                    (read.reference_name == 'attR_amplicon' and 'CAS' in row['id'] and alignment_start < 10 and alignment_end > 57 and p_reg_sequence in seq, 
                     ['complete_attR', 'cargo_attR']),

                    (read.reference_name == 'beacon_amplicon' and 'CAS' in row['id'] and ((alignment_start <= 21 and alignment_end < 66) or (alignment_start > 21 and alignment_end >= 66)), 
                     ['partial_beacon']),
                    (read.reference_name == 'beacon_amplicon' and 'CAS' in row['id'] and alignment_start <= 21 and alignment_end >= 66, 
                     ['partial_beacon', 'complete_beacon']),
                    
                    # slightly different rule compared to CAS sites because cryptic B based on 46 bp attB and the beacon written is 38 bp
                    # note the difference in alignment ends (15 bp versus 11 bp)
                    (read.reference_name == 'attL_amplicon' and 'AA' in row['id'] and alignment_start < 12 and alignment_end > 34, 
                     ['ambiguous_attL']),
                    (read.reference_name == 'attL_amplicon' and 'AA' in row['id'] and alignment_start < 12 and alignment_end > 45 and p_prime_sequence in seq, 
                     ['complete_attL']),
                    (read.reference_name == 'attL_amplicon' and 'AA' in row['id'] and alignment_start < 12 and alignment_end > 55 and p_prime_sequence in seq, 
                     ['complete_attL', 'cargo_attL']),

                    (read.reference_name == 'attR_amplicon' and 'AA' in row['id'] and alignment_start < 30 and alignment_end > 57, 
                     ['ambiguous_attR']),
                    (read.reference_name == 'attR_amplicon' and 'AA' in row['id'] and alignment_start < 20 and alignment_end > 57 and p_reg_sequence in seq, 
                     ['complete_attR']),
                    (read.reference_name == 'attR_amplicon' and 'AA' in row['id'] and alignment_start < 10 and alignment_end > 57 and p_reg_sequence in seq, 
                     ['complete_attR', 'cargo_attR']),

                    (read.reference_name == 'beacon_amplicon' and 'AA' in row['id'] and alignment_start > 21 or alignment_end < (21 + len(full_attb_sequence)), 
                     ['ambiguous_beacon']),
                    (read.reference_name == 'beacon_amplicon' and 'AA' in row['id'] and alignment_start <= 21 and alignment_end >= (21 + len(full_attb_sequence)) and bool(re.search(full_attb_sequence, seq)) == False, 
                     ['partial_beacon']),
                    (read.reference_name == 'beacon_amplicon' and 'AA' in row['id'] and alignment_start <= 21 and alignment_end >= (21 + len(full_attb_sequence)) and bool(re.search(full_attb_sequence, seq)) == True, 
                     ['complete_beacon']),

                    (read.reference_name == 'wt_amplicon', ['wt'])
                ]

                # if row['id'] == 'AA876':
                #     if read.reference_name == 'beacon_amplicon':
                #         print(read.seq)
                #         print(alignment_start)
                #         print(alignment_end)
                #         print(rules)
                #         print(bool(re.search(full_attb_sequence, seq)))
                #         print('-----')
                        
                        

                for condition, cat_keys in rules:
                    if condition:
                        for key in cat_keys:
                            categories[key].add(umi)
                            if 'attL' in key:
                                attL_records[umi] = seq_record
                            elif 'attR' in key:
                                attR_records[umi] = seq_record 
                            # only care about written beacon alignments, not cryptic beacons (aka regular genomic)
                            elif 'beacon' in key and 'AA' in row['id']:
                                beacon_records[umi] = seq_record                
                    
        for prefix, records in [(f'{sample_name}_attL', attL_records), (f'{sample_name}_attR', attR_records),(f'{sample_name}_beacon', beacon_records)]:
            if records:
                integration_side = prefix.split('_')[-1]
                path = f"{prefix}_extracted_reads/{row['id']}_{integration_side}.fastq"
                SeqIO.write(records.values(), path, "fastq")

        counts = {k: len(v) for k, v in categories.items()}
        beacon_total = sum([counts[key] for key in ['partial_beacon', 'complete_beacon']])
        beacon_placement_denominator = sum([counts[key] for key in ['complete_attL', 'complete_attR', 'partial_beacon', 'complete_beacon','wt']])
        partial_total = sum([counts[key] for key in ['partial_attL', 'partial_attR', 'partial_beacon', 'wt']])
        complete_P_total = sum([counts[key] for key in ['complete_attL', 'complete_attR', 'complete_beacon', 'wt']])
        cargo_P_total = sum([counts[key] for key in ['cargo_attL', 'cargo_attR', 'complete_beacon', 'wt']])
        
        def calc_percentage(numerator, denominator):
            return (numerator / denominator * 100) if denominator > 0 else 0.0

        integration_dict[row['id']] = (
            counts['wt'], 
            counts['partial_attL'], 
            counts['complete_attL'], 
            counts['cargo_attL'],
            counts['partial_attR'], 
            counts['complete_attR'], 
            counts['cargo_attR'], 
            counts['partial_beacon'],
            counts['partial_beacon'] if 'CAS' in row['id'] else counts['complete_beacon'],
            calc_percentage(counts['partial_attL'] + counts['partial_attR'], partial_total),
            calc_percentage(counts['complete_attL'] + counts['complete_attR'], complete_P_total),
            calc_percentage(counts['cargo_attL'] + counts['cargo_attR'], cargo_P_total),
            100 if 'CAS' in row['id'] else calc_percentage(counts['partial_beacon'] + counts['complete_beacon'] + counts['complete_attL'] + counts['complete_attR'], beacon_placement_denominator),
            100 if 'CAS' in row['id'] else calc_percentage(counts['complete_beacon'] + counts['complete_attL'] + counts['complete_attR'], beacon_placement_denominator),
            100 if 'CAS' in row['id'] else calc_percentage(counts['complete_beacon'], beacon_total),
        ) + tuple(row[key] for key in ['gene_name', 'gene_strand', 'gene_distance', 'same_strand', 'overlapping_feature', 'threat_tier'])

    return integration_dict


def write_integration_percentage(integration_dict, sample_name):
    # Open the CSV file in write mode
    filename = f'{sample_name}_integration_stats.csv'
    with open(filename, 'w') as file:
        # Write the header row
        file.write('Target,Number of WT Reads,Number of AttL Complete Reads,Number of AttL Cargo Reads,Number of AttR Partial Reads,Number of AttR Complete Reads,Number of AttR Cargo Reads,Number of Partial Beacon Reads,Number of Complete Beacon Reads,Partial P Integration Percentage,Complete P Integration Percentage,Cargo and P Integration Percentage,Partial Beacon Placement,Complete Beacon Placement,Beacon Fidelity,Closest Gene Name,Gene Strand,Distance from Gene,Same Strand as Cryptic,Overlapping Feature,Threat Tier\n')

        # Write the data rows
        for key, value in integration_dict.items():
            file.write(f'{key},{value[1]},{value[2]},{value[3]},{value[4]},{value[5]},{value[6]},{value[7]},{value[8]},{value[9]},{value[10]},{value[11]},{value[12]},{value[13]},{value[14]},{value[15]},{value[16]},{value[17]},{value[18]},{value[19]},{value[20]}\n')


if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(description="Extract reads overlapping a certain genomic position.")

    # Add the arguments
    parser.add_argument("--target_info", required=True, type=str, help="Metadata file")
    parser.add_argument("--alignment_dir", required=True, type=str, help="Alignment directory")
    parser.add_argument("--sample_name", required=True, type=str, help="Sample name")
    # Parse the arguments
    args = parser.parse_args()

    integration_dict = compute_integration_percentage(args.target_info, args.alignment_dir,args.sample_name)
    write_integration_percentage(integration_dict,args.sample_name)