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


def compute_integration_percentage(target_info, alignment_dir):
    
    os.mkdir("attL_extracted_reads")
    os.mkdir("attR_extracted_reads")
    target_info_df = pd.read_csv(target_info)
    integration_dict = {}

    for index, row in target_info_df.iterrows():
        alignment_file = alignment_dir + '/' + row['id'] + '_alignment.sam'
        alignment = pysam.AlignmentFile(alignment_file, "r")

        attL_umis = set()
        attR_umis = set()
        beacon_umis = set()
        attL_records = {}
        attR_records = {}

        for read in alignment:
            # Ignore if there's supplementary alignment
            if read.has_tag('XA'):
                continue
            
            if int(read.get_tag('AS')) < 30:
                continue

            # Determine if the selected read maps to the cargo or wildtype amplicon
            read_mapped_to_attL = read.reference_name == 'attL_amplicon'
            read_mapped_to_attR = read.reference_name == 'attR_amplicon'
            read_mapped_to_beacon = read.reference_name == 'beacon_amplicon'

            # Extract UMI from the read name
            umi = read.query_name.split(":")[-1]

            # Calculate start and end using CIGAR string
            start, end = 0, len(read.seq)
            start_soft_clip = 0
            end_soft_clip = 0

            # measure softclipping in order to cleave off ends for preparation into cs2
            for i, (operation, length) in enumerate(read.cigartuples):
                #beginning of sequence is soft clip
                if i == 0 and operation == 4:
                    start_soft_clip = length
                    start = start + start_soft_clip
                # end of sequence is softclipped
                elif i != 0 and operation == 4:
                    end_soft_clip = length
                    end = end - end_soft_clip
            

            # Now start and end should be the indices slicing the aligned part of the sequence
            seq = read.seq[start:end]
            qual = read.query_qualities[start:end]

            # Create SeqRecord
            seq_record = SeqRecord(Seq(seq), id=read.qname, description="",
                                   letter_annotations={"phred_quality": qual})
            
            print(attL_umis)
            print(attR_umis)
            print(beacon_umis)

            if read_mapped_to_attL:
                attL_umis.add(umi)
                # Only keep the record if this UMI hasn't been seen before for attL
                if umi not in attL_records:
                    attL_records[umi] = seq_record
            elif read_mapped_to_attR:
                attR_umis.add(umi)
                # Only keep the record if this UMI hasn't been seen before for attR
                if umi not in attR_records:
                    attR_records[umi] = seq_record
            elif read_mapped_to_beacon:
                beacon_umis.add(umi)

        alignment.close()

        

        # Write records to a FASTQ file
        if len(attL_records.values()) > 0:
            with open("attL_extracted_reads" + '/' + row['id'] + '_attL.fastq', "w") as attL_output_handle:
                SeqIO.write(attL_records.values(), attL_output_handle, "fastq")
        if len(attR_records.values()) > 0:
            with open("attR_extracted_reads" + '/' + row['id'] + '_attR.fastq', "w") as attR_output_handle:
                SeqIO.write(attR_records.values(), attR_output_handle, "fastq")

        total_count = len(attL_umis) + len(attR_umis) + len(beacon_umis)

        if total_count > 0:
            integration_percentage = (len(attL_umis) + len(attR_umis)) / total_count * 100
        else:
            integration_percentage = 0.0
        
        gene_name = row['gene_name']
        gene_strand = row['gene_strand']
        gene_distance = row['gene_distance']
        same_strand = row['same_strand']
        threat_tier = row['threat_tier']
        num_attL_reads = str(len(attL_umis))
        num_attR_reads = str(len(attR_umis))
        num_beacon_reads = str(len(beacon_umis))
        
        integration_dict[row['id']] = (num_attL_reads, num_attR_reads, num_beacon_reads, integration_percentage,gene_name,gene_strand,gene_distance,same_strand,threat_tier)

        print(row['id'] + ' integration_percentage: ' + str(integration_percentage))
    return(integration_dict)


def write_integration_percentage(integration_dict):
    # Open the CSV file in write mode
    filename = 'integration_stats.csv'
    with open(filename, 'w') as file:
        # Write the header row
        file.write('Target,Number of AttL reads,Number of AttR Reads,Number of Beacon reads,Integration Percentage,Closest Gene Name,Gene Strand,Distance from Gene,Same Strand as Cryptic,Threat Tier\n')

        # Write the data rows
        for key, value in integration_dict.items():
            file.write(f'{key},{value[0]},{value[1]},{value[2]},{value[3]},{value[4]},{value[5]},{value[6]},{value[7]},{value[8]}\n')


if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(description="Extract reads overlapping a certain genomic position.")

    # Add the arguments
    parser.add_argument("--target_info", required=True, type=str, help="Metadata file")
    parser.add_argument("--alignment_dir", required=True, type=str, help="Metadata file")
    # Parse the arguments
    args = parser.parse_args()

    integration_dict = compute_integration_percentage(args.target_info, args.alignment_dir)
    write_integration_percentage(integration_dict)