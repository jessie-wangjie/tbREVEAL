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

def extract_reads(target_info, bam_file, window, sample_name):

    # Specify the directory path
    fastq_dir = f'{sample_name}_extracted_reads'
    os.makedirs(fastq_dir, exist_ok=True)

    # Open the BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")

    target_info_df = pd.read_csv(target_info)

    # Prepare a dictionary to hold ID and read counts
    read_counts = {}

    for index, row in target_info_df.iterrows():
        start_position = row['start']
        end_position = row['end']
        chromosome = row['chromosome']
        id = row['id']
        start = start_position - window
        end = end_position + window

        print('Extracting reads for ' + id + '...')
        
        command = f'samtools view {bam_file} {chromosome}:{start}-{end} -b | samtools fastq -N - > {fastq_dir}/{id}.fastq'
        subprocess.run(command, shell=True)

        # Count the number of reads in the fastq file using Biopython's SeqIO
        count = sum(1 for _ in SeqIO.parse(f'{fastq_dir}/{id}.fastq', 'fastq'))
        read_counts[id] = count

    bam.close()

    # Save the read counts to a CSV file
    read_counts_df = pd.DataFrame(list(read_counts.items()), columns=['id', 'read_count'])
    read_counts_df.to_csv(f'{sample_name}_read_counts_per_site.csv', index=False)


if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(description="Extract reads overlapping a certain genomic position.")

    # Add the arguments
    parser.add_argument("--target_info", required=True, type=str, help="Metadata file")
    parser.add_argument("--window_size", type=int, default=1000, help="The size of the window around the position of interest")
    parser.add_argument("--bam_file", required=True, type=str, help="The name of the BAM file")
    parser.add_argument("--sample_name", required=True, type=str, help="Sample name")
    
    # Parse the arguments
    args = parser.parse_args()

    extract_reads(args.target_info, args.bam_file, args.window_size, args.sample_name)
