#!/usr/bin/env python

import argparse
import pandas as pd
import pysam
import subprocess
from Bio import SeqIO

def extract_reads(target_info, bam_file, window, sample_name):

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
        start = max(1, start_position - window)
        end = end_position + window

        print('Extracting reads for ' + id + '...')

        command = f'samtools view {bam_file} {chromosome}:{start}-{end} -b -P | samtools fastq -N - -F256 > {sample_name}_{id}.fastq'
        subprocess.run(command, shell=True)

        # Count the number of reads in the fastq file using Biopython's SeqIO
        count = sum(1 for _ in SeqIO.parse(f'{sample_name}_{id}.fastq', 'fastq'))
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