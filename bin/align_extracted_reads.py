#!/usr/bin/env python

import pandas as pd
import subprocess
import argparse
import os

def align(target_info,fastq_dir,amplicon_dir):
    target_info_df = pd.read_csv(target_info)
    # Specify the directory path
    alignments_dir = 'alignments'
    os.makedirs(alignments_dir)
    for index, row in target_info_df.iterrows():
        fastq_file = fastq_dir + '/' + row['id'] + '.fastq'
        fasta_file = amplicon_dir + '/' + row['id'] + '_amplicon.fasta'
        output_file = alignments_dir + '/' + row['id'] + '_alignment.sam'
        #subprocess.run(["minimap2", "-ax", "sr", fasta_file, fastq1_file, fastq2_file, "-o", output_file],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)

        subprocess.run(["bwa", "index", fasta_file])
        subprocess.run(["bwa", "mem", "-t", "28", "-r", "0.1", fasta_file, fastq_file, "-o", output_file])

if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(description="Extract reads overlapping a certain genomic position.")

    # Add the arguments
    parser.add_argument("--target_info", required=True, type=str, help="Metadata file")
    parser.add_argument("--fastq_dir", required=True, type=str, help="fastq dir")
    parser.add_argument("--amplicon_dir", required=True, type=str, help="amplicon dir")
    # Parse the arguments
    args = parser.parse_args()

    align(args.target_info, args.fastq_dir, args.amplicon_dir)