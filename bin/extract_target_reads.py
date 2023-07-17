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

def extract_reads(target_info, bam_file, window):

    # Specify the directory path
    fastq_dir = 'extracted_reads'
    os.makedirs(fastq_dir, exist_ok=True)

    # Open the BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")


    target_info_df = pd.read_csv(target_info)

    for index, row in target_info_df.iterrows():
        position = row['start']
        chromosome = row['chromosome']
        id = row['id']
        start = position - window
        end = position + window

        print('Extracting reads for ' + id + '...')
        reads_to_extract = set()
        
        command = f'samtools view {bam_file} {chromosome}:{start}-{end} -b | samtools fastq -N - -o {fastq_dir}/{id}.fastq'
        subprocess.run(command, shell=True)

    bam.close()


if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(description="Extract reads overlapping a certain genomic position.")

    # Add the arguments
    parser.add_argument("--target_info", required=True, type=str, help="Metadata file")
    parser.add_argument("--window_size", type=int, default=1000, help="The size of the window around the position of interest")
    parser.add_argument("--bam_file", required=True, type=str, help="The name of the BAM file")
    
    # Parse the arguments
    args = parser.parse_args()

    extract_reads(args.target_info, args.bam_file, args.window_size)
