#!/usr/bin/env python

import pandas as pd
import subprocess
import argparse
import os
import multiprocessing

def align_row(args):
    """
    Worker function to align each row.
    """
    fastq_dir, amplicon_dir, alignments_dir, row = args
    fastq_file = os.path.join(fastq_dir, row['id'] + '.fastq')
    fasta_file = os.path.join(amplicon_dir, row['id'] + '_amplicon.fasta')
    output_file = os.path.join(alignments_dir, row['id'] + '_alignment.sam')

    subprocess.run(["bwa", "index", fasta_file])
    subprocess.run(["bwa", "mem", "-t", "28", "-r", "1", "-k", "11", "-A", "8", "-E", "1", fasta_file, fastq_file, "-o", output_file])

def align(target_info, fastq_dir, amplicon_dir, sample_name):
    target_info_df = pd.read_csv(target_info)

    # Specify the directory path
    alignments_dir = f'{sample_name}_alignments'
    os.makedirs(alignments_dir, exist_ok=True)

    # Prepare jobs for multiprocessing
    jobs = [(fastq_dir, amplicon_dir, alignments_dir, row) for index, row in target_info_df.iterrows()]

    # Create a multiprocessing pool
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())

    # Process each row in parallel
    pool.map(align_row, jobs)

    # Close the pool and wait for the work to finish
    pool.close()
    pool.join()

if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(description="Extract reads overlapping a certain genomic position.")

    # Add the arguments
    parser.add_argument("--target_info", required=True, type=str, help="Metadata file")
    parser.add_argument("--fastq_dir", required=True, type=str, help="fastq dir")
    parser.add_argument("--amplicon_dir", required=True, type=str, help="amplicon dir")
    parser.add_argument("--sample_name", required=True, type=str, help="Sample name")
    # Parse the arguments
    args = parser.parse_args()

    align(args.target_info, args.fastq_dir, args.amplicon_dir, args.sample_name)