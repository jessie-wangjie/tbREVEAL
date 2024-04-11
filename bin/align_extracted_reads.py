#!/usr/bin/env python

import pandas as pd
import subprocess
import argparse
import os
import multiprocessing

import os
import subprocess

import os
import subprocess

def align_row(args):
    """
    Worker function to align each row.
    """
    fastq_dir, amplicon_dir, row = args
    fastq_file = f"{fastq_dir}/{row['id']}.fastq"
    fasta_file = f"{amplicon_dir}/{row['id']}_amplicon.fasta"
    output_file = f"{row['id']}_alignment.sam"
    output_bam = f"{row['id']}_alignment.bam"  # This will now directly receive the sorted BAM output

    subprocess.run(["bwa", "index", fasta_file])
    subprocess.run(["bwa", "mem", "-t", "28", "-r", "1", "-k", "11", "-A", "8", "-E", "1", fasta_file, fastq_file, "-o", output_file])

    # Pipe samtools view output to samtools sort and output directly to the final BAM file
    subprocess.run(f"samtools view -b {output_file} | samtools sort -o {output_bam}", shell=True)

    subprocess.run(f"samtools index {output_bam}", shell=True)

def align(target_info, fastq_dir, amplicon_dir):
    target_info_df = pd.read_csv(target_info)

    # Prepare jobs for multiprocessing
    jobs = [(fastq_dir, amplicon_dir, row) for index, row in target_info_df.iterrows()]

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
    # Parse the arguments
    args = parser.parse_args()

    align(args.target_info, args.fastq_dir, args.amplicon_dir)