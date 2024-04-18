#!/usr/bin/env python

import argparse
import os
import json
import pysam
from subprocess import Popen, PIPE
from functools import reduce

def parse_data(json_file,original_bam_fn,deduped_bam_fn,sample_name):

    output_fn = f'{sample_name}_qc_summary.csv'
    # Load JSON data
    with open(json_file, 'r') as file:
        data = json.load(file)

    def bam_read_count(bamfile):
        """Return a tuple of the number of mapped and unmapped reads in a bam file"""
        # Corrected command using subprocess.Popen
        cmd = ['samtools', 'flagstat', bamfile, '-O', 'tsv']
        p = Popen(cmd, stdout=PIPE)

        # Use another Popen for the cut and head commands or process the output in Python
        # It's more efficient to do the processing directly in Python
        mapped = 0
        unmapped = 0

        for i, line in enumerate(p.stdout):
            if i == 0:  # Only the first line contains the total number of reads
                total_reads = int(line.decode().split()[0])
                print(total_reads)
                # In this case, samtools flagstat already provides mapped and unmapped counts
                # so parsing as initially intended might not be directly applicable
                # You may want to adjust this logic based on the actual output format
                continue  # Assuming further processing or skipping as needed
            # Example processing, adjust according to the actual flagstat output format
            elif i == 6:
                mapped = int(line.decode().split()[0])
                print(mapped)

        return (mapped,total_reads)

    original_mapped_read_count,original_read_count = bam_read_count(original_bam_fn)
    deduped_mapped_read_count, deduped_read_count = bam_read_count(deduped_bam_fn)


    # Get before_filtering and after_filtering as separate dictionaries
    before_filtering = data["summary"]["before_filtering"]
    after_filtering = data["summary"]["after_filtering"]

    # Calculate the reads and bases filtered
    reads_filtered = before_filtering["total_reads"] - after_filtering["total_reads"]
    bases_filtered = before_filtering["total_bases"] - after_filtering["total_bases"]

    # Prepare data for CSV
    data_for_csv = [
        ["before filter reads", before_filtering['total_reads']],
        ["before filter bases", before_filtering['total_bases']],
        ["before filter q20 bases", before_filtering['q20_bases']],
        ["before filter q30 bases", before_filtering['q30_bases']],
        ["after filter reads", after_filtering['total_reads']],
        ["after filter bases", after_filtering['total_bases']],
        ["after filter q20 bases", after_filtering['q20_bases']],
        ["after filter q30 bases", after_filtering['q30_bases']],
        ["reads filtered", reads_filtered],
        ["bases filtered", bases_filtered],
        ["total_reads",original_read_count],
        ["total_reads_aligned",original_mapped_read_count],
        ["total_reads_align%",100*original_mapped_read_count/original_read_count],
        ["deduped_reads",deduped_read_count],
    ]
    # Write data to CSV file
    with open(output_fn, 'w') as file:
        for row in data_for_csv:
            file.write(f'{row[0]},{row[1]}\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse JSON file and count lines in FASTQ files.')
    parser.add_argument('--json_file', help='The JSON file to be parsed.')
    parser.add_argument('--original_bam_file', help='The original BAM file to be parsed.')
    parser.add_argument('--deduped_bam_file', help='The deduped BAM file to be parsed.')
    parser.add_argument('--sample_name', help='Sample name')
    args = parser.parse_args()

    parse_data(args.json_file, args.original_bam_file, args.deduped_bam_file, args.sample_name)