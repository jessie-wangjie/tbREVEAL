#!/usr/bin/env python

import argparse
import os
import json

def parse_data(json_file,sample_name):

    output_fn = f'{sample_name}_qc_summary.csv'
    # Load JSON data
    with open(json_file, 'r') as file:
        data = json.load(file)

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
        ["bases filtered", bases_filtered]
    ]

    # Write data to CSV file
    with open(output_fn, 'w') as file:
        for row in data_for_csv:
            file.write(f'{row[0]},{row[1]}\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse JSON file and count lines in FASTQ files.')
    parser.add_argument('--json_file', help='The JSON file to be parsed.')
    parser.add_argument('--fastq_dir', help='The directory containing the FASTQ files.')
    parser.add_argument('--sample_name', help='The directory containing the FASTQ files.')
    args = parser.parse_args()
    
    parse_data(args.json_file, args.sample_name)