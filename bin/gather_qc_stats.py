#!/usr/bin/env python

import argparse
import os
import json

def count_lines(directory):
    # Initialize total line count to 0
    total_line_count = 0
    # Create a list to store the output data
    data = []

    output_fn = 'probe_read_counts.csv'
    # Iterate over all files in the directory
    for filename in os.listdir(directory):
        # Only process .fastq files
        if filename.endswith(".fastq"):
            # Remove the file extension to get the sample name
            sample_name = filename[:-6]

            # Count the number of lines in the file
            with open(os.path.join(directory, filename), 'r') as file:
                # divide by 4 bc its fastq
                line_count = sum(1 for _ in file) / 4

            # Add this file's line count to the total line count
            total_line_count += line_count

            # Append the sample name and line count to the data list
            data.append((sample_name, line_count))

    # Write the data to a CSV file
    with open(output_fn, 'w') as file:
        file.write("site,reads\n")
        for sample_name, line_count in data:
            file.write(f"{sample_name},{line_count}\n")

    # divide by 4 bc its fastq
    return total_line_count/4

def parse_data(json_file, reads_near_probe):

    output_fn = 'qc_summary.csv'
    # Load JSON data
    with open(json_file, 'r') as file:
        data = json.load(file)

    # Get before_filtering and after_filtering as separate dictionaries
    before_filtering = data["summary"]["before_filtering"]
    after_filtering = data["summary"]["after_filtering"]

    # Calculate the reads and bases filtered
    reads_filtered = before_filtering["total_reads"] - after_filtering["total_reads"]
    bases_filtered = before_filtering["total_bases"] - after_filtering["total_bases"]

    # Calculate reads near probe %
    reads_near_probe_percent = reads_near_probe / after_filtering["total_reads"]

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
        ["reads near probe", reads_near_probe],
        ["reads near probe %", f'{reads_near_probe_percent:.2%}']
    ]

    # Write data to CSV file
    with open(output_fn, 'w') as file:
        for row in data_for_csv:
            file.write(f'{row[0]},{row[1]}\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse JSON file and count lines in FASTQ files.')
    parser.add_argument('--json_file', help='The JSON file to be parsed.')
    parser.add_argument('--fastq_dir', help='The directory containing the FASTQ files.')
    args = parser.parse_args()
    
    reads_near_probe = count_lines(args.fastq_dir)
    parse_data(args.json_file, reads_near_probe)