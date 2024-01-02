#!/usr/bin/env python

# import argparse
# import pysam
# import subprocess
# from Bio import SeqIO
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from collections import defaultdict
# import pandas as pd 
# from utils.base import *
# import subprocess
# import sys
# import psycopg2
# import multiprocessing

# def get_hardclipped_reads_in_region(bam_file, chromosome, start, end):
#     # Dictionary to store supplementary read ID and primary alignment sequence
#     results = {}
#     # Set to keep track of processed supplementary reads
#     processed_reads = set()

#     # Open the BAM file
#     with pysam.AlignmentFile(bam_file, "rb") as bam:
#         # Find all supplementary hard clipped reads in the region
#         for read in bam.fetch(chromosome, start, end):
#             if read.cigarstring and 'H' in read.cigarstring and read.is_supplementary and read.query_name not in processed_reads:
#                 processed_reads.add(read.query_name)
#                 # Check for the SA tag
#                 if read.has_tag('SA'):
#                     sa_tag = read.get_tag('SA')
#                     # Parse the SA tag to get the primary alignment location
#                     sa_info = sa_tag.split(';')[0].split(',')
#                     primary_chromosome, primary_position = sa_info[0], int(sa_info[1])

#                     # Fetch the primary alignment using the position from the SA tag
#                     if primary_position <= 100:
#                         primary_position = 100
#                     for primary_read in bam.fetch(primary_chromosome, primary_position-100, primary_position + 100):
#                         # Check if the positions match and the read is not supplementary
#                         if primary_read.reference_start == primary_position - 1 and not primary_read.is_supplementary:
#                             quality = ''.join(map(lambda x: chr( x+33 ), primary_read.query_qualities))
#                             results[read.query_name] = (primary_read.query_sequence,quality)
#                             break

#     return results

# def process_region(args):
#     """
#     Worker function to process each region.
#     """
#     bam_file, chromosome, start, end, id, fastq_dir = args

#     hardclipped_reads = get_hardclipped_reads_in_region(bam_file, chromosome, start, end)

#     command = f'samtools view {bam_file} {chromosome}:{start}-{end} -b | samtools fastq -N -F0 - > {fastq_dir}/{id}.fastq'
#     subprocess.run(command, shell=True)

#     # Append the hard-clipped supplementary read sequences to the FASTQ file
#     with open(f'{fastq_dir}/{id}.fastq', 'a') as fastq_file:
#         for key in hardclipped_reads:
#             sequence, quality_score = hardclipped_reads[key]
#             fastq_file.write(f'@{key}\n{sequence}\n+\n{quality_score}\n')

#     # Count the number of reads
#     count = sum(1 for _ in SeqIO.parse(f'{fastq_dir}/{id}.fastq', 'fastq'))
#     return id, count

# def extract_reads(target_info, bam_file, window, sample_name):
#     fastq_dir = f'{sample_name}_extracted_reads'
#     os.makedirs(fastq_dir, exist_ok=True)

#     target_info_df = pd.read_csv(target_info)

#     # Prepare the multiprocessing pool
#     pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())

#     # Prepare jobs for multiprocessing
#     jobs = []
#     for index, row in target_info_df.iterrows():
#         start = row['start'] - window
#         end = row['end'] + window
#         jobs.append((bam_file, row['chromosome'], start, end, row['id'], fastq_dir))

#     # Execute jobs in parallel
#     results = pool.map(process_region, jobs)

#     # Close the pool and wait for the work to finish
#     pool.close()
#     pool.join()

#     # Collect results
#     read_counts = dict(results)

#     # Save the read counts to a CSV file
#     read_counts_df = pd.DataFrame(list(read_counts.items()), columns=['id', 'read_count'])
#     read_counts_df.to_csv(f'{sample_name}_read_counts_per_site.csv', index=False)


# if __name__ == "__main__":
#     # Create the parser
#     parser = argparse.ArgumentParser(description="Extract reads overlapping a certain genomic position.")

#     # Add the arguments
#     parser.add_argument("--target_info", required=True, type=str, help="Metadata file")
#     parser.add_argument("--window_size", type=int, default=1000, help="The size of the window around the position of interest")
#     parser.add_argument("--bam_file", required=True, type=str, help="The name of the BAM file")
#     parser.add_argument("--sample_name", required=True, type=str, help="Sample name")
    
#     # Parse the arguments
#     args = parser.parse_args()

#     extract_reads(args.target_info, args.bam_file, args.window_size, args.sample_name)

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