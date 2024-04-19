#!/usr/bin/env python

import os
import argparse
import pandas as pd
import subprocess
import shutil
import glob
import sys

def read_target_info(target_info):
    return pd.read_csv(target_info)

def get_amplicon_sequence(amplicon_path, chr_name):
    amplicon_sequence_command = f'samtools faidx {amplicon_path} {chr_name}'
    amplicon_sequence = ''.join(subprocess.check_output(amplicon_sequence_command, shell=True).decode(sys.stdout.encoding).split('\n')[1:]).upper()
    return amplicon_sequence

def alignment_visualization(bam_dir, target_info_df, junction_type, sample_name):

    for _, row in target_info_df.iterrows():
        id = row['id']
        quant_window = row.get(f'{junction_type}_quant_window', None)
        if quant_window == '':
            continue
        quant_window_range = None
        if quant_window:
            quant_window_range = quant_window.split(':')[2]
            quant_window_lower_bound = int(quant_window_range.split('-')[0])
            quant_window_upper_bound = int(quant_window_range.split('-')[1])
            quant_window_range = f"{quant_window_lower_bound}-{quant_window_upper_bound}"

        fastq_fn = f"{id}_{junction_type}.fastq"

        if os.path.exists(fastq_fn):
            amplicon_path = f"{id}_amplicon.fasta"
            bam_path = f"{id}_alignment.bam"
            amplicon_sequence = get_amplicon_sequence(amplicon_path, f"{junction_type}_amplicon")

            allele2html_command = f"/data/tbHCA/bin/utils/bam2html.py -s {bam_path} -f {amplicon_path} -r {junction_type}_amplicon -o {id}_{junction_type}_alignment.html"
            print(allele2html_command)
            if quant_window:
                allele2html_command += f" -b {quant_window_range}"
            subprocess.call(allele2html_command, shell=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process reads and amplicons.')
    parser.add_argument('--bam', type=str, nargs = '+', help='Path to the directory containing alignments in BAM format.')
    parser.add_argument('--target_info', type=str, help='Path to the target information file.')
    parser.add_argument('--sample_name', type=str, help='Sample name')

    args = parser.parse_args()

    target_info_df = read_target_info(args.target_info)
    for junction_type in ['attL', 'attR', 'beacon', 'wt']:
        alignment_visualization(args.bam, target_info_df, junction_type, args.sample_name)