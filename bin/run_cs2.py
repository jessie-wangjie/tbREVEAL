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

def run_cs2_for_type(amplicon_dir, target_info_df, junction_type, sample_name):

    reads_dir = f'{sample_name}_{junction_type}_extracted_reads'
    
    for _, row in target_info_df.iterrows():
        id = row['id']
        quant_window = row.get(f'{junction_type}_quant_window', None)
        quant_window_range = None
        if quant_window:
            quant_window_range = quant_window.split(':')[2]
            quant_window_lower_bound = int(quant_window_range.split('-')[0]) - 1
            quant_window_upper_bound = int(quant_window_range.split('-')[1]) - 1
            quant_window_range = f"{quant_window_lower_bound}-{quant_window_upper_bound}"

        fastq_fn = f"{reads_dir}/{id}_{junction_type}.fastq"
        
        if os.path.exists(fastq_fn):
            amplicon_path = f"{amplicon_dir}/{id}_amplicon.fasta"
            amplicon_sequence = get_amplicon_sequence(amplicon_path, f"{junction_type}_amplicon")
            
            crispresso_command = [
                "CRISPResso", "--fastq_r1", fastq_fn, 
                "--amplicon_seq", amplicon_sequence, 
                "--amplicon_name", id, 
                "--name", f"{id}_{junction_type}", 
                "--write_detailed_allele_table", 
                "--bam_output", 
                "--exclude_bp_from_left", "0", 
                "--exclude_bp_from_right", "0", 
                "--amplicon_min_alignment_score", "5"
            ]
            if quant_window_range:
                crispresso_command.extend(["--quantification_window_coordinates", quant_window_range])
            
            subprocess.run(crispresso_command)
            
            os.makedirs(f"CRISPResso_on_{id}_{junction_type}/cs2_alignment_html", exist_ok=True)
            
            allele2html_command = f"/data/tbHCA/bin/utils/allele2html.py -f CRISPResso_on_{id}_{junction_type}/ -r {id}"
            if quant_window:
                allele2html_command += f" -b {quant_window}"
            subprocess.call(allele2html_command, shell=True)

    pattern = f"CRISPResso_on_*_{junction_type}"
    destination_folder = f"{sample_name}_cs2_{junction_type}"

    matching_files = glob.glob(pattern)
    for file_path in matching_files:
        file_name = os.path.basename(file_path)
        destination_path = os.path.join(destination_folder, file_name)
        shutil.move(file_path, destination_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process reads and amplicons.')
    parser.add_argument('--amplicon_dir', type=str, help='Path to the directory containing amplicons.')
    parser.add_argument('--target_info', type=str, help='Path to the target information file.')
    parser.add_argument('--sample_name', type=str, help='Sample name')

    args = parser.parse_args()

    os.makedirs(f'{args.sample_name}_cs2_attL', exist_ok=True)
    os.makedirs(f'{args.sample_name}_cs2_attR', exist_ok=True)
    os.makedirs(f'{args.sample_name}_cs2_beacon', exist_ok=True)

    target_info_df = read_target_info(args.target_info)
    for junction_type in ['attL', 'attR', 'beacon']:
        run_cs2_for_type(args.amplicon_dir, target_info_df, junction_type, args.sample_name)
