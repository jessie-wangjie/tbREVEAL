#!/usr/bin/env python

import subprocess
import pandas as pd
import argparse
import os
import sys
import glob
import shutil
from utils.base import *

def run_cs2(amplicon_dir,target_info, junction_type):

    target_info_df = pd.read_csv(target_info)

    if junction_type == 'attL':
        reads_dir = 'attL_extracted_reads'
        for index, row in target_info_df.iterrows():
            id = row['id']
            quant_window = row['attL_quant_window']
            quant_window_range = quant_window.split(':')[2]
            quant_window_lower_bound_adjusted_for_cs2 = int(quant_window_range.split('-')[0]) - 1
            quant_window_upper_bound_adjusted_for_cs2 = int(quant_window_range.split('-')[1]) - 1
            quant_window_range = str(quant_window_lower_bound_adjusted_for_cs2) + '-' + str(quant_window_upper_bound_adjusted_for_cs2)
            fastq_fn = reads_dir + '/' + id + '_attL.fastq'
            if os.path.exists(fastq_fn):
                # get amplicon
                amplicon_path = amplicon_dir +'/'+id+'_amplicon.fasta'
                chr = "attL_amplicon"
                amplicon_sequence_command = f'samtools faidx {amplicon_path} {chr}'
                amplicon_sequence = ''.join(subprocess.check_output(amplicon_sequence_command, shell=True).decode(sys.stdout.encoding).split('\n')[1:]).upper()
                subprocess.run(["CRISPResso", "--fastq_r1", fastq_fn, "--amplicon_seq", amplicon_sequence, "--amplicon_name",id,"--name", id + "_attL","--write_detailed_allele_table","--bam_output","--quantification_window_coordinates",quant_window_range,"--exclude_bp_from_left","0","--exclude_bp_from_right","0"])
                os.makedirs(os.path.join("CRISPResso_on_" + id + "_attL", "cs2_alignment_html"), exist_ok=True)
                path_to_cs2_output = "CRISPResso_on_" + id + "_attL"
                allele2html_command = "allele2html.py -f %s -r %s -b %s" % (path_to_cs2_output+'/', id, quant_window)
                subprocess.call(allele2html_command, shell=True)
        
        pattern = "CRISPResso_on_*_attL"
        destination_folder = "cs2_attL"

        # Find all files matching the pattern
        matching_files = glob.glob(pattern)
        # Move each file to the destination folder
        for file_path in matching_files:
            # Extract the file name from the path
            file_name = os.path.basename(file_path)
            
            # Create the destination path
            destination_path = os.path.join(destination_folder, file_name)
            
            # Move the file to the destination folder
            shutil.move(file_path, destination_path)


    elif junction_type == 'attR':
        reads_dir = 'attR_extracted_reads'
        for index, row in target_info_df.iterrows():
            id = row['id']
            quant_window = row['attR_quant_window']
            quant_window_range = quant_window.split(':')[2]
            quant_window_lower_bound_adjusted_for_cs2 = int(quant_window_range.split('-')[0]) - 1
            quant_window_upper_bound_adjusted_for_cs2 = int(quant_window_range.split('-')[1]) - 1
            quant_window_range = str(quant_window_lower_bound_adjusted_for_cs2) + '-' + str(quant_window_upper_bound_adjusted_for_cs2)
            fastq_fn = reads_dir + '/' + id + '_attR.fastq'
            if os.path.exists(fastq_fn):
                # get amplicon
                amplicon_path = amplicon_dir +'/'+id+'_amplicon.fasta'
                chr = "attR_amplicon"
                amplicon_sequence_command = f'samtools faidx {amplicon_path} {chr}'
                amplicon_sequence = ''.join(subprocess.check_output(amplicon_sequence_command, shell=True).decode(sys.stdout.encoding).split('\n')[1:]).upper()
                subprocess.run(["CRISPResso", "--fastq_r1", fastq_fn, "--amplicon_seq", amplicon_sequence, "--amplicon_name",id,"--name",id + "_attR","--write_detailed_allele_table","--bam_output","--quantification_window_coordinates",quant_window_range,"--exclude_bp_from_left","0","--exclude_bp_from_right","0"])
                os.makedirs(os.path.join("CRISPResso_on_" + id + "_attR", "cs2_alignment_html"), exist_ok=True)
                path_to_cs2_output = "CRISPResso_on_" + id + "_attR"
                allele2html_command = "allele2html.py -f %s -r %s -b %s" % (path_to_cs2_output+'/', id, quant_window)
                subprocess.call(allele2html_command, shell=True)
        pattern = "CRISPResso_on_*_attR"
        destination_folder = "cs2_attR"

        # Find all files matching the pattern
        matching_files = glob.glob(pattern)
        # Move each file to the destination folder
        for file_path in matching_files:
            # Extract the file name from the path
            file_name = os.path.basename(file_path)
            
            # Create the destination path
            destination_path = os.path.join(destination_folder, file_name)
            
            # Move the file to the destination folder
            shutil.move(file_path, destination_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process reads and amplicons.')
    #parser.add_argument('--reads_dir', type=str, help='Path to the directory containing reads.')
    parser.add_argument('--amplicon_dir', type=str, help='Path to the directory containing amplicons.')
    parser.add_argument('--target_info', type=str, help='Path to the target information file.')

    args = parser.parse_args()

    os.mkdir("cs2_attL")
    os.mkdir("cs2_attR")

    run_cs2(args.amplicon_dir, args.target_info, 'attL')
    run_cs2(args.amplicon_dir, args.target_info, 'attR')


