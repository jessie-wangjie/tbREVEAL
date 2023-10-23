#!/usr/bin/env python

import subprocess
import pandas as pd
import argparse
import os
import sys
import glob
import shutil
from utils.base import *

def run_cs2(amplicon_dir,target_info, junction_type, sample_name):

    target_info_df = pd.read_csv(target_info)

    if junction_type == 'attL':
        reads_dir = f'{sample_name}_attL_extracted_reads'
        for index, row in target_info_df.iterrows():
            id = row['id']
            quant_window = row['attL_quant_window']
            quant_window_range = quant_window.split(':')[2]
            quant_window_lower_bound_adjusted_for_cs2 = int(quant_window_range.split('-')[0]) - 1
            quant_window_upper_bound_adjusted_for_cs2 = int(quant_window_range.split('-')[1]) - 1
            quant_window_range = str(quant_window_lower_bound_adjusted_for_cs2) + '-' + str(quant_window_upper_bound_adjusted_for_cs2)
            fastq_fn = f"{reads_dir}/{id}_attL.fastq"
            if os.path.exists(fastq_fn):
                print('hi')
                # get amplicon
                amplicon_path = f"{amplicon_dir}/{id}_amplicon.fasta"
                chr = "attL_amplicon"
                amplicon_sequence_command = f'samtools faidx {amplicon_path} {chr}'
                amplicon_sequence = ''.join(subprocess.check_output(amplicon_sequence_command, shell=True).decode(sys.stdout.encoding).split('\n')[1:]).upper()
                subprocess.run(["CRISPResso", "--fastq_r1", fastq_fn, "--amplicon_seq", amplicon_sequence, "--amplicon_name",id,"--name", id + "_attL","--write_detailed_allele_table","--bam_output","--quantification_window_coordinates",quant_window_range,"--exclude_bp_from_left","0","--exclude_bp_from_right","0","--amplicon_min_alignment_score","5"])
                os.makedirs(os.path.join("CRISPResso_on_" + id + "_attL", "cs2_alignment_html"), exist_ok=True)
                path_to_cs2_output = "CRISPResso_on_" + id + "_attL"
                allele2html_command = "/data/tbHCA/bin/utils/allele2html.py -f %s -r %s -b %s" % (path_to_cs2_output+'/', id, quant_window)
                subprocess.call(allele2html_command, shell=True)
        
        pattern = "CRISPResso_on_*_attL"
        destination_folder = f"{sample_name}_cs2_attL"

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
        reads_dir = f'{sample_name}_attR_extracted_reads'
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
                subprocess.run(["CRISPResso", "--fastq_r1", fastq_fn, "--amplicon_seq", amplicon_sequence, "--amplicon_name",id,"--name",id + "_attR","--write_detailed_allele_table","--bam_output","--quantification_window_coordinates",quant_window_range,"--exclude_bp_from_left","0","--exclude_bp_from_right","0","--amplicon_min_alignment_score","5"])
                os.makedirs(os.path.join("CRISPResso_on_" + id + "_attR", "cs2_alignment_html"), exist_ok=True)
                path_to_cs2_output = "CRISPResso_on_" + id + "_attR"
                allele2html_command = "/data/tbHCA/bin/utils/allele2html.py -f %s -r %s -b %s" % (path_to_cs2_output+'/', id, quant_window)
                subprocess.call(allele2html_command, shell=True)
        pattern = "CRISPResso_on_*_attR"
        destination_folder = f"{sample_name}_cs2_attR"

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

    elif junction_type == 'beacon':
        reads_dir = f'{sample_name}_beacon_extracted_reads'
        for index, row in target_info_df.iterrows():
            id = row['id']
            #quant_window = row['beacon_quant_window']
            #quant_window_range = quant_window.split(':')[2]
            #quant_window_lower_bound_adjusted_for_cs2 = int(quant_window_range.split('-')[0]) - 1
            #quant_window_upper_bound_adjusted_for_cs2 = int(quant_window_range.split('-')[1]) - 1
            #quant_window_range = str(quant_window_lower_bound_adjusted_for_cs2) + '-' + str(quant_window_upper_bound_adjusted_for_cs2)
            fastq_fn = reads_dir + '/' + id + '_beacon.fastq'
            if os.path.exists(fastq_fn):
                # get amplicon
                amplicon_path = amplicon_dir +'/'+id+'_amplicon.fasta'
                chr = "beacon_amplicon"
                amplicon_sequence_command = f'samtools faidx {amplicon_path} {chr}'
                amplicon_sequence = ''.join(subprocess.check_output(amplicon_sequence_command, shell=True).decode(sys.stdout.encoding).split('\n')[1:]).upper()
                # subprocess.run(["CRISPResso", "--fastq_r1", fastq_fn, "--amplicon_seq", amplicon_sequence, "--amplicon_name",id,"--name",id + "_beacon","--write_detailed_allele_table","--bam_output","--quantification_window_coordinates",quant_window_range,"--exclude_bp_from_left","0","--exclude_bp_from_right","0","--amplicon_min_alignment_score","5"])
                subprocess.run(["CRISPResso", "--fastq_r1", fastq_fn, "--amplicon_seq", amplicon_sequence, "--amplicon_name",id,"--name",id + "_beacon","--write_detailed_allele_table","--bam_output","--exclude_bp_from_left","0","--exclude_bp_from_right","0","--amplicon_min_alignment_score","5"])
                os.makedirs(os.path.join("CRISPResso_on_" + id + "_beacon", "cs2_alignment_html"), exist_ok=True)
                path_to_cs2_output = "CRISPResso_on_" + id + "_beacon"
                allele2html_command = "/data/tbHCA/bin/utils/allele2html.py -f %s -r %s" % (path_to_cs2_output+'/', id)
                subprocess.call(allele2html_command, shell=True)
        pattern = "CRISPResso_on_*_beacon"
        destination_folder = f"{sample_name}_cs2_beacon"

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
    parser.add_argument('--sample_name', type=str, help='Sample name')

    args = parser.parse_args()

    os.makedirs(f'{args.sample_name}_cs2_attL',exist_ok=True)
    os.makedirs(f'{args.sample_name}_cs2_attR',exist_ok=True)
    os.makedirs(f'{args.sample_name}_cs2_beacon',exist_ok=True)

    run_cs2(args.amplicon_dir, args.target_info, 'attL', args.sample_name)
    run_cs2(args.amplicon_dir, args.target_info, 'attR', args.sample_name)
    run_cs2(args.amplicon_dir, args.target_info, 'beacon', args.sample_name)


