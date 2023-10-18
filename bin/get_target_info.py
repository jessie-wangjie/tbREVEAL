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

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement[base] for base in reversed(seq))

def get_target_info(metadata_fn, attp_reg_seq, attp_prime_seq):

    ids = []
    chrs = []
    starts = []
    ends = []
    wts = []
    beacons = []
    attLs = []
    attRs = []
    attL_quant_windows = []
    attR_quant_windows = []
    threat_tiers = []
    overlapping_features = []
    gene_names = []
    gene_strands = []
    gene_distances = []
    same_strands = []

    df = pd.read_csv(metadata_fn)

    for index, row in df.iterrows():
        id = row['Target']
        if 'AA' in id:
            target = df[df.Target == id]['Target'].iloc[0]
            chr = df[df.Target == id]['Chromosome'].iloc[0]
            start = df[df.Target == id]['Start'].iloc[0]
            end = df[df.Target == id]['Stop'].iloc[0]
            strand = df[df.Target == id]['Strand'].iloc[0]

            query = '''
            SELECT
                bases, left_half, right_half, spacer_table.jmin, target_gene_name, direction_of_transcription
            FROM
                atg_atg AS atg_table
            JOIN 
                attachment_sequence as att_seq ON att_seq.id = atg_table.expected_beacon
            JOIN 
                dna_sequence AS dna_sequence ON dna_sequence.id = atg_table.expected_beacon
            JOIN
               spacer_pair as spacer_table ON atg_table.spacer_pair = spacer_table.id
            JOIN
                target_gene as gene_table ON gene_table.id = spacer_table.target_gene
            WHERE
                atg_table.file_registry_id$ = %s
            '''
            cur.execute(query, [id])
            result = cur.fetchone()

            len_cargo_to_include = 20
            len_genomic_seq_to_include = 20

            if result:
                bases, left_half, right_half, spacer_jmin_cutsite, closest_gene, direction_of_transcription = result
                cargo_reference_path = f'/data/references/AAVG097.fa'
                reference_path = f'/data/references/hg38.fa'

            chr = 'chr' + str(chr)
            
            left_half_list = [int(x) for x in left_half.strip("[]").split(",")]
            right_half_list = [int(x) for x in right_half.strip("[]").split(",")]

            b_reg_sequence = bases[left_half_list[0]-1:left_half_list[1]+2]
            b_prime_sequence = bases[right_half_list[0]-1:right_half_list[1]]

            wt_command = f'samtools faidx {reference_path} {chr}:{int(spacer_jmin_cutsite)}-{int(spacer_jmin_cutsite)+len(bases)}'

            wt_sequence = ''.join(subprocess.check_output(wt_command, shell=True).decode(sys.stdout.encoding).split('\n')[1:]).upper()
            attL = (b_reg_sequence + attp_prime_seq).upper()
            attR = (attp_reg_seq + b_prime_sequence).upper()

            # store cargo sequence so we can search for subsequence locations
            with open(cargo_reference_path, 'r') as f:
                cargo_sequence = ''.join(line.strip() for line in f if not line.startswith('>')).upper()
            
            cargo_pprime_end_loc = cargo_sequence.find(attp_prime_seq) + len(attp_prime_seq)
            cargo_preg_start_loc = cargo_sequence.find(attp_reg_seq)

            cargo_20bp_before_preg = cargo_sequence[cargo_preg_start_loc-len_cargo_to_include:cargo_preg_start_loc]
            cargo_20bp_after_pprime = cargo_sequence[cargo_pprime_end_loc:cargo_pprime_end_loc+len_cargo_to_include]

            attL = attL + cargo_20bp_after_pprime
            attR = cargo_20bp_before_preg + attR

            length_of_dinucleotide = 2

            position_of_attL_dinucleotide_lower = len(b_reg_sequence) - length_of_dinucleotide + 1
            position_of_attL_dinucleotide_upper = len(b_reg_sequence)

            position_of_attR_dinucleotide_lower = len_cargo_to_include + len(attp_reg_seq) - length_of_dinucleotide + 1
            position_of_attR_dinucleotide_upper = len_cargo_to_include + len(attp_reg_seq)

            attL_quant_lower_bound = position_of_attL_dinucleotide_lower - 3
            attL_quant_upper_bound = position_of_attL_dinucleotide_upper + 3

            attR_quant_lower_bound = position_of_attR_dinucleotide_lower - 3
            attR_quant_upper_bound = position_of_attR_dinucleotide_upper + 3

            attL_quant_window_string = id + ':dinuc:' + str(attL_quant_lower_bound) + '-' + str(attL_quant_upper_bound) + ':0'
            attR_quant_window_string = id + ':dinuc:' + str(attR_quant_lower_bound) + '-' + str(attR_quant_upper_bound) + ':0'

            threat_tier = 'N/A'
            overlapping_feature = 'N/A'
            
            ids.append(id)
            chrs.append(chr)
            starts.append(start)
            ends.append(end)
            wts.append(wt_sequence)
            beacons.append(bases)
            attLs.append(attL)
            attRs.append(attR)
            attL_quant_windows.append(attL_quant_window_string)
            attR_quant_windows.append(attR_quant_window_string)
            threat_tiers.append(threat_tier)
            overlapping_features.append(overlapping_feature)
            gene_names.append(closest_gene)
            gene_strands.append(direction_of_transcription)
            gene_distances.append('0')
            same_strands.append('True')

        elif 'CAS' in id:
            target = df[df.Target == id]['Target'].iloc[0]
            chr = df[df.Target == id]['Chromosome'].iloc[0]
            start = df[df.Target == id]['Start'].iloc[0]
            end = df[df.Target == id]['Stop'].iloc[0]
            strand = df[df.Target == id]['Strand'].iloc[0]

            query = '''
            SELECT
                genome_build,
                chr,
                start,
                att_site_table.end,
                central_nucleotides_start,
                central_nucleotides_seq,
                strand,
                closest_gene,
                threat_tier, 
                overlapping_feature
            FROM
                integrase_cryptic_attachment_site AS att_site_table
            WHERE
                file_registry_id$ = %s
            '''

            cur.execute(query, [id])
            result = cur.fetchone()

            len_cargo_to_include = 20
            len_genomic_seq_to_include = 20

            if result:
                genome_build, chr, start, end, central_dinucleotide_start, central_nucleotides_seq, strand,closest_gene,threat_tier,overlapping_feature = result
                reference_path = f'/data/references/hg38.fa'
                cargo_reference_path = f'/data/references/AAVG097.fa'
                if strand == '+':
                    cryptic_b_reg_command = f'samtools faidx {reference_path} {chr}:{start-len_genomic_seq_to_include}-{central_dinucleotide_start + len(central_nucleotides_seq) - 1}'
                    cryptic_b_prime_command = f'samtools faidx {reference_path} {chr}:{central_dinucleotide_start + len(central_nucleotides_seq)}-{end+len_genomic_seq_to_include}'
                else:
                    cryptic_b_reg_command = f'samtools faidx {reference_path} {chr}:{start-len_genomic_seq_to_include}-{central_dinucleotide_start - 1}'
                    cryptic_b_prime_command = f'samtools faidx {reference_path} {chr}:{central_dinucleotide_start}-{end+len_genomic_seq_to_include}'
            
            cryptic_b_reg_sequence = ''.join(subprocess.check_output(cryptic_b_reg_command, shell=True).decode(sys.stdout.encoding).split('\n')[1:]).upper()
            cryptic_b_prime_sequence = ''.join(subprocess.check_output(cryptic_b_prime_command, shell=True).decode(sys.stdout.encoding).split('\n')[1:]).upper()

            if strand == '-':
                cryptic_b_reg_sequence, cryptic_b_prime_sequence = cryptic_b_prime_sequence, cryptic_b_reg_sequence
                cryptic_b_reg_sequence = reverse_complement(cryptic_b_reg_sequence)
                cryptic_b_prime_sequence = reverse_complement(cryptic_b_prime_sequence)

            cryptic_beacon = (cryptic_b_reg_sequence + cryptic_b_prime_sequence).upper()
            wt_sequence = cryptic_beacon
            cryptic_AttL_sequence = (cryptic_b_reg_sequence + attp_prime_seq).upper()
            cryptic_AttR_sequence = (attp_reg_seq + cryptic_b_prime_sequence).upper()

            # store cargo sequence so we can search for subsequence locations
            with open(cargo_reference_path, 'r') as f:
                cargo_sequence = ''.join(line.strip() for line in f if not line.startswith('>')).upper()
            
            cargo_pprime_end_loc = cargo_sequence.find(attp_prime_seq) + len(attp_prime_seq)
            cargo_preg_start_loc = cargo_sequence.find(attp_reg_seq)

            cargo_20bp_before_preg = cargo_sequence[cargo_preg_start_loc-len_cargo_to_include:cargo_preg_start_loc]
            cargo_20bp_after_pprime = cargo_sequence[cargo_pprime_end_loc:cargo_pprime_end_loc+len_cargo_to_include]

            cryptic_AttL_sequence = cryptic_AttL_sequence + cargo_20bp_after_pprime
            cryptic_AttR_sequence = cargo_20bp_before_preg + cryptic_AttR_sequence

            length_of_dinucleotide = len(central_nucleotides_seq)

            position_of_attL_dinucleotide_lower = len_genomic_seq_to_include + (central_dinucleotide_start - start) + 1 
            position_of_attL_dinucleotide_upper = len_genomic_seq_to_include + (central_dinucleotide_start - start) + length_of_dinucleotide

            position_of_attR_dinucleotide_lower = len_cargo_to_include + len(attp_reg_seq) - length_of_dinucleotide + 1
            position_of_attR_dinucleotide_upper = len_cargo_to_include + len(attp_reg_seq)

            attL_quant_lower_bound = position_of_attL_dinucleotide_lower - 3
            attL_quant_upper_bound = position_of_attL_dinucleotide_upper + 3

            attR_quant_lower_bound = position_of_attR_dinucleotide_lower - 3
            attR_quant_upper_bound = position_of_attR_dinucleotide_upper + 3

            attL_quant_window_string = id + ':dinuc:' + str(attL_quant_lower_bound) + '-' + str(attL_quant_upper_bound) + ':0'
            attR_quant_window_string = id + ':dinuc:' + str(attR_quant_lower_bound) + '-' + str(attR_quant_upper_bound) + ':0'

            gene_list = closest_gene.split(";")

            gene_name = [el.split(",")[0].strip() for el in gene_list]
            gene_strand = [el.split(",")[1].strip() for el in gene_list]
            gene_distance = [el.split(",")[2].strip() for el in gene_list]


            if strand in gene_strand:
                same_strand="True"
            else:
                same_strand="False"

            gene_name = "/".join(gene_name)
            gene_strand = "/".join(gene_strand)
            gene_distance = "/".join(gene_distance)

            print(gene_list)
            print(gene_name, gene_strand,gene_distance)

            ids.append(id)
            chrs.append(chr)
            starts.append(start)
            ends.append(end)
            wts.append(wt_sequence)
            beacons.append(cryptic_beacon)
            attLs.append(cryptic_AttL_sequence)
            attRs.append(cryptic_AttR_sequence)
            attL_quant_windows.append(attL_quant_window_string)
            attR_quant_windows.append(attR_quant_window_string)
            threat_tiers.append(threat_tier)
            overlapping_features.append(overlapping_feature)
            gene_names.append(gene_name)
            gene_strands.append(gene_strand)
            gene_distances.append(gene_distance)
            same_strands.append(same_strand)
    
    
    data = {
        'id': ids,
        'chromosome': chrs,
        'start': starts,
        'end': ends,
        'wt': wts,
        'beacon': beacons,
        'attL': attLs,
        'attR': attRs,
        'attL_quant_window': attL_quant_windows,
        'attR_quant_window': attR_quant_windows,
        'threat_tier': threat_tiers,
        'same_strand': same_strands,
        'overlapping_feature': overlapping_features,
        'gene_name': gene_names,
        'gene_strand': gene_strands,
        'gene_distance': gene_distances
    }

    df = pd.DataFrame(data)

    df.to_csv('target_info.csv', index=False)
    return(df)

if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(description="Extract reads overlapping a certain genomic position.")

    # Add the arguments
    parser.add_argument("--metadata", required=True, type=str, help="Metadata file")
    parser.add_argument("--attp_reg", required=True, type=str, help="P sequence")
    parser.add_argument("--attp_prime", required=True, type=str, help="P prime sequence")
    # Parse the arguments
    args = parser.parse_args()

    get_target_info(args.metadata, args.attp_reg, args.attp_prime)