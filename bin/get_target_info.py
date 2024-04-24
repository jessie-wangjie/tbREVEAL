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
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','N':'N'}
    return "".join(complement[base] for base in reversed(seq))

def get_target_info(metadata_fn, cosmic_info,attp_reg_seq, attp_prime_seq,reference_path,cargo_reference_path, sample_name):

    ids = []
    chrs = []
    starts = []
    ends = []
    wts = []
    beacons = []
    attLs = []
    attRs = []
    beacon_quant_windows = []
    attL_quant_windows = []
    attR_quant_windows = []
    wt_quant_windows = []
    threat_tiers = []
    overlapping_features = []
    gene_names = []
    gene_strands = []
    gene_distances = []
    same_strands = []

    df = pd.read_csv(metadata_fn)

    updated_records = []

    for index, row in df.iterrows():

        id = row['Target']
        pair_ids = row['Target'].split(';')
        for pair_id in pair_ids:
            if 'CAS' in pair_id:
                query = '''
                SELECT
                    strand
                FROM
                    integrase_cryptic_attachment_site AS att_site_table
                WHERE
                    file_registry_id$ = %s
                '''
                cur.execute(query, [pair_id])
                result = cur.fetchone()
                if result:
                    strand = result
                    final_record = [pair_id,row['Chromosome'],row['Start'],row['Stop'],strand]
                    updated_records.append(final_record)
            elif 'OT' in pair_id:
                query = '''
                SELECT
                    strand
                FROM
                    spacer_off_target AS spacer_off_target
                WHERE
                    file_registry_id$ = %s
                '''
                cur.execute(query, [pair_id])
                result = cur.fetchone()
                if result:
                    strand = result
                    final_record = [pair_id,row['Chromosome'],row['Start'],row['Stop'],strand]
                    updated_records.append(final_record)
            else:
                strand = row['Strand']
                final_record = [pair_id,row['Chromosome'],row['Start'],row['Stop'],strand]
                updated_records.append(final_record)

    updated_df = pd.DataFrame(updated_records, columns=['Target', 'Chromosome', 'Start', 'Stop', 'Strand'])

    for index, row in updated_df.iterrows():
        id = row['Target']
        # OT are Cas9 mediated off-targets
        if 'OT' in id:
            target = updated_df[updated_df.Target == id]['Target'].iloc[0]
            chr = updated_df[updated_df.Target == id]['Chromosome'].iloc[0]
            start = updated_df[updated_df.Target == id]['Start'].iloc[0]
            end = updated_df[updated_df.Target == id]['Stop'].iloc[0]
            strand = updated_df[updated_df.Target == id]['Strand'].iloc[0]

            if '_' in chr:
                continue

            if 'chr' not in chr:
                chr = 'chr' + chr

            query = '''
            SELECT
                chr, start, spacer_off_target.end, strand, overlap_gene, overlap_gene_biotype, threat_tier
            FROM
                spacer_off_target AS spacer_off_target
            WHERE
                file_registry_id$ = %s
            '''
            cur.execute(query, [id])
            result = cur.fetchone()

            len_genomic_seq_to_include = 30

            if result:
                chr, start, end, strand, overlap_gene, overlap_gene_biotype, threat_tier = result

            if strand == '-':
                cut_pos = len_genomic_seq_to_include + 4
            elif strand == '+':
                cut_pos = len_genomic_seq_to_include + (int(end)-int(start)) - 1

            flank_size = 3
            ot_quant_lower_bound = cut_pos - flank_size
            ot_quant_upper_bound = cut_pos + flank_size + 1

            amplicon_sequence_command = f'samtools faidx {reference_path} {chr}:{int(start)-len_genomic_seq_to_include}-{int(end)+len_genomic_seq_to_include}'
            amplicon_sequence = ''.join(subprocess.check_output(amplicon_sequence_command, shell=True).decode(sys.stdout.encoding).split('\n')[1:]).upper()
            ot_quant_window_string = id + ':dinuc:' + str(ot_quant_lower_bound) + '-' + str(ot_quant_upper_bound) + ':0'

            ids.append(id)
            chrs.append(chr)
            starts.append(start)
            ends.append(end)
            wts.append(amplicon_sequence)
            beacons.append('')
            attLs.append('')
            attRs.append('')
            wt_quant_windows.append(ot_quant_window_string)
            beacon_quant_windows.append(ot_quant_window_string)
            attL_quant_windows.append(ot_quant_window_string)
            attR_quant_windows.append(ot_quant_window_string)
            threat_tiers.append(threat_tier)
            overlapping_features.append(overlap_gene_biotype)
            gene_names.append(overlap_gene)
            gene_strands.append(strand)
            gene_distances.append('0')
            same_strands.append('True')

        elif 'AA' in id:
            target = updated_df[updated_df.Target == id]['Target'].iloc[0]
            chr = updated_df[updated_df.Target == id]['Chromosome'].iloc[0]
            start = updated_df[updated_df.Target == id]['Start'].iloc[0]
            end = updated_df[updated_df.Target == id]['Stop'].iloc[0]
            strand = updated_df[updated_df.Target == id]['Strand'].iloc[0]

            if '_' in chr:
                continue

            if 'chr' not in chr:
                chr = 'chr' + chr

            query = '''
            SELECT
                bases, left_half, right_half, spacer_table.jmin, spacer_table.jmax, target_gene_name, direction_of_transcription
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
                bases, left_half, right_half, spacer_jmin_cutsite, spacer_jmax_cutsite, closest_gene, direction_of_transcription = result

            left_half_list = [int(x) for x in left_half.strip("[]").split(",")]
            right_half_list = [int(x) for x in right_half.strip("[]").split(",")]

            b_reg_sequence = bases[left_half_list[0]-1:left_half_list[1]+2]
            b_prime_sequence = bases[right_half_list[0]-1:right_half_list[1]]

            beacon_beginning_sequence_command = f'samtools faidx {reference_path} {chr}:{int(spacer_jmin_cutsite)-len_genomic_seq_to_include}-{int(spacer_jmin_cutsite)-1}'
            beacon_end_sequence_command = f'samtools faidx {reference_path} {chr}:{int(spacer_jmax_cutsite)}-{int(spacer_jmax_cutsite)+len_genomic_seq_to_include}'

            beacon_beginning_sequence = ''.join(subprocess.check_output(beacon_beginning_sequence_command, shell=True).decode(sys.stdout.encoding).split('\n')[1:]).upper()
            beacon_end_sequence = ''.join(subprocess.check_output(beacon_end_sequence_command, shell=True).decode(sys.stdout.encoding).split('\n')[1:]).upper()

            beacon_sequence = beacon_beginning_sequence + bases + beacon_end_sequence

            wt_command = f'samtools faidx {reference_path} {chr}:{int(spacer_jmin_cutsite)}-{int(spacer_jmax_cutsite)}'

            wt_sequence = ''.join(subprocess.check_output(wt_command, shell=True).decode(sys.stdout.encoding).split('\n')[1:]).upper()
            attL = (beacon_beginning_sequence + b_reg_sequence + attp_prime_seq).upper()
            attR = (attp_reg_seq + b_prime_sequence + beacon_end_sequence).upper()

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

            position_of_beacon_dinucleotide_lower = len_genomic_seq_to_include + len(b_reg_sequence) - length_of_dinucleotide + 1
            position_of_beacon_dinucleotide_upper = len_genomic_seq_to_include + len(b_reg_sequence)

            position_of_attL_dinucleotide_lower = len_genomic_seq_to_include + len(b_reg_sequence) - length_of_dinucleotide + 1
            position_of_attL_dinucleotide_upper = len_genomic_seq_to_include + len(b_reg_sequence)

            position_of_attR_dinucleotide_lower = len_cargo_to_include + len(attp_reg_seq) - length_of_dinucleotide + 1
            position_of_attR_dinucleotide_upper = len_cargo_to_include + len(attp_reg_seq)

            beacon_quant_lower_bound = position_of_beacon_dinucleotide_lower - 3
            beacon_quant_upper_bound = position_of_beacon_dinucleotide_upper + 3

            attL_quant_lower_bound = position_of_attL_dinucleotide_lower - 3
            attL_quant_upper_bound = position_of_attL_dinucleotide_upper + 3

            attR_quant_lower_bound = position_of_attR_dinucleotide_lower - 3
            attR_quant_upper_bound = position_of_attR_dinucleotide_upper + 3

            wt_quant_lower_bound = 1
            wt_quant_upper_bound = len(wt_sequence)

            beacon_quant_window_string = id + ':dinuc:' + str(beacon_quant_lower_bound) + '-' + str(beacon_quant_upper_bound) + ':0'
            attL_quant_window_string = id + ':dinuc:' + str(attL_quant_lower_bound) + '-' + str(attL_quant_upper_bound) + ':0'
            attR_quant_window_string = id + ':dinuc:' + str(attR_quant_lower_bound) + '-' + str(attR_quant_upper_bound) + ':0'
            wt_quant_window_string = id + ':dinuc:' + str(wt_quant_lower_bound) + '-' + str(wt_quant_upper_bound) + ':0'

            threat_tier = 'N/A'
            overlapping_feature = 'N/A'

            ids.append(id)
            chrs.append(chr)
            starts.append(start)
            ends.append(end)
            wts.append(wt_sequence)
            beacons.append(beacon_sequence)
            attLs.append(attL)
            attRs.append(attR)
            beacon_quant_windows.append(beacon_quant_window_string)
            attL_quant_windows.append(attL_quant_window_string)
            attR_quant_windows.append(attR_quant_window_string)
            wt_quant_windows.append(wt_quant_window_string)
            threat_tiers.append(threat_tier)
            overlapping_features.append(overlapping_feature)
            gene_names.append(closest_gene)
            gene_strands.append(direction_of_transcription)
            gene_distances.append('0')
            same_strands.append('True')

        elif 'CAS' in id:
            target = updated_df[updated_df.Target == id]['Target'].iloc[0]
            chr = updated_df[updated_df.Target == id]['Chromosome'].iloc[0]
            start = updated_df[updated_df.Target == id]['Start'].iloc[0]
            end = updated_df[updated_df.Target == id]['Stop'].iloc[0]
            strand = updated_df[updated_df.Target == id]['Strand'].iloc[0]

            if '_' in chr:
                continue

            if 'chr' not in chr:
                chr = 'chr' + chr

            query = '''
            SELECT
                genome_build,
                chr,
                start,
                att_site_table.end,
                central_nucleotides_start,
                central_nucleotides_seq,
                strand,
                overlap_gene,
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

            position_of_beacon_dinucleotide_lower = len_genomic_seq_to_include + (central_dinucleotide_start - start) + 1
            position_of_beacon_dinucleotide_upper = len_genomic_seq_to_include + (central_dinucleotide_start - start) + length_of_dinucleotide

            position_of_attL_dinucleotide_lower = len_genomic_seq_to_include + (central_dinucleotide_start - start) + 1
            position_of_attL_dinucleotide_upper = len_genomic_seq_to_include + (central_dinucleotide_start - start) + length_of_dinucleotide

            position_of_attR_dinucleotide_lower = len_cargo_to_include + len(attp_reg_seq) - length_of_dinucleotide + 1
            position_of_attR_dinucleotide_upper = len_cargo_to_include + len(attp_reg_seq)

            beacon_quant_lower_bound = position_of_beacon_dinucleotide_lower - 3
            beacon_quant_upper_bound = position_of_beacon_dinucleotide_upper + 3

            attL_quant_lower_bound = position_of_attL_dinucleotide_lower - 3
            attL_quant_upper_bound = position_of_attL_dinucleotide_upper + 3

            attR_quant_lower_bound = position_of_attR_dinucleotide_lower - 3
            attR_quant_upper_bound = position_of_attR_dinucleotide_upper + 3

            beacon_quant_window_string = id + ':dinuc:' + str(beacon_quant_lower_bound) + '-' + str(beacon_quant_upper_bound) + ':0'
            attL_quant_window_string = id + ':dinuc:' + str(attL_quant_lower_bound) + '-' + str(attL_quant_upper_bound) + ':0'
            attR_quant_window_string = id + ':dinuc:' + str(attR_quant_lower_bound) + '-' + str(attR_quant_upper_bound) + ':0'
            wt_quant_window_string = beacon_quant_window_string

            if closest_gene is None:
                gene_list = []
            else:
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

            ids.append(id)
            chrs.append(chr)
            starts.append(start)
            ends.append(end)
            wts.append(wt_sequence)
            beacons.append(cryptic_beacon)
            attLs.append(cryptic_AttL_sequence)
            attRs.append(cryptic_AttR_sequence)
            beacon_quant_windows.append(beacon_quant_window_string)
            attL_quant_windows.append(attL_quant_window_string)
            attR_quant_windows.append(attR_quant_window_string)
            wt_quant_windows.append(wt_quant_window_string)
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
        'beacon_quant_window': beacon_quant_windows,
        'attL_quant_window': attL_quant_windows,
        'attR_quant_window': attR_quant_windows,
        'wt_quant_window': wt_quant_windows,
        'threat_tier': threat_tiers,
        'same_strand': same_strands,
        'overlapping_feature': overlapping_features,
        'gene_name': gene_names,
        'gene_strand': gene_strands,
        'gene_distance': gene_distances
    }

    df = pd.DataFrame(data)
    cosmic_info_df = pd.read_csv(cosmic_info,sep=',')

    cosmic_info_df = cosmic_info_df[['Gene Symbol','Role in Cancer']]

    final_df = pd.merge(df,cosmic_info_df,left_on='gene_name',right_on='Gene Symbol',how='left')

    final_df['Cosmic Gene'] = final_df['Role in Cancer'].apply(lambda x: 'No' if pd.isna(x) else 'Yes')

    final_df = final_df.drop('Gene Symbol', axis=1)

    final_df.to_csv(f'{sample_name}_target_info.csv', index=False)
    return(final_df)

if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(description="Extract reads overlapping a certain genomic position.")

    # Add the arguments
    parser.add_argument("--metadata", required=True, type=str, help="Metadata file")
    parser.add_argument("--cosmic_info", required=True, type=str, help="COSMIC info")
    parser.add_argument("--attp_reg", required=True, type=str, help="P sequence")
    parser.add_argument("--attp_prime", required=True, type=str, help="P prime sequence")
    parser.add_argument("--reference", required=True, type=str, help="P prime sequence")
    parser.add_argument("--cargo", required=True, type=str, help="P prime sequence")
    parser.add_argument("--sample_name", required=True, type=str, help="P prime sequence")
    # Parse the arguments
    args = parser.parse_args()

    get_target_info(args.metadata, args.cosmic_info, args.attp_reg, args.attp_prime, args.reference, args.cargo, args.sample_name)