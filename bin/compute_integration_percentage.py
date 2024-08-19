#!/usr/bin/env python

import argparse
import os
from collections import defaultdict
import pandas as pd
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import numpy as np


def has_indel_in_range(alignment, start_range, end_range):
    """
    Check if the given alignment has an indel (insertion or deletion) within a specific range on the reference.

    Parameters:
    alignment (pysam.AlignedSegment): An alignment from a BAM/SAM file.
    start_range (int): The start position of the range on the reference (0-based).
    end_range (int): The end position of the range on the reference (0-based, inclusive).

    Returns:
    bool: True if the alignment has an indel within the specified range, False otherwise.
    """

    if not alignment.cigartuples:
        return False

    ref_pos = alignment.reference_start
    for op, length in alignment.cigartuples:
        if op == 0 or op == 2 or op == 3:  # M, D, or N (skip or reference skip)
            # Before checking for indels, see if we've reached the range of interest
            if ref_pos > end_range:
                break  # Past the range of interest, stop checking
            if (op == 2 and start_range <= ref_pos <= end_range) or (
                op == 2 and ref_pos + length - 1 >= start_range and ref_pos <= end_range
            ):
                # Deletion falls within the range
                return True
            ref_pos += length
        elif op == 1:  # I
            if start_range <= ref_pos <= end_range:
                # Insertion at the start of the range or before moving into the range
                return True
        # For simplicity, this doesn't handle S (soft clip), H (hard clip), P (padding), and X/= (sequence match/mismatch)
        # as they don't affect reference position in the same way. Adjust if those are relevant for your use case.

    return False


def compute_integration_percentage(target_info, alignment_dir, sample_name):
    target_info_df = pd.read_csv(target_info)
    integration_dict = {}

    def get_umi(read):
        return read.query_name

    def calculate_soft_clipped_indices(read):
        start, end = 0, len(read.seq)
        for i, (operation, length) in enumerate(read.cigartuples):
            if i == 0 and operation == 4:
                start += length
            elif operation == 4:
                end -= length
        return start, end

    target_info_df["Role in Cancer"] = (
        target_info_df["Role in Cancer"].fillna("").astype(str)
    )

    for index, row in target_info_df.iterrows():
        alignment_file = f"{row['id']}_alignment.bam"

        attL_records = {}
        attR_records = {}
        beacon_records = {}
        wt_records = {}

        if not isinstance(row["beacon"], float):
            full_attb_sequence = row["beacon"][20 : len(row["beacon"]) - 21]
        else:
            full_attb_sequence = ""

        categories = {
            "complete_attL": set(),
            "complete_attR": set(),
            "partial_beacon": set(),
            "complete_beacon": set(),
            "wt": set(),
        }

        indel_counter = 0
        read_counter = 0

        beacon_quant_lower = int(row["beacon_quant_window"].split(":")[2].split("-")[0])
        beacon_quant_upper = int(row["beacon_quant_window"].split(":")[2].split("-")[1])

        wt_quant_lower = int(row["wt_quant_window"].split(":")[2].split("-")[0])
        wt_quant_upper = int(row["wt_quant_window"].split(":")[2].split("-")[1])

        if os.path.exists(alignment_file):
            with pysam.AlignmentFile(alignment_file, "rb") as alignment:
                for read in alignment:
                    if read.flag in [
                        2048,
                        2064,
                        4,
                    ]:  # Supplementary, Secondary, or Unaligned
                        continue
                    num_mismatch = read.get_tag("NM")
                    umi = get_umi(read)
                    start, end = calculate_soft_clipped_indices(read)
                    seq = read.seq[start:end]
                    qual = read.query_qualities[start:end]

                    seq_record = SeqRecord(
                        Seq(seq),
                        id=read.qname,
                        description="",
                        letter_annotations={"phred_quality": qual},
                    )
                    alignment_start, alignment_end = (
                        read.reference_start + 1,
                        read.reference_end,
                    )

                    rules = [
                        ## CAS
                        (
                            read.reference_name == "attL_amplicon"
                            and "CAS" in row["id"]
                            and alignment_start <= 21
                            and alignment_end >= 69
                            and num_mismatch <= 5,
                            ["complete_attL"],
                        ),
                        (
                            read.reference_name == "attR_amplicon"
                            and "CAS" in row["id"]
                            and alignment_start <= 21
                            and alignment_end >= 69
                            and num_mismatch <= 5,
                            ["complete_attR"],
                        ),
                        (
                            read.reference_name == "beacon_amplicon"
                            and "CAS" in row["id"]
                            and alignment_start <= 21
                            and alignment_end >= 66,
                            ["complete_beacon"],
                        ),
                        # AA
                        (
                            read.reference_name == "attL_amplicon"
                            and "AA" in row["id"]
                            and (alignment_start <= 21 and alignment_end >= 65)
                            and num_mismatch <= 5,
                            ["complete_attL"],
                        ),
                        (
                            read.reference_name == "attR_amplicon"
                            and "AA" in row["id"]
                            and (alignment_start <= 21 and alignment_end >= 65)
                            and num_mismatch <= 5,
                            ["complete_attR"],
                        ),
                        (
                            read.reference_name == "beacon_amplicon"
                            and "AA" in row["id"]
                            and (
                                alignment_start <= 21
                                and alignment_end >= (21 + len(full_attb_sequence) - 1)
                            )
                            and (full_attb_sequence not in seq),
                            ["partial_beacon"],
                        ),
                        (
                            read.reference_name == "beacon_amplicon"
                            and "AA" in row["id"]
                            and (
                                alignment_start <= 21
                                and alignment_end >= (21 + len(full_attb_sequence) - 1)
                            )
                            and (full_attb_sequence in seq),
                            ["complete_beacon"],
                        ),
                        (read.reference_name == "wt_amplicon", ["wt"]),
                    ]

                    read_counter += 1

                    for condition, cat_keys in rules:
                        if condition:
                            for key in cat_keys:
                                categories[key].add(umi)
                                if "complete_attL" in key:
                                    attL_records[umi] = seq_record
                                    if has_indel_in_range(
                                        read, beacon_quant_lower, beacon_quant_upper
                                    ):
                                        indel_counter += 1
                                elif "complete_attR" in key:
                                    attR_records[umi] = seq_record
                                    if has_indel_in_range(
                                        read, beacon_quant_lower, beacon_quant_upper
                                    ):
                                        indel_counter += 1
                                elif "complete_beacon" in key:
                                    beacon_records[umi] = seq_record
                                    if has_indel_in_range(
                                        read, beacon_quant_lower, beacon_quant_upper
                                    ):
                                        indel_counter += 1
                                elif "wt" in key:
                                    wt_records[umi] = seq_record
                                    if has_indel_in_range(
                                        read, wt_quant_lower, wt_quant_upper
                                    ):
                                        indel_counter += 1

            if read_counter == 0:
                indel_percent = 0
            else:
                indel_percent = 100 * indel_counter / read_counter

            for prefix, records in [
                (f"{sample_name}_attL", attL_records),
                (f"{sample_name}_attR", attR_records),
                (f"{sample_name}_beacon", beacon_records),
                (f"{sample_name}_wt", wt_records),
            ]:
                if records:
                    integration_side = prefix.split("_")[-1]
                    path = f"{row['id']}_{integration_side}.fastq"
                    SeqIO.write(records.values(), path, "fastq")

            counts = {k: len(v) for k, v in categories.items()}

            complete_att_max = max(counts["complete_attL"], counts["complete_attR"])
            complete_beacon_total = counts["complete_beacon"]
            partial_beacon_total = counts["partial_beacon"]
            wt_total = counts["wt"]

            total_conversion_percentage = 0
            total_complete_PGI_percentage = 0
            total_beacon_percentage = 0
            complete_beacon_percentage = 0
            beacon_fidelity_percentage = 0

            if complete_att_max + complete_beacon_total != 0:
                total_conversion_percentage = (
                    100 * (complete_att_max) / (complete_att_max + complete_beacon_total)
                )
            if (
                complete_att_max + complete_beacon_total + partial_beacon_total + wt_total
                != 0
            ):
                total_complete_PGI_percentage = (
                    100
                    * (complete_att_max)
                    / (
                        complete_att_max
                        + complete_beacon_total
                        + partial_beacon_total
                        + wt_total
                    )
                )
            if (
                complete_att_max + partial_beacon_total + complete_beacon_total + wt_total
                != 0
            ):
                total_beacon_percentage = (
                    100
                    * (partial_beacon_total + complete_beacon_total + complete_att_max)
                    / (
                        complete_att_max
                        + partial_beacon_total
                        + complete_beacon_total
                        + wt_total
                    )
                )
            if (
                complete_att_max + complete_beacon_total + partial_beacon_total + wt_total
                != 0
            ):
                complete_beacon_percentage = (
                    100
                    * (complete_beacon_total + complete_att_max)
                    / (
                        complete_att_max
                        + complete_beacon_total
                        + partial_beacon_total
                        + wt_total
                    )
                )
            if total_beacon_percentage != 0:
                beacon_fidelity_percentage = (
                    100 * (complete_beacon_percentage) / (total_beacon_percentage)
                )

            informative_reads = (
                counts["complete_attL"]
                + counts["complete_attR"]
                + complete_beacon_total
                + partial_beacon_total
                + wt_total
            )
            # indel_percentage = 100 * (indel_total / (attL_total + attR_total + wt_total + complete_beacon_total + partial_beacon_total))

            def calc_percentage(numerator, denominator):
                return (numerator / denominator * 100) if denominator > 0 else 0.0

            row["Role in Cancer"] = row["Role in Cancer"].replace(", ", ";")

            integration_dict[row["id"]] = (
                counts["wt"],
                counts["complete_attL"],
                counts["complete_attR"],
                complete_att_max,
                counts["partial_beacon"],
                counts["complete_beacon"],
                indel_counter,
                informative_reads,
                total_complete_PGI_percentage,
                total_beacon_percentage,
                complete_beacon_percentage,
                beacon_fidelity_percentage,
                indel_percent,
            ) + tuple(
                row[key]
                for key in [
                    "gene_name",
                    "gene_strand",
                    "gene_distance",
                    "same_strand",
                    "overlapping_feature",
                    "threat_tier",
                    "Cosmic Gene",
                    "Role in Cancer",
                    "Liver Expression (TPM)",
                ]
            )
        else:
            integration_dict[row["id"]] = (
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ) + tuple(
                row[key]
                for key in [
                    "gene_name",
                    "gene_strand",
                    "gene_distance",
                    "same_strand",
                    "overlapping_feature",
                    "threat_tier",
                    "Cosmic Gene",
                    "Role in Cancer",
                    "Liver Expression (TPM)",
                ]
            )
    return integration_dict


def write_integration_percentage(integration_dict, sample_name):
    # Open the CSV file in write mode
    filename = f"{sample_name}_integration_stats.csv"
    with open(filename, "w") as file:
        # Write the header row
        file.write(
            "Target,Number of WT Reads,Number of AttL Reads,Number of AttR Reads,Max Recombined Reads,Number of Partial Beacon Reads,Number of Complete Beacon Reads,Number of Indel Reads,Number of Informative Reads,Integration Percentage,Partial Beacon Placement,Complete Beacon Placement,Beacon Fidelity,Indels,Closest Gene Name,Gene Strand,Distance from Gene,Same Strand as Cryptic,Overlapping Feature,Threat Tier,Cosmic Gene,Role in Cancer,Liver Expression (TPM)\n"
        )

        # Write the data rows
        for key, value in integration_dict.items():
            file.write(
                f"{key},{value[0]},{value[1]},{value[2]},{value[3]},{value[4]},{value[5]},{value[6]},{value[7]},{value[8]},{value[9]},{value[10]},{value[11]},{value[12]},{value[13]},{value[14]},{value[15]},{value[16]},{value[17]},{value[18]},{value[19]},{value[20]},{value[21]}\n"
            )


if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(
        description="Extract reads overlapping a certain genomic position."
    )

    # Add the arguments
    parser.add_argument("--target_info", required=True, type=str, help="Metadata file")
    parser.add_argument("--bam", nargs="+", required=True, type=str, help="BAM files")
    parser.add_argument("--sample_name", required=True, type=str, help="Sample name")
    # Parse the arguments
    args = parser.parse_args()

    integration_dict = compute_integration_percentage(
        args.target_info, args.bam, args.sample_name
    )
    write_integration_percentage(integration_dict, args.sample_name)
