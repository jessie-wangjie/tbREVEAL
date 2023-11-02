#!/usr/bin/env python

import os
import glob
import pandas as pd
import argparse

def generate_indel_table(directory, output_fn):
    # Define the wildcard pattern
    file_pattern = "CRISPResso_quantification_of_editing_frequency.txt"

    # Initialize the final table as an empty DataFrame
    final_table = pd.DataFrame(columns=["Amplicon","Input Reads","Unmodified Reads", "Modified Reads","Unmodified%", "Modified%", "Insertions", "Deletions",
                                        "Substitutions", "Insertion%", "Deletion%", "Substitution%"])


    # Find all matching files in the directory
    file_paths = glob.glob(os.path.join(directory, file_pattern))

    # Loop through each file
    for file_path in file_paths:
        # Read the file and extract the required data
        with open(file_path, "r") as file:
            lines = file.readlines()

            # Extract the data from the second row (assuming it's always the second row)
            data = lines[1].strip().split("\t")
            amplicon_name = str(data[0])
            unmodified_percent = float(data[1])
            modified_percent = float(data[2])
            input_reads = int(data[3])
            unmodified_reads = int(data[6])
            modified_reads = int(data[7])
            insertions = int(data[9])
            deletions = int(data[10])
            substitutions = int(data[11])
            insertion_percent = float(data[12]) / input_reads
            deletion_percent = float(data[13]) / input_reads
            substitution_percent = float(data[14]) / input_reads

            # New line of code
            final_table = pd.concat([final_table, pd.DataFrame([{
                "Amplicon": amplicon_name,
                "Input Reads": input_reads,
                "Unmodified Reads": unmodified_reads,
                "Modified Reads": modified_reads,
                "Unmodified%": unmodified_percent,
                "Modified%": modified_percent,
                "Insertions": insertions,
                "Deletions": deletions,
                "Substitutions": substitutions,
                "Insertion%": insertion_percent,
                "Deletion%": deletion_percent,
                "Substitution%": substitution_percent
            }])], ignore_index=True)

    # Save the final table to a CSV file
    final_table.to_csv(output_fn, index=False)

    # Print the final table
    print(final_table)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a final table from files in a directory.")
    parser.add_argument("--cs2_directory", type=str, help="Path to the directory containing the files.")
    parser.add_argument("--output_fn", type=str, help="Path to the directory containing the files.")
    args = parser.parse_args()

    generate_indel_table(args.cs2_directory, args.output_fn)
