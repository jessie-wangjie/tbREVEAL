#!/usr/bin/env python

import os
import pandas as pd
import argparse
import glob
from openpyxl.styles import Border, Side
import numpy as np

def count_fastq_reads(directory):
    read_headers = set()  # Use a set to store read headers because it automatically handles duplicates
    total_reads = 0

    for file_path in directory:
        if file_path.endswith('.fastq'):  # Ensures that we're only processing FASTQ files
            with open(file_path, 'r') as file:
                for line in file:
                    if line.startswith('@'):  # FASTQ read headers typically start with '@'
                        if line not in read_headers:
                            read_headers.add(line)
                            total_reads += 1
                        # Skip the next three lines as they are part of the current read
                        next(file)
                        next(file)
                        next(file)
    return total_reads

def collate_integration_files(project_info, integration_stats_filenames, project_name, collapse_condition):

    project_info_df  = pd.read_csv(project_info)

    # Dictionary initialization
    dfs_integration = {}
    dfs_reads = {}
    dfs_short_integration = {}

    # Merge dataframes for each group
    merged_integration_dfs = []
    merged_reads_dfs = []
    merged_short_integration_dfs = {}

    # Iterate through each sample directory and read the CSV file
    for index,row in project_info_df.iterrows():
        sample_name = row['sample_name']
        group = row['group']

        file_path = [i for i in integration_stats_filenames if str(sample_name) in i][0]

        # Check if the file exists
        if os.path.exists(file_path):
            df = pd.read_csv(file_path)
            df['Fidelity Percentage'] = 100*df[f'Complete Beacon Placement']/df[f'Partial Beacon Placement']
            df['Max Complete Recombination'] = df[['Number of AttL Complete Reads','Number of AttR Complete Reads']].max(axis=1)
            df['Max Cargo Recombination'] = df[['Number of AttL Cargo Reads','Number of AttR Cargo Reads']].max(axis=1)

            # Define the base columns
            base_columns = ['Target', 'Closest Gene Name', 'Gene Strand', 'Distance from Gene', 'Same Strand as Cryptic', 'Overlapping Feature', 'Threat Tier','Cosmic Gene','Role in Cancer','Liver Expression (TPM)']

            base_columns_df = df[base_columns]
            # Filter the columns based on requirements
            integration_cols = [col for col in df.columns if "Integration" in col or "Percentage" in col or "Placement" in col or 'Indels' in col]
            reads_cols = [col for col in df.columns if "Reads" in col or 'reads' in col]
            if collapse_condition == 'Complete':
                integration_short_cols = [col for col in df.columns if 'Complete P Integration Percentage' in col or 'Complete Reads' in col or 'Complete Beacon Reads' or 'Partial Beacon Reads' in col or 'WT Reads' in col or 'Max Complete Recombination' in col]
            elif collapse_condition == 'Partial':
                integration_short_cols = [col for col in df.columns if 'Partial P Integration Percentage' in col or 'Partial Reads' in col or 'Partial Beacon Reads' in col or 'WT Reads' in col or 'Max Partial Recombination' in col]
            elif collapse_condition == 'Cargo':
                integration_short_cols = [col for col in df.columns if 'Cargo and P Integration Percentage' in col or 'Cargo Reads' in col or 'Complete Beacon Reads' in col or 'Partial Beacon Reads' in col or 'WT Reads' in col or 'Max Cargo Recombination' in col]


            # Rename columns to include sample name
            renamed_integration_cols = [f"{sample_name}_{col}" for col in integration_cols]
            renamed_reads_cols = [f"{sample_name}_{col}" for col in reads_cols]
            renamed_integration_short_cols = [f"{sample_name}_{col}" for col in integration_short_cols]

            df_integration = df[integration_cols]

            df_reads = df[reads_cols]

            df_integration_short = df[integration_short_cols]

            df_integration.columns = renamed_integration_cols
            df_reads.columns = renamed_reads_cols
            df_integration_short.columns = renamed_integration_short_cols

            # Update dictionary with group as key
            if group not in dfs_integration:
                dfs_integration[group] = [df_integration]
                dfs_reads[group] = [df_reads]
                dfs_short_integration[group] = [df_integration_short]
            else:
                dfs_integration[group].append(df_integration)
                dfs_reads[group].append(df_reads)
                dfs_short_integration[group].append(df_integration_short)

    for key in dfs_integration:
        temp_df = pd.concat(dfs_integration[key],axis=1)
        temp_df = temp_df.loc[:, ~temp_df.columns.duplicated()]
        merged_integration_dfs.append(temp_df)

    merged_integrations_dfs_combined = pd.concat(merged_integration_dfs,axis=1)
    merged_integrations_dfs_combined = pd.concat([base_columns_df,merged_integrations_dfs_combined],axis=1)

    for key in dfs_reads:
        temp_df = pd.concat(dfs_reads[key],axis=1)
        temp_df = temp_df.loc[:, ~temp_df.columns.duplicated()]
        merged_reads_dfs.append(temp_df)

    merged_reads_dfs_combined = pd.concat(merged_reads_dfs,axis=1)
    merged_reads_dfs_combined = pd.concat([base_columns_df,merged_reads_dfs_combined],axis=1)

    tmp_lst = []
    for key in dfs_short_integration:
        temp_df = pd.concat(dfs_short_integration[key],axis=1)
        temp_df = temp_df.loc[:, ~temp_df.columns.duplicated()]
        if key not in merged_short_integration_dfs:
            merged_short_integration_dfs[key] = [temp_df]
        else:
            merged_short_integration_dfs[key].append(temp_df)
        tmp_lst.append(temp_df)

    merged_short_integration_dfs_combined = pd.concat(tmp_lst,axis=1)
    merged_short_integration_dfs_combined = pd.concat([base_columns_df,merged_short_integration_dfs_combined],axis=1)

    dfs_by_condition = []

    for key in merged_short_integration_dfs:

        temp_df = pd.concat(merged_short_integration_dfs[key],axis=1)
        temp_df = pd.concat([base_columns_df,temp_df],axis=1)
        attL_total = temp_df.filter(like=f'AttL {collapse_condition}').sum(axis=1)
        attR_total = temp_df.filter(like=f'AttR {collapse_condition}').sum(axis=1)
        indel_total = temp_df.filter(like='Number of Indel Reads').sum(axis=1)

        att_max_total = temp_df.filter(like=f'Max {collapse_condition} Recombination').sum(axis=1)

        partial_beacon_total = temp_df.filter(like='Partial Beacon Reads').sum(axis=1)
        complete_beacon_total = temp_df.filter(like='Complete Beacon Reads').sum(axis=1)
        wt_total = temp_df.filter(like='WT Reads').sum(axis=1)
        total_conversion_percentage = 100*(att_max_total) / (att_max_total + complete_beacon_total)
        total_PGI_percentage = 100*(att_max_total) / (att_max_total + complete_beacon_total + partial_beacon_total + wt_total)
        total_beacon_percentage = 100*(partial_beacon_total + complete_beacon_total + att_max_total) / (att_max_total + partial_beacon_total + complete_beacon_total + wt_total)
        complete_beacon_percentage = 100*(complete_beacon_total + att_max_total) / (att_max_total + complete_beacon_total + partial_beacon_total + wt_total)
        beacon_fidelity_percentage = 100*(complete_beacon_percentage)/(total_beacon_percentage)
        indel_percentage = 100 * (indel_total / (attL_total + attR_total + wt_total + complete_beacon_total + partial_beacon_total))
        df_by_condition = pd.DataFrame({f'{key} AttL Total': attL_total,
                                    f'{key} AttR Total': attR_total,
                                    f'{key} Indel Total': indel_total,
                                    f'{key} Recombined Total': att_max_total,
                                    f'{key} Partial Beacon Total': partial_beacon_total,
                                    f'{key} Complete Beacon Total': complete_beacon_total,
                                    f'{key} WT Total': wt_total,
                                    f'{key} Total Beacon %': total_beacon_percentage,
                                    f'{key} Complete Beacon %': complete_beacon_percentage,
                                    f'{key} Beacon Fidelity %': beacon_fidelity_percentage,
                                    f'{key} PGI %': total_PGI_percentage,
                                    f'{key} Conversion %':total_conversion_percentage,
                                    f'{key} Indel %':indel_percentage})
        df_by_condition.reset_index(drop=True,inplace=True)
        dfs_by_condition.append(df_by_condition)

    dfs_by_condition_concat = pd.concat(dfs_by_condition,axis=1)
    dfs_by_condition_concat = pd.concat([base_columns_df,dfs_by_condition_concat],axis=1)
    dfs_by_condition_concat = dfs_by_condition_concat.fillna(0)

    output_fn = f'{project_name}_results.xlsx'
    # Save the two sheets into an Excel file
    with pd.ExcelWriter(output_fn, engine='openpyxl') as writer:
        merged_integrations_dfs_combined.to_excel(writer, sheet_name='Integration Percent', index=False)
        merged_reads_dfs_combined.to_excel(writer, sheet_name='Integration Reads', index=False)
        merged_short_integration_dfs_combined.to_excel(writer, sheet_name='Condensed Results', index=False)
        dfs_by_condition_concat.to_excel(writer, sheet_name='Results by Condition', index=False)

    return(output_fn)


def collate_reads_per_site_files(project_info,read_counts_per_site_filenames, excel_filename):
    project_info_df  = pd.read_csv(project_info)

    unchanged_column = "id"

    # Initialize empty dataframe
    collated_indel_df = pd.DataFrame()

    for index,row in project_info_df.iterrows():
        sample_name = row['sample_name']
        group = row['group']

        # Paths for attL and attR indel tables

        filepath = [i for i in read_counts_per_site_filenames if str(sample_name) in i][0]
        df = pd.read_csv(filepath)
        # Change column names
        for column in df.columns:
            if column != unchanged_column:
                df.rename(columns={column: f"{sample_name} {column}"}, inplace=True)

        # If this is the first file, assign it directly. Else, merge with existing collated_indel_df
        if collated_indel_df.empty:
            collated_indel_df = df
        else:
            collated_indel_df = pd.merge(collated_indel_df, df, how='outer',on=unchanged_column)

    # Append the collated dataframe as a new sheet to the existing Excel file
    with pd.ExcelWriter(excel_filename, engine='openpyxl', mode='a') as writer:
        collated_indel_df.to_excel(writer, sheet_name="Read counts", index=False)


def collate_qc_files(project_info,qc_filenames, extracted_reads_dirs, excel_filename):
    project_info_df  = pd.read_csv(project_info)

    # Initialize empty dataframe
    collated_indel_df = pd.DataFrame()

    for index,row in project_info_df.iterrows():
        sample_name = row['sample_name']
        group = row['group']

        # Load the csv file
        qc_filepath = [i for i in qc_filenames if str(sample_name) in i][0]
        extracted_reads_filepath = [i for i in extracted_reads_dirs if str(sample_name) in i]
        df = pd.read_csv(qc_filepath, header=None)

        # Rename columns
        df.columns = ['', sample_name]

        # Pivot dataframe and combine with collated_indel_df
        df = df.set_index('')

        if collated_indel_df.empty:
            collated_indel_df = df
        else:
            collated_indel_df = pd.concat([collated_indel_df, df], axis=1)

        # Add "reads near probe" row
        reads_near_probe = count_fastq_reads(extracted_reads_filepath)

        collated_indel_df.loc["reads near probe", sample_name] = reads_near_probe

        # Add "reads near probe %" row
        deduped_reads = float(collated_indel_df.loc["deduped reads", sample_name])
        collated_indel_df.loc["reads near probe %", sample_name] = round(reads_near_probe / deduped_reads * 100,2)  # as a percentage

    # Append the collated dataframe as a new sheet to the existing Excel file
    with pd.ExcelWriter(excel_filename, engine='openpyxl', mode='a') as writer:
        collated_indel_df.to_excel(writer, sheet_name=f"QC Summary")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine results from all samples")
    # parser.add_argument("--results_dir", help="Parent directory containing the subfolders")
    parser.add_argument("--project_name", help="Name of the project -- only affects output filenames")
    parser.add_argument("--project_config_file", required=True, help="Project config file")
    parser.add_argument("--integration_stats_files",nargs='+', required=True, help="Integration stats files")
    parser.add_argument("--read_counts_per_site_files",nargs='+', required=True, help="Read counts per site files")
    parser.add_argument("--qc_summary_files",nargs='+', required=True, help="QC Summary files")
    parser.add_argument("--extracted_reads_dirs",nargs='+', required=True, help="Extracted reads dirs")
    parser.add_argument("--collapse_condition_by", help="Collapse by PARTIAL P, COMPLETE P, or CARGO?",default='Complete', choices=['Partial','Complete','Cargo'])

    args = parser.parse_args()

    results_output_fn = f"{args.project_name}_results.xlsx"

    integration_fn = collate_integration_files(args.project_config_file, args.integration_stats_files, args.project_name, args.collapse_condition_by)
    collate_reads_per_site_files(args.project_config_file,args.read_counts_per_site_files, integration_fn)
    collate_qc_files(args.project_config_file,args.qc_summary_files, args.extracted_reads_dirs,integration_fn)
