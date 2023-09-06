import os
import pandas as pd
import argparse
import glob 

def count_fastq_reads(directory):
    file_paths = glob.glob(os.path.join(directory, "*.fastq"))
    total_lines = 0
    for file_path in file_paths:
        with open(file_path, 'r') as file:
            total_lines += sum(1 for line in file)
    return total_lines // 4  # divide by 4 as each read in a fastq file is represented by 4 lines

def collate_integration_files(parent_dir, project_name):
    subfolders = [f.path for f in os.scandir(parent_dir) if f.is_dir()]
    
    # Specify columns that don't change from folder to folder
    unchanged_columns = ["Target", "Closest Gene Name", "Gene Strand", 
                         "Distance from Gene", "Overlapping Feature", "Same Strand as Cryptic", "Threat Tier"]

    # Initialize empty dataframes
    collated_indel_df = pd.DataFrame()
    collated_integration_percentage_df = pd.DataFrame()

    for subfolder in subfolders:
        # Obtain sample name from subfolder name
        sample_name = os.path.basename(subfolder).split("_L1_ds")[0]

        # Load the csv file
        filepath = os.path.join(subfolder, "integration_and_indel_stats.csv")
        df = pd.read_csv(filepath)

        df.columns = [column.strip() for column in df.columns]

        # Change column names and create additional dataframe for "Integration Percentage" columns
        for column in df.columns:
            if column not in unchanged_columns:
                df.rename(columns={column: f"{sample_name} {column}"}, inplace=True)
                if "Integration Percentage" in column or "Number of AttL" in column or "Number of AttR" in column or "Number of Beacon" in column:
                    if collated_integration_percentage_df.empty:
                        collated_integration_percentage_df = df[unchanged_columns + [f"{sample_name} {column}"]]
                    else:
                        collated_integration_percentage_df = pd.merge(collated_integration_percentage_df, df[unchanged_columns + [f"{sample_name} {column}"]], on=unchanged_columns)

        # If this is the first file, assign it directly. Else, merge with existing collated_indel_df
        if collated_indel_df.empty:
            collated_indel_df = df
        else:
            collated_indel_df = pd.merge(collated_indel_df, df, on=unchanged_columns)

    # Remove the integration columns from the collated_indel_df
    for column in collated_indel_df.columns:
        if "Integration Percentage" in column:
            collated_indel_df = collated_indel_df.drop(columns=[column])

    # Reorder columns to put 'Target' at the beginning and the rest of unchanged_columns at the end
    indel_column_order = ["Target"] + [col for col in collated_indel_df.columns if col not in unchanged_columns and col != "Target"] + [col for col in unchanged_columns if col != "Target"]
    collated_indel_df = collated_indel_df[indel_column_order]

    integration_column_order = ["Target"] + [col for col in collated_integration_percentage_df.columns if col not in unchanged_columns and col != "Target" and "Integration" not in col] + [col for col in collated_integration_percentage_df.columns if col not in unchanged_columns and col != "Target" and "Integration" in col] + [col for col in unchanged_columns if col != "Target"]
    collated_integration_percentage_df = collated_integration_percentage_df[integration_column_order]


    # Save the collated dataframes to new csv files
    collated_indel_df.to_csv(f"{project_name}_indel_results.csv", index=False)
    collated_integration_percentage_df.to_csv(f"{project_name}_integration_results.csv", index=False)

    return(f"{project_name}_indel_results.csv", f"{project_name}_integration_results.csv")


def collate_reads_per_site_files(parent_dir, project_name):
    subfolders = [f.path for f in os.scandir(parent_dir) if f.is_dir()]

    unchanged_column = "id"

    # Initialize empty dataframe
    collated_indel_df = pd.DataFrame()

    for subfolder in subfolders:
        # Obtain sample name from subfolder name
        sample_name = os.path.basename(subfolder).split("_L1_ds")[0]

        # Load the csv file
        filepath = os.path.join(subfolder, "read_counts_per_site.csv")
        df = pd.read_csv(filepath)

        # Change column names
        for column in df.columns:
            if column != unchanged_column:
                df.rename(columns={column: f"{sample_name} {column}"}, inplace=True)

        # If this is the first file, assign it directly. Else, merge with existing collated_indel_df
        if collated_indel_df.empty:
            collated_indel_df = df
        else:
            collated_indel_df = pd.merge(collated_indel_df, df, on=unchanged_column)

    # Save the collated dataframe to a new csv file
    collated_indel_df.to_csv(f"{project_name}_read_counts_per_site.csv", index=False)

    return(f"{project_name}_read_counts_per_site.csv")

def collate_qc_files(parent_dir, project_name):
    subfolders = [f.path for f in os.scandir(parent_dir) if f.is_dir()]

    # Initialize empty dataframe
    collated_indel_df = pd.DataFrame()

    for subfolder in subfolders:
        # Obtain sample name from subfolder name
        sample_name = os.path.basename(subfolder).split("_L1_ds")[0]

        # Load the csv file
        filepath = os.path.join(subfolder, "qc", "qc_summary.csv")
        df = pd.read_csv(filepath, header=None)

        # Rename columns
        df.columns = ['', sample_name]

        # Pivot dataframe and combine with collated_indel_df
        df = df.set_index('')

        if collated_indel_df.empty:
            collated_indel_df = df
        else:
            collated_indel_df = pd.concat([collated_indel_df, df], axis=1)

        # Add "reads near probe" row
        reads_near_probe = count_fastq_reads(os.path.join(subfolder, "extracted_reads"))
        collated_indel_df.loc["reads near probe", sample_name] = reads_near_probe

        # Add "reads near probe %" row
        after_filter_reads = float(collated_indel_df.loc["after filter reads", sample_name])
        collated_indel_df.loc["reads near probe %", sample_name] = round(reads_near_probe / after_filter_reads * 100,2)  # as a percentage

    # Save the collated dataframe to a new csv file
    collated_indel_df.to_csv(f"{project_name}_qc_summary.csv")

    return(f"{project_name}_qc_summary.csv")

def combine_csv_files_into_excel(integration_stats_filename, indel_stats_filename, read_counts_filename, qc_summary_filename, output_filename):
    # Load data from CSV files into pandas dataframes
    df_integration = pd.read_csv(integration_stats_filename)
    df_indels = pd.read_csv(indel_stats_filename)
    df_reads = pd.read_csv(read_counts_filename)
    df_qc = pd.read_csv(qc_summary_filename)

    # Create a Pandas Excel writer using XlsxWriter as the engine
    writer = pd.ExcelWriter(output_filename, engine='xlsxwriter')

    # Write each dataframe to a different worksheet
    df_integration.to_excel(writer, sheet_name='Integration Stats', index=False)
    df_indels.to_excel(writer, sheet_name='Indel Stats', index=False)
    df_reads.to_excel(writer, sheet_name='Read Counts per Probe', index=False)
    df_qc.to_excel(writer, sheet_name='QC Summary', index=False)

    # Close the Pandas Excel writer and output the Excel file
    writer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine results from all samples")
    parser.add_argument("--results_dir", help="Parent directory containing the subfolders")
    parser.add_argument("--project_name", help="Name of the project -- only affects output filenames")

    args = parser.parse_args()

    results_output_fn = f"{args.project_name}_results.xlsx"

    indel_fn, integration_fn = collate_integration_files(args.results_dir, args.project_name)
    reads_per_site_fn = collate_reads_per_site_files(args.results_dir,args.project_name)
    qc_fn = collate_qc_files(args.results_dir, args.project_name)
    combine_csv_files_into_excel(integration_fn,indel_fn,reads_per_site_fn,qc_fn,results_output_fn)
