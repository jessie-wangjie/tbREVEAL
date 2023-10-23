#!/usr/bin/env python
import seaborn as sns
from matplotlib import pyplot as plt
import pandas as pd
import argparse

def create_plot(excel_fn):
    condensed_results_df = pd.read_excel(excel_fn,sheet_name='Condensed Results')
    reads_cols = [col for col in condensed_results_df if 'Reads' in col or 'reads' in col]
    read_counts_df = condensed_results_df[reads_cols]

    # Get a list of unique sample names
    samples = set(col.split('_')[0] for col in read_counts_df.columns)

    # Store recombination data for all samples
    recombination_data = {}

    for sample in samples:
        attL_col = f'{sample}_Number of AttL Cargo Reads'
        attR_col = f'{sample}_Number of AttR Cargo Reads'
        WT_col = f'{sample}_Number of WT Reads'
        beacon_col = f'{sample}_Number of Complete Beacon Reads'
        
        if attL_col in read_counts_df and attR_col in read_counts_df:
            # Get the max of attL and attR
            max_attL_attR = read_counts_df[[attL_col, attR_col]].max(axis=1)
            
            # Calculate recombination
            recombination = 100 * max_attL_attR / (max_attL_attR + read_counts_df[WT_col] + read_counts_df[beacon_col])
            
            # Store in the dictionary
            recombination_data[sample] = recombination

    # Convert the recombination data to a new dataframe
    recombination_df = pd.DataFrame(recombination_data)

    recombination_df = recombination_df.sort_index(axis=1)

    recombination_df['Target'] = condensed_results_df['Target']

    recombination_df_nonzero = recombination_df[recombination_df.drop('Target',axis=1).sum(axis=1) > 0]

    numeric_cols = recombination_df_nonzero.select_dtypes(include=['number'])
    recombination_df_nonzero['Average'] = numeric_cols.mean(axis=1)

    recombination_df_nonzero = recombination_df_nonzero.sort_values(by='Average', ascending=False)

    recombination_df_nonzero.drop(columns='Average', inplace=True)

    recombination_df_nonzero = recombination_df_nonzero.fillna(0)

    # Draw the heatmap with the mask and correct aspect ratio
    plt.figure(figsize=(30,20))
    ax = sns.heatmap(recombination_df_nonzero.drop('Target',axis=1).transpose(),
                
                xticklabels=recombination_df_nonzero['Target'],
                square=True,
                linewidths=0.2,
                linecolor='black',
                cmap='magma',
                vmin=0,
                cbar_kws={"shrink": 0.5,'location':'top','label':'Recombination %'}
                )

    ax.hlines([len(recombination_df_nonzero.columns)], *ax.get_xlim(), colors=['white'],linewidth=3)
    plt.savefig('recombination_heatmap.png')

if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(description="Extract reads overlapping a certain genomic position.")

    # Add the arguments
    parser.add_argument("--excel_report", required=True, type=str, help="Excel report file")
    
    # Parse the arguments
    args = parser.parse_args()

    create_plot(args.excel_report)