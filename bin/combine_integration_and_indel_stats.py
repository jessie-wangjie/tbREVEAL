#!/usr/bin/env python

import argparse
import pandas as pd

def merge_tables(attL_indel_table, attR_indel_table, integration_stats, output_fn):
    # read the csv files into dataframes
    df1 = pd.read_csv(attL_indel_table)
    df2 = pd.read_csv(attR_indel_table)
    df3 = pd.read_csv(integration_stats)

    df1['Indels'] = df1['Insertions'] + df1['Deletions']
    df1['Indel%'] = df1['Indels'] / df1['Input Reads']

    df2['Indels'] = df2['Insertions'] + df2['Deletions']
    df2['Indel%'] = df2['Indels'] / df2['Input Reads']

    # concatenate df1 and df2 vertically
    df = pd.concat([df1, df2])

    df = df.sort_values('Input Reads', ascending=False)
    df = df.drop_duplicates(subset='Amplicon', keep='first')

    df = df.reset_index(drop=True)

    print(df)
    print(df3)

    # merge the new dataframe with df3
    merged_df = pd.merge(df, df3, how='right', left_on='Amplicon', right_on='Target')

    # drop the 'Target' column
    merged_df.drop('Target', axis=1, inplace=True)
    
    # rename the 'Amplicon' column to 'Target'
    merged_df.rename(columns={'Amplicon': 'Target'}, inplace=True)

    merged_df['Target'] = df3['Target']
    
    merged_df.fillna({'Target':0, 'Input Reads':0, 'Unmodified Reads':0, 'Modified Reads':0, 'Unmodified%':0, 
                      'Modified%':0, 'Insertions':0, 'Deletions':0, 'Substitutions':0, 'Insertion%':0, 
                      'Deletion%':0, 'Substitution%':0, 'Indels':0,'Indel%': 0,'Number of WT Reads':0, 'Number of AttL reads':0,
                      'Number of AttR Reads':0,'Number of Beacon reads':0,'Beacon Integration':0,'Integration Percentage':0}, inplace=True)

    merged_df.to_csv(output_fn,index=False)
    
    # return the merged dataframe
    return merged_df    

if __name__ == '__main__':
    # initialize ArgumentParser
    parser = argparse.ArgumentParser(description='Process some CSV files.')

    # add arguments
    parser.add_argument('--attL_indel_table', help='The name of the attL indel table CSV file')
    parser.add_argument('--attR_indel_table', help='The name of the attR indel table CSV file')
    parser.add_argument('--integration_table', help='The name of the integration stats CSV file')
    parser.add_argument('--output_fn', help='The name of the integration stats CSV file')
    

    # parse arguments
    args = parser.parse_args()

    # call merge_tables function with parsed arguments
    merged_df = merge_tables(args.attL_indel_table, args.attR_indel_table, args.integration_table, args.output_fn)
    
    # print the merged dataframe
    print(merged_df)
