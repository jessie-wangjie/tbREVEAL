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

def create_amplicon_files(target_info):
    # Specify the directory path
    amplicons_dir = 'amplicons'
    os.makedirs(amplicons_dir)
    target_info_df = pd.read_csv(target_info)
    for index, row in target_info_df.iterrows():
        output_fasta_fn = amplicons_dir + '/' + row['id'] + '_amplicon.fasta'
        attL_amplicon = row['attL']
        attR_amplicon = row['attR']
        beacon_amplicon = row['beacon']
        wt_amplicon = row['wt']
        print(row['id'], beacon_amplicon)
        if 'CAS' in row['id']:
            records = [
                SeqRecord(Seq(attL_amplicon), id="attL_amplicon", description=""),
                SeqRecord(Seq(attR_amplicon), id="attR_amplicon", description=""),
                SeqRecord(Seq(beacon_amplicon), id="beacon_amplicon", description=""),
            ]
        elif 'AA' in row['id']:
            records = [
                SeqRecord(Seq(attL_amplicon), id="attL_amplicon", description=""),
                SeqRecord(Seq(attR_amplicon), id="attR_amplicon", description=""),
                SeqRecord(Seq(beacon_amplicon), id="beacon_amplicon", description=""),
                SeqRecord(Seq(wt_amplicon), id="wt_amplicon", description="")
            ]

        SeqIO.write(records, output_fasta_fn, "fasta")

if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(description="Extract reads overlapping a certain genomic position.")

    # Add the arguments
    parser.add_argument("--target_info", required=True, type=str, help="Metadata file")
    # Parse the arguments
    args = parser.parse_args()

    create_amplicon_files(args.target_info)