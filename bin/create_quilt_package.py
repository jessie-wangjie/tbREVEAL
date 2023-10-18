import argparse
import glob
import sys
from io import StringIO
import altair as alt
import pandas as pd
import quilt3
import re
import os
import json
from utils.base import *

class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self

    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio  # free up some memory
        sys.stdout = self._stdout


if __name__ == "__main__":


    p = quilt3.Package()

    parser = argparse.ArgumentParser(description='Process run_id argument')

    parser.add_argument('--run_id', type=str, help='run id')
    parser.add_argument('--parent_dir', type=str, help='folder where all results are contained (usually after running collate script)')

    args = parser.parse_args()

    run_id = args.run_id
    parent_folder = args.parent_dir
 

    # for f in glob.glob(os.path.join(f"{parent_folder}/*L1_ds*", "integration_and_indel_stats.csv")):
    #     sample_name = f.split('/')[1].split("L1_ds")[0].strip('_')
    #     quilt_path = 'integration_stats' + '/' + sample_name + '_integration_stats.csv'
    #     p.set(quilt_path, f)
    
    # for f in glob.glob(f"{parent_folder}/*L1_ds*/att*_alignments_html/*.html"):
    #     sample_name = f.split('/')[1].split("L1_ds")[0].strip('_')
    #     junction_type = f.split('/')[-1].split('_')[1].split('.')[0]
    #     site = f.split('/')[-1].split('_')[1].split('.')[1]
    #     quilt_path = 'attL_attR_alignment_files/'  + sample_name + '/' + junction_type  + '_alignments/' + site + '_' + junction_type + '.html'
    #     p.set(quilt_path, f)
    
    # for f in glob.glob(f"{parent_folder}/*L1_ds*/qc/*.csv"):
    #     sample_name = f.split('/')[1].split("L1_ds")[0].strip('_')
    #     qc_filename = f.split('/')[-1].split('.')[0]
    #     quilt_path = 'qc' + '/' + sample_name + '/' + qc_filename + '.csv'
    #     p.set(quilt_path, f)

    for f in glob.glob(f"{parent_folder}/*/cs2_output/cs2_*/CRISPResso_on*/cs2_alignment_html/*.html"):
        sample_name = f.split('/')[1].split("L1_ds")[0].strip('_')
        junction_type = f.split('/')[-1].split('_')[1].split('.')[0]
        site = f.split('/')[-1].split('_')[1].split('.')[1]
        quilt_path = 'attL_attR_beacon_alignment_files/'  + sample_name + '/' + junction_type  + '_alignments/' + site + '_' + junction_type + '.html'
        p.set(quilt_path, f)
    
    # for f in glob.glob(f"{parent_folder}/{run_id}*.csv"):
    #     quilt_path = f.split('/')[1]
    #     p.set(quilt_path, f)
        
    for f in glob.glob(f"{parent_folder}/*.xlsx"):
        quilt_path = f.split('/')[1]
        p.set(quilt_path, f)
    
    for f in glob.glob(f"{parent_folder}/*.bam*"):
        quilt_path = f.split('/')[1]
        p.set(quilt_path, f)
    
    for f in glob.glob(f"{parent_folder}/*.fa"):
        quilt_path = f.split('/')[1]
        p.set(quilt_path, f)

    # Pushing a package to a remote registry
    with Capturing() as output:
        p.push("HybridCapture/" + run_id, "s3://tb-ngs-genomics-quilt/", force=True)
    base_url = output[1].split()[-1]
    full_url = f"{base_url}/tree/{p.top_hash}"
    print(full_url)
