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

    args = parser.parse_args()

    run_id = args.run_id
 

    for f in glob.glob(os.path.join("*L1_ds*", "integration_and_indel_stats.csv")):
        # preview.append(os.path.basename(f))
        sample_name = f.split('/')[0].split("_L1_ds")[0]
        quilt_path = 'integration_stats' + '/' + sample_name + '_integration_stats.csv'
        p.set(quilt_path, f)
    
    for f in glob.glob("*L1_ds*/att*_alignments_html/*.html"):
        sample_name = f.split('/')[0].split("_L1_ds")[0]
        junction_type = f.split('/')[-1].split('_')[1].split('.')[0]
        site = f.split('/')[-1].split('_')[1].split('.')[1]
        quilt_path = 'attL_attR_alignment_files/'  + sample_name + '/' + junction_type  + '_alignments/' + site + '_' + junction_type + '.html'
        p.set(quilt_path, f)
    
    for f in glob.glob("*L1_ds*/qc/*.csv"):
        sample_name = f.split('/')[0].split("_L1_ds")[0]
        qc_filename = f.split('/')[-1].split('.')[0]
        quilt_path = 'qc' + '/' + sample_name + '/' + qc_filename + '.csv'
        p.set(quilt_path, f)

    for f in glob.glob("*L1_ds*/cs2_output/cs2_att*/CRISPResso_on*/cs2_alignment_html/*.html"):
        print(f)
        sample_name = f.split('/')[0].split("_L1_ds")[0]
        junction_type = f.split('/')[-1].split('_')[1].split('.')[0]
        site = f.split('/')[-1].split('_')[1].split('.')[1]
        quilt_path = 'attL_attR_alignment_files/'  + sample_name + '/' + junction_type  + '_alignments/' + site + '_' + junction_type + '.html'
        p.set(quilt_path, f)
        
    # p.set("integration_stats/S4-down-350-rep1_integration_and_indel_stats.csv","S4-down-350-rep1_L1_ds.d22b471c136648e3abae0dd35a0ccb19/integration_and_indel_stats.csv")
    # p.set_dir("fastq/" + pipeline_run_id[:-1], pipeline_run_id[:-1])

    # # output package
    # preview = ["platemap.json", "alignment_stats.json"]
    # for f in glob.glob(os.path.join(input, "stats.*.csv")):
    #     preview.append(os.path.basename(f))
    #     p.set(os.path.join(pipeline_run_id, os.path.basename(f)), f)
    
    

    # p.set(pipeline_run_id + "/qw_stats.csv", input + "/" + "qw_stats.csv")
    # p.set(pipeline_run_id + "/" + pipeline_run_id + ".stats.xlsx", input + "/" + input + ".stats.xlsx")
    # p.set(pipeline_run_id + "/platemap.json", input + "/platemap.json")
    # p.set(pipeline_run_id + "/alignment_stats.json", input + "/alignment_stats.json")
    # p.set(pipeline_run_id + "/status.txt", input + "/" + "status.txt")
    # p.set_dir(pipeline_run_id + "/cs2_alignment_html", input + "/cs2_alignment_html/")
    # p.set_meta({"Benchling Entry": entry_name, "Benchling URL": entry_url})
    # pd.Series(preview).to_json(input + "/quilt_summarize.json", orient="records")
    # p.set(pipeline_run_id + "/quilt_summarize.json", input + "/quilt_summarize.json")

    # Pushing a package to a remote registry
    with Capturing() as output:
        p.push("HybridCapture/" + run_id, "s3://tb-ngs-genomics-quilt/", force=True)
    base_url = output[1].split()[-1]
    full_url = f"{base_url}/tree/{p.top_hash}"
    print(full_url)
