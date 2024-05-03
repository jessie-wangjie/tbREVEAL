#!/usr/bin/env python

from utils.base import *
import argparse
from benchling_sdk.services.v2.stable import blob_service
from benchling_sdk.benchling import Benchling
from benchling_sdk.auth.api_key_auth import ApiKeyAuth
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd

def get_project_info(project_id,reads_parent_dir):
    project_query = '''
    SELECT assay,sequencing_run_id,sequencing_project_name
    FROM custom_tracking
    WHERE file_registry_id$ = %s
    '''
    cur.execute(project_query, [project_id])
    project_query_result = cur.fetchone()
    assay, sequencing_run_id, sequencing_project_name = project_query_result

    metadata_query = '''
    SELECT sample_name,attb,attp,primers,species,probes_or_barcodes,donor,"group"
    FROM tomebiosciences.genomic_assays_metadata$raw
    WHERE ctb_id = %s
    '''
    cur.execute(metadata_query, [project_id])
    metadata_query_result = cur.fetchall()
    sample_info = {}
    for result in metadata_query_result:
        sample_name,attb,attp,umi_type,species,probes,donor,group = result



        for entry in os.listdir(reads_parent_dir):
            entry_path = os.path.join(reads_parent_dir, entry)
            if os.path.isdir(entry_path) and sample_name in entry:
                reads_dir = reads_parent_dir + entry

        sample_info[sample_name] = (reads_dir,species,attb,attp,umi_type,probes,donor,group)

    samplesheet_df = pd.DataFrame.from_dict(sample_info, orient="index").reset_index()
    samplesheet_df.columns = ['sample_name','reads_dir','species','attb','attp','umi_type','probes_name','cargo_name','group']

    samplesheet_df.to_csv('samplesheet.csv',index=False)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="get project info")
    # parser.add_argument("--results_dir", help="Parent directory containing the subfolders")
    parser.add_argument("--project_id", help="project id (CTBxxx)")
    parser.add_argument("--reads_parent_dir", help="Parent directory where the reads are")

    args = parser.parse_args()

    get_project_info(args.project_id,args.reads_parent_dir)





