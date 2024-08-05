#!/usr/bin/env python

from utils.base import *
import argparse
from benchling_sdk.services.v2.stable import blob_service
from benchling_sdk.benchling import Benchling
from benchling_sdk.auth.api_key_auth import ApiKeyAuth
import pandas as pd

def get_project_info(project_id):

    metadata_query = '''
    SELECT sample_name,sequencing_sample_name,attb,attp,primers,hybrid_capture_panel,probes_or_barcodes,umi_type,species,probes_or_barcodes,donor,"group"
    FROM tomebiosciences.genomic_assays_metadata$raw
    WHERE genomics_project_queue = %s and archived$ = false
    '''
    cur.execute(metadata_query, [project_id])
    metadata_query_result = cur.fetchall()
    sample_info = {}
    for result in metadata_query_result:
        sample_name, sequencing_sample_name,attb, attp, primers, hybrid_capture_panel, probes_or_barcodes, umi_type, species, probes, donor, group = result

        R1, R2 = None, None  # Initialize R1 and R2 for each sample_name
        for entry in os.listdir('.'):
            if sequencing_sample_name in entry:
                entry_path = os.path.join(os.getcwd(), entry)  # Get the full path to the directory
                if "R1" in entry_path:
                    R1 = entry_path  # Store the absolute path
                elif "R2" in entry_path:
                    R2 = entry_path  # Store the absolute path

        if R1 and R2:  # Ensure both R1 and R2 are found before adding to the dictionary
            if primers and not umi_type:
                if umi_type != 'LMPCR':
                    sample_info[sample_name] = (R1, R2, species, attb, attp, primers, hybrid_capture_panel, donor, group)
                else:
                    sample_info[sample_name] = (R1, R2, species, attb, attp, primers, probes_or_barcodes, donor, group)
            else:
                if umi_type != 'LMPCR':
                    sample_info[sample_name] = (R1, R2, species, attb, attp, umi_type, hybrid_capture_panel, donor, group)
                else:
                    sample_info[sample_name] = (R1, R2, species, attb, attp, umi_type, probes_or_barcodes, donor, group)


    samplesheet_df = pd.DataFrame.from_dict(sample_info, orient="index").reset_index()
    samplesheet_df.columns = ['sample_name', 'read1', 'read2', 'species', 'attb', 'attp', 'umi_type', 'probes_name', 'cargo_name', 'group']

    samplesheet_df.to_csv('samplesheet.csv', index=False)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="get project info")
    # parser.add_argument("--results_dir", help="Parent directory containing the subfolders")
    parser.add_argument("--project_id", help="project id (CTBxxx)")

    args = parser.parse_args()

    get_project_info(args.project_id)





