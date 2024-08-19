#!/usr/bin/env python
from benchling_sdk.auth.api_key_auth import ApiKeyAuth
from benchling_sdk.benchling import Benchling
from benchling_sdk.models import DnaSequenceUpdate
from benchling_sdk.helpers.serialization_helpers import fields
from utils.base import *
import pandas as pd
import argparse

def add_sites(project_id, integration_csv_fn):

    integration_df = pd.read_csv(integration_csv_fn,sep=',')
    validated_sites = list(integration_df[integration_df['Integration Percentage'] > 0].Target)
    benchling = Benchling(url="https://tome.benchling.com", auth_method=ApiKeyAuth(api_key))
    for cas_id in validated_sites:
        if 'CAS' not in cas_id:
            continue
        entity = benchling.dna_sequences.list(name=cas_id)
        try:
            cas_info = entity.first()
        except:
            continue

        curr_site_validated_project_ids = cas_info.fields['validation tracking id list'].value
        cas_id_raw = cas_info.id

        if curr_site_validated_project_ids is not None:
            curr_site_validated_project_ids = curr_site_validated_project_ids.split(',')
        else:
            curr_site_validated_project_ids = []

        if project_id not in curr_site_validated_project_ids:
            curr_site_validated_project_ids.append(project_id)
            curr_site_validated_project_ids = list(sorted(curr_site_validated_project_ids))
            update_value_string = ','.join(curr_site_validated_project_ids)
            entityToUpdate = DnaSequenceUpdate(
                fields=fields({
                # Set this field to a new value
                "validation tracking id list": {"value": update_value_string}
                })
            )
            try:
                benchling.dna_sequences.update(dna_sequence_id=cas_id_raw, dna_sequence=entityToUpdate)
            except:
                continue

if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(description="Add sites to Benchling field")

    # Add the arguments
    parser.add_argument("--project_id", required=True, type=str, help="Project id")
    parser.add_argument("--integration_csv", required=True, type=str, help="Integration CSV")

    # Parse the arguments
    args = parser.parse_args()

    add_sites(args.project_id, args.integration_csv)