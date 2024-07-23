#!/usr/bin/env python

import argparse
import glob
import sys
from io import StringIO
import pandas as pd
import quilt3
import re
import os
import json
import tqdm
from benchling_sdk.auth.api_key_auth import ApiKeyAuth
from benchling_sdk.benchling import Benchling
from benchling_sdk.models import CustomEntityUpdate
from benchling_sdk.helpers.serialization_helpers import fields
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

    parser.add_argument('--output_folder', type=str, help='output_folder')
    parser.add_argument('--bucket_name', type=str, help='bucket name')
    parser.add_argument('--package_name', type=str, help='package name')
    parser.add_argument('--project_id', type=str, help='folder where all results are contained (usually after running collate script)')

    args = parser.parse_args()

    project_id = args.project_id
    project_id_base = project_id.split('_')[0]
    output_folder = args.output_folder
    package_name = args.package_name
    bucket_name = args.bucket_name

    package_name = package_name.replace(' ','_')

    benchling = Benchling(url="https://tome.benchling.com", auth_method=ApiKeyAuth(api_key))
    entity = benchling.custom_entities.list(name=project_id_base)
    print(entity.first())

    def print_absolute_paths(root_dir):
        paths = []
        for root, dirs, files in os.walk(root_dir):
            for file in files:
                paths.append(os.path.join(root, file))
        return(paths)

    paths = print_absolute_paths(output_folder)
    for path in paths:
        quilt_path = path.split(output_folder+'/')[-1]
        try:
            p.set(quilt_path, path)
        except:
            continue

    # Pushing a package to a remote registry
    with Capturing() as output:
        print(output)
        p.push(package_name, bucket_name, force=True)
    base_url = output[0]
    hash_string = f"/tree/{p.top_hash}"
    url_prefix = 'https://quilt.tome.bio/b/'
    package_name_top_name = package_name.split('/')[0]
    package_name_bottom_name = package_name.split('/')[1]

    print(package_name)

    full_url = url_prefix + "tb-ngs-genomics-quilt" + '/packages/' + package_name  + hash_string
    print(full_url)

    entityToUpdate = CustomEntityUpdate(
        fields=fields({
        # Set this field to a new value
        "Results link": {"value": full_url},
        "Status": {"value": "sfso_sDL7Wjdy"}
        })
    )
    benchling.custom_entities.update(entity_id=entity.first().id, entity=entityToUpdate)
