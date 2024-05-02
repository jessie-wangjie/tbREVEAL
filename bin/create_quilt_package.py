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
    parser.add_argument('--project_name', type=str, help='folder where all results are contained (usually after running collate script)')

    args = parser.parse_args()

    project_name = args.project_name
    output_folder = args.output_folder
    package_name = args.package_name
    bucket_name = args.bucket_name

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
        p.push(package_name, bucket_name, force=True)
    base_url = output[1].split()[-1]
    full_url = f"{base_url}/tree/{p.top_hash}"
    print(full_url)
