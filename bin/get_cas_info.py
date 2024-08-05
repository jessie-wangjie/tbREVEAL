#!/usr/bin/env python

import argparse
from utils.base import *

if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(description="Get CAS cut inforamtion from Benchling")

    # Add the arguments
    parser.add_argument("--cas", required=False, type=str, help="List of CAS id", default="")
    # Parse the arguments
    args = parser.parse_args()

    with open(args.cas) as fn:
        cas_list = [f"file_registry_id$ = '{row.rstrip()}'" for row in fn]

    query = '''
        SELECT chr, central_nucleotides_start - 1, central_nucleotides_start, file_registry_id$ 
        FROM integrase_cryptic_attachment_site
        WHERE
        '''
    query += " or ".join(cas_list)
    cur.execute(query)

    for row in cur.fetchall():
        print("\t".join(str(i) for i in row))
