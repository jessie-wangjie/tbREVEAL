#!/usr/bin/env python

import argparse
import pandas as pd
from utils.base import *
import subprocess
import sys


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement[base] for base in reversed(seq))


def get_cas_info(dinucleotides, species):
    query = '''
            SELECT
                chr,
                central_nucleotides_start - 1,
                central_nucleotides_start,
                file_registry_id$,
                strand,
                threat_tier
            FROM
                integrase_cryptic_attachment_site
            where species = %s
            '''

    if dinucleotides != "":
        di = set()
        for i in dinucleotides.split(","):
            di.add("central_nucleotides_seq = '" + i + "'")
            di.add("central_nucleotides_seq = '" + reverse_complement(i) + "'")

        query += '''
            and 
                ''' + " or ".join(di)
        cur.execute(query, [species.lower(), dinucleotides])
    else:
        cur.execute(query, [species.lower()])

    for row in cur.fetchall():
        print("\t".join(str(i) for i in row))


if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(description="Get CAS cut inforamtion from Benchling")

    # Add the arguments
    parser.add_argument("--dinucleotides", required=False, type=str, help="Metadata file", default="")
    parser.add_argument("--species", required=False, type=str, help="species", default="human")
    # Parse the arguments
    args = parser.parse_args()
    get_cas_info(args.dinucleotides, args.species)

