#!/usr/bin/env python

import argparse
from utils.base import *
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(description="Get cargo sequences from Benchling")

    # Add the arguments
    parser.add_argument("--cargo", required=False, type=str, help="List of cargo id", default="")
    # Parse the arguments
    args = parser.parse_args()

    cur.execute("SELECT registry_entity.file_registry_id, dna_sequence.bases FROM registry_entity "
                "JOIN dna_sequence on dna_sequence.id = registry_entity.id "
                "WHERE registry_entity.file_registry_id = %s", [args.cargo])
    cargo_name, cargo_bases = cur.fetchone()
    cargo_bases = cargo_bases.upper()

    seq_record = SeqRecord(Seq(cargo_bases), id=cargo_name, description='')

    with open("cargo.fasta", "w") as output_handle:
        SeqIO.write(seq_record, output_handle, 'fasta-2line')