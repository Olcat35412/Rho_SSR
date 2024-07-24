#!/usr/bin/env python3
import re
from Bio import SeqIO
import argparse
def filter_fasta(input_file, output_file):
    with open(input_file, "r") as f:
        records = SeqIO.parse(f, "fasta")
        valid_records = []
        illegal_chars_pattern = re.compile(r"[^A-Za-z0-9\s]")
        for record in records:
            if not illegal_chars_pattern.search(str(record.seq)):
                valid_records.append(record)
    with open(output_file, "w") as f:
        SeqIO.write(valid_records, f, "fasta")
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter FASTA file by removing sequences with illegal characters")
    parser.add_argument('-i',"--input_file", help="Input FASTA file")
    parser.add_argument('-o',"--output_file", help="Output FASTA file")
    args = parser.parse_args()
    filter_fasta(args.input_file, args.output_file)
