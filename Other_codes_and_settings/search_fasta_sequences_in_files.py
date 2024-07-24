#!/usr/bin/env python3
from Bio import SeqIO
import os
import csv
import argparse
import concurrent.futures
import pandas as pd

all_data = {}

def search_sequence_in_fasta(fasta_file, sequence_to_search):
    """
    Search for a specific sequence in a given FASTA file.
    """
    with open(fasta_file, 'r') as file:
        content = file.read()
        return sequence_to_search in content

def search_sequence_in_file(fasta_file, sequence_to_search):
    """
    Search for a single sequence in a single file.
    """
    if search_sequence_in_fasta(fasta_file, sequence_to_search):
        return 1
    else:
        return 0

def read_fasta_sequences(file_path):
    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequences[str(record.seq)] = record.id
    print(sequences)
    return sequences

def search_sequence_in_files(folder_path, sequence_to_search):
    search_results = []
    target_files = [f for f in os.listdir(folder_path) if f.endswith(".fa") or f.endswith(".fasta")]
    for filename in target_files:
        file_path = os.path.join(folder_path, filename)
        result = search_sequence_in_file(file_path, sequence_to_search)
        search_results.append(result)
    return search_results, target_files

def search_sequences_in_files(folder_path, sequences):
    all_data['filename'] = [f for f in os.listdir(folder_path) if f.endswith(".fa") or f.endswith(".fasta")]
    for sequence, _ in sequences.items():
        search_results, _ = search_sequence_in_files(folder_path, sequence)
        all_data[sequence] = search_results
        print(all_data[sequence])

def main():
    parser = argparse.ArgumentParser(description='Search for sequences in FASTA files.')
    parser.add_argument('input_fasta', type=str, help='Path to the input FASTA file')
    parser.add_argument('folder_path', type=str, help='Path to the folder containing target FASTA files')
    parser.add_argument('output_file', type=str, help='Path to the output CSV file')
    args = parser.parse_args()

    # Read the input FASTA file
    search_sequences = read_fasta_sequences(args.input_fasta)

    # Search for sequences in the target files
    search_sequences_in_files(args.folder_path, search_sequences)

    # Save the results to a CSV file
    pd.DataFrame(all_data).to_csv(args.output_file, index=False, sep='\t')

if __name__ == "__main__":
    main()

