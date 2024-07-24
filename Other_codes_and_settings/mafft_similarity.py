#!/usr/bin/env python3
from Bio import AlignIO
import numpy as np
import pandas as pd

def calculate_similarity_matrix(alignment):
    num_sequences = len(alignment)
    similarity_matrix = np.zeros((num_sequences, num_sequences))

    for i in range(num_sequences):
        for j in range(num_sequences):
            if i <= j:  # Matrix is symmetrical, so only compute half
                similarity = calculate_sequence_similarity(alignment[i].seq, alignment[j].seq)
                similarity_matrix[i][j] = similarity
                similarity_matrix[j][i] = similarity

    return similarity_matrix

def calculate_sequence_similarity(seq1, seq2):
    match_count = sum(1 for a, b in zip(seq1, seq2) if a == b)
    return match_count / len(seq1)

#put_your_mafft_combined_file_here
alignment_file = "/SSR_gene_data/tongyong_data/Phylogenetic_tree/OrthoFinder/Results_Jun03_1/all_mafft.notermal"
alignment = AlignIO.read(alignment_file, "fasta")
similarity_matrix = calculate_similarity_matrix(alignment)
sequence_names = [record.id for record in alignment]


df = pd.DataFrame(similarity_matrix, columns=sequence_names, index=sequence_names)

print("Similarity Matrix:")
print(df)
#put your output_file here
output_file = ""
df.to_csv(output_file, sep='\t')

print(f"Similarity matrix saved to {output_file}")
