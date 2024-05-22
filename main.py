from typing import Dict

from Bio import SeqIO



def read_sequence(file_path: str, format: str) -> str:
    for record in SeqIO.parse(file_path + "." + format, format):
        return record.seq

seq = read_sequence("r", "fastq")

print(len(seq))


def get_distribution(seq, k):
    k_mer = {}

    for i in range(len(seq) - k + 1):
        if seq[i:i + k] not in k_mer:
            k_mer[seq[i:i + k]] = 0
        k_mer[seq[i:i + k]] += 1

    k_mer_distribution = {}
    n = sum(k_mer.values())
    for x, y in k_mer.items():
        k_mer_distribution[x] = y / n

    return k_mer_distribution


seq2 = read_sequence("klebsiella_pneumoniae_reference-1", "fasta")

seq2 = seq2.replace("N", "T")

distr1 = get_distribution(seq, 3)
distr2 = get_distribution(seq2, 3)

print(set(distr2.keys()))

print(len(distr2))
import pandas as pd
from scipy.spatial.distance import cosine

line_s = pd.Series(distr1)
query_s = pd.Series(distr2)

print(len(line_s))
print(len(query_s))

print(1 - cosine(line_s, query_s))
