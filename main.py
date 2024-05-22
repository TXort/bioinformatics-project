from typing import Dict

import Bio
from Bio import SeqIO
import pandas as pd
from Bio.Seq import Seq
from scipy.spatial.distance import cosine


def read_sequence(file_path: str, format: str) -> Seq:
    for record in SeqIO.parse(file_path + "." + format, format):
        return record.seq


def get_distribution(seq: str, k: int) -> Dict[str, float]:
    k_mer: Dict[str, int] = {}

    for i in range(len(seq) - k + 1):
        if seq[i:i + k] not in k_mer:
            k_mer[seq[i:i + k]] = 0
        k_mer[seq[i:i + k]] += 1

    k_mer_distribution: Dict[str, float] = {}

    div: int = sum(k_mer.values())

    for x, y in k_mer.items():
        k_mer_distribution[x] = y / div

    return k_mer_distribution


seq1: Seq = read_sequence("r", "fastq")
seq2: Seq = read_sequence("klebsiella_pneumoniae_reference-1", "fasta")


seq2 = seq2.replace("N", "T")

distr1: Dict[str, float] = get_distribution(seq1, 3)
distr2: Dict[str, float] = get_distribution(seq2, 3)

print(set(distr2.keys()))

print(len(distr2))

line_s = pd.Series(distr1)
query_s = pd.Series(distr2)

print(len(line_s))
print(len(query_s))

print(1 - cosine(line_s, query_s))
