from typing import Dict

from Bio import SeqIO

seq = ""

for seq_record in SeqIO.parse("r.fastq", "fastq"):
    print(seq_record.id)
    print(seq_record.seq)
    print(type(seq_record.seq))
    print(len(seq_record))
    seq = seq_record.seq
    break


print(len(seq))

def get_distribution(seq, k):
    k_mer = {}

    for i in range(len(seq) - k + 1):
        if seq[i:i+k] not in k_mer:
            k_mer[seq[i:i+k]] = 0
        k_mer[seq[i:i + k]] += 1

    k_mer_distribution = {}
    n = sum(k_mer.values())
    for x, y in k_mer.items():
        k_mer_distribution[x] = y / n

    return k_mer_distribution

seq2 = ""

for seq_record in SeqIO.parse("klebsiella_pneumoniae_reference-1.fasta", "fasta"):
    seq2 = seq_record.seq
    break

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

