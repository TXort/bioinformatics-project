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



