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

