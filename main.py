from typing import Dict, Tuple
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


def calc_cosine_similarity(distr1: Dict[str, float], distr2: Dict[str, float]) -> float:
    line_s: pd.Series = pd.Series(distr1)
    query_s: pd.Series = pd.Series(distr2)

    return 1 - cosine(line_s, query_s)


def construct_matrix(d: Dict[Tuple[str, str], float]) -> pd.DataFrame:
    df: pd.DataFrame = pd.DataFrame(d.values(), index=pd.MultiIndex.from_tuples(d.keys())).unstack().fillna(1)
    return df


seq1: Seq = read_sequence("r", "fastq")
seq2: Seq = read_sequence("klebsiella_pneumoniae_reference-1", "fasta")

seq2 = seq2.replace("N", "T")

distr1: Dict[str, float] = get_distribution(seq1, 3)
distr2: Dict[str, float] = get_distribution(seq2, 3)

similarity: float = calc_cosine_similarity(distr1, distr2)

d: Dict[Tuple[str, str], float] = {}

d[("seq1", "seq2")] = similarity
d[("seq2", "seq1")] = similarity

df: pd.DataFrame = construct_matrix(d)

print(df)
