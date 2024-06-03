import os
from collections import defaultdict
from itertools import product
from random import randint
from typing import Dict, List, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import multiprocessing as mp

from scipy import spatial

CPU_COUNT: int = mp.cpu_count()
CONSTANT_K: int = 5
DEFAULT_DICT: Dict[str, int] = {"".join(x): 0 for x in product("ACGT", repeat=CONSTANT_K)}

def get_distribution(seq: Seq) -> Dict[str, float]:
    k_mer: Dict[str, int] = defaultdict(int, {k: 0 for k in DEFAULT_DICT.keys()})

    for i in range(len(seq) - CONSTANT_K + 1):
        k_mer[seq[i:i + CONSTANT_K]] += 1

    total = sum(k_mer.values())

    return {k: v / total for k, v in k_mer.items()}

def preprocess_distribution(distr: Dict[str, float]) -> List[float]:
    return [distr[x] for x in sorted(distr.keys())]

def calc_cosine_similarity(seq1: Seq, seq2: Seq) -> float:
    return 1 - spatial.distance.cosine(
        preprocess_distribution(get_distribution(seq1)),
        preprocess_distribution(get_distribution(seq2))
    )

def create_confusion_matrix(ref: Dict[str, Seq]) -> pd.DataFrame:

    similarity_matrix: defaultdict[str, Dict[str, float]] = defaultdict(lambda: defaultdict(float))

    l: List[Tuple[str, str, Seq, Seq]] = []

    for i, (k1, v1) in enumerate(ref.items()):
        for j, (k2, v2) in enumerate(ref.items()):
            if j >= i:
                continue
            l.append((k1, k2, v1, v2))
            l.append((k2, k1, v1.complement(), v2))
            l.append((k1, k2, v1.reverse_complement(), v2))

    with mp.Pool(CPU_COUNT) as pool:
        results = pool.starmap(calc_cosine_similarity, [(x[2], x[3]) for x in l])

    for i, (k1, k2, v1, v2) in enumerate(l):
        similarity_matrix[k1][k2] = max(similarity_matrix[k1][k2], results[i])
        similarity_matrix[k2][k1] = max(similarity_matrix[k2][k1], results[i])

    df = pd.DataFrame(index=ref.keys(), columns=ref.keys())

    for k1 in similarity_matrix.keys():
        for k2 in similarity_matrix[k1].keys():
            df.at[k1, k2] = similarity_matrix[k1][k2]

    for key in ref.keys():
        df.at[key, key] = 1.0

    return df

def replace_with_random(seq: Seq, rlist: List[str] = ["A", "C", "G", "T"]) -> Seq:
    return Seq("".join([rlist[randint(0, 3)] if x == "N" else x for x in seq]))

def read_references(file_path: str, format: str) -> Dict[str, Seq]:
    references: Dict[str, Seq] = {}
    for file in os.listdir(file_path):
        for record in SeqIO.parse(file_path + file, format):
            references[file] = replace_with_random(record.seq)
            break
    return references

references_dir: str = "references/"
sequences: Dict[str, Seq] = read_references(references_dir, "fasta")

cf_matrix = create_confusion_matrix(sequences)

print(cf_matrix.to_string())

cf_matrix.to_csv("confusion_matrix_k_" + str(CONSTANT_K) + ".csv")