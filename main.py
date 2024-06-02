import os
from random import randint
from typing import List, Dict, Tuple

from Bio import SeqIO
import pandas as pd
from Bio.Seq import Seq
from scipy import spatial
from itertools import product
from collections import defaultdict
from time import time
import multiprocessing as mp


CPU_COUNT = mp.cpu_count()
CONSTANT_K = 3
TOLERANCE = 0.7
SAMPLES_PER_BACTERIA = 500
UNCLASSIFIED = "unclassified"
DEFAULT_DICT = {"".join(x): 0 for x in product("ACGT", repeat=CONSTANT_K)}

def replace_with_random(seq: Seq, rlist: List[str] = ["A", "C", "G", "T"]) -> Seq:
    return Seq("".join([rlist[randint(0, 3)] if x == "N" else x for x in seq]))

def preprocess_distribution(distr: Dict[str, float]) -> List[float]:
    return [distr[x] for x in sorted(distr.keys())]

def calc_cosine_similarity(distr1: List[float], distr2: List[float]) -> float:
    return 1 - spatial.distance.cosine(distr1, distr2)

def get_distribution(seq: Seq, k: int = 3) -> Dict[str, float]:
    k_mer: Dict[str, int] = defaultdict(int, {k: 0 for k in DEFAULT_DICT.keys()})

    for i in range(len(seq) - k + 1):
        k_mer[seq[i:i + k]] += 1

    total: int = sum(k_mer.values())

    return {k: v / total for k, v in k_mer.items()}

def get_preprocessed_distribution(seq: Seq, k: int = 3) -> List[float]:
    return preprocess_distribution(get_distribution(seq, k))

def read_references(file_path: str, format: str) -> Dict[str, Seq]:
    references: Dict[str, Seq] = {}
    for file in os.listdir(file_path):
        for record in SeqIO.parse(file_path + file, format):
            references[file] = replace_with_random(record.seq)
            break
    return references

def read_sequences(file_path: str, format: str, limit: int) -> List[Seq]:
    sequences: List[Seq] = []
    for record in SeqIO.parse(file_path, format):
        if len(record.seq) < CONSTANT_K:
            continue
        sequences.append(replace_with_random(record.seq))
        if len(sequences) >= limit:
            break
    return sequences

def compare_and_assign(ref_d: Dict[str, List[float]], seqs: List[Seq]) -> Dict[str, int]:
    results: Dict[str, int] = defaultdict(int)

    for seq in seqs:

        max_similarity: float = 0
        max_key: str = UNCLASSIFIED

        distributions: List[List[float]] = [
            preprocess_distribution(get_distribution(seq, CONSTANT_K)),
            preprocess_distribution(get_distribution(seq.reverse_complement(), CONSTANT_K)),
            preprocess_distribution(get_distribution(seq.complement(), CONSTANT_K))
        ]

        for key, value in ref_d.items():

            similarity: float = max(
                calc_cosine_similarity(value, distributions[0]),
                calc_cosine_similarity(value, distributions[1]),
                calc_cosine_similarity(value, distributions[2])
            )

            if similarity < TOLERANCE:
                continue

            if similarity > max_similarity:
                max_similarity = similarity
                max_key = key

        results[max_key] += 1

    return results

if __name__ == "__main__":
    references_dir: str = "references/"
    tests_dir: str = "tests/"