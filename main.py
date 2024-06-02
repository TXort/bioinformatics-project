import argparse
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


def parse_arguments() -> argparse.Namespace:
    parser: argparse.ArgumentParser = argparse.ArgumentParser(
        description='Bacterial classification using cosine similarity')

    parser.add_argument('--references_dir', default='references/', type=str, help='Path to the references directory')
    parser.add_argument('--tests_dir', default='tests/', type=str, help='Path to the tests directory')
    parser.add_argument('--constant_k', default=5, type=int, help='Constant K value')
    parser.add_argument('--tolerance', default=0.7, type=float, help='Tolerance value')
    parser.add_argument('--samples_per_bacteria', default=500, type=int, help='Number of samples per bacteria')
    parser.add_argument('--save_to_file', default=False, type=bool, help='Save results to file')
    parser.add_argument('--cpu_count', default=mp.cpu_count(), type=int, help='Number of CPUs to use')

    if not os.path.exists(parser.parse_args().references_dir):
        raise FileNotFoundError("References directory not found")
    if not os.path.exists(parser.parse_args().tests_dir):
        raise FileNotFoundError("Tests directory not found")

    if not 3 <= parser.parse_args().constant_k <= 10:
        raise ValueError("Constant K must be between 3 and 10")

    if not 0 <= parser.parse_args().tolerance <= 1:
        raise ValueError("Tolerance must be between 0 and 1")

    if not 1 <= parser.parse_args().cpu_count <= mp.cpu_count():
        raise ValueError("CPU count must be between 1 and the number of CPUs available")

    return parser.parse_args()


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


def create_dataframe(ref: Dict[str, Seq], samples: Dict[str, List[Seq]]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    start: float = time()

    with mp.Pool(CPU_COUNT) as pool:
        items: List[List[float]] = pool.starmap(get_preprocessed_distribution,
                                                [(value, CONSTANT_K) for key, value in ref.items()])

    red_d_prep: Dict[str, List[float]] = dict(zip(ref.keys(), items))

    end: float = time()
    print("Time for reference distribution processing: ", end - start)

    df: pd.DataFrame = pd.DataFrame(index=ref.keys(), columns=samples.keys())

    df.loc[UNCLASSIFIED] = 0

    result_dict: Dict[Tuple[str, str], int] = defaultdict(int)

    time_start: float = time()

    with mp.Pool(CPU_COUNT) as pool:
        global_list: List[Dict[str, int]] = pool.starmap(compare_and_assign,
                                                         [(red_d_prep, seqs) for seqs in samples.values()])

    for data, sample_name in zip(global_list, samples.keys()):
        for key, value in data.items():
            result_dict[(key, sample_name)] += value

    time_end: float = time()

    print("Time for comparison: ", time_end - time_start)

    for (key, sample_name), value in result_dict.items():
        df.loc[key, sample_name] = value

    df.fillna(0, inplace=True)

    return df


if __name__ == "__main__":
    args = parse_arguments()

    CONSTANT_K = args.constant_k
    TOLERANCE = args.tolerance
    SAMPLES_PER_BACTERIA = args.samples_per_bacteria
    CPU_COUNT = args.cpu_count
    DEFAULT_DICT = {"".join(x): 0 for x in product("ACGT", repeat=CONSTANT_K)}

    references_dir: str = "references/"
    tests_dir: str = "tests/"

    start: float = time()
    sequences: Dict[str, Seq] = read_references(references_dir, "fasta")

    samples: Dict[str, List[Seq]] = {
        x: read_sequences(tests_dir + x, "fastq", SAMPLES_PER_BACTERIA) for x in os.listdir(tests_dir)
    }
    end: float = time()

    print("Time for reading: ", end - start)

    df: pd.DataFrame = create_dataframe(sequences, samples)

    if args.save_to_file:
        df.to_csv(f"results_k_{CONSTANT_K}_t_{TOLERANCE}.csv")

    print(df.to_string())

    total_end: float = time()
    print("Total time: ", total_end - start)