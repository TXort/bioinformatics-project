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


def parse_arguments() -> argparse.Namespace | ValueError | FileNotFoundError:
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
    """
    This function replaces any 'N' nucleotides in a sequence with a random nucleotide.

    Parameters:
    seq (Seq): The sequence to be processed.
    rlist (List[str]): A list of nucleotides to choose from when replacing 'N'. Default value is ["A", "C", "G", "T"].

    Returns:
    Seq: The processed sequence, where any 'N' nucleotides have been replaced with a random nucleotide.

    Note:
    The function iterates over the sequence, replacing any 'N' nucleotides it encounters with a random nucleotide from 'rlist'. The random nucleotide is chosen using the randint function from the random module.
    """
    return Seq("".join([rlist[randint(0, 3)] if x == "N" else x for x in seq]))


def preprocess_distribution(distr: Dict[str, float]) -> List[float]:
    """
    This function preprocesses a distribution of k-mers.

    Parameters:
    distr (Dict[str, float]): The distribution to be preprocessed. This is a dictionary where the keys are the k-mers and the values are the frequencies of the k-mers.

    Returns:
    List[float]: A list of floats representing the preprocessed distribution. The k-mers are sorted in lexicographical order, and the corresponding frequencies are returned in the same order.

    Note:
    The function sorts the keys of the distribution (i.e., the k-mers) in lexicographical order, and then returns a list of the corresponding values (i.e., the frequencies of the k-mers).
    """
    return [distr[x] for x in sorted(distr.keys())]


def calc_cosine_similarity(distr1: List[float], distr2: List[float]) -> float:
    """
    This function calculates the cosine similarity between two distributions.

    Parameters:
    distr1 (List[float]): The first distribution, represented as a list of floats.
    distr2 (List[float]): The second distribution, represented as a list of floats.

    Returns:
    float: The cosine similarity between the two distributions. This is a value between 0 and 1, where 1 means the distributions are identical and 0 means they are completely dissimilar.

    Note:
    The function uses the scipy.spatial.distance.cosine function to calculate the cosine distance between the two distributions, and then subtracts this from 1 to get the cosine similarity.
    """
    return 1 - spatial.distance.cosine(distr1, distr2)


def get_distribution(seq: Seq, k: int = 3) -> Dict[str, float]:
    """
    This function calculates the distribution of k-mers in a given sequence.

    Parameters:
    seq (Seq): The sequence for which the k-mer distribution is to be calculated.
    k (int): The length of the k-mers. Default value is 3.

    Returns:
    Dict[str, float]: A dictionary where the keys are the k-mers and the values are the frequencies of the k-mers in the sequence.

    Note:
    The function first initializes a dictionary with all possible k-mers as keys and 0 as values. It then iterates over the sequence, incrementing the count of each k-mer it encounters. Finally, it calculates the total count of all k-mers and returns a dictionary where the values are the frequencies of the k-mers (i.e., the count of each k-mer divided by the total count).
    """
    k_mer: Dict[str, int] = defaultdict(int, {k: 0 for k in DEFAULT_DICT.keys()})

    for i in range(len(seq) - k + 1):
        k_mer[seq[i:i + k]] += 1

    total: int = sum(k_mer.values())

    return {k: v / total for k, v in k_mer.items()}


def get_preprocessed_distribution(seq: Seq, k: int = 3) -> List[float]:
    """
    This function takes a sequence and an integer k, and returns a preprocessed distribution of the sequence.

    Parameters:
    seq (Seq): The sequence to be processed.
    k (int): The length of the k-mers to be used in the distribution. Default value is 3.

    Returns:
    List[float]: A list of floats representing the preprocessed distribution of the sequence.

    Note:
    The function first gets the distribution of the sequence using the get_distribution function, and then preprocesses this distribution using the preprocess_distribution function.
    """
    return preprocess_distribution(get_distribution(seq, k))


def read_references(file_path: str, format: str) -> Dict[str, Seq]:
    """
    This function reads the first sequence from each file in a directory and returns a dictionary mapping file names to sequences.

    Parameters:
    file_path (str): The path to the directory containing the files.
    format (str): The format of the files (e.g., 'fasta', 'fastq').

    Returns:
    Dict[str, Seq]: A dictionary where the keys are the names of the files in the directory and the values are the first sequence from each file. The sequences are processed to replace any 'N' nucleotides with a random nucleotide.

    Note:
    The function will only read the first sequence from each file. If a file contains more than one sequence, all sequences after the first are ignored.
    """
    references: Dict[str, Seq] = {}
    for file in os.listdir(file_path):
        for record in SeqIO.parse(file_path + file, format):
            references[file] = replace_with_random(record.seq)
            break
    return references


def read_sequences(file_path: str, format: str, limit: int) -> List[Seq]:
    """
    This function reads sequences from a file and returns a list of sequences.

    Parameters:
    file_path (str): The path to the file containing the sequences.
    format (str): The format of the file (e.g., 'fasta', 'fastq').
    limit (int): The maximum number of sequences to read from the file.

    Returns:
    List[Seq]: A list of sequences read from the file. The sequences are processed to replace any 'N' nucleotides with a random nucleotide.

    Note:
    The function will stop reading the file once it has read 'limit' number of sequences. Also, any sequence that is shorter than the global constant 'CONSTANT_K' is skipped.
    """
    sequences: List[Seq] = []
    for record in SeqIO.parse(file_path, format):
        if len(record.seq) < CONSTANT_K:
            continue
        sequences.append(replace_with_random(record.seq))
        if len(sequences) >= limit:
            break
    return sequences


def compare_and_assign(ref_d: Dict[str, List[float]], seqs: List[Seq]) -> Dict[str, int]:
    """
    This function compares each sequence in a list with a reference distribution and assigns the sequence to the most similar reference.

    Parameters:
    ref_d (Dict[str, List[float]]): The reference distribution. This is a dictionary where the keys are the names of the references and the values are the preprocessed distributions of the references.
    seqs (List[Seq]): The list of sequences to be compared with the reference distribution.

    Returns:
    Dict[str, int]: A dictionary where the keys are the names of the references and the values are the number of sequences assigned to each reference. If a sequence is not similar enough to any reference (i.e., the cosine similarity is less than the global constant 'TOLERANCE'), it is assigned to 'UNCLASSIFIED'.

    Note:
    The function first preprocesses the distribution of each sequence (including its reverse complement and complement). It then calculates the cosine similarity between each preprocessed distribution and the reference distribution. If the cosine similarity is greater than the current maximum similarity and greater than or equal to 'TOLERANCE', the sequence is assigned to the current reference. If the sequence is not similar enough to any reference, it is assigned to 'UNCLASSIFIED'.
    """
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
    """
    This function creates a pandas DataFrame that represents the classification results of the sequences.

    Parameters:
    ref (Dict[str, Seq]): The reference sequences. This is a dictionary where the keys are the names of the references and the values are the sequences of the references.
    samples (Dict[str, List[Seq]]): The sample sequences. This is a dictionary where the keys are the names of the samples and the values are lists of sequences in the samples.

    Returns:
    Tuple[pd.DataFrame, pd.DataFrame]: A pandas DataFrame where the index is the names of the references and the columns are the names of the samples. The values in the DataFrame are the number of sequences in each sample that were classified as each reference. If a sequence in a sample was not classified as any reference (i.e., it was 'UNCLASSIFIED'), it is counted in the 'UNCLASSIFIED' row.

    Note:
    The function first preprocesses the distribution of each reference sequence using multiprocessing. It then initializes a DataFrame with the names of the references as the index and the names of the samples as the columns, and sets the initial value of all cells to 0. It then compares each sequence in each sample with the preprocessed distributions of the references using multiprocessing, and increments the count in the corresponding cell of the DataFrame. Finally, it fills any remaining NaN values in the DataFrame with 0 and returns the DataFrame.
    """
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
