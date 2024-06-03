import csv
import os
import re
from typing import Tuple

def read_csv_files() -> list[str]:
    csv_files = []
    for file in os.listdir():
        if file.endswith(".csv"):
            csv_files.append(file)
    return csv_files

files: list[str] = read_csv_files()

for file in files:

    print(f"Analyzing {file}")

    dict2d: dict[str, dict[str, str]] = {}

    with open(file, 'r') as data:
        for line in csv.DictReader(data):
            key_raw: str = line['']
            if "_" in key_raw:
                key: str = line[''][0:line[''].find('_reference.fasta')]
            else:
                key: str = line['']
            dict2d[key] = {key: value for key, value in line.items() if key != ''}

    dict2d_int: dict[str, dict[str, int]] = {}

    for key in dict2d.keys():
        dict2d_int[key] = {key: int(value) for key, value in dict2d[key].items()}

    unclassified: dict[str, int] = dict2d_int['unclassified']
    del dict2d_int['unclassified']

    dict2d_max: dict[str, int] = {}

    for key in dict2d_int.keys():
        for key2 in dict2d_int[key].keys():
            if key2 not in dict2d_max.keys():
                dict2d_max[key2] = dict2d_int[key][key2]
            else:
                if dict2d_int[key][key2] > dict2d_max[key2]:
                    dict2d_max[key2] = dict2d_int[key][key2]

    for key_ref in dict2d_int.keys():
        max_key: str = ''
        max_val: int = 0

        for key_gen in dict2d_int[key_ref].keys():
            if dict2d_int[key_ref][key_gen] > max_val:
                max_val = dict2d_int[key_ref][key_gen]
                max_key = key_gen

        same: re.Match[str] | None = re.search(key_ref, max_key)
        if same:
            print(f"\tOk, max value for {key_ref} is {max_val} and it is {max_key}")
        else:
            print(f"\tError, max value for {key_ref} is {max_val} and it is {max_key}")

    correct_list: list[Tuple[int, int]] = []

    for k1, k2 in zip(sorted(dict2d_max.keys()), sorted(unclassified.keys())):
        if k1 != k2:
            print(f"\tError, {k1} is not equal to {k2}")
            continue
        sum_col: int = sum(dict2d_int[x][k1] for x in dict2d_int)
        calc: float = dict2d_max[k1] / (sum_col)
        print("\t" + k1 + "  correct: ",  calc )
        correct_list.append((dict2d_max[k1], sum_col))

    avg: float = 0
    for c in correct_list:
        avg += c[0] / c[1]
    avg /= len(correct_list)

    print(f"\tAverage: {avg}")
    print()

