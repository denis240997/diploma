import argparse
import os
import sys
import time
from typing import Tuple


def load_dataset(sequences_path: str) -> Tuple[dict, int]:
    data_list = []
    with open(sequences_path) as fh:
        for line in fh:
            first = line.split()[0]
            data_list.append(first)
    data_list.sort()
    return {('',): data_list}, len(data_list[0])


def add_combinations(diagonal: tuple, sequences: list, result: dict,
                     rest_start: int, group_size: int):
    res_start_str = str(rest_start)
    for i in range(1, group_size):
        for j in range(0, i):
            new_diagonal = diagonal + (sys.intern(sequences[j][rest_start] +
                                                  res_start_str +
                                                  sequences[i][rest_start]),)
            if new_diagonal in result:
                result[new_diagonal].append(sequences[i])
            else:
                result[new_diagonal] = [sequences[i]]


def process_recursive(diagonal: tuple, sequences: list, result: dict,
                      start: int, end: int, rest_start: int, rest_end: int):
    groups_dict = {}
    for sequence in sequences:
        subsequence = sequence[start:end]
        if subsequence in groups_dict:
            groups_dict[subsequence].append(sequence)
        else:
            groups_dict[subsequence] = [sequence]

    half = (rest_start + rest_end) // 2
    left_ready = rest_start == half and rest_end - half == 1
    left_deeper = half != rest_end
    right_ready = rest_end == half and half - rest_start == 1
    right_deeper = half != rest_start
    for _, group in groups_dict.items():
        group_size = len(group)
        if group_size > 1:
            if left_ready:
                add_combinations(diagonal, group, result, half, group_size)
            elif left_deeper:
                process_recursive(diagonal, group, result,
                                  rest_start, half, half, rest_end)
            if right_ready:
                add_combinations(diagonal, group, result, rest_start, group_size)
            elif right_deeper:
                process_recursive(diagonal, group, result,
                                  half, rest_end, rest_start, half)


def process_diagonals(diagonals: dict, seq_len: int) -> dict:
    result = {}
    for diag, sequences in diagonals.items():
        diag_end = int(diag[-1][1:-1]) + 1 if diag[-1] else 0
        process_recursive(diag, sequences, result, 0, diag_end, diag_end, seq_len)
    return result


def process_dataset(sequences_path: str, dst_path: str):
    dataset, seq_len = load_dataset(sequences_path)

    dimension = 1
    while True:
        dataset = process_diagonals(dataset, seq_len)
        number_diagonals = len(dataset)
        print('Number of %2d dimensional diagonals: %10d'
              % (dimension, number_diagonals))
        number_hypercubes = 0

        with open(os.path.join(dst_path, 'hypercubes_' + str(dimension) + '.txt'),
                  'w') as f:
            for diag, sequences in sorted(dataset.items(), key=lambda x: x[0]):
                sequences.sort()  # crucial for next iterations
                number_hypercubes += len(sequences)
                diag_name = ':'.join(diag[1:])
                f.write(''.join(diag_name + '\t' + seq + '\n' for seq in sequences))

        print('Number of %2d dimensional hypercubes: %9d'
              % (dimension, number_hypercubes))
        dimension += 1
        if number_diagonals == number_hypercubes:
            break


if __name__ == '__main__':
    start_time = time.time()

    p = argparse.ArgumentParser()
    p.add_argument('-s', '--sequences',
                   help='the filename with the list of measured genotypes',
                   required=True)
    p.add_argument('-d', '--dst',
                   help='name of destination folder for hypercubes'
                        '"hypercubes" by default',
                   default='hypercubes')
    args = p.parse_args()

    if not os.path.exists(args.dst):
        os.makedirs(args.dst)

    process_dataset(args.sequences, args.dst)

    print(time.time() - start_time)
