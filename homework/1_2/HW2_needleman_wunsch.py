import numpy as np
from random import choice, randint
from Bio import pairwise2
from Bio.Align import substitution_matrices

global gap
global matrix
gap = -2
matrix = substitution_matrices.load('BLOSUM62')


def needleman_wunsch(seq1, seq2):
    dim_i = len(seq1) + 1
    dim_j = len(seq2) + 1

    D = np.zeros((dim_i, dim_j), dtype=int)

    D[1:, 0] = np.arange(1, dim_i) * gap
    D[0, 1:] = np.arange(1, dim_j) * gap

    for i in range(1, dim_i):
        for j in range(1, dim_j):
            D[i][j] = max(D[i - 1][j - 1] + matrix[seq1[i - 1]][seq2[j - 1]], D[i - 1][j] + gap, D[i][j - 1] + gap)

    align_seq1 = ''
    align_seq2 = ''
    i = len(seq1)
    j = len(seq2)
    while i > 0 or j > 0:
        if i > 0 and j > 0 and D[i][j] == D[i - 1][j - 1] + matrix[seq1[i - 1]][seq2[j - 1]]:
            align_seq1 += seq1[i - 1]
            align_seq2 += seq2[j - 1]
            i -= 1
            j -= 1
        elif D[i][j] == D[i - 1][j] + gap:
            align_seq1 += seq1[i - 1]
            align_seq2 += '-'
            i -= 1
        else:
            align_seq1 += '-'
            align_seq2 += seq2[j - 1]
            j -= 1

    return align_seq1[::-1], align_seq2[::-1], D[-1][-1]


def generate_sequences(min_len=2, max_len=10):
    len = randint(min_len, max_len)
    return ''.join(choice('ATGC') for _ in range(len))


def test():
    k_error = 0
    sequences = [generate_sequences() for _ in range(10)]

    for first_seq in range(len(sequences)):
        for second_seq in range(first_seq + 1, len(sequences)):
            s1, s2, score = needleman_wunsch(sequences[first_seq], sequences[second_seq])

            alignments = pairwise2.align.globalds(
                sequences[first_seq],
                sequences[second_seq],
                matrix,
                gap,
                gap
            )

            if alignments:
                pairwise_alignment = alignments[0]
                pairwise_score = pairwise_alignment.score
                pairwise_seq1 = pairwise_alignment.seqA
                pairwise_seq2 = pairwise_alignment.seqB

                if score != pairwise_score:
                    print(f"Исходные последовательности: {sequences[first_seq]}, {sequences[second_seq]}")
                    print(f"Score: {score}")
                    print(f"{s1}\n{s2}")
                    print(f"Bio.pairwise2:\n{pairwise_seq1}\n{pairwise_seq2}")
                    print(f"Bio.pairwise2 score: {pairwise_score}")
                    print('__________')
                    k_error += 1
    return k_error


print(test())
