import numpy as np
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

    return align_seq1, align_seq2, D[-1][-1]
