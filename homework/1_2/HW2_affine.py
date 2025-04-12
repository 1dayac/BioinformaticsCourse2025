import numpy as np
from Bio.Align import substitution_matrices
from Bio import pairwise2

global MIN
global gap_open
global gap_extend
global matrix
MIN = -np.inf
gap_open = -10
gap_extend = -1
matrix = substitution_matrices.load('BLOSUM62')


def needleman_wunsch_affine(seq1, seq2):
    dim_i = len(seq2) + 1
    dim_j = len(seq1) + 1

    I_x = np.full((dim_i, dim_j), MIN)
    I_x[0, :] = gap_open + gap_extend * np.arange(dim_j)
    I_x[0, 0] = 0
    I_x[1:, :] = I_x[0, :]
    I_x[1:, 0] = MIN

    I_y = np.full((dim_i, dim_j), MIN)
    I_y[:, 0] = gap_open + gap_extend * np.arange(dim_i)
    I_y[0, 0] = 0
    I_y[:, 1:] = I_y[:, 0:1]
    I_y[0, 1:] = MIN

    D = np.full((dim_i, dim_j), MIN)
    D[0, 0] = 0
    D[1:, 1:] = 0

    for j in range(1, dim_j):
        for i in range(1, dim_i):
            I_x[i][j] = max(gap_open + D[i][j - 1], gap_extend + I_x[i][j - 1])
            I_y[i][j] = max(gap_open + D[i - 1][j], gap_extend + I_y[i - 1][j])
            D[i][j] = max(matrix[seq1[j - 1]][seq2[i - 1]] + D[i - 1][j - 1], I_x[i][j], I_y[i][j])

    align_seq1 = ''
    align_seq2 = ''
    i = len(seq2)
    j = len(seq1)
    while i > 0 or j > 0:
        if i > 0 and j > 0 and D[i][j] == D[i - 1][j - 1] + matrix[seq1[j - 1]][seq2[i - 1]]:
            align_seq1 += seq1[j - 1]
            align_seq2 += seq2[i - 1]
            i -= 1
            j -= 1
        elif i > 0 and D[i][j] == I_y[i][j]:
            align_seq1 += '-'
            align_seq2 += seq2[i - 1]
            i -= 1
        elif j > 0 and D[i][j] == I_x[i][j]:
            align_seq1 += seq1[j - 1]
            align_seq2 += '-'
            j -= 1

    return align_seq1[::-1], align_seq2[::-1], D[-1][-1]
