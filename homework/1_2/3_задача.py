
import itertools
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from Bio.Data import CodonTable
from scipy.stats import pearsonr

blosum62 = {
    'A': {'A': 4,  'C': 0,  'D': -2, 'E': -1, 'F': -2, 'G': 0,  'H': -2, 'I': -1, 'K': -1, 'L': -1,
          'M': -1, 'N': -2, 'P': -1, 'Q': -1, 'R': -1, 'S': 1,  'T': 0,  'V': 0,  'W': -3, 'Y': -2},
    'C': {'A': 0,  'C': 9,  'D': -3, 'E': -4, 'F': -2, 'G': -3, 'H': -3, 'I': -1, 'K': -3, 'L': -1,
          'M': -1, 'N': -3, 'P': -3, 'Q': -3, 'R': -3, 'S': -1, 'T': -1, 'V': -1, 'W': -2, 'Y': -2},
    'D': {'A': -2, 'C': -3, 'D': 6,  'E': 2,  'F': -3, 'G': -1, 'H': -1, 'I': -2, 'K': -1, 'L': -4,
          'M': -3, 'N': 1,  'P': -1, 'Q': 0,  'R': -2, 'S': 0,  'T': -1, 'V': -2, 'W': -4, 'Y': -3},
    'E': {'A': -1, 'C': -4, 'D': 2,  'E': 5,  'F': -3, 'G': -2, 'H': 0,  'I': -2, 'K': 1,  'L': -3,
          'M': -2, 'N': 0,  'P': -1, 'Q': 2,  'R': 0,  'S': 0,  'T': -1, 'V': -2, 'W': -3, 'Y': -2},
    'F': {'A': -2, 'C': -2, 'D': -3, 'E': -3, 'F': 6,  'G': -3, 'H': -1, 'I': 0,  'K': -3, 'L': 0,
          'M': 0,  'N': -3, 'P': -4, 'Q': -3, 'R': -3, 'S': -2, 'T': -2, 'V': -1, 'W': 1,  'Y': 3},
    'G': {'A': 0,  'C': -3, 'D': -1, 'E': -2, 'F': -3, 'G': 6,  'H': -2, 'I': -4, 'K': -2, 'L': -4,
          'M': -3, 'N': 0,  'P': -2, 'Q': -2, 'R': -2, 'S': 1,  'T': 0,  'V': -2, 'W': -2, 'Y': -3},
    'H': {'A': -2, 'C': -3, 'D': -1, 'E': 0,  'F': -1, 'G': -2, 'H': 8,  'I': -3, 'K': -1, 'L': -3,
          'M': -2, 'N': 1,  'P': -2, 'Q': -1, 'R': 0,  'S': -1, 'T': -2, 'V': -3, 'W': 2,  'Y': 2},
    'I': {'A': -1, 'C': -1, 'D': -2, 'E': -2, 'F': 0,  'G': -4, 'H': -3, 'I': 4,  'K': -3, 'L': 2,
          'M': 1,  'N': -2, 'P': -3, 'Q': -2, 'R': -3, 'S': -2, 'T': -1, 'V': 3,  'W': -3, 'Y': -1},
    'K': {'A': -1, 'C': -3, 'D': -1, 'E': 1,  'F': -3, 'G': -2, 'H': -1, 'I': -3, 'K': 5,  'L': -2,
          'M': -1, 'N': 0,  'P': -1, 'Q': 1,  'R': 2,  'S': 0,  'T': -1, 'V': -2, 'W': -3, 'Y': -2},
    'L': {'A': -1, 'C': -1, 'D': -4, 'E': -3, 'F': 0,  'G': -4, 'H': -3, 'I': 2,  'K': -2, 'L': 4,
          'M': 2,  'N': -2, 'P': -3, 'Q': -2, 'R': -2, 'S': -2, 'T': -1, 'V': 1,  'W': -2, 'Y': -1},
    'M': {'A': -1, 'C': -1, 'D': -3, 'E': -2, 'F': 0,  'G': -3, 'H': -2, 'I': 1,  'K': -1, 'L': 2,
          'M': 5,  'N': -2, 'P': -2, 'Q': 0,  'R': -1, 'S': -1, 'T': -1, 'V': 1,  'W': -1, 'Y': -1},
    'N': {'A': -2, 'C': -3, 'D': 1,  'E': 0,  'F': -3, 'G': 0,  'H': 1,  'I': -2, 'K': 0,  'L': -2,
          'M': -2, 'N': 6,  'P': -2, 'Q': 0,  'R': 0,  'S': 1,  'T': 0,  'V': -2, 'W': -4, 'Y': -2},
    'P': {'A': -1, 'C': -3, 'D': -1, 'E': -1, 'F': -4, 'G': -2, 'H': -2, 'I': -3, 'K': -1, 'L': -3,
          'M': -2, 'N': -2, 'P': 7,  'Q': -1, 'R': -2, 'S': -1, 'T': -1, 'V': -2, 'W': -4, 'Y': -3},
    'Q': {'A': -1, 'C': -3, 'D': 0,  'E': 2,  'F': -3, 'G': -2, 'H': -1, 'I': -2, 'K': 1,  'L': -2,
          'M': 0,  'N': 0,  'P': -1, 'Q': 5,  'R': 1,  'S': 0,  'T': -1, 'V': -2, 'W': -2, 'Y': -1},
    'R': {'A': -1, 'C': -3, 'D': -2, 'E': 0,  'F': -3, 'G': -2, 'H': 0,  'I': -3, 'K': 2,  'L': -2,
          'M': -1, 'N': 0,  'P': -2, 'Q': 1,  'R': 5,  'S': -1, 'T': -1, 'V': -2, 'W': -3, 'Y': -2},
    'S': {'A': 1,  'C': -1, 'D': 0,  'E': 0,  'F': -2, 'G': 1,  'H': -1, 'I': -2, 'K': 0,  'L': -2,
          'M': -1, 'N': 1,  'P': -1, 'Q': 0,  'R': -1, 'S': 4,  'T': 1,  'V': -2, 'W': -3, 'Y': -2},
    'T': {'A': 0,  'C': -1, 'D': -1, 'E': -1, 'F': -2, 'G': 0,  'H': -2, 'I': -1, 'K': -1, 'L': -1,
          'M': -1, 'N': 0,  'P': -1, 'Q': -1, 'R': -1, 'S': 1,  'T': 5,  'V': 0,  'W': -2, 'Y': -2},
    'V': {'A': 0,  'C': -1, 'D': -2, 'E': -2, 'F': -1, 'G': -2, 'H': -3, 'I': 3,  'K': -2, 'L': 1,
          'M': 1,  'N': -2, 'P': -2, 'Q': -2, 'R': -2, 'S': -2, 'T': 0,  'V': 4,  'W': -3, 'Y': -1},
    'W': {'A': -3, 'C': -2, 'D': -4, 'E': -3, 'F': 1,  'G': -2, 'H': 2,  'I': -3, 'K': -3, 'L': -2,
          'M': -1, 'N': -4, 'P': -4, 'Q': -2, 'R': -3, 'S': -3, 'T': -2, 'V': -3, 'W': 11, 'Y': 2},
    'Y': {'A': -2, 'C': -2, 'D': -3, 'E': -2, 'F': 3,  'G': -3, 'H': 2,  'I': -1, 'K': -2, 'L': -1,
          'M': -1, 'N': -2, 'P': -3, 'Q': -1, 'R': -2, 'S': -2, 'T': -2, 'V': -1, 'W': 2,  'Y': 7} }


table = CodonTable.unambiguous_dna_by_name["Standard"]

aa_to_codons = {}
for codon, aa in table.forward_table.items():
    aa_to_codons.setdefault(aa, []).append(codon)
aa_list = sorted(set(aa_to_codons.keys()))




A = pd.DataFrame(index=aa_list, columns=aa_list, dtype=float)

for aa1, aa2 in itertools.product(aa_list, repeat=2):
    cod_1 = aa_to_codons[aa1]
    cod_2 = aa_to_codons[aa2]

    dists = []
    for c1 in cod_1:
        for c2 in cod_2:
            diff = 0
            for a, b in zip(c1, c2):
                if a != b:
                    diff += 1
            dists.append(diff)
    A.loc[aa1, aa2] = np.mean(dists)

print("Матрица средних расстояний между кодонами (A):")
print(A.round(2))

blosum = pd.DataFrame(index=aa_list, columns=aa_list, dtype=float)

for aa1 in aa_list:
    for aa2 in aa_list:
        score = None
        if aa1 in blosum62 and aa2 in blosum62[aa1]:
            score = blosum62[aa1][aa2]
        elif aa2 in blosum62 and aa1 in blosum62[aa2]:
            score = blosum62[aa2][aa1]
        blosum.loc[aa1, aa2] = score

f_A = A.values[np.triu_indices(len(aa_list), k=1)]
f_B = blosum.values[np.triu_indices(len(aa_list), k=1)]

mask = ~np.isnan(f_B)
corr, pval = spearmanr(f_A[mask], f_B[mask])
print(f"\nКорреляция Спирмена с BLOSUM62: r = {corr:.3f}, p-value = {pval:.3g}")

corr_p, pval_p = pearsonr(f_A[mask], f_B[mask])
print(f"Корреляция с BLOSUM62: r = {corr_p:.3f}, p-value = {pval_p:.3g}")