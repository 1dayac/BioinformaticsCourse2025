import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
# from adjustText import adjust_text    #not recommended because plot is dirtier, for labels.



###1  count distance between acids

genetic_code = {
    'A': ['GCU', 'GCC', 'GCA', 'GCG'],
    'C': ['UGU', 'UGC'],
    'D': ['GAU', 'GAC'],
    'E': ['GAA', 'GAG'],
    'F': ['UUU', 'UUC'],
    'G': ['GGU', 'GGC', 'GGA', 'GGG'],
    'H': ['CAU', 'CAC'],
    'I': ['AUU', 'AUC', 'AUA'],
    'K': ['AAA', 'AAG'],
    'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
    'M': ['AUG'],
    'N': ['AAU', 'AAC'],
    'P': ['CCU', 'CCC', 'CCA', 'CCG'],
    'Q': ['CAA', 'CAG'],
    'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
    'T': ['ACU', 'ACC', 'ACA', 'ACG'],
    'V': ['GUU', 'GUC', 'GUA', 'GUG'],
    'W': ['UGG'],
    'Y': ['UAU', 'UAC']}
# print(len(genetic_code))

def triplet_distance(codon1, codon2):
    return sum(n1 != n2 for n1, n2 in zip(codon1, codon2))

#calculation of HW metric

distances = defaultdict(dict)
for A1 in genetic_code:
    for A2 in genetic_code:
        
        avg_dist, cnt = 0, 0
        if A1!=A2:
            for tr1 in genetic_code[A1]:
                for tr2 in genetic_code[A2]:
                    avg_dist += triplet_distance(tr1,tr2);  cnt+=1
        else:
            cnt = 1

        distances[A1][A2] = avg_dist / cnt












##2 figure out how to compare to blosum (flatten + correlate)


# BLOSUM62 matrix as a dictionary of dictionaries
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


L1,L2 = [],[]
for A1 in genetic_code:
    for A2 in genetic_code:
        L1.append(distances[A1][A2]); L2.append(blosum62[A1][A2])


upper_indices = np.triu_indices(20, k=1) # keep comparisons between diagonals, but no twice counting


L1np = np.array(L1)
L2np = np.array(L2)
L1_trig = np.array(L1).reshape(20, 20)[upper_indices].flatten()
L2_trig = np.array(L2).reshape(20, 20)[upper_indices].flatten()



pearson_corr = np.corrcoef(L1_trig, L2_trig)[0, 1]
print(f"Pearson correlation: {pearson_corr:.3f}")

def plot():    

    plt.figure(figsize=(8, 6))
    plt.scatter(L1_trig, L2_trig, alpha=0.6, color='blue', edgecolor='w')



 ###LABELS DIRTY the pic, SO I removed them.
    # amino_acids = sorted(genetic_code.keys())  # Ensure consistent order
    # n = len(amino_acids)
    # labels = []
    # for i in range(n):
    #     for j in range(i+1, n):  
    #         labels.append(f"{amino_acids[i]}-{amino_acids[j]}")
    # texts = [
    #     plt.text(x, y, label, fontsize=8, ha='center', va='center') 
    #     for x, y, label in zip(L1_trig, L2_trig, labels)
    # ]
    # adjust_text(texts)
    

    
    m, b = np.polyfit(L1_trig, L2_trig, 1)
    plt.plot(L1_trig, m * L1_trig + b, 'r--', label=f'Fit: y = {m:.2f}x + {b:.2f}')

    plt.xlabel("Codon Distance (my HW metric)", fontsize=12)
    plt.ylabel("BLOSUM62 Score", fontsize=12)
    plt.title(f"Correlation: r = {pearson_corr:.3f}", fontsize=14)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.show()
plot()

