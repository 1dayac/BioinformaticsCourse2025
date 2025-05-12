import os

from Bio import SeqIO as Sq
from Bio.Align import PairwiseAligner

from hw1 import find_hamming_dist, find_closest_substring, find_levenshtein_dist

# Task 1
script_dir = os.path.dirname(os.path.abspath(__file__))

path_to_fasta = os.path.join(script_dir, "data", "gattaca.fasta")
records = list(Sq.parse(path_to_fasta, "fasta"))

s1 = records[0]
s2 = records[1]

print(find_hamming_dist(s1, s2))

# Task 2
path_to_fasta = os.path.join(script_dir, "data", "f8.fasta")
records = list(Sq.parse(path_to_fasta, "fasta"))

s1 = str(records[0].seq)
s2 = str(records[1].seq)

index, substring, distance = find_closest_substring(s1, s2)
print(index, distance)

# Task 3
print(find_levenshtein_dist(s1, s2))
aligner = PairwiseAligner()
aligner.mode = 'global'

s1, s2 = aligner.align(s1, s2)[0]
distance = find_hamming_dist(s1, s2)

assert find_levenshtein_dist(s1, s2) == distance
