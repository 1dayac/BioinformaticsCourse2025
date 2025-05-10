from pathlib import Path

from Bio import SeqIO as Sq
from Bio.Align import PairwiseAligner

from hw1 import hamming_distance, find_closest_substring, levenshtein_distance

# Task 1
script_dir = Path(__file__).parent
FF = script_dir / "data" / "gattaca.fasta"
records = list(Sq.parse(FF, "fasta"))

seq1 = records[0]
seq2 = records[1]

print(hamming_distance(seq1, seq2))

# Task 2
FF = script_dir / "data" / "f8.fasta"
records = list(Sq.parse(FF, "fasta"))
seq1 = str(records[0].seq)
seq2 = str(records[1].seq)

if len(seq1) >= len(seq2):
    long_seq = seq1
    short_seq = seq2
else:
    long_seq = seq2
    short_seq = seq1

index, substring, distance = find_closest_substring(long_seq, short_seq)
print(index, distance)

# Task 3
print(levenshtein_distance(seq1, seq2))
aligner = PairwiseAligner()
aligner.mode = 'global'

s1, s2 = aligner.align(seq1, seq2)[0]
distance = hamming_distance(s1, s2)

assert levenshtein_distance(seq1, seq2) == distance
