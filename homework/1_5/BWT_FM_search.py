from pysuffixarray.core import SuffixArray
from Bio import SeqIO


class BWT_searcher():
    def __init__(self, reference):
        self.reference = reference
        self.suffix_array = SuffixArray(reference).suffix_array()
        self.bwt = ''
        self.alphabet = sorted('$ACGT')
        self.occurence = {}
        for char in self.alphabet:
            self.occurence[char] = []

        for i in range(len(reference) + 1):
            suff_i = self.suffix_array[i]
            char = reference[suff_i - 1] if suff_i else '$'
            self.bwt += char

            for key in self.alphabet:
                last_cnt = self.occurence[key][-1] if self.occurence[key] else 0
                self.occurence[key].append(last_cnt)
            self.occurence[char][-1] += 1

        self.cnt = {}
        sum_ = 0
        for char in self.alphabet:
            self.cnt[char] = sum_
            sum_ += self.occurence[char][-1]

    def bwt_pattern_search(self, pattern):
        R_L = 0
        R_H = len(self.bwt) - 1
        for char in reversed(pattern):
            R_L = self.cnt[char] + (0 if R_L == 0 else self.occurence[char][R_L - 1])
            R_H = self.cnt[char] + self.occurence[char][R_H] - 1
            if R_L > R_H:
                return []
        return self.suffix_array[R_L:R_H + 1]

    def bwt_pattern_search_with_mutation(self, pattern, n_mutations):
        # разбиваю на n_mutations + 1 частей
        # если мутаций ровно n_mutations, значит как минимум в одной из частей не будет мутаций,
        # её можно будет найти обычным поиском bwt_pattern_search,
        # затем относительно неё пройтись по участку генома и сравнить с целым ридом

        k = len(pattern) // (n_mutations + 1)
        segments = [pattern[i * k:i * k + k] for i in range(n_mutations + 1)]

        for i, segment in enumerate(segments):
            list_suff_i = self.bwt_pattern_search(segment)
            if not list_suff_i:
                continue

            len_pattern = len(pattern)
            len_reference = len(self.reference)

            for suff_i in list_suff_i:
                start = suff_i - i * k
                if start < 0 or start + len_pattern > len_reference:
                    continue
                mismatches = 0
                for j in range(len_pattern):
                    if self.reference[start + j] != pattern[j]:
                        mismatches += 1
                        if mismatches > n_mutations:
                            break
                if mismatches == n_mutations:
                    return True

        return False


reference = str(SeqIO.read("./BWT_folder/genome.fa", "fasta").seq)
genome_bwt = BWT_searcher(reference)

k0, k1, k5 = 0, 0, 0
reads0, reads1, reads5 = [], [], []
with open("./BWT_folder/sample_reads.fasta", "r") as file:
    for record in SeqIO.parse(file, "fasta"):
        read_sequence = str(record.seq)
        if len(genome_bwt.bwt_pattern_search(read_sequence)) > 0:
            k0 += 1
            reads0.append(read_sequence)
        elif genome_bwt.bwt_pattern_search_with_mutation(read_sequence, 1):
            k1 += 1
            reads1.append(read_sequence)
        elif genome_bwt.bwt_pattern_search_with_mutation(read_sequence, 5):
            k5 += 1
            reads5.append(read_sequence)

print('Нет мутаций, количество ридов:', k0)
'''
for j in reads0:
    print(j)
'''
print('1 мутация, количество ридов:', k1)
'''
for j in reads1:
    print(j)
'''
print('5 мутаций, количество ридов:', k5)
'''
for j in reads5:
    print(j)
'''