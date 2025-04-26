from Bio import SeqIO, pairwise2
from Bio.Align import substitution_matrices
from code import needleman_wunsch

def read_two_seqs_from_fasta(path):
    records = list(SeqIO.parse(path, "fasta"))
    assert len(records) >= 2, "В файле должно быть как минимум две последовательности"
    return str(records[0].seq), str(records[1].seq)

sequences = [
    "AAA", 
    "AAA",
    "GATTACA",
    "GCATGCT",
    "MEANLY",
    "PLEASANTLY",
    "PAWHEAE",
    "HEAGAWGHEE"
]

fasta_paths = [
    "C:\\Users\\Redmi\\Desktop\\Биоинформатика\\BioinformaticsCourse2025\\homework\\1_1\\data\\gattaca.fasta",
    "C:\\Users\\Redmi\\Desktop\\Биоинформатика\\BioinformaticsCourse2025\\homework\\1_1\\data\\f8.fasta"
]

for path in fasta_paths:
    seq1, seq2 = read_two_seqs_from_fasta(path)
    sequences.append(seq1)
    sequences.append(seq2)

matrix = substitution_matrices.load("BLOSUM62")
gap_open = -10
gap_extend = -10 

def test_score_matches(seq1, seq2):
    bio_alignments = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)
    bio_score = bio_alignments[0][2]

    try:
        my_aln1, my_aln2, my_score = needleman_wunsch(seq1, seq2, matrix, gap_open)
    except Exception as e:
        print("Ошибка в реализации:")
        raise e

    assert abs(my_score - bio_score) < 1e-6, f"Несовпадение весов: мой {my_score}, Biopython {bio_score}"
    
    print(f"Мой вес = {my_score}, Biopython = {bio_score}")

if __name__ == "__main__":
    print("Тест: Нидлман–Вунш на всех последовательностях из списка")

    for i in range(0, len(sequences), 2):
        seq1 = sequences[i]
        seq2 = sequences[i+1]
        test_score_matches(seq1, seq2)
