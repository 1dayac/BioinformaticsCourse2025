import unittest

from pathlib import Path
from Bio import SeqIO as Sq
from Bio.Align import substitution_matrices
from Bio.Align import PairwiseAligner

from needleman_wunsch import needleman_wunsch

class TestNeedlemanWunsch(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.matrix = substitution_matrices.load("BLOSUM62")
        cls.gap_cost = -2
        cls.align_aligner = PairwiseAligner()
        cls.align_aligner.mode = 'global'
        cls.align_aligner.substitution_matrix = cls.matrix
        cls.align_aligner.open_gap_score = cls.gap_cost
        cls.align_aligner.extend_gap_score = 0

    def test_example_sequences(self):
        seq1 = "PLEASANT"
        seq2 = "ELEPHANT"
        
        _, score = needleman_wunsch(seq1, seq2, self.matrix, self.gap_cost)
        alignments = self.align_aligner.align(seq1, seq2)
        pa_alignment = next(iter(alignments))
        pa_score = pa_alignment.score
        
        self.assertAlmostEqual(score, pa_score, places=5,
                               msg=f"Оценка функции {score} не совпадает с оценкой PairwiseAligner {pa_score}")

    def test_known_alignment(self):
        seq1 = "MEEPQSDPSV"
        seq2 = "MEESQSDPSV"
        
        _, score = needleman_wunsch(seq1, seq2, self.matrix, self.gap_cost)
        alignments = self.align_aligner.align(seq1, seq2)
        pa_alignment = next(iter(alignments))
        pa_score = pa_alignment.score

        self.assertAlmostEqual(score, pa_score, places=5,
                               msg="Оценка выравнивания для известных последовательностей не совпадает с PairwiseAligner")
    
    def test_sequnce_from_hw1(self):
        script_dir = Path(__file__).parent.parent

        with open(script_dir / "1_1" / "data" / "gattaca.fasta", "rt") as FF:
            records = list(Sq.parse(FF, "fasta"))
        
        seq1 = records[0]
        seq2 = records[1]

        _, score = needleman_wunsch(seq1, seq2, self.matrix, self.gap_cost)
        alignments = self.align_aligner.align(seq1, seq2)
        pa_alignment = next(iter(alignments))
        pa_score = pa_alignment.score

        self.assertAlmostEqual(score, pa_score, places=5,
                               msg="Оценка выравнивания для последовательности из ДЗ1 не совпадает с PairwiseAligner")

if __name__ == '__main__':
    unittest.main()