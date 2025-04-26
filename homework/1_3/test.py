import sys
import os
import unittest
from Bio import SeqIO

sys.path.insert(0, os.path.dirname(__file__))

from hirschberg_recursion_tree.algo import hirschberg_tree, Node

class TestHirschberg(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if len(sys.argv) != 2:
            raise RuntimeError("Usage: python test.py <input.fasta>")
        fasta_file = sys.argv[1]

        with open(fasta_file, "r") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
        if len(records) < 2:
            raise RuntimeError(f"В файле {fasta_file} менее двух записей FASTA.")

        cls.headerA = records[0].id
        cls.A = str(records[0].seq)
        cls.headerB = records[1].id
        cls.B = str(records[1].seq)

        cls.del_cost = 2
        cls.ins_cost = 2
        cls.mismatch_cost = 1
        cls.match_score = 3

        cls.tree = hirschberg_tree(
            cls.A, cls.B,
            del_cost=cls.del_cost,
            ins_cost=cls.ins_cost,
            mismatch_cost=cls.mismatch_cost,
            match_score=cls.match_score
        )

    def test_root_sequences(self):
        """Корневой узел содержит исходные последовательности"""
        self.assertEqual(self.tree.A, self.A)
        self.assertEqual(self.tree.B, self.B)

    def test_children_on_nontrivial(self):
        """Если длины >1, у корня должны быть два потомка"""
        if len(self.A) > 1 and len(self.B) > 1:
            self.assertEqual(len(self.tree.children), 2)

    def test_leaves_are_base(self):
        """Во всех листовых узлах хотя бы одна из строк пустая или длина 1"""
        def check(node: Node):
            if not node.children:
                self.assertTrue(
                    len(node.A) <= 1 or len(node.B) <= 1,
                    f"Leaf node has lengths A={len(node.A)}, B={len(node.B)}"
                )
            for child in node.children:
                check(child)
        check(self.tree)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python test.py <input.fasta>")
        sys.exit(1)
    fasta_file = sys.argv[1]
    TestHirschberg.fasta_file = fasta_file
    unittest.main(argv=[sys.argv[0]])
