import unittest
from solution import build_hirschberg_tree


class SimpleHirschbergTests(unittest.TestCase):
    def test_basic_example(self):
        """Test with the basic example provided in the assignment."""
        seq1 = "AGTACGCA"
        seq2 = "TATGC"

        # Build tree
        tree = build_hirschberg_tree(
            seq1, seq2,
            del_cost=2, ins_cost=2,
            mismatch_cost=1, match_score=3
        )

        # Verify root content
        self.assertEqual(tree.seq1, seq1)
        self.assertEqual(tree.seq2, seq2)

        # Check that tree has the correct structure
        self.assertIsNotNone(tree.left)
        self.assertIsNotNone(tree.right)

        # First sequence should be split at index 4
        self.assertEqual(tree.left.seq1, "AGTA")
        self.assertEqual(tree.right.seq1, "CGCA")

    def test_identical_sequences(self):
        """Test with identical sequences."""
        seq = "ACGT"

        tree = build_hirschberg_tree(
            seq, seq,
            del_cost=2, ins_cost=2,
            mismatch_cost=1, match_score=3
        )

        self.assertEqual(tree.seq1, seq)
        self.assertEqual(tree.seq2, seq)

        # Should be split in half
        self.assertEqual(tree.left.seq1, "AC")
        self.assertEqual(tree.left.seq2, "AC")
        self.assertEqual(tree.right.seq1, "GT")
        self.assertEqual(tree.right.seq2, "GT")

    def test_empty_sequence(self):
        """Test with an empty sequence."""
        seq = "ACGT"

        # Empty first sequence
        tree1 = build_hirschberg_tree(
            "", seq,
            del_cost=2, ins_cost=2,
            mismatch_cost=1, match_score=3
        )

        self.assertEqual(tree1.seq1, "")
        self.assertEqual(tree1.seq2, seq)
        self.assertIsNone(tree1.left)
        self.assertIsNone(tree1.right)

        # Empty second sequence
        tree2 = build_hirschberg_tree(
            seq, "",
            del_cost=2, ins_cost=2,
            mismatch_cost=1, match_score=3
        )

        self.assertEqual(tree2.seq1, seq)
        self.assertEqual(tree2.seq2, "")
        self.assertIsNone(tree2.left)
        self.assertIsNone(tree2.right)

    def test_single_character(self):
        """Test with single character sequences."""
        tree = build_hirschberg_tree(
            "A", "G",
            del_cost=2, ins_cost=2,
            mismatch_cost=1, match_score=3
        )

        self.assertEqual(tree.seq1, "A")
        self.assertEqual(tree.seq2, "G")
        self.assertIsNone(tree.left)
        self.assertIsNone(tree.right)

    def test_different_length_sequences(self):
        """Test with sequences of very different lengths."""
        seq1 = "ACGTACGTACGT"  # 12 characters
        seq2 = "ACG"  # 3 characters

        tree = build_hirschberg_tree(
            seq1, seq2,
            del_cost=2, ins_cost=2,
            mismatch_cost=1, match_score=3
        )

        # Root should contain original sequences
        self.assertEqual(tree.seq1, seq1)
        self.assertEqual(tree.seq2, seq2)

        # First sequence should be split at index 6
        self.assertEqual(tree.left.seq1, "ACGTAC")
        self.assertEqual(tree.right.seq1, "GTACGT")

    def test_base_case_handling(self):
        """Test that base cases are handled correctly."""

        def check_leaf_condition(node):
            """Check that all leaf nodes satisfy the base case condition."""
            if node.left is None and node.right is None:  # This is a leaf node
                # At least one sequence should have length 0 or 1
                self.assertTrue(len(node.seq1) <= 1 or len(node.seq2) <= 1)
            else:
                # Check children recursively
                if node.left:
                    check_leaf_condition(node.left)
                if node.right:
                    check_leaf_condition(node.right)

        # Test with longer sequences to ensure multiple recursive calls
        seq1 = "ACGTACGTACGTACGT"
        seq2 = "TGCATGCATGCA"

        tree = build_hirschberg_tree(
            seq1, seq2,
            del_cost=2, ins_cost=2,
            mismatch_cost=1, match_score=3
        )

        # Check that all leaf nodes satisfy the base case condition
        check_leaf_condition(tree)

    def test_different_scoring_parameters(self):
        """Test with different scoring parameters."""
        seq1 = "ACGT"
        seq2 = "ACTT"

        # Higher match score should prefer matching A and C
        tree1 = build_hirschberg_tree(
            seq1, seq2,
            del_cost=2, ins_cost=2,
            mismatch_cost=1, match_score=5  # High match score
        )

        # Higher deletion cost should prefer mismatches
        tree2 = build_hirschberg_tree(
            seq1, seq2,
            del_cost=10, ins_cost=2,  # High deletion cost
            mismatch_cost=1, match_score=3
        )

        # Trees should be different due to different optimal alignments
        # Just check they're constructed without errors
        self.assertIsNotNone(tree1)
        self.assertIsNotNone(tree2)


if __name__ == "__main__":
    unittest.main()
