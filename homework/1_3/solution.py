class RecursiveCallNode:
    """Node in the Hirschberg algorithm recursive call tree."""

    def __init__(self, seq1, seq2):
        # Sequences being aligned at this node
        self.seq1 = seq1
        self.seq2 = seq2
        # Child nodes (left and right subproblems)
        self.left = None
        self.right = None

    def __str__(self):
        """String representation of the node."""
        return f"('{self.seq1}', '{self.seq2}')"


def calculate_score_row(seq1, seq2, del_cost, ins_cost, mismatch_cost, match_score):
    """
    Calculate the final row of Needleman-Wunsch scoring matrix using linear space.
    """
    previous_row = [j * -ins_cost for j in range(len(seq2) + 1)]

    for i in range(1, len(seq1) + 1):
        current_row = [i * -del_cost] + [0] * len(seq2)

        for j in range(1, len(seq2) + 1):
            if seq1[i - 1] == seq2[j - 1]:
                diag_score = previous_row[j - 1] + match_score
            else:
                diag_score = previous_row[j - 1] - mismatch_cost

            del_score = previous_row[j] - del_cost
            ins_score = current_row[j - 1] - ins_cost

            current_row[j] = max(diag_score, del_score, ins_score)

        previous_row = current_row

    return previous_row


def build_hirschberg_tree(seq1, seq2, del_cost, ins_cost, mismatch_cost, match_score):
    """
    Build a recursive call tree for the Hirschberg algorithm.
    """

    current_node = RecursiveCallNode(seq1, seq2)

    if not seq1 or not seq2 or len(seq1) == 1 or len(seq2) == 1:
        return current_node

    mid_point = len(seq1) // 2
    left_seq1 = seq1[:mid_point]
    right_seq1 = seq1[mid_point:]

    forward_scores = calculate_score_row(
        left_seq1, seq2, del_cost, ins_cost, mismatch_cost, match_score
    )

    backward_scores = calculate_score_row(
        right_seq1[::-1], seq2[::-1], del_cost, ins_cost, mismatch_cost, match_score
    )

    optimal_split = 0
    best_score = float('-inf')

    for k in range(len(seq2) + 1):
        combined_score = forward_scores[k] + backward_scores[len(seq2) - k]

        if combined_score > best_score:
            best_score = combined_score
            optimal_split = k

    left_seq2 = seq2[:optimal_split]
    right_seq2 = seq2[optimal_split:]

    current_node.left = build_hirschberg_tree(
        left_seq1, left_seq2, del_cost, ins_cost, mismatch_cost, match_score
    )

    current_node.right = build_hirschberg_tree(
        right_seq1, right_seq2, del_cost, ins_cost, mismatch_cost, match_score
    )

    return current_node


def display_recursive_tree(node, indent=0):
    """
    Print the recursive call tree in a readable format.
    """
    print(" " * indent + f"|-- {node}")

    if node.left:
        display_recursive_tree(node.left, indent + 4)
    if node.right:
        display_recursive_tree(node.right, indent + 4)
