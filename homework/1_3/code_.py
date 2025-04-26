from anytree import Node as ATNode, RenderTree

class CallNode:
    def __init__(self, seq1, seq2, left=None, right=None):
        self.seq1 = seq1
        self.seq2 = seq2
        self.left = left
        self.right = right

    def __repr__(self):
        return f"({self.seq1},{self.seq2})"

def score_row(a, b, DEL, INS, MATCH, MISMATCH):

    prev = [j * INS for j in range(len(b) + 1)]
    for i in range(1, len(a) + 1):
        curr = [0] * (len(b) + 1)
        curr[0] = i * DEL
        for j in range(1, len(b) + 1):
            score_diag = prev[j-1] + (MATCH if a[i-1] == b[j-1] else MISMATCH)
            score_del = prev[j] + DEL
            score_ins = curr[j-1] + INS
            curr[j] = max(score_diag, score_del, score_ins)
        prev = curr
    return prev

def nw_align(a, b, DEL, INS, MATCH, MISMATCH):
    m, n = len(a), len(b)
    prev = [j * INS for j in range(n + 1)]
    path = [[None] * (n + 1) for _ in range(m + 1)]

    for i in range(1, m + 1):
        curr = [0] * (n + 1)
        curr[0] = i * DEL
        path[i][0] = (i - 1, 0)  # delete from above
        for j in range(1, n + 1):
            match_mismatch = prev[j-1] + (MATCH if a[i-1] == b[j-1] else MISMATCH)
            delete = prev[j] + DEL
            insert = curr[j-1] + INS
            best = max(match_mismatch, delete, insert)
            curr[j] = best
            if best == match_mismatch:
                path[i][j] = (i-1, j-1)
            elif best == delete:
                path[i][j] = (i-1, j)
            else:
                path[i][j] = (i, j-1)
        prev = curr

    align_a, align_b = "", ""
    i, j = m, n
    while i > 0 or j > 0:
        pi, pj = path[i][j]
        if pi == i-1 and pj == j-1:
            align_a = a[i-1] + align_a
            align_b = b[j-1] + align_b
        elif pi == i-1 and pj == j:
            align_a = a[i-1] + align_a
            align_b = '-' + align_b
        else:  # pi == i and pj == j-1
            align_a = '-' + align_a
            align_b = b[j-1] + align_b
        i, j = pi, pj

    return align_a, align_b

def algo_hirschberg(a, b,DEL, INS,MATCH, MISMATCH):
    node = CallNode(a, b)
    if len(a) == 0:
        return '-' * len(b), b, node
    if len(b) == 0:
        return a, '-' * len(a), node
    if len(a) <= 1 or len(b) <= 1:
        align_a, align_b = nw_align(a, b, DEL, INS, MATCH, MISMATCH)
        return align_a, align_b, node
    mid = len(a) // 2
    left_a = a[:mid]
    right_a = a[mid:]
    left_scores = score_row(left_a, b, DEL, INS, MATCH, MISMATCH)
    right_scores = score_row(right_a[::-1], b[::-1], DEL, INS, MATCH, MISMATCH)[::-1]
    combined = [l + r for l, r in zip(left_scores, right_scores)]
    split_idx = combined.index(max(combined))
    align_left_1, align_left_2, node_left = algo_hirschberg(
        left_a, b[:split_idx],
        DEL, INS, MATCH, MISMATCH
    )
    align_right_1, align_right_2, node_right = algo_hirschberg(
        right_a, b[split_idx:],
        DEL, INS, MATCH, MISMATCH
    )

    node.left = node_left
    node.right = node_right

    return align_left_1 + align_right_1, align_left_2 + align_right_2, node


def tree_terminal(root):
    def to_anytree(n, parent=None):
        label = f"({n.seq1},{n.seq2})"
        at = ATNode(label, parent=parent)
        if n.left:
            to_anytree(n.left, at)
        if n.right:
            to_anytree(n.right, at)
        return at

    at_root = to_anytree(root)
    for pre, _, nd in RenderTree(at_root):
        print(f"{pre}{nd.name}")


# def score_alignment(aln1: str, aln2: str, MATCH: int, MISMATCH: int, inDEL: int) -> int:
#     score = 0
#     for a, b in zip(aln1, aln2):
#         if a == '-' or b == '-':
#             score += inDEL
#         elif a == b:
#             score += MATCH
#         else:
#             score += MISMATCH
#     return score


if __name__ == "__main__":
    a = "AGTACGCA"
    b = "TATGC"
    DEL = -2
    INS = -2
    MATCH = 2
    MISMATCH = -1

    align1, align2, root = algo_hirschberg(a, b,
                                                       DEL, INS,
                                                       MATCH, MISMATCH)
    print("Дерево рекурсивных вызовов алгоритма Хиршберга:")
    tree_terminal(root)

    print("\nВыравнивание:")
    print(align1)
    print(align2)
