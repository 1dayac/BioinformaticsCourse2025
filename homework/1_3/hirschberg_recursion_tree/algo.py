from typing import List, Optional

class Node:
    def __init__(self, A: str, B: str):
        self.A: str = A
        self.B: str = B
        self.children: List['Node'] = []

    def __repr__(self) -> str:
        return f"Node(A={self.A!r}, B={self.B!r}, children={len(self.children)})"


def nw_score(
    A: str,
    B: str,
    del_cost: int,
    ins_cost: int,
    mismatch_cost: int,
    match_score: int
) -> List[int]:
    m, n = len(A), len(B)
    prev = [j * -ins_cost for j in range(n + 1)]
    for i in range(1, m + 1):
        curr = [i * -del_cost] + [0] * n
        for j in range(1, n + 1):
            match_mis = prev[j - 1] + (match_score if A[i-1] == B[j-1] else -mismatch_cost)
            delete = prev[j] - del_cost
            insert = curr[j-1] - ins_cost
            curr[j] = max(match_mis, delete, insert)
        prev = curr
    return prev


def hirschberg_tree(
    A: str,
    B: str,
    del_cost: int,
    ins_cost: int,
    mismatch_cost: int,
    match_score: int
) -> Node:
    node = Node(A, B)
    n, m = len(A), len(B)

    if n == 0 or m == 0 or n == 1 or m == 1:
        return node

    mid = n // 2
    scoreL = nw_score(A[:mid], B, del_cost, ins_cost, mismatch_cost, match_score)
    scoreR = nw_score(A[mid:][::-1], B[::-1], del_cost, ins_cost, mismatch_cost, match_score)

    best_k: int = 0
    best_score: Optional[int] = None
    for i in range(m + 1):
        s = scoreL[i] + scoreR[m - i]
        if best_score is None or s > best_score:
            best_score = s
            best_k = i

    left = hirschberg_tree(A[:mid], B[:best_k], del_cost, ins_cost, mismatch_cost, match_score)
    right = hirschberg_tree(A[mid:], B[best_k:], del_cost, ins_cost, mismatch_cost, match_score)
    node.children = [left, right]

    return node


if __name__ == "__main__":
    A = "AGTACGCA"
    B = "TATGC"
    tree = hirschberg_tree(A, B, del_cost=2, ins_cost=2, mismatch_cost=1, match_score=3)

    def print_tree(node: Node, indent: int = 0) -> None:
        print(" " * indent + f"|-('{node.A}', '{node.B}')")
        for child in node.children:
            print_tree(child, indent + 1)

    print_tree(tree)
