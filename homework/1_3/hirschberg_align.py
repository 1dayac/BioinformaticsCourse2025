from random import choice, randint
from Bio import pairwise2


class Node:
    def __init__(self, seq1, seq2, left=None, right=None):
        self.value = (seq1, seq2)
        self.left = left
        self.right = right

def dynamic_two_lines(seq1, seq2):
    len_seq1 = len(seq1)
    len_seq2 = len(seq2)
    prev_sum_row = [j * INS for j in range(len_seq2 + 1)]
    curr_sum_row = [0] * (len_seq2 + 1)

    for i in range(1, len_seq1 + 1):
        curr_sum_row[0] = prev_sum_row[0] + DEL
        for j in range(1, len_seq2 + 1):
            score_match = MATCH if seq1[i - 1] == seq2[j - 1] else MISMATCH
            score_cmp = prev_sum_row[j - 1] + score_match
            score_del = prev_sum_row[j] + DEL
            score_ins = curr_sum_row[j - 1] + INS
            curr_sum_row[j] = max(score_cmp, score_del, score_ins)

        prev_sum_row = curr_sum_row
        curr_sum_row = [0] * (len_seq2 + 1)

    return prev_sum_row


def hirschberg_global_align(seq1, seq2):
    len_seq1, len_seq2 = len(seq1), len(seq2)
    node = Node(seq1, seq2)

    if len_seq1 == 0:
        return '-' * len(seq2), seq2, node
    elif len_seq2 == 0:
        return seq1, '-' * len(seq1), node
    elif len_seq1 == 1:
        first_ind = seq2.find(seq1)
        if first_ind == -1:
            return seq1 + (len_seq2 - 1) * '-', seq2, node
        else:
            return first_ind * '-' + seq1 + (len_seq2 - first_ind - 1) * '-', seq2, node
    else:
        mid_seq1 = int(len_seq1 / 2)
        prefix = seq1[:mid_seq1]
        suffix = seq1[mid_seq1:]

        final_sum_row_prefix = dynamic_two_lines(prefix, seq2)
        final_sum_row_suffix = dynamic_two_lines(suffix[::-1], seq2[::-1])[::-1]

        sum_row = [l + r for l, r in zip(final_sum_row_prefix, final_sum_row_suffix)]
        mid_b, _ = max(enumerate(sum_row), key=lambda x: x[1])

        align_seq1_left, align_seq2_left, node_left = hirschberg_global_align(prefix, seq2[:mid_b])
        align_seq1_right, align_seq2_right, node_right = hirschberg_global_align(suffix, seq2[mid_b:])

        node.left = node_left
        node.right = node_right

    return align_seq1_left + align_seq1_right, align_seq2_left + align_seq2_right, node


def score(align_s1, align_s2):
    score = 0
    for a, b in zip(align_s1, align_s2):
        if a == b:
            score += MATCH
        else:
            if a == '-':
                score += INS
            elif b == '-':
                score += DEL
            else:
                score += MISMATCH
    return score


def build_tree(root, curr_ind=0, ind=False, delimiter='_'):
    if root is None:
        return [], 0, 0, 0

    line1 = []
    line2 = []

    node_text = f'({root.value[0]},{root.value[1]})'
    new_root_width = gap_size = len(node_text)

    l_box, l_box_width, l_root_start, l_root_end = build_tree(root.left, 2 * curr_ind + 1, ind, delimiter)
    r_box, r_box_width, r_root_start, r_root_end = build_tree(root.right, 2 * curr_ind + 2, ind, delimiter)

    if l_box_width > 0:
        l_root = (l_root_start + l_root_end) // 2 + 1
        line1.append(' ' * (l_root + 1))
        line1.append('_' * (l_box_width - l_root))
        line2.append(' ' * l_root + '/')
        line2.append(' ' * (l_box_width - l_root))
        gap_size += 1
    else:
        line1.append('')
        line2.append('')

    line1.append(node_text)
    line2.append(' ' * new_root_width)

    if r_box_width > 0:
        r_root = (r_root_start + r_root_end) // 2
        line1.append('_' * r_root)
        line1.append(' ' * (r_box_width - r_root))
        line2.append(' ' * r_root + '\\')
        line2.append(' ' * (r_box_width - r_root))
        gap_size += 1
    else:
        line1.append('')
        line2.append('')

    gap = ' ' * gap_size
    new_box = [''.join(line1), ''.join(line2)]

    for i in range(max(len(l_box), len(r_box))):
        l_line = l_box[i] if i < len(l_box) else ' ' * l_box_width
        r_line = r_box[i] if i < len(r_box) else ' ' * r_box_width
        new_box.append(l_line + gap + r_line)

    return new_box, len(new_box[0]), l_box_width + gap_size // 2, l_box_width + gap_size // 2 + new_root_width


def print_tree(root):
    lines, *_ = build_tree(root)
    for line in lines:
        print(line)
        
def generate_amino_seq(min_len=2, max_len=10):
    len = randint(min_len, max_len)
    return ''.join(choice('ATGC') for _ in range(len))


def test():
    n_error = 0
    n_seq = 10
    all_seq = [generate_amino_seq() for _ in range(n_seq)]

    for ind_s1 in range(n_seq):
        for ind_s2 in range(ind_s1 + 1, n_seq):
            align_s1, align_s2, root = hirschberg_global_align(all_seq[ind_s1], all_seq[ind_s2])
            hirschberg_score = score(align_s1, align_s2)

            alignments = pairwise2.align.globalmd(
                all_seq[ind_s1], all_seq[ind_s2],
                MATCH, MISMATCH,
                INS, INS,
                DEL, DEL,
                one_alignment_only=True)

            if alignments:
                pairwise_alignment = alignments[0]
                pairwise_score = pairwise_alignment.score
                pairwise_seq1 = pairwise_alignment.seqA
                pairwise_seq2 = pairwise_alignment.seqB

                if hirschberg_score != pairwise_score:
                    print(f"Исходные последовательности: {all_seq[ind_s1]}, {all_seq[ind_s2]}")
                    print(f"Hirschberg score: {hirschberg_score}")
                    print(f"{align_s1}\n{align_s2}")
                    print(f"Bio.pairwise2:\n{pairwise_seq1}\n{pairwise_seq2}")
                    print(f"Bio.pairwise2 score: {pairwise_score}")
                    print('_____________')
                    n_error += 1
    return n_error


# INS — штраф за вставку (gap в первой последовательности seq1)
INS = -2
# DEL — штраф за удаление (gap во второй последовательности seq2)
DEL = -2
MATCH = 2
MISMATCH = -1

print('Количество неверных выравниваний:', test())
print('Пример дерева рекурсивных вызовов алгоритма Хиршберга:')
s1 = generate_amino_seq()
s2 = generate_amino_seq()

align_s1, align_s2, root = hirschberg_global_align(s1, s2)
print_tree(root)
