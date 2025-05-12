from collections import deque

def Hirschberg_tree(a, b, dl = -2, ins = -2, match = 2, mismatch = -1):
    '''Для обеспечения линейности по памяти не делаю рекурсивные вызовы в глубину,
    а использую очередь для обхода в ширину'''

    def nw(s1, s2):
        n, m = len(s1), len(s2)
        prev = [i * dl for i in range(m + 1)]
        for j in range(1, n + 1):
            curr = [j * ins] + [0] * m
            for i in range(1, m + 1):
                if s2[i - 1] == s1[j - 1]:
                    cost_match = prev[i - 1] + match
                else:
                    cost_match = prev[i - 1] + mismatch

                cost_del = prev[i] + dl
                cost_ins = curr[i - 1] + ins

                curr[i] = max(cost_match, cost_del, cost_ins)

            prev = curr
        return prev

    k, len_l, len_r = 0, 0, (max(len(a), len(b)) + 1) * 5
    flg = 0
    s, p  = '', ''
    queue = deque([(a, b, k, len_l, len_r, '_')])
    while queue:
        curr_a, curr_b, kl, len_l, len_r, nap = queue.popleft()
        if kl != k:
            if flg: print(p)
            print(s)
            k, s, flg, p = kl,  '', 1, ''
        n, m = len(curr_a), len(curr_b)
        tmp = (len_r-len_l-m-n-3)//2
        l_cur = len(s)
        s += f'{' ' * (tmp + len_l - l_cur)}({curr_a},{curr_b})'
        if nap == '/':
            n_cur = len(s)-3-len(p)
        else: n_cur =len(s)-1-n-m-len(p)
        p += ' '*n_cur+nap

        if n == 0 or m == 0 or n == 1 or m == 1:
            continue

        mid = n // 2
        score_left = nw(curr_a[:mid], curr_b)
        score_right = nw(curr_a[mid:][::-1], curr_b[::-1])

        j = max(range(m + 1), key=lambda j: score_left[j] + score_right[m - j])

        s1_left, s1_right = curr_a[:mid], curr_a[mid:]
        s2_left, s2_right = curr_b[:j], curr_b[j:]
        kl_new, c_med = kl+1, (len_l+len_r)//2
        queue.append((s1_left, s2_left, kl_new, len_l, c_med, '/'))
        queue.append((s1_right, s2_right, kl_new, c_med, len_r, '\\'))
    if flg: print(p)
    print(s)

    return

a = 'AGTACGCA'
b = 'TATGC'

#a, b = 'CG', 'T'

#a = 'AGTACGCA'
#b = 'TAGAGC'

#a = 'ACGTACGTACGT'
#b = 'AGTACCTACCGT'
#a = 'CCCCCCCCCCCCCCCCCCCCCCCC'
#b = 'TTTTTTTTTTTTTTTTTTTTTTT'

#a = 'ACGTACGTACGTAGTACGCAAGTACGCAACGTACGTACGTAGTACGCAAGTACGCAACGTACGTACGTAGTACGCAAGTACGCAACGTACGTACGTAGTACGCAAGTACGCA'
#b = 'AGTACCTACCGTAGTACCTACCGTAGTAAGTACCTACCGTTAGAGCTAGAGCAGTACGCAAGTACCTACCGTCGTAGTAAGTACCTACCGTTAGAGCTAGAGCAGTATACGA'
dl = -2
ins = -2
match = 2
mismatch = -3

Hirschberg_tree(a, b, dl, ins,match, mismatch)








