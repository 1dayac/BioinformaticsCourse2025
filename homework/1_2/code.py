# Task1

def needleman_wunsch(seq1, seq2, matrix, cost):
    n, m = len(seq1), len(seq2)
    m_scor = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        m_scor[i][0] = i * cost
    for j in range(1, m + 1):
        m_scor[0][j] = j * cost
    for i in range(n):
        for j in range(m):
            x, y = seq1[i], seq2[j]
            match = matrix.get((x, y)) or matrix.get((y, x))
            d = m_scor[i][j] + match
            u = m_scor[i][j + 1] + cost
            l = m_scor[i + 1][j] + cost
            m_scor[i + 1][j + 1] = max(d, u, l)

    ans1, ans2 = [], []
    i, j = n, m

    while i > 0 or j > 0:
        a, b = seq1[i - 1], seq2[j - 1]

        match = matrix.get((a, b)) or matrix.get((b, a))
        if i > 0 and m_scor[i][j] == m_scor[i - 1][j] + cost:
            ans1.append(a)
            ans2.append("-")
            i -= 1
        elif m_scor[i][j] == m_scor[i - 1][j - 1] + match:
            ans1.append(a)
            ans2.append(b)
            i -= 1
            j -= 1
        else:
            ans1.append("-")
            ans2.append(b)
            j -= 1

    ans1 = ''.join(ans1)
    ans2 = ''.join(ans2)

    return ans1, ans2, m_scor[n][m]


# Task2
def needleman_wunsch_af(seq1, seq2, matrix, g, h):
    n, m = len(seq1), len(seq2)

    m_scor = [[float('-inf')] * (m + 1) for _ in range(n + 1)]
    m_gap_v = [[float('-inf')] * (m + 1) for _ in range(n + 1)]
    m_gap_h = [[float('-inf')] * (m + 1) for _ in range(n + 1)] 

    m_scor[0][0] = 0
    for i in range(1, n + 1):
        m_gap_v[i][0] = g + (i - 1) * h
        m_scor[i][0] = m_gap_v[i][0]
    for j in range(1, m + 1):
        m_gap_h[0][j] = g + (j - 1) * h
        m_scor[0][j] = m_gap_h[0][j]

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            x, y = seq1[i - 1], seq2[j - 1]
            match = matrix.get((x, y)) or matrix.get((y, x))

            m_gap_v[i][j] = max(m_scor[i - 1][j] + g, m_gap_v[i - 1][j] + h)
            m_gap_h[i][j] = max(m_scor[i][j - 1] + g, m_gap_h[i][j - 1] + h)
            m_scor[i][j] = max(m_scor[i - 1][j - 1] + match, m_gap_v[i][j], m_gap_h[i][j])

    ans1, ans2 = [], []
    i, j = n, m
    curr = 'm_scor'

    while i > 0 or j > 0:
        if curr == 'm_scor':
            if i > 0 and j > 0:
                x, y = seq1[i - 1], seq2[j - 1]
                match = matrix.get((x, y)) or matrix.get((y, x))
                if m_scor[i][j] == m_scor[i - 1][j - 1] + match:
                    ans1.append(x)
                    ans2.append(y)
                    i -= 1
                    j -= 1
                    continue
            if m_scor[i][j] == m_gap_v[i][j]:
                curr = 'gap_v'
            elif m_scor[i][j] == m_gap_h[i][j]:
                curr = 'gap_h'

        elif curr == 'gap_v':
            x = seq1[i - 1]
            ans1.append(x)
            ans2.append("-")
            if m_gap_v[i][j] == m_scor[i - 1][j] + g:
                curr = 'm_scor'
            i -= 1

        elif curr == 'gap_h':
            y = seq2[j - 1]
            ans1.append("-")
            ans2.append(y)
            if m_gap_h[i][j] == m_scor[i][j - 1] + g:
                curr = 'm_scor'
            j -= 1

    return ''.join(reversed(ans1)), ''.join(reversed(ans2)), m_scor[n][m]

