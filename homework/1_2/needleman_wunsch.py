def needleman_wunsch(seq1, seq2, matrix, gap_cost):    
    m = len(seq1)
    n = len(seq2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    traceback = [[''] * (n + 1) for _ in range(m + 1)]
    
    for i in range(1, m + 1):
        dp[i][0] = dp[i - 1][0] + gap_cost
        traceback[i][0] = 'U'
    for j in range(1, n + 1):
        dp[0][j] = dp[0][j - 1] + gap_cost
        traceback[0][j] = 'L'

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            a = seq1[i - 1]
            b = seq2[j - 1]
            
            if (a, b) in matrix:
                score_sub = matrix[(a, b)]
            elif (b, a) in matrix:
                score_sub = matrix[(b, a)]
            else:
                score_sub = -1
            
            diag = dp[i - 1][j - 1] + score_sub  
            up = dp[i - 1][j] + gap_cost           
            left = dp[i][j - 1] + gap_cost         
            
            dp[i][j] = max(diag, up, left)
            if dp[i][j] == diag:
                traceback[i][j] = 'D'
            elif dp[i][j] == up:
                traceback[i][j] = 'U'
            else:
                traceback[i][j] = 'L'
    
    align1, align2 = [], []
    i, j = m, n
    
    while i > 0 or j > 0:
        if i > 0 and j > 0 and traceback[i][j] == 'D':
            align1.append(seq1[i - 1])
            align2.append(seq2[j - 1])
            i -= 1
            j -= 1

        elif i > 0 and traceback[i][j] == 'U':
            align1.append(seq1[i - 1])
            align2.append('-')
            i -= 1

        elif j > 0 and traceback[i][j] == 'L':
            align1.append('-')
            align2.append(seq2[j - 1])
            j -= 1
    
    align1.reverse()
    align2.reverse()
    
    alignment = ("".join(align1), "".join(align2))
    score = dp[m][n]
    return alignment, score
