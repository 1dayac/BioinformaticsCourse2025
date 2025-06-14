#old HW file, I only moved it to proper folder.

# from Bio.SubsMat import MatrixInfo as matlist
# from Bio import Seq


#basis from
def main(A, B, compare, gap_penalty):
    def get_score(a, b, matrix):
        """Получает оценку замены из матрицы"""
        return matrix.get((a, b), matrix.get((b, a), -1))  # Пробуем оба варианта
    def sequence_alignment(A, B, matrix, gap_payment): # Needleman Wunsch
        m, n = len(A), len(B)
        F = [[0] * (n + 1) for _ in range(m + 1)]
        prev = [[0] * (n + 1) for _ in range(m + 1)]
        for i in range(m + 1):
            F[i][0] = gap_payment * i
        for j in range(n + 1):
            F[0][j] = gap_payment * j
        
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                match = F[i-1][j-1] + get_score(A[i-1], B[j-1], matrix)  #matrix
                delete = F[i-1][j] + gap_payment #gap in top B
                insert = F[i][j-1] + gap_payment #gap in column A
                options = [ delete, insert, match]
                F[i][j] = max(options)
                prev[i][j] = options.index(F[i][j])
                

        
        return F[m][n], F, prev

    #how to restore? backtrack.

    def restore_path(A, B, prev):
        i, j = len(A), len(B)
        path = [(i, j)]
        
        aligned_A = []
        aligned_B = []

        while i > 0 or j > 0:
            if i == 0:  # Only move left (gap in A)
                j -= 1
            elif j == 0:  # Only move up (gap in B)
                i -= 1
            else:
                direction = prev[i][j]
                if direction == 0:  # Came from top 
                    i -= 1
                elif direction == 1:  # Came from left 
                    j -= 1
                else:  #  (match/mismatch)
                    i -= 1
                    j -= 1
            path.append((i, j))
        
        path.reverse()
        print("Alignment path:", path)
        a_prev, b_prev = path[0]
        
        for i in range(1,len(path)):
            a,b = path[i]
            fst = (a>a_prev)
            snd = (b> b_prev)
            a_prev, b_prev = a, b
            aligned_A.append(A[a-1] if fst else "_")
            aligned_B.append(B[b-1] if snd else "_")
        # print(*aligned_A)
        # print(*aligned_B)
        return aligned_A, aligned_B



    # A = "ACGTA"
    # B = "ATA"
    # A = "AAA"
    # B="BBBBBBAAA"
    # gap_penalty = -2

    score, dp_table, prev = sequence_alignment(A, B, compare, gap_penalty)
    # print("Optimal alignment score:", score)
    # print("DP table:")
    # for row in dp_table:
    #     print(row)

    aligned_A, aligned_B = restore_path(A,B,prev)
    return aligned_A, aligned_B, score
