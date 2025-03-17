
#Task 1

def hamming_distance(s1, s2):
   if (len(s1) != len(s2)):
       return 'Error. Different line lengths.'
   counter = 0
   
   for i in range(len(s1)):
       if s1[i] != s2[i]:
           counter += 1
   return counter


#Task 2

def substring(s1, s2):
    n_s1, n_s2 = len(s1), len(s2)
    if n_s1 < n_s2:
        long_s, short_s = s2, s1
        n_long_s, n_short_s = n_s2, n_s1
    else:
        long_s, short_s = s1, s2
        n_long_s, n_short_s = n_s1, n_s2

    h_d, pos, answer = float('inf'), float('inf'), '' 

    for i in range(n_long_s - n_short_s + 1):
        sub_s = long_s[i:i + n_short_s]
        c_d = hamming_distance(sub_s, short_s)
        if c_d < h_d:
            h_d, pos, answer = c_d, i, sub_s    
    return pos, answer, h_d


# Task 3

def levenshtein_distance(s1, s2):
    n, m = len(s1), len(s2)
    if n < m:
        s1, s2 = s2, s1
        n, m = m, n
    curr = [0] * (m + 1)
    prew = list(range(m + 1))
    for i in range(1, n + 1):
        curr[0] = i
        for j in range(1, m + 1):
            cost = 0 if s1[i - 1] == s2[j - 1] else 1
            D, I, C = prew[j] + 1,prew[j] + 1,prew[j - 1] + cost
            curr[j] = min(D, I, C)
        prew, curr = curr, prew
    return prew[m]

# s1 = 'dsds'
# s2 = 'ekjfejfjefie' 
# s1 = ''
# s2 = 'fs' 
# print(levenshtein_distance(s1, s2))