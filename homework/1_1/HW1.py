def hamming_dist(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Строки должны быть одинаковой длины.")
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def closest_substr(s1, s2):
    long, short = (s1, s2) if len(s1) > len(s2) else (s2, s1)

    min_hamming_dist = float('inf')
    index = -1

    for i in range(len(long) - len(short) + 1):
        dist = 0
        for j in range(len(short)):
            if long[i + j] != short[j]:
                dist += 1
            if dist >= min_hamming_dist:
                break

        if dist < min_hamming_dist:
            min_hamming_dist = dist
            index = i

    result = long[index:index + len(short)]

    if len(result) == len(short) and result != short:
        return 0, '', 0
    return index, result, min_hamming_dist


def levenshtein_dist(s1, s2):
    long, short = (s1, s2) if len(s1) > len(s2) else (s2, s1)
    len_long = len(long)
    len_short = len(short)

    prev = list(range(len_long + 1))

    for j in range(1, len_short + 1):
        curr = [0] * (len_long + 1)
        curr[0] = j
        for i in range(1, len_long + 1):
            curr[i] = min(curr[i - 1] + 1,
                          prev[i] + 1,
                          prev[i - 1] + (0 if long[i - 1] == short[j - 1] else 1))
        prev = curr
    return prev[len_long]
