def hamming_distance(aligned_seq1, aligned_seq2):
    return sum(ch1 != ch2 for ch1, ch2 in zip(aligned_seq1, aligned_seq2))


def find_closest_substring(long_seq, short_seq):
    best_distance = float('inf')
    best_index = None
    best_substring = None
    window_length = len(short_seq)

    for i in range(len(long_seq) - window_length + 1):
        current_substring = long_seq[i:i + window_length]
        distance = hamming_distance(current_substring, short_seq)
        if distance < best_distance:
            best_distance = distance
            best_index = i
            best_substring = current_substring
    return best_index, best_substring, best_distance


def levenshtein_distance(s, t):
    if len(s) > len(t):
        s, t = t, s
    m = len(s)
    n = len(t)

    previous_row = list(range(m + 1))

    for j in range(1, n + 1):
        current_row = [j] + [0] * m
        for i in range(1, m + 1):
            cost = 0 if s[i - 1] == t[j - 1] else 1
            current_row[i] = min(
                current_row[i - 1] + 1,
                previous_row[i] + 1,
                previous_row[i - 1] + cost
            )
        previous_row = current_row
    return previous_row[m]
