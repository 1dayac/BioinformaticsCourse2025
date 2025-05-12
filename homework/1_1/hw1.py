def find_hamming_dist(s1, s2):
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def find_closest_substring(s1, s2):
    min_dist = float('inf')
    best_idx = -1
    best_subs = ""

    long_s, short_s = (s1, s2) if len(s1) > len(s2) else (s2, s1)

    for i in range(len(long_s) - len(short_s) + 1):
        curr_subs = long_s[i:i + len(short_s)]
        curr_dist = find_hamming_dist(short_s, curr_subs)

        if curr_dist < min_dist:
            min_dist = curr_dist
            best_idx = i
            best_subs = curr_subs

        if min_dist == 0:
            break

    return best_idx, best_subs, min_dist


def find_levenshtein_dist(s1, s2):
    long_s, short_s = (s1, s2) if len(s1) > len(s2) else (s2, s1)
    prev_row = list(range(len(long_s) + 1))
    curr_row = [0] * (len(long_s) + 1)

    for i in range(1, len(short_s) + 1):
        curr_row[0] = i

        for j in range(1, len(long_s) + 1):
            if short_s[i - 1] == long_s[j - 1]:
                curr_row[j] = prev_row[j - 1]
            else:
                curr_row[j] = 1 + min(
                    prev_row[j],
                    curr_row[j - 1],
                    prev_row[j - 1]
                )
        prev_row, curr_row = curr_row, prev_row
    return prev_row[len(long_s)]
