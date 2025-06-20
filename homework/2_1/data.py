import random
random.seed(12)

length = 1000
read_len = 150
n = 2000
genome = ''.join(random.choices('ACGT', k=length))

reads = []
max_start = len(genome) - read_len
for _ in range(n):
    start = random.randint(0, max_start)
    reads.append(genome[start:start + read_len])

with open('syn_reads.fastq', 'w') as f:
    for i, seq in enumerate(reads):
        f.write(f"@read{i}\n")
        f.write(f"{seq}\n")
        f.write("+\n")
        f.write(f"{'I' * len(seq)}\n")


















