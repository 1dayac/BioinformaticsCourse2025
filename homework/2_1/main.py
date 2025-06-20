import sys
from debruijn.io import read_fastq_reads
from debruijn.graph import DeBruijnGraph


def main():
    if len(sys.argv) != 3:
        print("Usage: python main.py <fastq_file> <k>")
        sys.exit(1)

    fastq_file = sys.argv[1]
    k = int(sys.argv[2])

    reads = read_fastq_reads(fastq_file)
    dbg = DeBruijnGraph(k)
    dbg.build_from_reads(reads)

    print(f"Number of nodes: {len(dbg.graph)}")
    print(f"Number of edges: {sum(len(v) for v in dbg.graph.values())}")
    print(f"Total unique kmers: {len(dbg.coverage)}")


if __name__ == "__main__":
    main()
