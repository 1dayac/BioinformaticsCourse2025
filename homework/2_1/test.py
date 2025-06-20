from debruijn.graph import DeBruijnGraph
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random
import time


def test_simple_graph():
    reads = [
        SeqRecord(Seq("ACGTAC")),
        SeqRecord(Seq("CGTACG")),
    ]
    k = 4
    dbg = DeBruijnGraph(k)
    dbg.build_from_reads(reads)

    graph = dbg.get_graph()
    coverage = dbg.get_coverage()

    assert len(graph) > 0
    assert all(isinstance(v, list) for v in graph.values())
    assert all(c > 0 for c in coverage.values())


def test_kmer_counts():
    reads = [SeqRecord(Seq("AAAAG"))]
    k = 3
    dbg = DeBruijnGraph(k)
    dbg.build_from_reads(reads)
    cov = dbg.get_coverage()
    assert cov["AAA"] == 2
    assert cov["AAG"] == 1


def test_random_reads():
    genome = ''.join(random.choices("ACGT", k=1000))
    reads = []
    for _ in range(30):
        start = random.randint(0, 850)
        read = genome[start:start + 150]
        reads.append(
            SeqRecord(Seq(read), id=f"read_{start}", description="", letter_annotations={"phred_quality": [40] * 150}))

    k = 25
    dbg = DeBruijnGraph(k)
    dbg.build_from_reads(reads)
    cov = dbg.get_coverage()
    assert all(len(kmer) == k for kmer in cov)
    assert all(c > 0 for c in cov.values())


def test_compression():
    reads = [SeqRecord(Seq("AAATTTGGGCCC"))]
    k = 3
    dbg = DeBruijnGraph(k)
    dbg.build_from_reads(reads)
    before_edges = sum(len(v) for v in dbg.get_graph().values())
    dbg.compress()
    after_edges = sum(len(v) for v in dbg.get_graph().values())

    assert after_edges < before_edges, f"{after_edges} !< {before_edges}"
    for label in dbg.get_edge_labels().values():
        assert len(label) >= k


def test_remove_tips():
    k = 3
    dbg = DeBruijnGraph(k)
    main_path = "AAATTTGGGCCC"
    tips = ["AAAGG", "GGGTT", "TTTAA"]
    reads = [SeqRecord(Seq(main_path))] + [SeqRecord(Seq(t)) for t in tips]

    dbg.build_from_reads(reads)
    before = sum(len(v) for v in dbg.get_graph().values())
    dbg.remove_tips(quantile=1.0)
    after = sum(len(v) for v in dbg.get_graph().values())
    assert after < before


def test_on_real_trimmed_data():
    from debruijn.io import read_fastq_reads
    import os

    path = "/Users/n.mishanin/Trimmomatic-0.39/trimmed_1.fastq"
    if not os.path.exists(path):
        print("[SKIP] trimmed_1.fastq not found")
        return

    k = 31
    t0 = time.perf_counter()
    reads = read_fastq_reads(path, 10000)
    t1 = time.perf_counter()
    print(f"Loaded {len(reads)} reads in {t1 - t0:.2f} sec")

    dbg = DeBruijnGraph(k)
    dbg.build_from_reads(reads)
    t2 = time.perf_counter()
    print(f"Graph built: {len(dbg.graph)} nodes, {sum(len(v) for v in dbg.graph.values())} edges in {t2 - t1:.2f} sec")

    dbg.compress()
    t3 = time.perf_counter()
    print(
        f"After compression: {len(dbg.graph)} nodes, {sum(len(v) for v in dbg.graph.values())} edges in {t3 - t2:.2f} sec")

    dbg.remove_tips()
    t4 = time.perf_counter()
    print(
        f"After tip removal: {len(dbg.graph)} nodes, {sum(len(v) for v in dbg.graph.values())} edges in {t4 - t3:.2f} sec")

    edge_lens = [len(l) for l in dbg.get_edge_labels().values()]
    if edge_lens:
        print(f"Average edge label length: {sum(edge_lens) / len(edge_lens):.2f}")
    else:
        print("No edge labels remaining.")

    print(f"Total pipeline time: {t4 - t0:.2f} sec")


if __name__ == "__main__":
    test_simple_graph()
    test_kmer_counts()
    test_random_reads()
    test_compression()
    test_remove_tips()
    test_on_real_trimmed_data()
    print("All tests passed")
