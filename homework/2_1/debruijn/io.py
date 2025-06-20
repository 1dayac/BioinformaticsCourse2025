from Bio import SeqIO
from itertools import islice


def read_fastq_reads(file_path, limit=None):
    handle = SeqIO.parse(file_path, "fastq")
    if limit is not None:
        return list(islice(handle, limit))
    return list(handle)
