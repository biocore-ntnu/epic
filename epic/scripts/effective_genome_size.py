import time
from collections import defaultdict
from sys import stderr

from pyfaidx import Fasta

from joblib import Parallel, delayed


def effective_genome_size(fasta, read_length, nb_cores):
    """Compute effective genome size for genome."""

    idx = Fasta(fasta, read_ahead=int(10e5))

    genome_length = sum([len(c) for c in idx])

    results = Parallel(n_jobs=nb_cores)(
        delayed(compute_number_effective_chromosome_reads)(
            chromosome, read_length) for chromosome in idx)

    read_counts = sum(results)
    # total_reads = sum([r[1] for r in results])

    # effective_genome_size = read_counts / total_reads
    effective_genome_size_with_n = read_counts / genome_length

    return effective_genome_size, effective_genome_size_with_n


def compute_number_effective_chromosome_reads(chromosome, read_length):

    read_counts = defaultdict(int)

    for start in range(0, len(chromosome) - read_length):
        sequence = str(chromosome[start:start + read_length]).upper()

        if not "N" in sequence:
            read_counts[sequence] += 1

    nb_reads = len(read_counts)
    print(time.localtime.tm_hour, time.localtime.tm_min, file=stderr)
    print(chromosome + ":", nb_reads, file=stderr)
    return len(nb_reads)  #, sum(read_counts.values())

# def chromosome_chunks(chromosome, chunk_length, read_length):
#     """Split chromosome into chunks that overlap by read_length."""

#     for i in range(0, chunk_length):
#         yield str(chromosome[i:i + chunk_length - read_length])

#     for i in range(chunk_length, len(chromosome) - 1):
#         yield str(chromosome[i:i + chunk_length - read_length])

# for i in range(len(chromosome) - read_length, len(chromosome)):
#     yield str(chromosome[i:i + chunk_length])
