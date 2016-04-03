from pyfaidx import Fasta
from collections import defaultdict


def effective_genome_size(fasta, read_length=30):
    """Compute effective genome size for genome."""

    idx = Fasta(fasta, read_ahead=int(10e5))
    read_counts = defaultdict(int)

    genome_length = 0

    for chromosome in idx:

        genome_length += len(chromosome)

        for start in range(0, len(chromosome) - read_length):
            sequence = str(chromosome[start:start + read_length]).upper()

            if not "N" in sequence:
                read_counts[sequence] += 1

    effective_genome_size = len(read_counts) / sum(read_counts.values())
    effective_genome_size_with_n = len(read_counts) / genome_length

    return effective_genome_size, effective_genome_size_with_n
