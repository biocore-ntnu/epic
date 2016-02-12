from collections import OrderedDict


def find_chunks_start_and_ends(genome_chunks):
    """Return the starts and ends of the given genome chunks.

    Keyword Arguments:
    genome_chunks -- OrderedDict {chromosome: chunk_list}
    """

    chunks = OrderedDict()
    for chromosome in genome_chunks:

        chunk_starts_ends = []

        for chunk in genome_chunks[chromosome]:
            start, end = chunk.head(1), chunk.tail(1)

            start_bin = start.Bin.values[0]
            end_bin = end.Bin.values[0]

            chunk_starts_ends.append((start_bin, end_bin))

        chunks[chromosome] = chunk_starts_ends

    return chunks
