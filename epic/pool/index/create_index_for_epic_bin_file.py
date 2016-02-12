from os.path import basename
from collections import OrderedDict, defaultdict

import pandas as pd
import numpy as np


def create_index_for_epic_bin_file(epic_bin_file_path, chunk_start_end):
    """Return positions in the file containing specific genomic regions.

    The chunk_start_ends is a list of tuples indicating the start and end of a
    genomic region.

    Keyword Arguments:
    epic_bin_file_path -- file to find file positions for.
    chunk_start_ends   -- the genomic region to find file positons for.
    """

    df = pd.read_table(epic_bin_file_path, index_col=[0, 1])
    df.insert(0, "line_number", np.arange(1, len(df) + 1))
    line_numbers = OrderedDict()

    chromosomes_in_df = df.index.get_level_values(0)
    chromosomes_in_common = sorted(set(chromosomes_in_df).intersection(
        chunk_start_end))

    for chromosome in chromosomes_in_common:
        chromosome_df = df.xs(chromosome, level="Chromosome", drop_level=False)

        chromosome_chunks = chunk_start_end[chromosome]
        bins = pd.Series(
            chromosome_df.index.get_level_values(1),
            index=chromosome_df.line_number)

        chromosome_line_numbers = []
        for start, end in chromosome_chunks:
            bigger_than_start = bins[bins >= start]
            smaller_than_end = bins[bins <= end]

            if not (bigger_than_start.empty or smaller_than_end.empty):
                start_line_nb = bigger_than_start.index[0]
                end_line_nb = smaller_than_end.index[-1]
            else:
                start_line_nb = end_line_nb = 0

            chromosome_line_numbers.append((start_line_nb, end_line_nb))

        line_numbers[chromosome] = chromosome_line_numbers

    return line_numbers
