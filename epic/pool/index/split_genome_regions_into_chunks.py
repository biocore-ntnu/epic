from collections import OrderedDict
import numpy as np


def split_chromosome_df(df, nb_lines_per_chunk):

    if len(df) > nb_lines_per_chunk:
        number_splits = round(len(df) / nb_lines_per_chunk)
    else:
        number_splits = 1

    df_splits = []
    for d in np.array_split(df, number_splits):
        split = d.reset_index(drop=True)
        df_splits.append(split)

    return df_splits


def split_genome_regions_into_chunks(genome_df, nb_lines_per_chunk):
    """Split genome df into smaller parts.

    Keyword Arguments:
    genome_df          -- DataFrame w/cols Chromosome (str), Bin (int)
    nb_lines_per_chunk -- Size of chunks you want to split into.
    """
    chromosome_dfs = split_df_into_chromosome_dfs(genome_df)
    chromosome_df_chunks = OrderedDict()
    for chromosome_df in chromosome_dfs:
        chromosome = chromosome_df.head(1).Chromosome.values[0]
        chromosome_df_chunks[chromosome] = split_chromosome_df(
            chromosome_df, nb_lines_per_chunk)

    return chromosome_df_chunks


def split_df_into_chromosome_dfs(genome_df):
    """Split df by chromosome and return list of chromosome dfs.

    Keyword Arguments: genome_df -- DataFrame where first column is
    "Chromosome" (str)
    """
    chromosome_dfs = []
    genome_df = genome_df.set_index("Chromosome")
    for chromosome in genome_df.index.get_level_values(
            "Chromosome").drop_duplicates():
        chromosome_df = genome_df[genome_df.index.get_level_values(
            "Chromosome") == chromosome].reset_index()
        chromosome_dfs.append(chromosome_df)

    return chromosome_dfs
