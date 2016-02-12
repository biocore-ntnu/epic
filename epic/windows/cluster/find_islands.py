import logging
from functools import partial

import pandas as pd

from joblib import Parallel, delayed


def find_islands(dfs, window_size, allowed_gaps, score_threshold, nb_cpu):
    logging.info("Merging bins into islands.")
    parallel_find_islands = partial(_find_islands, window_size, allowed_gaps,
                                    score_threshold)
    clustered_islands = Parallel(n_jobs=nb_cpu)(
        delayed(parallel_find_islands)(df) for df in dfs)
    return clustered_islands


def _find_islands(window_size, allowed_gaps, score_threshold, df):

    if df.empty:
        return df

    distance_allowed = window_size * allowed_gaps
    chromosome = df.iloc[0].Chromosome

    # Code below could probably be optimized....
    # Idea: use one function to do everything
    island_starts = df.groupby(((df.Bin - df.Bin.shift()).abs(
    ) > distance_allowed).cumsum()).apply(lambda island: island.Bin.min())
    island_ends = df.groupby(((df.Bin - df.Bin.shift()).abs(
    ) > distance_allowed).cumsum()).apply(
        lambda island: island.Bin.max() + (window_size - 1))
    island_chip_counts = df.groupby(((df.Bin - df.Bin.shift()).abs(
    ) > distance_allowed).cumsum()).apply(lambda island: island.ChIP.sum())
    island_input_counts = df.groupby(((df.Bin - df.Bin.shift()).abs(
    ) > distance_allowed).cumsum()).apply(lambda island: island.Input.sum())
    island_scores = df.groupby(((df.Bin - df.Bin.shift()).abs(
    ) > distance_allowed).cumsum()).apply(lambda island: island.Score.sum())
    logging.debug("Concatenating for " + chromosome)
    df = pd.concat(
        [island_starts, island_ends, island_chip_counts, island_input_counts,
         island_scores],
        axis=1).reset_index(drop=True)
    logging.debug("Concatenation done for " + chromosome)

    df.columns = ["Start", "End", "ChIP", "Input", "Score"]
    df.insert(0, "Chromosome", chromosome)
    df = df[df["Score"] > score_threshold]

    return df.reset_index(drop=True)
