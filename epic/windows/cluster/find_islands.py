import logging
from functools import partial
import sys

import pandas as pd

from joblib import Parallel, delayed

from typing import Iterable
from argparse import Namespace

from epic.src.find_islands import _find_islands_cython

def find_islands(dfs, score_threshold, args):
    # type: (Iterable[pd.DataFrame], float, Namespace) -> Iterable[pd.DataFrame]
    logging.info("Clustering bins into islands.")
    parallel_find_islands = partial(_find_islands, args.window_size,
                                    args.gaps_allowed, score_threshold)
    clustered_islands = Parallel(n_jobs=args.number_cores)(
        delayed(parallel_find_islands)(df) for df in dfs)
    return clustered_islands


def _find_islands(window_size, gaps_allowed, score_threshold, df):
    # type: (int, int, float, pd.DataFrame) -> pd.DataFrame


    if df.empty:
        return df

    # adding one to allowed gaps, because if one bin starts at x and another at x + window_size
    # there is no distance between them
    distance_allowed = window_size * (gaps_allowed + 1)
    chromosome = df.iloc[0].Chromosome

    df = df[df.ChIP > 0]

    df = df.sort_values("Bin")

    # print(chromosome)
    # print(df.head())

    df.loc[:, ["ChIP", "Input"]] = df[["ChIP", "Input"]].astype(int)
    # rowdicts = []
    # for group, gdf in df.groupby(((df.Bin - df.Bin.shift()).abs() > distance_allowed).cumsum()):

    #     start = gdf.Bin.iloc[0]
    #     end = gdf.Bin.iloc[-1] + (window_size - 1)
    #     chip_count = gdf.ChIP.sum()
    #     input_count = gdf.Input.sum()
    #     island_score = gdf.Score.sum()

    #     rowdicts.append({"Start": start, "End": end, "ChIP": chip_count,
    #                      "Input": input_count, "Score": island_score})


    # df = pd.DataFrame.from_dict(rowdicts)[["Start", "End", "ChIP", "Input", "Score"]]
    # print(df.dtypes, file=sys.stderr)
    starts, ends, chip, background, scores = _find_islands_cython(df.Bin.values, df.ChIP.values, df.Input.values, df.Score.values, distance_allowed, score_threshold, window_size)

    df = pd.DataFrame({"Start": starts, "End": ends, "ChIP": chip, "Input": background, "Score": scores})[["Start", "End", "ChIP", "Input", "Score"]]

    # island_starts = df.groupby(((df.Bin - df.Bin.shift()).abs(
    # ) > distance_allowed).cumsum()).apply(lambda island: island.Bin.min())
    # island_ends = df.groupby(((df.Bin - df.Bin.shift()).abs(
    # ) > distance_allowed).cumsum()).apply(
    #     lambda island: island.Bin.max() + (window_size - 1))
    # island_chip_counts = df.groupby(((df.Bin - df.Bin.shift()).abs(
    # ) > distance_allowed).cumsum()).apply(lambda island: island.ChIP.sum())
    # island_input_counts = df.groupby(((df.Bin - df.Bin.shift()).abs(
    # ) > distance_allowed).cumsum()).apply(lambda island: island.Input.sum())
    # island_scores = df.groupby(((df.Bin - df.Bin.shift()).abs(
    # ) > distance_allowed).cumsum()).apply(lambda island: island.Score.sum())

    # df = pd.concat(
    #     [island_starts, island_ends, island_chip_counts, island_input_counts,
    #      island_scores],
    #     axis=1).reset_index(drop=True)

    # df.columns = ["Start", "End", "ChIP", "Input", "Score"]
    df.insert(0, "Chromosome", chromosome)
    # df = df[df["Score"] > score_threshold]

    return df.reset_index(drop=True)
