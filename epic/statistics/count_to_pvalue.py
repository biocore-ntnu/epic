from operator import neg as minus
from math import log
from epic.statistics.compute_poisson import _poisson
import pandas as pd
from functools import partial
from joblib import Parallel, delayed
import logging
from typing import List, Sequence

def count_to_pvalue(merged_dfs, island_enriched_threshold,
                    average_window_readcount, nb_cpu):
    # type: (Sequence[pd.DataFrame], int, float, int) -> List[pd.DataFrame]

    logging.info("Giving bins poisson score.")
    parallel_count_to_pvalue = partial(_count_to_pvalue, island_enriched_threshold,
                             average_window_readcount)
    return Parallel(n_jobs=nb_cpu)(delayed(parallel_count_to_pvalue)(df)
                                   for df in merged_dfs)


def _count_to_pvalue(island_enriched_threshold,
                     average_window_readcount, df):
    # type: (int, float, pd.DataFrame) -> pd.DataFrame

    df = df.loc[df["ChIP"] >= island_enriched_threshold]
    scores = df["ChIP"].apply(
        lambda count: _poisson(count, average_window_readcount) + 0.000001)
    df.insert(len(df.columns), "Score", scores)

    pd.options.mode.chained_assignment = None

    df["Score"] = df["Score"].apply(log)
    df["Score"] = df["Score"].apply(minus)

    pd.options.mode.chained_assignment = "warn"

    return df
