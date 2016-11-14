import pandas as pd
from numpy import int32
from typing import Any, Tuple, Sequence

def merge_chromosome_dfs(df_tuple):
    # type: (Tuple[pd.DataFrame, pd.DataFrame]) -> pd.DataFrame
    """Merges data from the two strands into strand-agnostic counts."""

    plus_df, minus_df = df_tuple
    index_cols = "Chromosome Bin".split()
    count_column = plus_df.columns[0]

    if plus_df.empty:
        return return_other(minus_df, count_column, index_cols)
    if minus_df.empty:
        return return_other(plus_df, count_column, index_cols)

    # sum duplicate bins
    # TODO: why are there duplicate bins here in the first place?
    plus_df = plus_df.groupby(index_cols).sum()
    minus_df = minus_df.groupby(index_cols).sum()

    # first sum the two bins from each strand
    df = pd.concat([plus_df, minus_df], axis=1).fillna(0).sum(axis=1)
    df = df.reset_index().sort_values(by="Bin")

    df.columns = ["Chromosome", "Bin", count_column]

    df = df.sort_values(["Chromosome", "Bin"])
    df[["Bin", count_column]] = df[["Bin", count_column]].astype(int32)
    df = df[[count_column, "Chromosome", "Bin"]]
    return df.reset_index(drop=True)


def return_other(df, count_column, index_cols):
    # type: (pd.DataFrame, Any, Sequence[Any]) -> pd.DataFrame

    df[[count_column, "Bin"]] = df[[count_column, "Bin"]].astype(int32)
    df = df.groupby(index_cols).sum().reset_index()
    return df[[count_column, "Chromosome", "Bin"]]

# multiprocessing.pool.RemoteTraceback:
# """
# Traceback (most recent call last):
#   File "/local/home/endrebak/anaconda3/lib/python3.5/site-packages/joblib-0.9.4-py3.5.egg/joblib/parallel.py", line 130, in __call__
#     return self.func(*args, **kwargs)
#   File "/local/home/endrebak/anaconda3/lib/python3.5/site-packages/joblib-0.9.4-py3.5.egg/joblib/parallel.py", line 72, in __call__
#     return [func(*args, **kwargs) for func, args, kwargs in self.items]
#   File "/local/home/endrebak/anaconda3/lib/python3.5/site-packages/joblib-0.9.4-py3.5.egg/joblib/parallel.py", line 72, in <listcomp>
#     return [func(*args, **kwargs) for func, args, kwargs in self.items]
#   File "/local/home/endrebak/anaconda3/lib/python3.5/site-packages/bioepic-0.0.9-py3.5.egg/epic/windows/count/merge_chromosome_dfs.py", line 21, in merge_chromosome_dfs
#     df = pd.concat([plus_df, minus_df], axis=1).fillna(0).sum(axis=1)
#   File "/local/home/endrebak/anaconda3/lib/python3.5/site-packages/pandas/tools/merge.py", line 846, in concat
#     return op.get_result()
#   File "/local/home/endrebak/anaconda3/lib/python3.5/site-packages/pandas/tools/merge.py", line 1031, in get_result
#     indexers[ax] = obj_labels.reindex(new_labels)[1]
#   File "/local/home/endrebak/anaconda3/lib/python3.5/site-packages/pandas/indexes/multi.py", line 1422, in reindex
#     raise Exception("cannot handle a non-unique multi-index!")
# Exception: cannot handle a non-unique multi-index!
