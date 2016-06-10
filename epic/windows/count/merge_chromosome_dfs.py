import pandas as pd
from numpy import int32


def merge_chromosome_dfs(df_tuple):
    """Merges data from the two strands into strand-agnostic counts."""

    plus_df, minus_df = df_tuple
    if plus_df.empty:
        return minus_df
    if minus_df.empty:
        return plus_df

    index_cols = "Chromosome Bin".split()
    plus_df, minus_df = plus_df.set_index(index_cols), minus_df.set_index(
        index_cols)

    count_column = plus_df.columns[0]

    # first sum the two bins from each strand
    df = pd.concat([plus_df, minus_df], axis=1).fillna(0).sum(axis=1)
    df = df.reset_index().sort_values(by="Bin")

    df.columns = ["Chromosome", "Bin", count_column]

    df = df.sort_values(["Chromosome", "Bin"])
    df[["Bin", count_column]] = df[["Bin", count_column]].astype(int32)
    df = df[[count_column, "Chromosome", "Bin"]]
    return df.reset_index(drop=True)
