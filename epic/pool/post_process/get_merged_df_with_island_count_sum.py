import pandas as pd


def get_merged_df_with_island_count_sum(df, index_cols=[0, 1]):
    """Count nb islands in each row and add sum as column.

    Keyword Arguments:
    df --
    """

    island_cols = [c for c in df.columns if c.endswith("_Island")]
    island_df = df[island_cols]

    new_df = df.drop(island_cols, axis=1)

    nb_islands_per_row = island_df.sum(axis=1)

    new_df.insert(0, "Islands", nb_islands_per_row)

    return new_df
