import pandas as pd


def merge_chromosome_dfs(df_tuple):

    plus_df, minus_df = df_tuple
    if plus_df.empty:
        return minus_df
    if minus_df.empty:
        return plus_df

    df = pd.concat([plus_df, minus_df]).sort_values(by="Bin").reset_index(
        drop=True)

    df["Count"] = df.groupby(["Bin"])["Count"].transform("sum")

    # Only need one entry per bin
    df = df.drop_duplicates(["Chromosome", "Bin"])

    df[["Bin", "Count"]] = df[["Bin", "Count"]]

    return df.reset_index(drop=True)
