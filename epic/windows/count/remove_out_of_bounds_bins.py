def remove_out_of_bounds_bins(df, chromosome_size):
    """Remove all reads that were shifted outside of the genome endpoints."""

    # The dataframe is empty and contains no bins out of bounds
    if "Bin" not in df:
        return df

    df = df.drop(df[df.Bin > chromosome_size].index)

    return df.drop(df[df.Bin < 0].index)
