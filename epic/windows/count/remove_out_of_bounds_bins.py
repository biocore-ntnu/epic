import pandas as pd

def remove_out_of_bounds_bins(df, chromosome_size):
    # type: (pd.DataFrame, int) -> pd.DataFrame
    """Remove all reads that were shifted outside of the genome endpoints."""

    # The dataframe is empty and contains no bins out of bounds
    if "Bin" not in df:
        return df

    df = df.drop(df[df.Bin > chromosome_size].index)

    return df.drop(df[df.Bin < 0].index)


def remove_bins_with_ends_out_of_bounds(df, chromosome_size,
                                        window_size):
    # type: (pd.DataFrame, int, int) -> pd.DataFrame
    """Remove all reads that were shifted outside of the genome endpoints."""

    # The dataframe is empty and contains no bins out of bounds

    # print(df.head(2))
    # print(chromosome_size)
    # print(window_size)
    out_of_bounds = df[df.index.get_level_values("Bin") + window_size >
                       chromosome_size].index
    # print(len(out_of_bounds))
    df = df.drop(out_of_bounds)

    return df

    # dfms = Parallel(n_jobs=args.number_cores)(
    #     delayed(_create_matrixes)(chromosome, chip, input, islands)
    #     for chromosome in all_chromosomes)
