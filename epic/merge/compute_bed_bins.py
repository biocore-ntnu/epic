import pandas as pd


def compute_bins(df, bin_size, name):
    bins = []

    for chromosome, start, end in zip(df.iloc[:, 0], df.iloc[:, 1], df.iloc[:, 2]):
        start, end = int(start), int(end)
        for chromosome_bin in range(start - (start % bin_size), end - (end % bin_size) + bin_size, bin_size):
            d = {"Chromosome": chromosome, "Bin": chromosome_bin}
            d[name] = 1
            bins.append(d)

    bins = pd.DataFrame.from_dict(bins).set_index(["Chromosome", "Bin"])

    return bins


def merge_bed_bins(dfs):
    from functools import reduce

    df = reduce(lambda l, r: l.join(r, how="outer"), dfs).fillna(0)
    df = df[~df.index.duplicated()]
    # s.name = "Enriched"

    return df
