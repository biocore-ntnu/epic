import os
from functools import reduce
import logging

try: # py3
    from math import gcd
except:
    from fractions import gcd

import pandas as pd

from epic.merge.compute_bed_bins import compute_bins, merge_bed_bins


def compute_bin_size(dfs):

    bin_sizes = []
    for df in dfs.values():
        bins = df.head(100000).index.get_level_values("Bin").astype(int)
        bin_size = reduce(gcd, bins)
        bin_sizes.append(bin_size)

    assert len(set(bin_sizes)) == 1, "Matrixes have different bin sizes: " + str(bin_sizes)
    bin_size = bin_sizes.pop()
    logging.info("Bin size: " + str(bin_size))

    return bin_size


def _remove_epic_enriched(dfs):

    new_dfs = {}
    for n, df in dfs.items():
        bad_cols = [c for c in df.columns if "Enriched_" in c]
        df = df.drop(bad_cols, axis=1)
        new_dfs[n] = df

    return new_dfs


def add_new_enriched_bins_matrixes(regions, dfs, bin_size):

    dfs = _remove_epic_enriched(dfs)

    # Create new enriched columns with the region file names

    names = ["Enriched_" + os.path.basename(r) for r in regions]
    regions = [(name, pd.read_table(r, sep="\s+", usecols=[0, 1, 2], header=None,
                                    names=["Chromosome", "Start", "End"])) for name, r in zip(names, regions)]
    bins = [compute_bins(df, bin_size, name) for name, df in regions]
    r = merge_bed_bins(bins)

    # intersect each df with the bins
    # feilen er her:
    new_dfs = {}

    for i, (n, df) in enumerate(dfs.items()):

        df = df.join(r, how="outer").fillna(0)

        if i != 0:
            df = df.drop(r.columns, axis=1)

        new_dfs[n] = df

    return new_dfs
