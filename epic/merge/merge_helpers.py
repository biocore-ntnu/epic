import os
from functools import reduce
import logging
from collections import OrderedDict

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

    new_dfs = OrderedDict()
    for n, df in dfs.items():
        bad_cols = [c for c in df.columns if "Enriched_" in c]
        df = df.drop(bad_cols, axis=1)
        new_dfs[n] = df

    return new_dfs

def region_files_to_bins(region_files, names, bin_size):

    region_files = [(name, pd.read_table(r, sep="\s+", usecols=[0, 1, 2], header=None,
                                    names=["Chromosome", "Start", "End"])) for name, r in zip(names, region_files)]
    bins = [compute_bins(df, bin_size, name) for name, df in region_files]
    regions = merge_bed_bins(bins)

    return regions


def add_new_enriched_bins_matrixes(region_files, dfs, bin_size):

    """Add enriched bins based on bed files.

    There is no way to find the correspondence between region file and matrix
    file, but it does not matter."""

    dfs = _remove_epic_enriched(dfs)

    names = ["Enriched_" + os.path.basename(r) for r in region_files]

    regions = region_files_to_bins(region_files, names, bin_size)

    new_dfs = OrderedDict()

    assert len(regions.columns) == len(dfs)
    for region, (n, df) in zip(regions, dfs.items()):

        region_col = regions[region]

        df = df.join(region_col, how="outer").fillna(0)

        new_dfs[n] = df

    return new_dfs
