import logging
import re
from fnmatch import fnmatch
from multiprocessing import Pool
from sys import exit as sys_exit

from epic.config.cache_settings import MEMORY
from joblib import Parallel, delayed

try:
    from functools import lru_cache
except ImportError:
    from functools32 import lru_cache

try:
    from itertools import izip_longest
except ImportError:
    from itertools import zip_longest as izip_longest


def _merge_chip_and_input(chip_df, input_df):

    chip_df_nb_bins = len(chip_df)
    merged_df = chip_df.merge(input_df,
                              how="left",
                              on=["Chromosome", "Bin"],
                              suffixes=[" ChIP", " Input"])
    merged_df = merged_df[["Chromosome", "Bin", "Count ChIP", "Count Input"]]
    merged_df.columns = ["Chromosome", "Bin", "ChIP", "Input"]

    merged_df = merged_df.fillna(0)

    assert len(merged_df) == chip_df_nb_bins

    return merged_df


def merge_chip_and_input(chip_dfs, input_dfs, nb_cpu):

    logging.info("Merging ChIP and Input data.")
    merged_chromosome_dfs = Parallel(n_jobs=nb_cpu)(
        delayed(_merge_chip_and_input)(chip_df, input_df)
        for chip_df, input_df in zip(chip_dfs, input_dfs))
    return merged_chromosome_dfs


def get_total_number_of_reads(dfs):
    return sum([df.Count.sum() for df in dfs])


def merge_same_files(sample1_dfs, sample2_dfs, nb_cpu):

    logging.info("Merging ChIP and Input data.")
    merged_chromosome_dfs = Parallel(n_jobs=nb_cpu)(
        delayed(_merge_same_files)(sample1_df, sample2_df)
        for sample1_df, sample2_df in zip(sample1_dfs, sample2_dfs))

    return merged_chromosome_dfs


def _merge_same_files(sample1_df, sample2_df):

    sample1_df_nb_bins = len(sample1_df)
    merged_df = sample1_df.merge(sample2_df,
                                 how="outer",
                                 on=["Chromosome", "Bin"])
    merged_df = merged_df.fillna(0)

    merged_df["Count"] = merged_df["Count_x"] + merged_df["Count_y"]

    merged_df = merged_df.drop(["Count_x", "Count_y"], axis=1)

    return merged_df


def _print_args(argv_zero, args):
    """Prints the cl-arguments used to stdout.

    (cannot just use sys.argv due to default values being implicit)
    """

    arg_list = ["#", argv_zero]
    for k, v in sorted(args.items()):
        if str(v) == "False":
            continue
        elif str(v) == "True":
            arg_list.append(k)
        else:
            arg_list.append(k + " " + str(v))

    print(" ".join(arg_list))
