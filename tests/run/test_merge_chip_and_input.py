import pytest

import pandas as pd
import numpy as np

import logging
from io import StringIO
from joblib import delayed, Parallel


@pytest.fixture
def input_data():
    pass


@pytest.fixture
def expected_result():
    pass


def merge_chip_and_input(windows, nb_cpu):
    """Merge lists of chromosome bin df chromosome-wise.


    Returns a list of dataframes, one per chromosome, with the collective count
    per bin for all files.


    Keyword Arguments:
    windows -- OrderedDict where the keys are files, the values are lists of
               dfs, one per chromosome.
    nb_cpu  -- cores to use
    """

    windows = iter(windows)
    merged = next(windows)

    for chromosome_dfs in windows:
        merged = merge_two_bin_dfs(merged, chromosome_dfs, nb_cpu)

    return merged

# @pytest.mark.unit
# def test_merge_two_bin_files(sample1_dfs, sample2_dfs):
#     """TODO: Need to test that the lists might not have the same/all chromosomes.

#     It might be possible that there are no sig islands on one chromosome in one
#     file, while there are in the others. Solve by taking in dict with chromos
#     instead of list with files?

#     You will probably be asked about a bug due to this some time.
#     """

#     print("Read run epic code. Begin there!\n" * 5)
#     result = merge_chip_and_input([sample2_dfs, sample2_dfs], 1)
#     print(result)
#     assert 1


def merge_two_bin_dfs(sample1_dfs, sample2_dfs, nb_cpu):
    merged_chromosome_dfs = Parallel(n_jobs=nb_cpu)(
        delayed(_merge_two_bin_dfs)(df1, df2)
        for df1, df2 in zip(sample1_dfs, sample2_dfs))

    return merged_chromosome_dfs


def _merge_two_bin_dfs(df1, df2):

    merged_df = df1.merge(df2,
                          how="outer",
                          on=["Chromosome", "Bin"],
                          suffixes=("_x", "_y"))
    merged_df = merged_df.fillna(0)

    merged_df["Count"] = merged_df["Count_x"] + merged_df["Count_y"]

    merged_df = merged_df.drop(["Count_x", "Count_y"], axis=1)

    return merged_df


@pytest.fixture
def sample1_dfs():
    return [pd.read_table(
        StringIO(u"""
Count Chromosome    Bin
1       chrM    400
1       chrM   2600
1       chrM   3600
1       chrM   3800
1       chrM  12800
1       chrM  14200"""),
        sep="\s+",
        header=0), pd.read_table(
            StringIO(u"""Count Chromosome        Bin
1       chrX    2820000
1       chrX    2854800
1       chrX    3001400
1       chrX    3354400
1       chrX    3489400
1       chrX    3560200
1       chrX    4011200
1       chrX    4644600
1       chrX    4653600
1       chrX    4793400
1       chrX    5136800
1       chrX    5572800
1       chrX    5589400
1       chrX    5792000
1       chrX    5961800
1       chrX    6951000
1       chrX    7125800
1       chrX    7199000
1       chrX    7443200
1       chrX    7606000
1       chrX    7627800
1       chrX    8035600
1       chrX    8073600
1       chrX    8367800
1       chrX    9021000
1       chrX    9472400
1       chrX    9620800
1       chrX    9652000
1       chrX    9801000
1       chrX    9953800"""),
            sep="\s+",
            header=0)]


@pytest.fixture
def sample2_dfs():
    return [pd.read_table(
        StringIO(u"""
Count Chromosome    Bin
1       chrM    400
1       chrM   2600
1       chrM   3600
1       chrM   3800
1       chrM  12800
1       chrM  14200"""),
        header=0,
        sep="\s+", ), pd.read_table(
            StringIO(u"""Count Chromosome        Bin
1       chrX    2820000
1       chrX    2854800
1       chrX    3001400
1       chrX    3354400
1       chrX    3489400
1       chrX    3560200
1       chrX    4011200
1       chrX    4644600
1       chrX    4653600
1       chrX    4793400
1       chrX    5136800
1       chrX    5572800
1       chrX    5589400
1       chrX    5792000
1       chrX    5961800
1       chrX    6951000
1       chrX    7125800
1       chrX    7199000
1       chrX    7443200
1       chrX    7606000
1       chrX    7627800
1       chrX    8035600
1       chrX    8073600
1       chrX    8367800
1       chrX    9021000
1       chrX    9472400
1       chrX    9620800
1       chrX    9652000
1       chrX    9801000
1       chrX    9953800"""),
            sep="\s+",
            header=0)]
