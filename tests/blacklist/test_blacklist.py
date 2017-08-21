import pytest

from io import StringIO
from collections import namedtuple

import pandas as pd

from epic.blacklist.compute_poisson import compute_poisson

@pytest.fixture()
def matrix():

    return pd.read_table(StringIO(u"""Chromosome Bin f1.bed f2.bed
chr1 0 1.0 1.0
chr1 200 1.0 1.0
chr1 400 1.0 10.0
chr1 600 1.0 1.0
chr1 800 10.0 500.0
chr1 1000 1.0 2.0
chr1 2000 1.0 2.0
chr1 2200 1.0 600.0
chr1 2400 1.0 2.0"""), sep=" ", header=0, index_col=[0, 1])


MockArgs = namedtuple("MockNamespace",
        ["number_cores", "genome", "keep_duplicates", "window_size",
         "fragment_size", "bonferroni", "effective_genome_fraction", "chromosome_sizes"])


@pytest.fixture()
def effective_genome_size_dict():
    return {"chr1": 2000}


@pytest.fixture()
def mock_args():

    return MockArgs(1, "hg38", True, 200, 150, 0.05, 2000, None)


@pytest.fixture()
def expected_result():
    return pd.read_table(StringIO(u"""Chromosome Bin End
0 chr1   800   999
1 chr1  2200   2399"""), index_col=0, header=0, sep="\s+")


def test_compute_poisson(matrix, mock_args, expected_result):

    result = compute_poisson(matrix, mock_args)

    print(result)
    print(expected_result)

    assert result.equals(expected_result)
