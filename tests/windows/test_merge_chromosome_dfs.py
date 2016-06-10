import pytest

from io import StringIO
import pandas as pd
from numpy import allclose, int32

from epic.windows.count.merge_chromosome_dfs import merge_chromosome_dfs


@pytest.mark.integration
def test_merge_chromosome_dfs(plus_df, minus_df, expected_df):

    actual_df = merge_chromosome_dfs((plus_df, minus_df))
    print(actual_df.to_csv(sep=" "), "actual")
    print(expected_df.to_csv(sep=" "), "expected")
    print(actual_df, "actual")
    print(expected_df, "expected")
    print(actual_df.dtypes, "actual")
    print(expected_df.dtypes, "expected")
    assert actual_df.equals(expected_df)


@pytest.fixture
def plus_df():

    return pd.read_table(
        StringIO(u"""
   Count Chromosome       Bin
2       chr1  39036800
1       chr1  73781000
"""),
        sep=r"\s+")


@pytest.fixture
def minus_df():

    return pd.read_table(
        StringIO(u"""
   Count Chromosome       Bin
5       chr1   39036800
1       chr1   90059600
"""),
        sep=r"\s+")


@pytest.fixture
def expected_df():
    return pd.read_table(
        StringIO(u"""
Count   Chromosome Bin
7      chr1       39036800
1      chr1           73781000
1     chr1        90059600"""),
        sep="\s+",
        dtype={"Count": int32,
               "Bin": int32})
