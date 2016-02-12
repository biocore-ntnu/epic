from __future__ import division
from sys import stdout

from collections import OrderedDict
from itertools import chain

import pytest

import pandas as pd
import numpy as np

from io import StringIO
from itertools import zip_longest

from epic.pool.index.split_genome_regions_into_chunks import split_chromosome_df, split_genome_regions_into_chunks


@pytest.fixture
def genome_df(chromosome_df1, chromosome_df2):
    return chromosome_df1.append(chromosome_df2)


@pytest.fixture
def chromosome_df1():

    return pd.read_table(
        StringIO("""Chromosome	Bin
chr1	10000
chr1	10050
chr1	10100
chr1	10150
chr1	10200
chr1	10250
chr1	10300
chr1	10350
chr1	10400
chr1	10450"""),
        header=0,
        sep="\s+")


@pytest.fixture
def chromosome_df2():

    return pd.read_table(
        StringIO("""Chromosome	Bin
chr2	10100
chr2	10250
chr2	10350
chr2	10450
chr2	10600
chr2	11350
chr2	11600
chr2	11650
chr2	11700
chr2	11750"""),
        header=0,
        sep="\s+")


# yapf: disable
@pytest.fixture
def expected_result1(chromosome_df1):
    return [pd.read_table(StringIO("""Chromosome Bin
chr1 10000
chr1 10050
chr1 10100
chr1 10150"""), sep=" "),
    pd.read_table(StringIO("""Chromosome Bin
chr1 10200
chr1 10250
chr1 10300"""), sep=" "),
    pd.read_table(StringIO("""Chromosome Bin
chr1 10350
chr1 10400
chr1 10450"""), sep=" ")]


@pytest.fixture
def expected_result2(chromosome_df2):
    return [pd.read_table(StringIO("""Chromosome Bin
chr2 10100
chr2 10250
chr2 10350
chr2 10450"""), sep=" "),
    pd.read_table(StringIO("""Chromosome Bin
chr2 10600
chr2 11350
chr2 11600"""), sep=" "),
    pd.read_table(StringIO("""Chromosome Bin
chr2 11650
chr2 11700
chr2 11750"""), sep=" ")]
# yapf: enable



def describe_split_genome_regions():
    @pytest.mark.mergepools
    def test_df1(chromosome_df1, expected_result1):

        desired_split_size = 3
        result = split_chromosome_df(chromosome_df1, desired_split_size)

        for d, x in zip_longest(result, expected_result1):
            print(d, "actual")
            print(x, "expected")

        assert all([r.equals(x) for r, x in zip_longest(result,
                                                        expected_result1)])

    @pytest.mark.mergepools
    def test_df2(chromosome_df2, expected_result2):

        desired_split_size = 3
        result = split_chromosome_df(chromosome_df2, desired_split_size)

        for d, x in zip_longest(result, expected_result2):
            print(d, "actual")
            print(x, "expected")

        assert all([r.equals(x) for r, x in zip_longest(result,
                                                        expected_result2)])

    @pytest.mark.mergepools
    def test_split_size_of_one(chromosome_df1):

        desired_split_size = 1
        result = split_chromosome_df(chromosome_df1, desired_split_size)

        assert len(result) == len(chromosome_df1)
        assert all([len(r) == 1 for r in result])


@pytest.fixture
def expected_result_split_genome_df(expected_result1, expected_result2):
    od = OrderedDict()
    od["chr1"] = expected_result1
    od["chr2"] = expected_result2
    return od


@pytest.mark.mergepools
def test_df_into_chunks(genome_df, expected_result_split_genome_df):
    """Check that dict of lists of dfs are equal."""
    result = split_genome_regions_into_chunks(genome_df, 3)

    res = list(result.values())
    exp = list(expected_result_split_genome_df.values())
    print(res)
    print(exp)

    for r, x in zip_longest(chain(*res), chain(*exp)):
        print(r, "actual")
        print(x, "expected")

        assert r.equals(x)
