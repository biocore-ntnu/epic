import pytest

from io import StringIO

import pandas as pd

from epic.merge.compute_bed_bins import compute_bins, merge_bed_bins

@pytest.fixture
def bed_dfs():

    c1 = """chr1    887771  887796  U0      0       -
chr1    994660  994685  U0      0       -
chr1    1041102 1041127 U0      0       +
chr1    1770383 1770408 U0      0       -
chr1    1995141 1995166 U0      0       -
chr1    2063984 2064009 U0      0       -
chr1    2129359 2129384 U0      0       +
chr1    2239108 2239133 U0      0       +
chr1    2318805 2318830 U0      0       +
chr1    3006132 3006157 U0      0       -"""

    c2 = """chr1    1325303 1325328 U0      0       -
chr1    1541598 1541623 U0      0       +
chr1    1599121 1599146 U0      0       +
chr1    1820285 1820310 U0      0       -
chr1    2448322 2448347 U0      0       -
chr1    3046141 3046166 U0      0       -
chr1    3437168 3437193 U0      0       -
chr1    3504032 3504057 U0      0       +
chr1    3637087 3637112 U0      0       -
chr1    3681903 3681928 U0      0       -"""

    return [pd.read_table(StringIO(c), sep="\s+", header=None) for c in [c1, c2]]


@pytest.fixture()
def expected_result():

    c = """Bin Chromosome
887600       chr1
994600       chr1
1041000       chr1
1770200       chr1
1770400       chr1
1995000       chr1
2063800       chr1
2064000       chr1
2129200       chr1
2239000       chr1
2318800       chr1
3006000       chr1"""

    return pd.read_table(StringIO(c), sep="\s+", index_col=[1, 0])

@pytest.fixture()
def expected_result2():
    c = """Bin Chromosome
1325200       chr1
1541400       chr1
1541600       chr1
1599000       chr1
1820200       chr1
2448200       chr1
3046000       chr1
3437000       chr1
3504000       chr1
3637000       chr1
3681800       chr1"""

    return pd.read_table(StringIO(c), sep="\s+", index_col=[1, 0])

@pytest.fixture
def matrixes_to_merge(expected_result, expected_result2):

    df1 = expected_result
    df2 = expected_result2

    df1.insert(0, "Enriched_f1", 1)
    df2.insert(0, "Enriched_f2", 1)

    return df1, df2

def test_compute_bins(bed_dfs, expected_result):

    result = compute_bins(bed_dfs[0], 200)
    print(result)
    print(expected_result)

    assert result.equals(expected_result)



@pytest.fixture()
def expected_result_merge_bed_bins():
    c = """Chromosome Bin Enriched_f1 Enriched_f2
chr1 887600 1.0 0.0
chr1 994600 1.0 0.0
chr1 1041000 1.0 0.0
chr1 1325200 0.0 1.0
chr1 1541400 0.0 1.0
chr1 1541600 0.0 1.0
chr1 1599000 0.0 1.0
chr1 1770200 1.0 0.0
chr1 1770400 1.0 0.0
chr1 1820200 0.0 1.0
chr1 1995000 1.0 0.0
chr1 2063800 1.0 0.0
chr1 2064000 1.0 0.0
chr1 2129200 1.0 0.0
chr1 2239000 1.0 0.0
chr1 2318800 1.0 0.0
chr1 2448200 0.0 1.0
chr1 3006000 1.0 0.0
chr1 3046000 0.0 1.0
chr1 3437000 0.0 1.0
chr1 3504000 0.0 1.0
chr1 3637000 0.0 1.0
chr1 3681800 0.0 1.0"""

    return pd.read_table(StringIO(c), sep=" ", header=0, index_col=[0, 1])


def test_merge_bed_bins(matrixes_to_merge, expected_result_merge_bed_bins):

    df = merge_bed_bins(matrixes_to_merge)

    print(df.head(2))
    print(expected_result_merge_bed_bins.head(2))

    assert df.equals(expected_result_merge_bed_bins)
