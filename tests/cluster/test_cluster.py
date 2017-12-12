import pytest

from io import StringIO

import pandas as pd

from epic.cluster.cluster import trunks_flanks_valleys


@pytest.fixture
def merged_matrix():

    c = u"""Chromosome Bin TotalEnriched A B
chr1 0 1.0 1.0 0.0
chr1 200 1.0 1.0 0.0
chr1 400 1.0 0.0 1.0
chr1 800 1.0 0.0 1.0
chr1 1000 1.0 1.0 0.0
chr1 1200 1.0 0.0 0.0
chr1 1400 1.0 1.0 0.0
chr1 1800 1.0 1.0 0.0
chr1 2000 1.0 1.0 0.0"""

    return pd.read_table(StringIO(c), sep=" ", header=0)


@pytest.fixture
def expected_result():

    c = u"""Chromosome IslandID RegionID Kind End Start MeanEnriched MinEnriched MaxEnrichedCluster A B
0 chr1 0 0:599 trunk 599 0 1.0 1.0 1.0 2.0 1.0
0 chr1 1 800:1599 trunk 1599 800 1.0 1.0 1.0 2.0 1.0
0 chr1 2 1800:2199 trunk 2199 1800 1.0 1.0 1.0 2.0 0.0"""

    return pd.read_table(StringIO(c), sep="\s+", header=0, index_col=0)


def test_trunks_flanks_valleys(merged_matrix, expected_result):

    result = trunks_flanks_valleys(merged_matrix)

    print(result.to_csv(sep=" "))
    print(result.dtypes)
    print(expected_result.to_csv(sep=" "))
    print(expected_result.dtypes)

    assert result.equals(expected_result)


@pytest.fixture
def merged_matrix2():

    c = u"""Chromosome Bin TotalEnriched A B
chr1 0 1.0 1.0 0.0
chr1 200 5.0 1.0 0.0
chr1 400 4.0 0.0 1.0
chr1 800 1.0 0.0 1.0
chr1 1000 1.0 1.0 0.0
chr1 1200 7.0 0.0 0.0
chr1 1400 2.0 1.0 0.0
chr1 1600 2.0 1.0 0.0
chr1 1800 6.0 1.0 0.0
chr1 2000 1.0 1.0 0.0"""

    return pd.read_table(StringIO(c), sep=" ", header=0)


@pytest.fixture
def expected_result2():

    c = u"""Chromosome IslandID RegionID Kind End Start MeanEnriched MinEnriched MaxEnrichedCluster A B
0 chr1 0 0:199 flank 199 0 1.0 1.0 5.0 1.0 0.0
0 chr1 0 200:599 trunk 599 200 4.5 4.0 5.0 1.0 1.0
0 chr1 1 800:1199 flank 1199 800 1.0 1.0 7.0 1.0 1.0
0 chr1 1 1200:1399 trunk 1399 1200 7.0 7.0 7.0 0.0 0.0
0 chr1 1 1400:1799 valley 1799 1400 2.0 2.0 7.0 2.0 0.0
0 chr1 1 1800:1999 trunk 1999 1800 6.0 6.0 7.0 1.0 0.0
0 chr1 1 2000:2199 flank 2199 2000 1.0 1.0 7.0 1.0 0.0"""

    return pd.read_table(StringIO(c), sep="\s+", header=0, index_col=0)


def test_trunks_flanks_valleys2(merged_matrix2, expected_result2):

    result = trunks_flanks_valleys(merged_matrix2)

    print(result.to_csv(sep=" "))
    print(result.dtypes)
    print(result.index)
    print(expected_result2.to_csv(sep=" "))
    print(expected_result2.dtypes)
    print(expected_result2.index)

    assert result.equals(expected_result2)


# @pytest.fixture
# def merged_matrix3():

#     c = u"""Chromosome Bin TotalEnriched A B
# chr1 0 1.0 1.0 0.0
# chr1 200 5.0 1.0 0.0
# chr1 400 4.0 0.0 1.0
# chr1 800 1.0 0.0 1.0
# chr1 1000 1.0 1.0 0.0
# chr1 1200 7.0 0.0 0.0
# chr1 1400 2.0 1.0 0.0
# chr1 1600 2.0 1.0 0.0
# chr1 1800 6.0 1.0 0.0
# chr1 2000 1.0 1.0 0.0"""

#     return pd.read_table(StringIO(c), sep=" ", header=0)


# @pytest.fixture
# def expected_result3():

#     c = u"""Chromosome IslandID RegionID Kind End Start MeanEnriched MinEnriched MaxEnrichedCluster A B
# 0 chr1 0 0:599 trunk 599 0 3.3333333333333335 1.0 5.0 2.0 1.0
# 0 chr1 1 800:2199 trunk 2199 800 2.857142857142857 1.0 7.0 5.0 1.0"""

#     return pd.read_table(StringIO(c), sep="\s+", header=0, index_col=0)


# def test_trunks_flanks_valleys3(merged_matrix3, expected_result3):

#     result = trunks_flanks_valleys(merged_matrix3, trunk_diff=6)

#     print(result.to_csv(sep=" "))
#     print(result.dtypes)
#     print(result.index)
#     print(expected_result3.to_csv(sep=" "))
#     print(expected_result3.dtypes)
#     print(expected_result3.index)

#     assert result.equals(expected_result3)



@pytest.fixture
def merged_matrix4():

    c = u"""Chromosome Bin TotalEnriched A B
chr1 0 1.0 1.0 0.0
chr1 200 5.0 1.0 0.0
chr1 400 4.0 0.0 1.0
chr1 800 1.0 0.0 1.0
chr1 1000 1.0 1.0 0.0
chr1 1200 7.0 0.0 0.0
chr1 1400 2.0 1.0 0.0
chr1 1600 2.0 1.0 0.0
chr1 1800 6.0 1.0 0.0
chr1 2000 1.0 1.0 0.0"""

    return pd.read_table(StringIO(c), sep=" ", header=0)


@pytest.fixture
def expected_result4():

    c = u"""Chromosome IslandID RegionID Kind End Start MeanEnriched MinEnriched MaxEnrichedCluster A B
0 chr1 0 0:199 flank 199 0 1.0 1.0 7.0 1.0 0.0
0 chr1 0 200:599 trunk 599 200 4.5 4.0 7.0 1.0 1.0
0 chr1 0 800:1199 valley 1199 800 1.0 1.0 7.0 1.0 1.0
0 chr1 0 1200:1399 trunk 1399 1200 7.0 7.0 7.0 0.0 0.0
0 chr1 0 1400:1799 valley 1799 1400 2.0 2.0 7.0 2.0 0.0
0 chr1 0 1800:1999 trunk 1999 1800 6.0 6.0 7.0 1.0 0.0
0 chr1 0 2000:2199 flank 2199 2000 1.0 1.0 7.0 1.0 0.0"""
    return pd.read_table(StringIO(c), sep="\s+", header=0, index_col=0)


def test_trunks_flanks_valleys4(merged_matrix4, expected_result4):

    result = trunks_flanks_valleys(merged_matrix4, trunk_diff=3, distance_allowed=400)

    print(result.to_csv(sep=" "))
    print(expected_result4.to_csv(sep=" "))

    assert result.equals(expected_result4)


@pytest.fixture
def merged_matrix5():

    c = u"""Chromosome Bin TotalEnriched A B
chr1 0 1.0 1.0 0.0
chr1 200 5.0 1.0 0.0
chr1 400 4.0 0.0 1.0
chr1 800 1.0 0.0 1.0
chr1 1000 1.0 1.0 0.0
chr1 1200 7.0 0.0 0.0
chr1 1400 2.0 1.0 0.0
chr1 1600 2.0 1.0 0.0
chr1 1800 6.0 1.0 0.0
chr1 2000 1.0 1.0 0.0"""

    return pd.read_table(StringIO(c), sep=" ", header=0)


@pytest.fixture
def expected_result5():

    c = u"""Chromosome IslandID RegionID Kind End Start MeanEnriched MinEnriched MaxEnrichedCluster A B
0 chr1 0 0:199 trunk 199 0 1.0 1.0 1.0 1.0 0.0
0 chr1 1 200:399 trunk 399 200 5.0 5.0 5.0 1.0 0.0
0 chr1 2 400:599 trunk 599 400 4.0 4.0 4.0 0.0 1.0
0 chr1 3 800:999 trunk 999 800 1.0 1.0 1.0 0.0 1.0
0 chr1 4 1000:1199 trunk 1199 1000 1.0 1.0 1.0 1.0 0.0
0 chr1 5 1200:1399 trunk 1399 1200 7.0 7.0 7.0 0.0 0.0
0 chr1 6 1400:1599 trunk 1599 1400 2.0 2.0 2.0 1.0 0.0
0 chr1 7 1600:1799 trunk 1799 1600 2.0 2.0 2.0 1.0 0.0
0 chr1 8 1800:1999 trunk 1999 1800 6.0 6.0 6.0 1.0 0.0
0 chr1 9 2000:2199 trunk 2199 2000 1.0 1.0 1.0 1.0 0.0"""

    return pd.read_table(StringIO(c), sep="\s+", header=0, index_col=0)


def test_trunks_flanks_valleys5(merged_matrix5, expected_result5):

    result = trunks_flanks_valleys(merged_matrix5, trunk_diff=1, distance_allowed=100)

    print(result.to_csv(sep=" "))
    print(expected_result5.to_csv(sep=" "))

    assert result.equals(expected_result5)
