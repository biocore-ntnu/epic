import pytest

from io import StringIO

import pandas as pd

from epic.cluster.cluster import trunks_flanks_valleys

@pytest.fixture
def merged_matrix():

    c = """Chromosome Bin TotalEnriched A B
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
    c = """Index A    B
chr1_0_0:599_trunk      2.0  1.0
chr1_1_800:1599_trunk   2.0  1.0
chr1_2_1800:2199_trunk  2.0  0.0"""

    return pd.read_table(StringIO(c), sep="\s+", header=0, index_col="Index")


def test_trunks_flanks_valleys(merged_matrix, expected_result):

    result = trunks_flanks_valleys(merged_matrix)

    print(result)
    print(expected_result)

    assert result.equals(expected_result)


@pytest.fixture
def merged_matrix2():

    c = """Chromosome Bin TotalEnriched A B
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

    c = """Index   A    B
chr1_0_0:199_flank       1.0  0.0
chr1_0_200:599_trunk     1.0  1.0
chr1_1_800:1199_flank    1.0  1.0
chr1_1_1200:1399_trunk   0.0  0.0
chr1_1_1400:1799_valley  2.0  0.0
chr1_1_1800:1999_trunk   1.0  0.0
chr1_1_2000:2199_flank   1.0  0.0"""

    return pd.read_table(StringIO(c), sep="\s+", header=0, index_col="Index")

def test_trunks_flanks_valleys2(merged_matrix2, expected_result2):

    result = trunks_flanks_valleys(merged_matrix2)

    assert result.equals(expected_result2)
