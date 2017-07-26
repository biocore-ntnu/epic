import pytest

import pandas as pd
from io import StringIO

from epic.merge.merge import merge_matrixes

@pytest.fixture()
def expected_result_no_enriched():

    c = u"""Chromosome Bin TotalEnriched data/align/AAG_KO_ChIP_1.bed data/align/AAG_KO_ChIP_2.bed data/align/AAG_KO_ChIP_3.bed data/align/AAG_KO_Input_1.bed data/align/AAG_KO_Input_2.bed data/align/AAG_KO_Input_3.bed data/align/ELP1_KO_ChIP_1.bed data/align/ELP1_KO_ChIP_2.bed data/align/ELP1_KO_ChIP_3.bed data/align/ELP1_KO_Input_1.bed data/align/ELP1_KO_Input_2.bed data/align/ELP1_KO_Input_3.bed data/align/WT_ChIP_1.bed data/align/WT_ChIP_2.bed data/align/WT_ChIP_3.bed data/align/WT_Input_1.bed data/align/WT_Input_2.bed data/align/WT_Input_3.bed
chr1 10000 0 0 2 0 5 1 3 0 3 0 3 2 1 0 0 1 2 2 0
chr1 10200 0 0 0 0 2 0 0 0 1 1 2 2 1 0 1 0 2 0 0
chr1 10400 0 0 0 0 1 0 1 0 0 0 2 1 1 0 0 0 2 3 4
chr1 10600 0 0 0 0 2 0 3 0 0 0 2 2 1 0 0 0 0 0 1
chr1 10800 0 0 0 0 6 8 5 0 1 0 0 11 4 0 0 0 2 1 15
chr1 11000 0 0 0 0 7 4 8 0 0 0 3 12 12 0 3 0 1 3 10
chr1 11200 0 0 0 0 4 0 1 0 0 0 1 0 3 0 0 0 0 1 0
chr1 11400 0 1 0 0 2 0 0 0 0 0 2 0 0 0 2 1 0 0 0
chr1 12200 0 0 0 0 5 0 1 0 0 0 3 6 8 0 0 0 3 0 2"""

    return pd.read_table(StringIO(c), sep=" ", index_col=[0, 1, 2])


@pytest.fixture()
def expected_result_one_enriched():

    c = u"""Chromosome Bin TotalEnriched data/align/AAG_KO_ChIP_1.bed data/align/AAG_KO_ChIP_2.bed data/align/AAG_KO_ChIP_3.bed data/align/AAG_KO_Input_1.bed data/align/AAG_KO_Input_2.bed data/align/AAG_KO_Input_3.bed data/align/ELP1_KO_ChIP_1.bed data/align/ELP1_KO_ChIP_2.bed data/align/ELP1_KO_ChIP_3.bed data/align/ELP1_KO_Input_1.bed data/align/ELP1_KO_Input_2.bed data/align/ELP1_KO_Input_3.bed data/align/WT_ChIP_1.bed data/align/WT_ChIP_2.bed data/align/WT_ChIP_3.bed data/align/WT_Input_1.bed data/align/WT_Input_2.bed data/align/WT_Input_3.bed
chr1 10200 1 0 0 0 2 0 0 0 1 1 2 2 1 0 1 0 2 0 0
chr1 10400 1 0 0 0 1 0 1 0 0 0 2 1 1 0 0 0 2 3 4
chr1 10600 1 0 0 0 2 0 3 0 0 0 2 2 1 0 0 0 0 0 1"""

    return pd.read_table(StringIO(c), sep=" ", index_col=[0, 1, 2])


@pytest.mark.unit
def test_merge_matrixes_no_enriched_raises(merge_matrixes_dfs_no_enriched, expected_result_no_enriched):

    dfs = merge_matrixes_dfs_no_enriched

    keep_nonenriched = False
    enriched_per_file = False
    nb_cpus = 2

    with pytest.raises(Exception):
        df = merge_matrixes(dfs, keep_nonenriched, enriched_per_file, False, nb_cpus)


@pytest.mark.unit
def test_merge_matrixes_one_enriched(merge_matrixes_dfs_one_enriched, expected_result_one_enriched):

    dfs = merge_matrixes_dfs_one_enriched

    keep_nonenriched = False
    enriched_per_file = False
    nb_cpus = 2

    df = merge_matrixes(dfs, keep_nonenriched, enriched_per_file, False, nb_cpus)

    print(df.head().to_csv(sep=" "))
    print(expected_result_one_enriched.head().to_csv(sep=" "))

    assert  df.equals(expected_result_one_enriched)
