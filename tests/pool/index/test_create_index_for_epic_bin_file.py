from __future__ import division
from sys import stdout

from collections import OrderedDict
import pytest

import pandas as pd
import numpy as np

from io import StringIO

from epic.pool.index.create_index_for_epic_bin_file import create_index_for_epic_bin_file  #, create_index_df_for_epic_bin_files


@pytest.fixture
def chunk_start_end():
    od = OrderedDict()
    od["chr1"] = [(10000, 10100), (10150, 10250), (10300, 10450)]
    od["chr2"] = [(10100, 10600), (11350, 11750)]
    return od


@pytest.fixture
def epic_bin_file_path(tmpdir):
    epic_bin_file = tmpdir.join("df1_file.txt")
    epic_bin_file.write(
        """Chromosome	Bin	Island	Exp1_3_Polymerase_II	Exp2_3_Polymerase_II	Exp2_3h_Input	Exp1_3h_Input
chr1	10000	False	3	0	4	4
chr1	10050	False	1	1	5	4
chr1	10100	False	3	1	1	3
chr1	10150	False	4	1	3	5
chr1	10200	False	3	0	6	1
chr1	10250	False	1	0	2	1
chr1	10300	False	2	0	3	3
chr1	10350	True	5	1	6	4
chr1	10400	True	3	1	1	1
chr1	10450	True	2	0	1	2
chr2	10100	False	1	0	0	1
chr2	10250	False	1	0	0	0
chr2	10350	False	1	1	0	0
chr2	10450	False	0	2	0	0
chr2	10600	False	0	2	1	0
chr2	11350	False	0	1	1	0
chr2	11600	False	1	0	0	1
chr2	11650	False	0	1	0	0
chr2	11700	False	0	1	1	0
chr2	11750	False	1	0	1	0""")
    return str(epic_bin_file)


@pytest.fixture
def epic_bin_file_path2(tmpdir):
    epic_bin_file = tmpdir.join("df2_file.txt")
    epic_bin_file.write(
        """Chromosome	Bin	Island	Exp1_21_Polymerase_II	Exp2_21_Polymerase_II	Exp2_21h_Input	Exp1_21h_Input
chr1	10000	False	3	0	4	4
chr1	10150	False	4	1	3	5
chr1	10200	False	3	0	6	1
chr1	10300	False	2	0	3	3
chr1	10350	True	5	1	6	4
chr1	10400	True	3	1	1	1
chr1	10450	True	2	0	1	2
chr2	10100	False	1	0	0	1
chr2	11600	False	1	0	0	1
chr2	11700	False	0	1	1	0
chr3	11750	False	1	0	1	0""")
    return str(epic_bin_file)


@pytest.fixture
def epic_bin_file_paths(epic_bin_file_path, epic_bin_file_path2):

    return [epic_bin_file_path, epic_bin_file_path2]


@pytest.mark.mergepools
def test_create_index_for_epic_bin_files(epic_bin_file_path,
                                         epic_bin_file_path2, chunk_start_end,
                                         expected_result):
    result = create_index_for_epic_bin_file(epic_bin_file_path,
                                            chunk_start_end)
    print(result, "result")
    print(expected_result, "expected_result")

    assert result == expected_result


@pytest.fixture
def expected_result():
    od = OrderedDict()
    od["chr1"] = [(1, 3), (4, 6), (7, 10)]
    od["chr2"] = [(11, 15), (16, 20)]
    return od

# yapf: disable
# @pytest.fixture
# def expected_result():

#     return pd.read_table(
#         StringIO("""Chromosome df_file1_start df_file1_end df_file2_start df_file2_end
# chr1 0 2 0 0
# chr1 3 5 1 2
# chr1 6 9 3 6
# chr2 10 14 7 7
# chr2 15 19 8 9
# chr3 -1 -1 10 10"""),
#         sep="\s+", header=0)
# # yapf: enable

# @pytest.mark.mergepools
# def test_create_index_df_for_epic_bin_files(
#         chunk_start_end, epic_bin_file_paths, expected_result):

#     result = create_index_df_for_epic_bin_files(epic_bin_file_paths,
#                                                 chunk_start_end)
#     print(result, "result")
#     print(type(result), "type(result)")

#     assert expected_result.equals(result)

@pytest.fixture
def chunk_start_end2():
    od = OrderedDict()
    od["chr1"] = [(9800, 10450)]
    return od


@pytest.mark.mergepools
def test_create_index_for_epic_bin_files2(epic_bin_file_path,
                                          epic_bin_file_path2, chunk_start_end2,
                                          expected_result2):
    result = create_index_for_epic_bin_file(epic_bin_file_path,
                                            chunk_start_end2)
    print(result, "result")
    print(expected_result2, "expected_result")

    assert result == expected_result2

@pytest.fixture
def expected_result2():
    od = OrderedDict()
    od["chr1"] = [(1, 10)]
    return od
