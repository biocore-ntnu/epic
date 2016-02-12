import pytest

import pandas as pd
import numpy as np

from io import StringIO

from collections import OrderedDict

from epic.pool.merge.get_chunks import get_chunks


@pytest.fixture
def file_indexes(epic_bin_file_path1, epic_bin_file_path2):
    return OrderedDict([(epic_bin_file_path1, (4, 6)), (epic_bin_file_path2, (
        2, 3))])


@pytest.fixture
def expected_result():
    return [
        pd.read_table(
            StringIO("""chr1 10150 0 4 1 3 5
chr1 10200 0 3 0 6 1
chr1 10250 0 1 0 2 1"""),
            sep="\s+",
            index_col=[0, 1],
            header=None), pd.read_table(
                StringIO("""
chr1 10150 0 4 1 3 5
chr1 10200 0 3 0 6 1
        """),
                sep="\s+",
                index_col=[0, 1],
                header=None)
    ]


@pytest.mark.mergepools
def test_get_chunks(file_indexes, expected_result):
    result = get_chunks(file_indexes)

    for a, x in zip(result, expected_result):

        assert a.equals(x)


@pytest.fixture
def epic_bin_file_path1(tmpdir):
    epic_bin_file = tmpdir.join("df1_file.txt")
    epic_bin_file.write(
        """Chromosome	Bin	Island	Exp1_3_Polymerase_II	Exp2_3_Polymerase_II	Exp2_3h_Input	Exp1_3h_Input
chr1	10000	0	3	0	4	4
chr1	10050	0	1	1	5	4
chr1	10100	0	3	1	1	3
chr1	10150	0	4	1	3	5
chr1	10200	0	3	0	6	1
chr1	10250	0	1	0	2	1
chr1	10300	0	2	0	3	3
chr1	10350	1	5	1	6	4
chr1	10400	1	3	1	1	1
chr1	10450	1	2	0	1	2
chr2	10100	0	1	0	0	1
chr2	10250	0	1	0	0	0
chr2	10350	0	1	1	0	0
chr2	10450	0	0	2	0	0
chr2	10600	0	0	2	1	0
chr2	11350	0	0	1	1	0
chr2	11600	0	1	0	0	1
chr2	11650	0	0	1	0	0
chr2	11700	0	0	1	1	0
chr2	11750	0	1	0	1	0""")
    return str(epic_bin_file)


@pytest.fixture
def epic_bin_file_path2(tmpdir):
    epic_bin_file = tmpdir.join("df2_file.txt")
    epic_bin_file.write(
        """Chromosome	Bin	Island	Exp1_21_Polymerase_II	Exp2_21_Polymerase_II	Exp2_21h_Input	Exp1_21h_Input
chr1	10000	0	3	0	4	4
chr1	10150	0	4	1	3	5
chr1	10200	0	3	0	6	1
chr1	10300	0	2	0	3	3
chr1	10350	1	5	1	6	4
chr1	10400	1	3	1	1	1
chr1	10450	1	2	0	1	2
chr2	10100	0	1	0	0	1
chr2	11600	0	1	0	0	1
chr2	11700	0	0	1	1	0
chr3	11750	0	1	0	1	0""")
    return str(epic_bin_file)
