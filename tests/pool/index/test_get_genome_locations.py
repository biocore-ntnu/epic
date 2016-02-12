from sys import stdout

import pytest

import pandas as pd
import numpy as np

from io import StringIO

from epic.pool.index.get_genome_locations import get_union_of_all_genome_locations


@pytest.mark.mergepools
def test_get_union_of_all_genome_locations(pool_df1, pool_df2,
                                           expected_result):

    indexes = get_union_of_all_genome_locations([pool_df1, pool_df2])
    print(indexes)
    print(expected_result)

    assert indexes.equals(expected_result)


@pytest.fixture
def expected_result():

    return pd.read_table(
        StringIO("""Chromosome Bin
chr1	10000
chr1	10050
chr1	10550
chr2	10450
chr2	10600"""),
        header=0,
        sep="\s+")


@pytest.fixture
def pool_df1(tmpdir):
    tmpfile1 = tmpdir.join("df1_file.txt")
    tmpfile1.write(
        """Chromosome	Bin	Island	Exp1_3_Polymerase_II	Exp2_3_Polymerase_II	Exp2_3h_Input	Exp1_3h_Input
chr1	10000	False	3	0	4	4
chr1	10550	False	1	1	5	4
chr2	10450	False	0	2	0	0
chr2	10600	False	0	2	1	0""")
    return str(tmpfile1)


@pytest.fixture
def pool_df2(tmpdir):
    tmpfile2 = tmpdir.join("df2_file.txt")
    tmpfile2.write(
        """Chromosome	Bin	Island	Exp1_3_Polymerase_II	Exp2_3_Polymerase_II	Exp2_3h_Input	Exp1_3h_Input
chr1	10000	False	3	0	4	4
chr1	10050	False	1	1	5	4
chr2	10450	False	0	2	0	0
chr2	10600	False	0	2	1	0""")

    return str(tmpfile2)

# @pytest.fixture
# def pool_df1():
#     return pd.read_table(StringIO("""Chromosome	Bin	Island Exp1_3_Polymerase_II	Exp2_3_Polymerase_II	Exp2_3h_Input	Exp1_3h_Input
# chr1	10000	False	3	0	4	4
# chr1	10050	False	1	1	5	4
# chr1	10100	False	3	1	1	3
# chr1	10150	False	4	1	3	5
# chr1	10200	False	3	0	6	1
# chr1	10250	False	1	0	2	1
# chr1	10300	False	2	0	3	3
# chr1	10350	True	5	1	6	4
# chr1	10400	True	3	1	1	1
# chr1	10450	True	2	0	1	2
# chr2	10100	False	1	0	0	1
# chr2	10250	False	1	0	0	0
# chr2	10350	False	1	1	0	0
# chr2	10450	False	0	2	0	0
# chr2	10600	False	0	2	1	0
# chr2	11350	False	0	1	1	0
# chr2	11600	False	1	0	0	1
# chr2	11650	False	0	1	0	0
# chr2	11700	False	0	1	1	0
# chr2	11750	False	1	0	1	0"""), index_col=[0, 1], header=0, sep="\s+").drop("Island", axis=1)

# @pytest.fixture
# def pool_df2():
#     return pd.read_table(StringIO("""Chromosome	Bin	Island	Exp1_24_Polymerase_II	Exp2_24_Polymerase_II	Exp1_24h_Input	Exp2_24h_Input
# chr1	10000	False	5	0	1	1
# chr1	10050	True	4	2	4	2
# chr1	10100	True	4	2	3	5
# chr1	10150	True	1	2	7	3
# chr1	10200	True	2	1	2	2
# chr1	10250	True	7	2	2	6
# chr1	10300	True	1	4	2	1
# chr1	10350	True	5	2	4	3
# chr1	10400	True	1	1	3	2
# chr1	10450	True	1	3	2	2
# chr2	10150	False	0	1	1	0
# chr2	10300	False	0	1	0	0
# chr2	10350	False	1	0	0	0
# chr2	11350	False	1	0	1	0
# chr2	11550	False	1	0	0	0
# chr2	11650	False	0	1	1	1
# chr2	11800	False	1	0	0	1
# chr2	11850	False	0	1	0	0
# chr2	11950	False	0	1	1	1
# chr2	12150	False	1	0	0	0"""), index_col=[0, 1], header=0, sep="\s+").drop("Island", axis=1)

# @pytest.fixture
# def expected_result():

#     return pd.read_table(StringIO("""Chromosome	Bin	Exp1_3_Polymerase_II	Exp2_3_Polymerase_II	Exp2_3h_Input	Exp1_3h_Input	Exp1_24_Polymerase_II	Exp2_24_Polymerase_II	Exp1_24h_Input	Exp2_24h_Input
# chr1	10000	3	0	4	4	5	0	1	1
# chr1	10050	1	1	5	4	4	2	4	2
# chr1	10100	3	1	1	3	4	2	3	5
# chr1	10150	4	1	3	5	1	2	7	3
# chr1	10200	3	0	6	1	2	1	2	2
# chr1	10250	1	0	2	1	7	2	2	6
# chr1	10300	2	0	3	3	1	4	2	1
# chr1	10350	5	1	6	4	5	2	4	3
# chr1	10400	3	1	1	1	1	1	3	2
# chr1	10450	2	0	1	2	1	3	2	2
# chr2	10100	1	0	0	1	0	0	0	0
# chr2	10150	0	0	0	0	0	1	1	0
# chr2	10250	1	0	0	0	0	0	0	0
# chr2	10300	0	0	0	0	0	1	0	0
# chr2	10350	1	1	0	0	1	0	0	0
# chr2	10450	0	2	0	0	0	0	0	0
# chr2	10600	0	2	1	0	0	0	0	0
# chr2	11350	0	1	1	0	1	0	1	0
# chr2	11550	0	0	0	0	1	0	0	0
# chr2	11600	1	0	0	1	0	0	0	0
# chr2	11650	0	1	0	0	0	1	1	1
# chr2	11700	0	1	1	0	0	0	0	0
# chr2	11750	1	0	1	0	0	0	0	0
# chr2	11800	0	0	0	0	1	0	0	1
# chr2	11850	0	0	0	0	0	1	0	0
# chr2	11950	0	0	0	0	0	1	1	1
# chr2	12150	0	0	0	0	1	0	0	0"""), index_col=[0, 1], header=0, sep="\s+")
