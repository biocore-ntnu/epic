import pytest
import pandas as pd

from io import StringIO

from epic.pool.post_process.get_number_of_islands_in_each_bin_per_sample import get_number_of_islands_in_each_bin_per_sample


@pytest.fixture
def df():

    contents = """Chromosome	Bin	18h_Island	Exp2_18_Polymerase_II	Exp1_18_Polymerase_II	Exp2_18h_Input	Exp1_18h_Input	21h_Island	Exp2_21_Polymerase_II	Exp1_21_Polymerase_II	Exp2_21h_Input	Exp1_21h_Input	24h_Island	Exp1_24_Polymerase_II	Exp2_24_Polymerase_II	Exp1_24h_Input	Exp2_24h_Input
chr1	10000	1	4	1	2	5	1	3	3	1	1	0	5	0	1	1
chr1	10050	1	2	1	2	2	1	3	1	1	2	1	4	2	4	2
chr1	10100	1	1	2	0	2	1	4	2	0	0	1	4	2	3	5
chr1	10150	1	5	6	2	2	1	6	6	0	8	1	1	2	7	3
chr1	10200	1	6	1	4	2	1	4	8	4	4	1	2	1	2	2
chr1	10250	1	2	2	1	3	1	0	3	2	2	1	7	2	2	6
chr1	10300	1	5	1	4	2	1	3	5	2	5	1	1	4	2	1
chr1	10350	1	5	4	6	4	1	6	2	2	2	1	5	2	4	3
chr1	10400	1	1	3	2	1	1	1	2	1	0	1	1	1	3	2
chrY	59362200	0	1	4	1	0	0	0	0	0	0	0	0	0	0	0
chrY	59362250	0	0	0	0	0	0	1	0	0	1	0	0	0	0	0
chrY	59362300	0	1	0	1	0	0	2	0	1	0	0	0	0	0	0
chrY	59362350	0	0	0	0	0	0	2	0	0	1	0	0	0	0	0
chrY	59362400	0	1	0	1	1	0	0	1	0	3	0	0	0	0	0
chrY	59362450	0	1	1	0	1	0	0	0	0	0	0	0	2	0	1
chrY	59362750	0	0	0	0	0	0	0	0	0	0	0	1	1	0	0
chrY	59362800	0	0	1	1	0	0	0	0	0	0	0	2	0	1	0
chrY	59362900	0	1	0	2	1	0	3	0	0	0	0	1	0	1	1
chrY	59362950	0	1	1	0	0	0	2	2	0	0	0	0	0	0	0
chrY	59363000	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1
chrY	59363100	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0
chrY	59363200	0	0	0	0	0	0	2	3	0	0	0	0	0	0	0
chrY	59363300	0	4	0	1	0	0	1	2	1	1	0	2	2	1	0
chrY	59363350	0	0	0	0	0	0	1	0	1	0	0	1	0	1	0
chrY	59363450	0	3	0	1	0	0	1	1	0	4	0	2	1	2	0
chrY	59363500	0	0	1	1	0	0	0	0	0	0	0	1	0	0	1"""

    return pd.read_table(StringIO(contents), sep="\s+", index_col=[0, 1])


@pytest.fixture
def expected_result():
    return pd.read_table(
        StringIO("""
Chromosome Bin 18h_Island 21h_Island 24h_Island
chr1 10000 1 1 0
chr1 10050 1 1 1
chr1 10100 1 1 1
chr1 10150 1 1 1
chr1 10200 1 1 1
chr1 10250 1 1 1
chr1 10300 1 1 1
chr1 10350 1 1 1
chr1 10400 1 1 1
chrY 59362200 0 0 0
chrY 59362250 0 0 0
chrY 59362300 0 0 0
chrY 59362350 0 0 0
chrY 59362400 0 0 0
chrY 59362450 0 0 0
chrY 59362750 0 0 0
chrY 59362800 0 0 0
chrY 59362900 0 0 0
chrY 59362950 0 0 0
chrY 59363000 0 0 0
chrY 59363100 0 0 0
chrY 59363200 0 0 0
chrY 59363300 0 0 0
chrY 59363350 0 0 0
chrY 59363450 0 0 0
chrY 59363500 0 0 0
    """),
        sep="\s+",
        index_col=[0, 1])


@pytest.mark.mergepools
def test_get_island_info(df, expected_result):
    result = get_number_of_islands_in_each_bin_per_sample(df)

    assert expected_result.equals(result)
