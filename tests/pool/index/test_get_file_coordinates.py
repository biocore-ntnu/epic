from collections import OrderedDict
import pytest

import pandas as pd
import numpy as np

from io import StringIO


@pytest.fixture
def f18(tmpdir):
    contents = """Chromosome	Bin	Island	Exp2_18_Polymerase_II	Exp1_18_Polymerase_II	Exp2_18h_Input	Exp1_18h_Input
1	chr1	10000	True	4	1	2	5
2	chr1	10050	True	2	1	2	2
3	chr1	10100	True	1	2	0	2
4	chr1	10150	True	5	6	2	2
5	chr1	10200	True	6	1	4	2
6	chr1	10250	True	2	2	1	3
7	chr1	10300	True	5	1	4	2
8	chr1	10350	True	5	4	6	4
9	chr1	10400	True	1	3	2	1
10	chrY	59362200	False	1	4	1	0
11	chrY	59362300	False	1	0	1	0
12	chrY	59362400	False	1	0	1	1
13	chrY	59362450	False	1	1	0	1
14	chrY	59362800	False	0	1	1	0
15	chrY	59362900	False	1	0	2	1
16	chrY	59362950	False	1	1	0	0
17	chrY	59363300	False	4	0	1	0
18	chrY	59363450	False	3	0	1	0
19	chrY	59363500	False	0	1	1	0"""

    f = tmpdir.join("f18")
    f.write(contents)

    return str(f)


@pytest.fixture
def f21(tmpdir):
    contents = """Chromosome	Bin	21h_Island	Exp2_21_Polymerase_II	Exp1_21_Polymerase_II	Exp2_21h_Input	Exp1_21h_Input
1	chr1	10000	True	3	3	1	1
2	chr1	10050	True	3	1	1	2
3	chr1	10100	True	4	2	0	0
4	chr1	10150	True	6	6	0	8
5	chr1	10200	True	4	8	4	4
6	chr1	10250	True	0	3	2	2
7	chr1	10300	True	3	5	2	5
8	chr1	10350	True	6	2	2	2
9	chr1	10400	True	1	2	1	0
10	chrY	59362250	False	1	0	0	1
11	chrY	59362300	False	2	0	1	0
12	chrY	59362350	False	2	0	0	1
13	chrY	59362400	False	0	1	0	3
14	chrY	59362900	False	3	0	0	0
15	chrY	59362950	False	2	2	0	0
16	chrY	59363200	False	2	3	0	0
17	chrY	59363300	False	1	2	1	1
18	chrY	59363350	False	1	0	1	0
19	chrY	59363450	False	1	1	0	4"""

    f = tmpdir.join("f21")
    f.write(contents)

    return str(f)


@pytest.fixture
def f24(tmpdir):

    contents = """Chromosome	Bin	24h_Island	Exp1_24_Polymerase_II	Exp2_24_Polymerase_II	Exp1_24h_Input	Exp2_24h_Input
1	chr1	10000	False	5	0	1	1
2	chr1	10050	True	4	2	4	2
3	chr1	10100	True	4	2	3	5
4	chr1	10150	True	1	2	7	3
5	chr1	10200	True	2	1	2	2
6	chr1	10250	True	7	2	2	6
7	chr1	10300	True	1	4	2	1
8	chr1	10350	True	5	2	4	3
9	chr1	10400	True	1	1	3	2
10	chrY	59362450	False	0	2	0	1
11	chrY	59362750	False	1	1	0	0
12	chrY	59362800	False	2	0	1	0
13	chrY	59362900	False	1	0	1	1
14	chrY	59363000	False	0	1	0	1
15	chrY	59363100	False	1	0	0	0
16	chrY	59363300	False	2	2	1	0
17	chrY	59363350	False	1	0	1	0
18	chrY	59363450	False	2	1	2	0
19	chrY	59363500	False	1	0	0	1"""

    f = tmpdir.join("f24")
    f.write(contents)

    return str(f)


@pytest.fixture
def files():
    return [f18, f21, f24]


@pytest.fixture
def starts_and_ends_of_chunks():
    od = OrderedDict()
    od["chr1"] = [(10000, 10200), (10250, 10400)]
    od["chrY"] = [(59362200, 59362800), (59362850, 59363500)]
    return od


@pytest.fixture
def expected_result(f18, f21, f24):

    all_file_coordinates = {}
    all_file_coordinates["chr1"] = [
        OrderedDict([(f18, (1, 5)), (f21, (1, 5)), (f24, (1, 5))]),
        OrderedDict([(f18, (6, 9)), (f21, (6, 9)), (f24, (6, 9))]),
    ]
    all_file_coordinates["chrY"] = [
        OrderedDict([(f18, (10, 14)), (f21, (10, 13)), (f24, (10, 12))]),
        OrderedDict([(f18, (14, 19)), (f21, (13, 19)), (f24, (12, 19))]),
    ]
    return all_file_coordinates


def get_file_coordinates(infiles, places_to_split):
    """Get the (line_start, line_end) of each chunk for each file.

    Keyword Arguments:
    infiles        --
    places_to_split --
    """


@pytest.mark.mergepools
def test_get_file_coordinates(files, expected_result):
    print(expected_result)

    # assert 0
    # assert result == expected_result
