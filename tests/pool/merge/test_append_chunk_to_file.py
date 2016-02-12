import tempfile
from os import close, devnull

import pytest

import pandas as pd
import numpy as np

from io import StringIO

from epic.pool.merge.append_chunk_to_file import append_chunk_to_file


@pytest.fixture
def existing_data(tmpdir):

    content = """chr1	10000	0	3	0	4	4	0	3	0	4	4
chr1	10050	0	1	1	5	4	0	0	0	0	0
chr1	10100	0	3	1	1	3	0	0	0	0	0
chr1	10150	0	4	1	3	5	0	4	1	3	5
chr1	10200	0	3	0	6	1	0	3	0	6	1
chr1	10250	0	1	0	2	1	0	0	0	0	0
"""

    f = tmpdir.join("existing_data.csv")
    f.write(content)

    return str(f)


@pytest.fixture
def chunk():
    return pd.read_table(
        StringIO("""chr1 10300 0 2 0 3 3 0 2 0 3 3
chr1 10350 1 5 1 6 4 1 5 1 6 4
chr1 10400 1 3 1 1 1 1 3 1 1 1
chr1 10450 1 2 0 1 2 1 2 0 1 2"""),
        sep="\s+",
        index_col=[0, 1],
        header=None)


@pytest.fixture
def expected_result(tmpdir):
    content = """chr1 10000 0 3 0 4 4 0 3 0 4 4
chr1 10050 0 1 1 5 4 0 0 0 0 0
chr1 10100 0 3 1 1 3 0 0 0 0 0
chr1 10150 0 4 1 3 5 0 4 1 3 5
chr1 10200 0 3 0 6 1 0 3 0 6 1
chr1 10250 0 1 0 2 1 0 0 0 0 0
chr1 10300 0 2 0 3 3 0 2 0 3 3
chr1 10350 1 5 1 6 4 1 5 1 6 4
chr1 10400 1 3 1 1 1 1 3 1 1 1
chr1 10450 1 2 0 1 2 1 2 0 1 2"""

    return pd.read_table(StringIO(content), sep="\s+", header=None)


@pytest.mark.mergepool
def test_append_chunk_to_file(chunk, existing_data, expected_result):

    append_chunk_to_file(chunk, existing_data)
    actual_result = pd.read_table(existing_data, sep="	", header=None)
    print(actual_result.to_csv(sep=" "), "a")
    print(actual_result.columns, "a")
    print(actual_result.shape, "a")
    print(expected_result.to_csv(sep=" "), "x")
    print(expected_result.columns, "x")
    print(expected_result.shape, "x")

    assert expected_result.equals(actual_result)
