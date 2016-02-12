from sys import stdout
from collections import OrderedDict

import pytest

import pandas as pd
import numpy as np

from io import StringIO

from epic.pool.index.find_chunks_start_and_ends import find_chunks_start_and_ends


# yapf: disable
@pytest.fixture
def expected_result():

    od = OrderedDict()
    od["chr1"] = [(10000, 10100), (10150, 10250), (10300, 10450)]
    od["chr2"] = [(10100, 10600), (11350, 11750)]
    return od
# yapf: disable

@pytest.fixture
def indata():

    od = OrderedDict()
    od["chr1"] = [pd.read_table(StringIO("""Chromosome Bin
chr1 10000
chr1 10050
chr1 10100"""), sep=" "),
    pd.read_table(StringIO("""Chromosome Bin
chr1 10150
chr1 10200
chr1 10250"""), sep=" "),
    pd.read_table(StringIO("""Chromosome Bin
chr1 10300
chr1 10350
chr1 10400
chr1 10450"""), sep=" ")]
    od["chr2"] = [pd.read_table(StringIO("""Chromosome Bin
chr2 10100
chr2 10250
chr2 10350
chr2 10450
chr2 10600"""), sep=" "),
pd.read_table(StringIO("""Chromosome Bin
chr2 11350
chr2 11600
chr2 11650
chr2 11700
chr2 11750"""), sep=" ")]

    return od



@pytest.mark.mergepools
def test_create_chunk_index(indata, expected_result):

    result = find_chunks_start_and_ends(indata)
    assert result == expected_result
