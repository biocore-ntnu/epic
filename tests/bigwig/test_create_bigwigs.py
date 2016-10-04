import pytest

import pandas as pd
import numpy as np

from os import stat
from io import StringIO

from epic.bigwig.create_bigwigs import _create_bigwig


@pytest.fixture
def input_data():
    contents = u"""Chromosome Bin End examples/test.bed
chr1 887600 887799 0
chr1 994600 994799 0
chr1 1041000 1041199 0
chr1 1325200 1325399 1
chr1 1541600 1541799 1
chr1 1599000 1599199 1
chr1 1770200 1770399 0
chr1 1820200 1820399 1
chr1 1995000 1995199 0
chr1 2063800 2063999 0
chr1 2129400 2129599 0
chr1 2239000 2239199 0
chr1 2318800 2318999 0
chr1 2448200 2448399 1
chr1 3006000 3006199 0
chr1 3046000 3046199 1
chr1 3089200 3089399 0
chr1 3093800 3093999 0
chr1 3096400 3096599 0"""

    return pd.read_table(StringIO(contents), sep="\s+", index_col=[0, 1, 2])


@pytest.fixture
def output_bigwig(tmpdir):
    p = tmpdir.mkdir("sub").join("outfile.bw")
    return str(p)


# not used due to travis/bioconda r problem
@pytest.mark.unit
def test_create_bigwigs(input_data, output_bigwig, args_200):

    _create_bigwig(input_data, output_bigwig, args_200)

    filesize = stat(output_bigwig).st_size

    print(filesize, "filesize")

    assert filesize > 0
