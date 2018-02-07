import pytest

from collections import namedtuple

import pandas as pd
import numpy as np

from os import stat
from io import StringIO

from epic.bigwig.create_bigwigs import _create_bigwig, create_log2fc_data
from epic.config.genomes import create_genome_size_dict

MockNamespace = namedtuple("MockNamespace", ["treatment", "control"])

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


@pytest.mark.unit
def test_create_bigwigs(input_data, output_bigwig, args_200_fast):

    print(input_data)
    print(output_bigwig)
    d = create_genome_size_dict(args_200_fast.genome)

    print(d)
    _create_bigwig(input_data, output_bigwig, d)

    filesize = stat(output_bigwig).st_size

    print(filesize, "filesize")

    assert filesize > 0


@pytest.fixture
def args():

    return MockNamespace(["/mnt/cargo/epigenome_roadmap/ChIP_1_keratinocyte.bed.gz", "/mnt/cargo/epigenome_roadmap/ChIP_2_keratinocyte.bed.gz", "/mnt/cargo/epigenome_roadmap/ChIP_3_keratinocyte.bed.gz"], ["/mnt/cargo/epigenome_roadmap/Input_1_keratinocyte.bed.gz", "/mnt/cargo/epigenome_roadmap/Input_2_keratinocyte.bed.gz", "/mnt/cargo/epigenome_roadmap/Input_3_keratinocyte.bed.gz"])


@pytest.fixture
def matrix():

    c = u"""Chromosome Bin Enriched /mnt/cargo/epigenome_roadmap/ChIP_1_keratinocyte.bed.gz /mnt/cargo/epigenome_roadmap/ChIP_2_keratinocyte.bed.gz /mnt/cargo/epigenome_roadmap/ChIP_3_keratinocyte.bed.gz /mnt/cargo/epigenome_roadmap/Input_1_keratinocyte.bed.gz /mnt/cargo/epigenome_roadmap/Input_2_keratinocyte.bed.gz /mnt/cargo/epigenome_roadmap/Input_3_keratinocyte.bed.gz
chr1 4163800 0 0 3 0 3 0 2
chr1 41630000 0 0 0 0 0 3 3
chr1 41630200 0 0 1 0 0 1 6
chr1 41630400 0 1 1 3 0 0 1
chr1 41630600 0 1 1 1 0 5 1
chr1 41630800 0 1 2 0 1 0 7
chr1 41631000 0 4 0 3 0 3 4
chr1 41631200 0 1 1 3 1 2 4
chr1 41631400 0 1 1 0 1 0 1
chr1 41631600 0 1 0 2 1 3 4
chr1 41631800 0 0 0 2 0 1 0
chr1 41632000 0 2 0 2 0 2 3
chr1 41632200 0 1 0 1 0 2 0
chr1 41632400 0 0 1 1 0 1 0
chr1 41632600 0 0 0 0 1 5 5
chr1 41632800 0 1 0 0 1 1 1"""

    return pd.read_table(StringIO(c), sep=" ", header=0)


@pytest.mark.unit
def test_create_log2fc_data(matrix, args):

    result = create_log2fc_data(matrix, args)

    print(result)
