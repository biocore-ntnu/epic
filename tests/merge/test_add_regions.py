import pytest
from collections import OrderedDict

import pandas as pd
from io import StringIO


from epic.merge.merge import merge_matrixes
from epic.merge.merge_helpers import add_new_enriched_bins_matrixes


@pytest.fixture
def regions(tmpdir):

    fs = []
    for i, c in enumerate(["""chr1    10000   10599   2.2761062711783457e-05  67.49046260339546       .
chr1    72000   91599   2.4770408905838545e-226 235.4664881299362       .""",
           """chr1    9800    15199   0.0048446172754557214   33.652547110032025      .
chr1    40000   41199   2.187570707966001e-08   1000.0  .""",
           """chr1    9800    10599   3.239383152206723e-79   204.30687218918862      .
chr1    38800   40799   2.4798100382025985e-11  1000.0  ."""]):

        name = str(i)
        f = tmpdir.mkdir(name).join(name)
        f.write(c)
        fs.append(str(f))

    return fs


@pytest.fixture
def dfs():

    od = OrderedDict()

    for n, c in [("fibroblast.matrix.gz","""Chromosome Bin Enriched_fibroblast chrX/ChIP_1_fibroblast.bed.gz chrX/ChIP_2_fibroblast.bed.gz chrX/ChIP_3_fibroblast.bed.gz chrX/Input_1_fibroblast.bed.gz chrX/Input_2_fibroblast.bed.gz chrX/Input_3_fibroblast.bed.gz
chr1 10000 1 4 14 13 2 4 14
chr1 10200 1 17 24 14 9 9 16
chr1 10400 1 3 1 1 1 0 2
chr1 11400 1 0 0 1 0 0 0
chr1 11600 1 0 0 0 0 0 0
chr1 11800 1 1 0 1 0 0 0
chr1 12000 1 0 0 0 0 0 0
chr1 12200 1 0 0 1 0 0 0
chr1 12400 1 0 0 0 0 0 0
chr1 12600 1 0 0 0 1 0 0
chr1 12800 1 0 0 0 0 0 0
chr1 13000 1 3 6 4 1 3 2
chr1 13200 1 0 1 0 1 1 0
chr1 13400 1 7 4 5 1 1 3
chr1 13600 1 0 0 0 0 0 0
chr1 13800 1 0 0 1 0 0 0
chr1 14000 1 0 0 0 0 0 0
chr1 14200 1 0 0 0 0 0 0
chr1 14400 1 0 0 0 0 0 0"""),
("keratinocyte.matrix.gz", """Chromosome Bin Enriched_keratinocyte chrX/ChIP_1_keratinocyte.bed.gz chrX/ChIP_2_keratinocyte.bed.gz chrX/ChIP_3_keratinocyte.bed.gz chrX/Input_1_keratinocyte.bed.gz chrX/Input_2_keratinocyte.bed.gz chrX/Input_3_keratinocyte.bed.gz
chr1 9800 1 1 0 0 2 0 0
chr1 10000 1 13 15 17 11 2 17
chr1 10200 1 13 25 23 16 2 24
chr1 10400 1 2 0 2 10 0 3
chr1 10600 1 0 0 0 0 0 0
chr1 10800 1 0 0 1 0 0 0
chr1 11000 1 0 0 0 0 0 0
chr1 11200 1 0 0 0 0 0 0
chr1 11400 1 0 0 0 0 0 0
chr1 11600 1 0 0 1 0 0 0
chr1 11800 1 0 0 0 0 0 0
chr1 12000 1 0 0 2 0 0 0
chr1 12200 1 0 0 0 0 0 0
chr1 12400 1 0 0 0 0 0 0
chr1 12600 1 0 0 0 0 0 1
chr1 12800 1 0 1 1 0 0 0
chr1 13000 1 0 0 6 7 1 3
chr1 13200 1 0 2 0 2 1 1
chr1 13400 1 1 1 6 4 0 2"""),
        ("melanocyte.matrix.gz", """Chromosome Bin Enriched_melanocyte chrX/ChIP_1_melanocyte.bed.gz chrX/ChIP_2_melanocyte.bed.gz chrX/ChIP_3_melanocyte.bed.gz chrX/Input_1_melanocyte.bed.gz chrX/Input_2_melanocyte.bed.gz chrX/Input_3_melanocyte.bed.gz
chr1 9800 1 0 0 2 0 0 0
chr1 10000 1 13 3 128 2 2 21
chr1 10200 1 15 8 96 5 3 23
chr1 10400 1 3 0 4 3 0 7
chr1 11800 0 1 0 0 0 0 0
chr1 12000 0 0 0 0 0 0 4
chr1 12200 0 0 0 0 0 0 1
chr1 12400 0 0 1 0 0 0 1
chr1 12600 0 1 0 0 0 0 0
chr1 13000 0 5 3 1 0 0 2
chr1 13200 0 1 3 0 0 0 1
chr1 13400 0 3 0 3 1 0 3
chr1 13800 0 1 0 0 0 0 0
chr1 14600 0 1 0 1 0 0 0
chr1 14800 0 1 0 4 1 0 4
chr1 15000 0 1 0 1 1 0 2
chr1 15600 0 0 1 0 0 0 1
chr1 15800 0 1 0 0 0 0 0
chr1 16000 0 1 1 2 1 2 2""")]:
        df = pd.read_table(StringIO(c), sep="\s+", header=0, index_col=[0, 1])
        od[n] = df

    return od


# def test_add_regions(regions, dfs):

#     result = add_new_enriched_bins_matrixes(regions, dfs, 200)

#     assert 0


@pytest.fixture
def expected_result():

    c = """Chromosome Bin TotalEnriched chrX/ChIP_1_fibroblast.bed.gz chrX/ChIP_1_keratinocyte.bed.gz chrX/ChIP_1_melanocyte.bed.gz chrX/ChIP_2_fibroblast.bed.gz chrX/ChIP_2_keratinocyte.bed.gz chrX/ChIP_2_melanocyte.bed.gz chrX/ChIP_3_fibroblast.bed.gz chrX/ChIP_3_keratinocyte.bed.gz chrX/ChIP_3_melanocyte.bed.gz chrX/Input_1_fibroblast.bed.gz chrX/Input_1_keratinocyte.bed.gz chrX/Input_1_melanocyte.bed.gz chrX/Input_2_fibroblast.bed.gz chrX/Input_2_keratinocyte.bed.gz chrX/Input_2_melanocyte.bed.gz chrX/Input_3_fibroblast.bed.gz chrX/Input_3_keratinocyte.bed.gz chrX/Input_3_melanocyte.bed.gz
chr1 9800 2.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 2.0 0.0 2.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
chr1 10000 3.0 4.0 13.0 13.0 14.0 15.0 3.0 13.0 17.0 128.0 2.0 11.0 2.0 4.0 2.0 2.0 14.0 17.0 21.0
chr1 10200 3.0 17.0 13.0 15.0 24.0 25.0 8.0 14.0 23.0 96.0 9.0 16.0 5.0 9.0 2.0 3.0 16.0 24.0 23.0
chr1 10400 3.0 3.0 2.0 3.0 1.0 0.0 0.0 1.0 2.0 4.0 1.0 10.0 3.0 0.0 0.0 0.0 2.0 3.0 7.0
chr1 10600 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
chr1 10800 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
chr1 11000 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
chr1 11200 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
chr1 11400 2.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
chr1 11600 2.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
chr1 11800 2.0 1.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
chr1 12000 2.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 4.0
chr1 12200 2.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0
chr1 12400 2.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0
chr1 12600 2.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0
chr1 12800 2.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
chr1 13000 2.0 3.0 0.0 5.0 6.0 0.0 3.0 4.0 6.0 1.0 1.0 7.0 0.0 3.0 1.0 0.0 2.0 3.0 2.0
chr1 13200 2.0 0.0 0.0 1.0 1.0 2.0 3.0 0.0 0.0 0.0 1.0 2.0 0.0 1.0 1.0 0.0 0.0 1.0 1.0
chr1 13400 2.0 7.0 1.0 3.0 4.0 1.0 0.0 5.0 6.0 3.0 1.0 4.0 1.0 1.0 0.0 0.0 3.0 2.0 3.0
chr1 13600 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
chr1 13800 1.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
chr1 14000 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
chr1 14200 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
chr1 14400 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
chr1 14600 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
chr1 14800 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 4.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 4.0
chr1 15000 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 2.0
chr1 15600 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0
chr1 15800 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
chr1 16000 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 2.0 0.0 0.0 1.0 0.0 0.0 2.0 0.0 0.0 2.0"""

    return pd.read_table(StringIO(c), sep=" ", index_col=[0, 1, 2])


@pytest.fixture
def simple_dfs():

    od = OrderedDict()

    for n, c in [("melanocyte.matrix", """Chromosome Bin Enriched chrX/ChIP_1_melanocyte.bed.gz chrX/ChIP_2_melanocyte.bed.gz chrX/Input_1_melanocyte.bed.gz chrX/Input_2_melanocyte.bed.gz
chr1 800 1 0 2 0 0
chr1 1200 1 13 128 2 2"""),
                 ("fibroblast.matrix", """Chromosome Bin Enriched chrX/ChIP_1_fibroblast.bed.gz chrX/ChIP_2_fibroblast.bed.gz chrX/Input_1_fibroblast.bed.gz chrX/Input_2_fibroblast.bed.gz
chr1 800 1 0 2 0 0
chr1 1200 1 13 128 2 2""")]:

        df = pd.read_table(StringIO(c), sep="\s+", header=0, index_col=[0, 1])
        od[n] = df

    return od
# Need to fix test data
# covered in test merge wiht without equal
# def test_merge_matrixes_with_regions(dfs, regions, expected_result):

#     enriched_per_file, keep_nonenriched = False, False

#     df = merge_matrixes(dfs, keep_nonenriched, regions, enriched_per_file, 1)

#     print(df.to_csv(sep=" "))
#     print(expected_result.to_csv(sep=" "))

#     assert df.equals(expected_result)




@pytest.fixture
def simple_regions(tmpdir):

    fs = []
    for i, c in enumerate([
            """chr1	600	1200	2.2761062711783457e-05	67.49046260339546	.""",
            """chr1	400	1600	0.0048446172754557214	33.652547110032025	."""],):

        name = str(i)
        f = tmpdir.mkdir(name).join(name)
        f.write(c)
        fs.append(str(f))

    return fs

@pytest.fixture
def simple_expected_result():

    pass
# melanocyte.matrix
# Chromosome Bin Enriched chrX/ChIP_1_melanocyte.bed.gz chrX/ChIP_2_melanocyte.bed.gz chrX/Input_1_melanocyte.bed.gz chrX/Input_2_melanocyte.bed.gz
# chr1 400 1.0 0.0 0.0 0.0 0.0
# chr1 600 2.0 0.0 0.0 0.0 0.0
# chr1 800 2.0 0.0 2.0 0.0 0.0
# chr1 1000 2.0 0.0 0.0 0.0 0.0
# chr1 1200 2.0 13.0 128.0 2.0 2.0
# chr1 1400 1.0 0.0 0.0 0.0 0.0
# chr1 1600 1.0 0.0 0.0 0.0 0.0

# fibroblast.matrix
# Chromosome Bin Enriched chrX/ChIP_1_fibroblast.bed.gz chrX/ChIP_2_fibroblast.bed.gz chrX/Input_1_fibroblast.bed.gz chrX/Input_2_fibroblast.bed.gz
# chr1 400 1.0 0.0 0.0 0.0 0.0
# chr1 600 2.0 0.0 0.0 0.0 0.0
# chr1 800 2.0 0.0 2.0 0.0 0.0
# chr1 1000 2.0 0.0 0.0 0.0 0.0
# chr1 1200 2.0 13.0 128.0 2.0 2.0
# chr1 1400 1.0 0.0 0.0 0.0 0.0
# chr1 1600 1.0 0.0 0.0 0.0 0.0

def test_simple_add_regions(simple_dfs, simple_regions):

    result = add_new_enriched_bins_matrixes(simple_regions, simple_dfs, 200)

    for k, v in result.items():
        print(k)
        print(v.to_csv(sep=" "))

    # print(result.to_csv(sep=" "))

    assert 0
