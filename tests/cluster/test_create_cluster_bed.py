import pytest


from io import StringIO

import pandas as pd
from numpy import allclose, int64, float64

from epic.cluster.cluster import trunks_flanks_valleys

@pytest.fixture
def input_bed():

    contents = u"""Chromosome	Start	End	ClusterID	Index	RegionKind	MaxEnrichedCluster	Bins	TotalEnriched
chr1	9800	10599	chr1_0	chr1_0_9800:10599_trunk	trunk	1	9800,10000,10200,10400	1,1,1,1
chr1	87000	88999	chr1_3	chr1_3_87000:88999_trunk	trunk	3	87000,87200,87400,87600,87800,88000,88200,88400,88600,88800	3,3,3,3,3,3,3,3,3,2
chr1	89000	89799	chr1_3	chr1_3_89000:89799_valley	valley	3	89000,89400,89600	1,1,1
chr1	89800	91199	chr1_3	chr1_3_89800:91199_trunk	trunk	3	89800,90000,90200,90600,90800,91000	2,2,2,2,2,2"""

    return pd.read_table(StringIO(contents), header=0)



@pytest.fixture
def input_bed_simple():

    contents = u"""Chromosome	Start	End	ClusterID	Index	RegionKind	MaxEnrichedCluster	Bins	TotalEnriched
chr1	89000	89799	chr1_3	chr1_3_89000:89799_trunk	trunk	3	89000,89400,89600	2,3,2"""

    return pd.read_table(StringIO(contents), header=0)


def create_cluster_bed(bed, bin_size):

    rowdicts = []
    for _, row in bed.iterrows():

        total_enriched = row.TotalEnriched.split(",")
        bins = row.Bins.split(",")
        cluster_id = row.Index
        chromosome = row.Chromosome

        for _bin, number in zip(bins, total_enriched):
            _bin, number = int(_bin), int(number)
            rowdict = {"Chromosome": chromosome, "Start": _bin, "End": _bin +
                       bin_size - 1, "Name": cluster_id, "Score": 0, "Strand": "."}

            rowdicts.extend([rowdict] * number)

    return pd.DataFrame.from_dict(rowdicts)["Chromosome Start End Name Score Strand".split()]
# ['1', '1', '1']
# ['89000', '89400', '89600']


def test_create_cluster_bed_simple(input_bed_simple, expected_result_simple):

    result = create_cluster_bed(input_bed_simple, 200)

    print(result.to_csv(sep=" "))
    print(expected_result_simple.to_csv(sep=" "))

    assert result.equals(expected_result_simple)

@pytest.fixture
def expected_result_simple():

    contents = u"""Chromosome	Start	End	Name	Score	Strand
chr1	89000	89199	chr1_3_89000:89799_trunk	0	.
chr1	89000	89199	chr1_3_89000:89799_trunk	0	.
chr1	89400	89599	chr1_3_89000:89799_trunk	0	.
chr1	89400	89599	chr1_3_89000:89799_trunk	0	.
chr1	89400	89599	chr1_3_89000:89799_trunk	0	.
chr1	89600	89799	chr1_3_89000:89799_trunk	0	.
chr1	89600	89799	chr1_3_89000:89799_trunk	0	."""

    return pd.read_table(StringIO(contents), header=0)

# chr1	89000	89799	chr1_3	chr1_3_89000:89799_valley	trunk	3	89000,89400,89600	2,3,2

def test_create_cluster_bed(): #input_bed, expected_result_simple):

    assert 0

    # result = create_cluster_bed(input_bed_simple, 200)

    # print(result.to_csv(sep=" "))
    # print(expected_result_simple.to_csv(sep=" "))

    # assert result.equals(expected_result_simple)
