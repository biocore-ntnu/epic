# import pytest


# from io import StringIO

# import pandas as pd
# from numpy import allclose, int64, float64

# from epic.cluster.cluster import create_cluster_bed

# @pytest.fixture
# def input_bed():

#     contents = u"""Chromosome	Start	End	ClusterID	Index	RegionKind	MaxEnrichedCluster	Bins	TotalEnriched
# chr1	9800	10599	chr1_0	chr1_0_9800:10599_trunk	trunk	1	9800,10000,10200,10400	1,1,1,1
# chr1	87000	88999	chr1_3	chr1_3_87000:88999_trunk	trunk	3	87000,87200,87400,87600,87800,88000,88200,88400,88600,88800	3,3,3,3,3,3,3,3,3,2
# chr1	89000	89799	chr1_3	chr1_3_89000:89799_valley	valley	3	89000,89400,89600	1,1,1
# chr1	89800	91199	chr1_3	chr1_3_89800:91199_trunk	trunk	3	89800,90000,90200,90600,90800,91000	2,2,2,2,2,2"""

#     return pd.read_table(StringIO(contents), header=0)



# @pytest.fixture
# def input_bed_simple():

#     contents = u"""Chromosome	Start	End	ClusterID	Index	RegionKind	MaxEnrichedCluster	Bins	TotalEnriched
# chr1	89000	89799	chr1_3	chr1_3_89000:89799_trunk	trunk	3	89000,89400,89600	2,3,2"""

#     return pd.read_table(StringIO(contents), header=0)


# def test_create_cluster_bed_simple(input_bed_simple, expected_result_simple):

#     result = create_cluster_bed(input_bed_simple, 200)

#     print(result.to_csv(sep=" "))
#     print(expected_result_simple.to_csv(sep=" "))

#     assert result.equals(expected_result_simple)


# @pytest.fixture
# def expected_result_simple():

#     contents = u"""Chromosome	Start	End	Name	Score	Strand
# chr1	89000	89199	chr1_3_89000:89799_trunk	0	.
# chr1	89000	89199	chr1_3_89000:89799_trunk	0	.
# chr1	89400	89599	chr1_3_89000:89799_trunk	0	.
# chr1	89400	89599	chr1_3_89000:89799_trunk	0	.
# chr1	89400	89599	chr1_3_89000:89799_trunk	0	.
# chr1	89600	89799	chr1_3_89000:89799_trunk	0	.
# chr1	89600	89799	chr1_3_89000:89799_trunk	0	."""

#     return pd.read_table(StringIO(contents), header=0)

# # chr1	89000	89799	chr1_3	chr1_3_89000:89799_valley	trunk	3	89000,89400,89600	2,3,2

# @pytest.fixture
# def expected_result():
#     contents = u"""Chromosome Start End Name Score Strand
# 0 chr1 9800 9999 chr1_0_9800:10599_trunk 0 .
# 1 chr1 10000 10199 chr1_0_9800:10599_trunk 0 .
# 2 chr1 10200 10399 chr1_0_9800:10599_trunk 0 .
# 3 chr1 10400 10599 chr1_0_9800:10599_trunk 0 .
# 4 chr1 87000 87199 chr1_3_87000:88999_trunk 0 .
# 5 chr1 87000 87199 chr1_3_87000:88999_trunk 0 .
# 6 chr1 87000 87199 chr1_3_87000:88999_trunk 0 .
# 7 chr1 87200 87399 chr1_3_87000:88999_trunk 0 .
# 8 chr1 87200 87399 chr1_3_87000:88999_trunk 0 .
# 9 chr1 87200 87399 chr1_3_87000:88999_trunk 0 .
# 10 chr1 87400 87599 chr1_3_87000:88999_trunk 0 .
# 11 chr1 87400 87599 chr1_3_87000:88999_trunk 0 .
# 12 chr1 87400 87599 chr1_3_87000:88999_trunk 0 .
# 13 chr1 87600 87799 chr1_3_87000:88999_trunk 0 .
# 14 chr1 87600 87799 chr1_3_87000:88999_trunk 0 .
# 15 chr1 87600 87799 chr1_3_87000:88999_trunk 0 .
# 16 chr1 87800 87999 chr1_3_87000:88999_trunk 0 .
# 17 chr1 87800 87999 chr1_3_87000:88999_trunk 0 .
# 18 chr1 87800 87999 chr1_3_87000:88999_trunk 0 .
# 19 chr1 88000 88199 chr1_3_87000:88999_trunk 0 .
# 20 chr1 88000 88199 chr1_3_87000:88999_trunk 0 .
# 21 chr1 88000 88199 chr1_3_87000:88999_trunk 0 .
# 22 chr1 88200 88399 chr1_3_87000:88999_trunk 0 .
# 23 chr1 88200 88399 chr1_3_87000:88999_trunk 0 .
# 24 chr1 88200 88399 chr1_3_87000:88999_trunk 0 .
# 25 chr1 88400 88599 chr1_3_87000:88999_trunk 0 .
# 26 chr1 88400 88599 chr1_3_87000:88999_trunk 0 .
# 27 chr1 88400 88599 chr1_3_87000:88999_trunk 0 .
# 28 chr1 88600 88799 chr1_3_87000:88999_trunk 0 .
# 29 chr1 88600 88799 chr1_3_87000:88999_trunk 0 .
# 30 chr1 88600 88799 chr1_3_87000:88999_trunk 0 .
# 31 chr1 88800 88999 chr1_3_87000:88999_trunk 0 .
# 32 chr1 88800 88999 chr1_3_87000:88999_trunk 0 .
# 33 chr1 89000 89199 chr1_3_89000:89799_valley 0 .
# 34 chr1 89400 89599 chr1_3_89000:89799_valley 0 .
# 35 chr1 89600 89799 chr1_3_89000:89799_valley 0 .
# 36 chr1 89800 89999 chr1_3_89800:91199_trunk 0 .
# 37 chr1 89800 89999 chr1_3_89800:91199_trunk 0 .
# 38 chr1 90000 90199 chr1_3_89800:91199_trunk 0 .
# 39 chr1 90000 90199 chr1_3_89800:91199_trunk 0 .
# 40 chr1 90200 90399 chr1_3_89800:91199_trunk 0 .
# 41 chr1 90200 90399 chr1_3_89800:91199_trunk 0 .
# 42 chr1 90600 90799 chr1_3_89800:91199_trunk 0 .
# 43 chr1 90600 90799 chr1_3_89800:91199_trunk 0 .
# 44 chr1 90800 90999 chr1_3_89800:91199_trunk 0 .
# 45 chr1 90800 90999 chr1_3_89800:91199_trunk 0 .
# 46 chr1 91000 91199 chr1_3_89800:91199_trunk 0 .
# 47 chr1 91000 91199 chr1_3_89800:91199_trunk 0 ."""

#     return pd.read_table(StringIO(contents), header=0, sep="\s+")

# def test_create_cluster_bed(input_bed, expected_result):

#     df = create_cluster_bed(input_bed, 200)

#     print(df.to_csv(sep=" "))
#     print(expected_result.to_csv(sep=" "))

#     assert df.equals(expected_result)

#     # result = create_cluster_bed(input_bed_simple, 200)

#     # print(result.to_csv(sep=" "))
#     # print(expected_result_simple.to_csv(sep=" "))

#     # assert result.equals(expected_result_simple)
