# import pytest

# from io import StringIO

# import pandas as pd

# from example_pipeline.scripts.add_color_trunks_flanks_valleys_bed import (add_color_trunks_flanks_valleys_bed)

# @pytest.fixture
# def df():
#     c = """Chromosome      Start   End     ClusterID       Index   RegionKind      MaxEnrichedCluster      Bins    TotalEnriched
# chr1    89000   89799   chr1_3  chr1_3_89000:89799_valley       valley  3       89000,89200,89400,89600 1,1,1,1
# chr1    89800   91199   chr1_3  chr1_3_89800:91199_trunk        trunk   3       89800,90000,90200,90400,90600,90800,91000       2,2,2,2,2,2,2
# chr1    123600  125399  chr1_4  chr1_4_123600:125399_trunk      trunk   1       123600,123800,124000,124200,124400,124600,124800,125000,125200       1,1,1,1,1,1,1,1,1
# chr1    234200  235999  chr1_5  chr1_5_234200:235999_trunk      trunk   3       234200,234400,234600,234800,235000,235200,235400,235600,235800       3,3,3,3,3,3,3,3,3
# chr1    236000  236199  chr1_5  chr1_5_236000:236199_flank      flank   3       236000  1"""

#     return pd.read_table(StringIO(c), sep="\s+", header=0)


# @pytest.fixture
# def expected_result():

#     c = """Chromosome Start End Index Score Strand ThickStart ThickEnd RegionKind
# chr1 89000 89799 chr1_3_89000:89799_valley 0 . 89000 89799 128,0,128
# chr1 89800 91199 chr1_3_89800:91199_trunk 0 . 89800 91199 255,215,0
# chr1 123600 125399 chr1_4_123600:125399_trunk 0 . 123600 125399 255,215,0
# chr1 234200 235999 chr1_5_234200:235999_trunk 0 . 234200 235999 255,215,0
# chr1 236000 236199 chr1_5_236000:236199_flank 0 . 236000 236199 255,140,0"""

#     return pd.read_table(StringIO(c), sep="\s+")


# def test_add_color_trunks_flanks_valleys_bed(df, expected_result):

#     result = add_color_trunks_flanks_valleys_bed(df)

#     print(result.to_csv(sep=" ", index=False))
#     print(expected_result)

#     assert result.equals(expected_result)
