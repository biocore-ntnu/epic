from sys import stdout

import pytest

import pandas as pd
import numpy as np

from io import StringIO

# def merge_epic_bins(genome_regions, epic_binfile_paths):

#     for genome_region in genome_regions:
#         epic_dataframes_to_merge = []
#         for epic_binfile_path in epic_binfile_paths:
#             df = pd.read_table(epic_binfile_path, index_col=[0, 1])

#         print(pd.concat(epic_dataframes_to_merge, axis=1).to_csv())

# @pytest.fixture
# def genome_regions(expected_result):
#     genome_regions = expected_result.reset_index()[["Chromosome", "Bin"]]

#     genome_regions
#     split_genome_regions = np.array_split(genome_regions, 3)
#     split_genome_regions_as_index = [genome_region.set_index(["Chromosome", "Bin"]) for genome_region in split_genome_regions]


@pytest.fixture
def expected_result():

    return pd.read_table(StringIO(
        """Chromosome	Bin	Island	Exp1_3_Polymerase_II	Exp2_3_Polymerase_II	Exp2_3h_Input	Exp1_3h_Input	Exp1_21_Polymerase_I	Exp2_21_Polymerase_II	Exp2_21h_Input	Exp1_21h_Input
chr1	10000	False	3	0	4	4	3	0	4	4
chr1	10050	False	0	0	0	0	1	1	5	4
chr1	10550	False	1	1	5	4	0	0	0	0
chr2	10450	False	0	2	0	0	0	2	0	0
chr2	10600	False	0	2	1	0	0	2	1	0"""))


@pytest.fixture
def epic_binfile_paths(pool_df1, pool_df2):
    return [pool_df1, pool_df2]


@pytest.fixture
def pool_df1(tmpdir):
    tmpfile1 = tmpdir.join("df1_file.txt")
    tmpfile1.write(
        """Chromosome	Bin	Island	Exp1_3_Polymerase_II	Exp2_3_Polymerase_II	Exp2_3h_Input	Exp1_3h_Input
chr1	10000	False	3	0	4	4
chr1	10550	False	1	1	5	4
chr2	10450	False	0	2	0	0
chr2	10600	False	0	2	1	0""")
    return str(tmpfile1)


@pytest.fixture
def pool_df2(tmpdir):
    tmpfile2 = tmpdir.join("df2_file.txt")
    tmpfile2.write(
        """Chromosome	Bin	Island	Exp1_3_Polymerase_II	Exp2_3_Polymerase_II	Exp2_3h_Input	Exp1_3h_Input
chr1	10000	False	3	0	4	4
chr1	10050	False	1	1	5	4
chr2	10450	False	0	2	0	0
chr2	10600	False	0	2	1	0""")

    return str(tmpfile2)

# from sys import stdout

# import pytest

# import pandas as pd
# from pandas import Index
# import numpy as np

# from io import StringIO

# @pytest.mark.mergepools
# def test_split_genome_location_df_into_chunks(genome_locations_df, expected_result):

#     split_chromosome_df = split_genome_location_df_into_chunks(expected_result, 1)
#     print(expected_result, "expected_result")
#     print(split_df, "split_df")

#     assert all([x.equals(s) for x, s in zip(expected_result, split_df)])

# def split_genome_location_df_into_chunks(genome_locations_df, nb_chunks):
#     return np.array_split(genome_locations_df, nb_chunks)

#     print(split_chromosome_df, "split_chromosome_df")
# def genome_locations_df():
#     assert all([x.equals(s) for x, s in zip(expected_result, split_chromosome_df)])
#     return pd.read_table(StringIO("""Chromosome Bin
# chr1	10000
# chr1	10050
# chr1	10550
# chr2	10450
# chr2	10600"""), index_col=[0, 1], header=0, sep="\s+").index

# @pytest.fixture
# def expected_result():

#     return [Index([('chr1', 10000)], dtype='object'), Index([('chr1', 10050)], dtype='object'), Index([('chr1', 10550)], dtype='object'), Index([('chr2', 10450)], dtype='object'), Index([('chr2', 10600)], dtype='object')]
