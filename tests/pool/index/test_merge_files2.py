from os.path import basename
from collections import OrderedDict
from natsort import natsorted
from itertools import chain
from sys import stdout

import pytest

import pandas as pd
import numpy as np

from io import StringIO


def merge_files(dict_of_file_locations):
    """Horizontally concatenate files according to the indexes provided.

    To horizontally concatenate (join) several files too large to work with in
    memory at the same time, each file (the keys of the dicts) is split
    according to the list of file positions given for that file (the items of
    the dicts). Index nb X for each file, corresponds to the same genomic
    region.

    BAD: include something about why split (how the method uses this split).

    Keyword Arguments: dict_of_file_locations -- dict of dicts. The key is the
    chromosome. The dicts within the dicts contain the file locations (value)
    for each file (key).

    """

    files = _get_all_files(dict_of_file_locations)
    header_dict = _get_headers_from_files(files)

    print(header_dict)
    for chromosome, file_locations_dict in dict_of_file_locations.items():

        pass
        # horizontally_aligned_chunk = pd.concat(df)
        # horizontally_aligned_chunks.append(df)
        # print(chunks_to_horizontally_align, "chunks_to_horizontally_align")
        # print(
        #     type(chunks_to_horizontally_align),
        #     "type(chunks_to_horizontally_align)")

        # cdf = pd.concat(horizontally_aligned_chunks)
        # print(cdf.to_csv(sep=" "), "cdf")
        # islands = df.Island
        # islands_sum = islands.sum()
        # df = df.drop("Island", axis=1)
        # print(df.to_csv(sep=" "))
        # print(islands_sum)
        # df.insert(0, "Nb_islands", islands_sum)

        # print(df.to_csv(sep=" "), "df")
        # print(df.columns, "df.columns")
        # print(df.Island)

        #     horizontally_aligned_chunk = pd.concat(chunks_to_horizontally_align,
        #                                            axis=1)

        #     horizontally_aligned_chunks.append(horizontally_aligned_chunk)

        # df = pd.concat(horizontally_aligned_chunks).reset_index()
        # df.columns = ["Chromosome", "Bin"] + list(df.columns[2:])

        # return df


@pytest.fixture
def chr1_indexes(epic_indexes):

    od = OrderedDict()
    od["chr1"] = OrderedDict([(epic_bin_file_path1, [(0, 3), (4, 6), (7, 9)]),
                              (epic_bin_file_path2, [(0, 1), (2, 3), (4, 6)])])
    return epic_indexes["chr1"]


@pytest.fixture
def expected_result_merge_chr1(expected_result):
    return expected_result.xs("chr1", level="Chromosome", drop_level=False)


# @pytest.mark.mergepool
def test_merge_chromosome(chr1_indexes, expected_result_merge_chr1,
                          header_dict):

    print(chr1_indexes, "chr1_indexes")
    print(
        expected_result_merge_chr1.to_csv(sep=" "),
        "expected_result_merge_chr1")
    result = merge_chromosome(chr1_indexes, header_dict)
    print(result.to_csv(sep=" "), "result")
    print(result.shape, expected_result_merge_chr1.shape)
    # assert 0
    # assert result == expected_result


def merge_chromosome(file_indexes, header_dict):

    total_nb_file_locations = _get_nb_of_file_locations_in_each_file(
        file_indexes)
    horizontally_aligned_chunks = []

    for i in range(total_nb_file_locations):

        print(i)
        for c in chunks_to_horizontally_align:
            print(c.to_csv(sep=" "))
        horizontally_aligned_chunk = pd.concat(chunks_to_horizontally_align,
                                               axis=1).fillna(False)
        horizontally_aligned_chunks.append(horizontally_aligned_chunk)

    return pd.concat(horizontally_aligned_chunks)


def _sum_island_columns(df):
    island_columns = [col for col in df.columns if col.endswith("_Island")]
    island_df = df[island_columns].sum(axis=1, numeric_only=True)
    df = df.drop(island_columns, axis=1)

    df.insert(0, "Nb_islands", island_df)
    return df


@pytest.fixture
def expected_result_sum_islands():
    return pd.read_table(
        StringIO(
            """Chromosome Bin Nb_islands Exp1_3_Polymerase_II Exp2_3_Polymerase_II Exp2_3h_Input Exp1_3h_Input Exp1_21_Polymerase_II Exp2_21_Polymerase_II Exp2_21h_Input Exp1_21h_Input
chr1 10300 0 2 0 3 3 2 0 3 3
chr1 10350 2 5 1 6 4 5 1 6 4
chr1 10400 2 3 1 1 1 3 1 1 1"""),
        sep="\s+",
        index_col=[0, 1])


@pytest.fixture
def merged_df_chunk():
    return pd.read_table(
        StringIO(
            """Chromosome Bin df1_file.txt_Island Exp1_3_Polymerase_II Exp2_3_Polymerase_II Exp2_3h_Input Exp1_3h_Input df2_file.txt_Island Exp1_21_Polymerase_II Exp2_21_Polymerase_II Exp2_21h_Input Exp1_21h_Input
chr1 10300 False 2 0 3 3 False 2 0 3 3
chr1 10350 True 5 1 6 4 True 5 1 6 4
chr1 10400 True 3 1 1 1 True 3 1 1 1"""),
        sep="\s+",
        index_col=[0, 1])


@pytest.mark.mergepool
def test__sum_island_columns(merged_df_chunk, expected_result_sum_islands):

    result = _sum_island_columns(merged_df_chunk)

    # assert result.equals(expected_result_sum_islands)


def _get_all_files(dict_of_file_locations):
    files = set()
    for d in dict_of_file_locations.values():
        files.update(set(d.keys()))
    return natsorted(list(files))


def _get_nb_of_file_locations_in_each_file(file_locations):

    nb_indexes = set()
    for chromosome, d in file_locations.items():
        nb_indexes.add(len(d))

    assert len(
        nb_indexes) == 1, "The lists of file locations do not all have the same length!"

    return list(nb_indexes)[0]


@pytest.mark.mergepool
def test_merge_multiple_files(epic_indexes, expected_result):

    result = merge_files(epic_indexes)
    print(result)
    print(expected_result)
    # assert 0  # need to create chunk index for file two first


@pytest.fixture
def expected_result():
    return pd.read_table(
        StringIO(
            """Chromosome Bin Island Exp1_3_Polymerase_II Exp2_3_Polymerase_II Exp2_3h_Input Exp1_3h_Input Exp1_21_Polymerase_II Exp2_21_Polymerase_II Exp2_21h_Input Exp1_21h_Input
chr1 10000 False 3.0 0.0 4.0 4.0 3.0 0.0 4.0 4.0
chr1 10050 False 1.0 1.0 5.0 4.0 0.0 0.0 0.0 0.0
chr1 10100 False 3.0 1.0 1.0 3.0 0.0 0.0 0.0 0.0
chr1 10150 False 4.0 1.0 3.0 5.0 4.0 1.0 3.0 5.0
chr1 10200 False 3.0 0.0 6.0 1.0 3.0 0.0 6.0 1.0
chr1 10250 False 1.0 0.0 2.0 1.0 0.0 0.0 0.0 0.0
chr1 10300 False 2.0 0.0 3.0 3.0 2.0 0.0 3.0 3.0
chr1 10350 True 5.0 1.0 6.0 4.0 5.0 1.0 6.0 4.0
chr1 10400 True 3.0 1.0 1.0 1.0 3.0 1.0 1.0 1.0
chr1 10450 True 2.0 0.0 1.0 2.0 2.0 0.0 1.0 2.0
chr2 10100 False 1.0 0.0 0.0 1.0 1.0 0.0 0.0 1.0
chr2 10250 False 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
chr2 10350 False 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0
chr2 10450 False 0.0 2.0 0.0 0.0 0.0 0.0 0.0 0.0
chr2 10600 False 0.0 2.0 1.0 0.0 0.0 0.0 0.0 0.0
chr2 11350 False 0.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0
chr2 11600 False 1.0 0.0 0.0 1.0 1.0 0.0 0.0 1.0
chr2 11650 False 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0
chr2 11700 False 0.0 1.0 1.0 0.0 0.0 1.0 1.0 0.0
chr2 11750 False 1.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0
chr3 11750 False 0.0 0.0 0.0 0.0 1.0 0.0 1.0 0.0"""),
        sep="\s+",
        index_col=[0, 1, 2])


@pytest.fixture
def expected_result_single(epic_bin_file_path1):
    return pd.read_table(epic_bin_file_path1, sep="\t")


@pytest.fixture
def file_locations(epic_bin_file_path1, epic_bin_file_path2):
    return [epic_bin_file_path1, epic_bin_file_path2]


@pytest.fixture
def epic_bin_file_path1(tmpdir):
    epic_bin_file = tmpdir.join("df1_file.txt")
    epic_bin_file.write(
        """Chromosome	Bin	Island	Exp1_3_Polymerase_II	Exp2_3_Polymerase_II	Exp2_3h_Input	Exp1_3h_Input
chr1	10000	False	3	0	4	4
chr1	10050	False	1	1	5	4
chr1	10100	False	3	1	1	3
chr1	10150	False	4	1	3	5
chr1	10200	False	3	0	6	1
chr1	10250	False	1	0	2	1
chr1	10300	False	2	0	3	3
chr1	10350	True	5	1	6	4
chr1	10400	True	3	1	1	1
chr1	10450	True	2	0	1	2
chr2	10100	False	1	0	0	1
chr2	10250	False	1	0	0	0
chr2	10350	False	1	1	0	0
chr2	10450	False	0	2	0	0
chr2	10600	False	0	2	1	0
chr2	11350	False	0	1	1	0
chr2	11600	False	1	0	0	1
chr2	11650	False	0	1	0	0
chr2	11700	False	0	1	1	0
chr2	11750	False	1	0	1	0""")
    return str(epic_bin_file)


@pytest.fixture
def epic_bin_file_path2(tmpdir):
    epic_bin_file = tmpdir.join("df2_file.txt")
    epic_bin_file.write(
        """Chromosome	Bin	Island	Exp1_21_Polymerase_II	Exp2_21_Polymerase_II	Exp2_21h_Input	Exp1_21h_Input
chr1	10000	False	3	0	4	4
chr1	10150	False	4	1	3	5
chr1	10200	False	3	0	6	1
chr1	10300	False	2	0	3	3
chr1	10350	True	5	1	6	4
chr1	10400	True	3	1	1	1
chr1	10450	True	2	0	1	2
chr2	10100	False	1	0	0	1
chr2	11600	False	1	0	0	1
chr2	11700	False	0	1	1	0
chr3	11750	False	1	0	1	0""")
    return str(epic_bin_file)


@pytest.fixture
def epic_indexes(epic_bin_file_path1, epic_bin_file_path2):
    od = OrderedDict()
    od["chr1"] = OrderedDict([(epic_bin_file_path1, [(0, 3), (4, 6), (7, 9)]),
                              (epic_bin_file_path2, [(0, 1), (2, 3), (4, 6)])])
    od["chr2"] = OrderedDict([(epic_bin_file_path1, [(10, 13), (14, 16), (
        17, 19)]), (epic_bin_file_path2, [(7, 7), (7, 8), (8, 9)])])
    od["chr3"] = OrderedDict([(epic_bin_file_path2, [(10, 10)])])
    return od
