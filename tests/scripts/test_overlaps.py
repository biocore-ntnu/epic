import pytest

import pandas as pd
import numpy as np

import pkg_resources, os

from io import StringIO

from epic.config.genomes import (create_genome_size_dict,
                                 get_effective_genome_length)

__author__ = "Endre Bakken Stovner https://github.com/endrebak/"
__license__ = "MIT"



@pytest.fixture
def region_matrixes(epic_overlap_intermediate_region_matrixes):

    return [pd.read_table(f, sep=" ", index_col=0) for f in epic_overlap_intermediate_region_matrixes]


@pytest.fixture
def expected_result_region_overlap():

    df = pd.read_table("examples/epic-overlaps/region_overlap_result.csv", sep=" ", index_col=0)

    return df


@pytest.mark.current
def test__compute_region_overlap(region_matrixes, expected_result_region_overlap):

    region_matrix = region_matrixes[0]

    df = _compute_region_overlap(region_matrix)
    print(df)

    # df.to_csv("examples/epic-overlaps/region_overlap_result.csv", sep=" ")
    # print(df)

    assert df.equals(expected_result_region_overlap)


@pytest.fixture
def expected_result_intermediate_region_matrix(epic_overlap_intermediate_region_matrixes):

    df = pd.read_table(epic_overlap_intermediate_region_matrixes[0], sep=" ", index_col=0)
    return df

@pytest.mark.unit
def test_create_region_matrix(epic_overlap_files, expected_result_intermediate_region_matrix):

    all_files = epic_overlap_files[:3]

    df = _create_overlap_matrix_regions(all_files[0], all_files)
    # for i, f in enumerate(all_files):
    #     df = _create_overlap_matrix_regions(f, all_files)
    #     b = basename(f).split(".")[0]
    #     df.to_csv("examples/epic-overlaps/" + b + "_region.matrix", sep=" ")
    print(df.head())
    print(expected_result_intermediate_region_matrix.head())

    assert df.equals(expected_result_intermediate_region_matrix)





@pytest.fixture
def expected_result_intermediate_nucleotide_matrix(epic_overlap_intermediate_nucleotide_matrixes):

    df = pd.read_table(epic_overlap_intermediate_nucleotide_matrixes[0], sep=" ", index_col=0)
    return df


@pytest.mark.unit
def test_create_nucleotide_matrix(epic_overlap_files, expected_result_intermediate_nucleotide_matrix):

    all_files = epic_overlap_files[:3]

    df = _create_overlap_matrix_nucleotides(all_files[0], all_files, 200)
    print(df, "actual")
    print(expected_result_intermediate_nucleotide_matrix, "expected")

    # print(all_files)
    # for i, f in enumerate(all_files):
    #     df = _create_overlap_matrix_nucleotides(f, all_files, 200)
    #     b = basename(f).split(".")[0]
    #     df.to_csv("examples/epic-overlaps/" + b + ".matrix", sep=" ")

    assert df.equals(expected_result_intermediate_nucleotide_matrix)



@pytest.fixture
def nucleotide_matrixes(epic_overlap_intermediate_nucleotide_matrixes):

    return [pd.read_table(f, sep=" ", index_col=0) for f in epic_overlap_intermediate_nucleotide_matrixes]


@pytest.fixture
def expected_result_nucleotide_overlap():

    df = pd.read_table("examples/epic-overlaps/nucleotide_overlap_result.csv", sep=" ", index_col=0)

    return df


@pytest.mark.unit
def test__compute_nucleotide_overlap(nucleotide_matrixes, expected_result_nucleotide_overlap):

    nucleotide_matrix = nucleotide_matrixes[0]

    df = _compute_nucleotide_overlap(nucleotide_matrix, 200)

    # df.to_csv("examples/epic-overlaps/nucleotide_overlap_result.csv", sep=" ")
    # print(df)

    assert df.equals(expected_result_nucleotide_overlap)
