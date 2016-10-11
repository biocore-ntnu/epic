import pytest

import pandas as pd
import numpy as np

import pkg_resources, os

from io import StringIO

from epic.scripts.overlaps.overlaps import (_compute_region_overlap,
                                            _create_overlap_matrix_regions)

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
