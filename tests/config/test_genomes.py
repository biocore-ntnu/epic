import pytest

import pandas as pd
import numpy as np

from io import StringIO

from epic.config.genomes import create_genome_size_dict


@pytest.fixture
def input_data():
    pass


@pytest.fixture
def expected_result():
    pass


@pytest.mark.unit
def test_create_genome_size_dict(input_data, expected_result):

    result = create_genome_size_dict("hg19")
    print(result)
    assert 0
    # assert result == expected_result
