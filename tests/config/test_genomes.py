import pytest

import pandas as pd
import numpy as np

import pkg_resources, os

from io import StringIO

from epic.config.genomes import create_genome_size_dict

@pytest.fixture
def input_data():
    pass


@pytest.fixture
def expected_result():
    return {'chr8': 146364022, 'chr12': 133851895, 'chr9': 141213431, 'chr21': 48129895, 'chr15': 102531392, 'chr14': 107349540, 'chrY': 59373566, 'chr10': 135534747, 'chr17': 81195210, 'chrX': 155270560, 'chr11': 135006516, 'chr18': 78077248, 'chr2': 243199373, 'chr6': 171115067, 'chr3': 198022430, 'chr20': 63025520, 'chr22': 51304566, 'chr13': 115169878, 'chr19': 59128983, 'chr5': 180915260, 'chrM': 16571, 'chr7': 159138663, 'chr16': 90354753, 'chr4': 191154276, 'chr1': 249250621} 


@pytest.mark.unit
def test_create_genome_size_dict(input_data, expected_result):

    result = create_genome_size_dict("hg19")
    assert result == expected_result
