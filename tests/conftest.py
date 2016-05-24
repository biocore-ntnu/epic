import pytest

from collections import namedtuple

import pandas as pd

__author__ = "Endre Bakken Stovner https://github.com/endrebak/"
__license__ = "MIT"

MockNamespace = namedtuple("MockNamespace",
                           ["number_cores", "genome", "keep_duplicates",
                            "window_size", "fragment_size", "paired_end",
                            "gaps_allowed", "false_discovery_rate_cutoff",
                            "effective_genome_length", "treatment", "control"])

egs = 2290813547.4  # this is the effective genome size used by the original sicer for hg19


@pytest.fixture(scope="session")
def args_200_fast():
    return MockNamespace(25, "hg19", False, 200, 150, False, 3, 1, egs,
                         ["examples/test.bed"], ["examples/control.bed"])


@pytest.fixture(scope="session")
def args_200():
    return MockNamespace(1, "hg19", False, 200, 150, False, 3, 0.05, egs,
                         ["examples/test.bed"], ["examples/control.bed"])


@pytest.fixture(scope="session")
def args_50():
    return MockNamespace(1, "hg19", False, 50, 150, False, 3, 0.05, egs,
                         ["examples/test.bed"], ["examples/control.bed"])


@pytest.fixture(scope="session")
def expected_result_example_input():
    return pd.read_table("examples/expected_results.csv", sep=" ", skiprows=1)
