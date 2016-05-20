import pytest

from collections import namedtuple

import pandas as pd

MockNamespace = namedtuple("MockNamespace",
                           ["number_cores", "genome", "keep_duplicates",
                            "window_size", "fragment_size", "paired_end",
                            "gaps_allowed", "false_discovery_rate_cutoff",
                            "treatment", "control"])


@pytest.fixture(scope="session", autouse=True)
def args_200_fast():
    return MockNamespace(25, "hg19", False, 200, 150, False, 3, 1,
                         ["examples/test.bed"], ["examples/control.bed"])


@pytest.fixture(scope="session", autouse=True)
def args_200():
    return MockNamespace(1, "hg19", False, 200, 150, False, 3, 0.05,
                         ["examples/test.bed"], ["examples/control.bed"])


@pytest.fixture(scope="session", autouse=True)
def args_50():
    return MockNamespace(1, "hg19", False, 50, 150, False, 3, 0.05,
                         ["examples/test.bed"], ["examples/control.bed"])


@pytest.fixture(scope="session", autouse=True)
def expected_result_example_input():
    return pd.read_table("examples/expected_results.csv", sep=" ", skiprows=1)
