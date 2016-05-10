import pytest

from collections import namedtuple

MockNamespace = namedtuple("MockNamespace",
                           ["number_cores", "genome", "keep_duplicates",
                            "window_size", "fragment_size", "paired_end"])


@pytest.fixture(scope="session", autouse=True)
def args_200():
    return MockNamespace(1, "hg19", False, 200, 150, False)


@pytest.fixture(scope="session", autouse=True)
def args_50():
    return MockNamespace(1, "hg19", False, 50, 150, False)

# @pytest.fixture(scope="session", autouse=True)
# def input_lmfit():
#     return pd.read_table("example_data/input_lmfit.csv", sep=" ")
