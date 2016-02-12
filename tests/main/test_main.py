import pytest

from epic.run.run_epic import run_epic


@pytest.fixture
def fun():
    pass


@pytest.mark.slow
def test_run_epic():

    print("Need to refactor and test main module!\n" * 5)
    assert 0
