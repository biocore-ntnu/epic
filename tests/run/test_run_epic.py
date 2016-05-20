import pytest

import pandas as pd
import numpy as np

import logging
from io import StringIO
from joblib import delayed, Parallel

from epic.run.run_epic import run_epic


def test_run_epic(expected_result_example_input, args_200_fast):
    df = run_epic(args_200_fast)

    int_cols = "Chromosome Start End ChIP Input".split()
    float_cols = "Score Fold_change P FDR".split()

    assert df[int_cols].equals(expected_result_example_input[int_cols])
    assert np.allclose(df[float_cols],
                       expected_result_example_input[float_cols])
