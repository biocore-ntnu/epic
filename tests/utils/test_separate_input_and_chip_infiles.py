import pytest

import pandas as pd
import numpy as np

from io import StringIO

from epic.utils.separate_input_and_chip_infiles import separate_input_and_chip_infiles


def describe_separate_input_and_chip_infiles():
    @pytest.mark.unit
    def with_a_pair_of_files_and_default_input_pattern():

        result = separate_input_and_chip_infiles(["bed/input1.bed",
                                                  "bed/test.bed"])
        assert result == (["bed/test.bed"], ["bed/input1.bed"])

    @pytest.mark.unit
    def with_a_pattern_that_matches_all_files():
        with pytest.raises(ValueError) as raise_info:
            separate_input_and_chip_infiles(
                ["bed/input1.bed", "bed/test.bed"], "bed")
        assert "All files considered input." in str(raise_info.value)

    @pytest.mark.unit
    def with_a_pattern_that_matches_no_files():
        with pytest.raises(ValueError) as raise_info:
            separate_input_and_chip_infiles(
                ["bed/input1.bed", "bed/test.bed"], "rural juror")
        assert "No input files found." in str(raise_info.value)
