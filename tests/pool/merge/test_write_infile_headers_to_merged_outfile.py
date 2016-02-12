import pytest

import pandas as pd
import numpy as np

from io import StringIO

from epic.pool.merge.write_infile_headers_to_merged_outfile import write_infile_headers_to_merged_outfile, get_index_and_headers


@pytest.fixture
def expected_result_get():

    return "\t".join("Chromosome	Bin".split(
    ) + """df1_file.txt_Island	Exp1_3_Polymerase_II	Exp2_3_Polymerase_II	Exp2_3h_Input	Exp1_3h_Input  df2_file.txt_Island	Exp1_21_Polymerase_II	Exp2_21_Polymerase_II	Exp2_21h_Input	Exp1_21h_Input""".split(
    ))


@pytest.mark.mergepools
def test_get_index_and_headers(infiles, expected_result_get):

    result = get_index_and_headers(infiles, [0, 1])
    assert result == expected_result_get


@pytest.fixture
def epic_bin_file_path1(tmpdir):
    epic_bin_file = tmpdir.join("df1_file.txt")
    epic_bin_file.write(
        """Chromosome	Bin	df1_file.txt_Island	Exp1_3_Polymerase_II	Exp2_3_Polymerase_II	Exp2_3h_Input	Exp1_3h_Input
chr1	10000	0	3	0	4	4
chr1	10050	0	1	1	5	4
chr2	11750	0	1	0	1	0""")
    return str(epic_bin_file)


@pytest.fixture
def epic_bin_file_path2(tmpdir):
    epic_bin_file = tmpdir.join("df2_file.txt")
    epic_bin_file.write(
        """Chromosome	Bin	df2_file.txt_Island	Exp1_21_Polymerase_II	Exp2_21_Polymerase_II	Exp2_21h_Input	Exp1_21h_Input
chr1	10000	0	3	0	4	4
chr1	10150	0	4	1	3	5
chr3	11750	0	1	0	1	0""")
    return str(epic_bin_file)


@pytest.fixture
def infiles(epic_bin_file_path1, epic_bin_file_path2):
    return [epic_bin_file_path1, epic_bin_file_path2]


@pytest.fixture
def expected_result_write(expected_result_get, tmpdir):
    f = tmpdir.join("expected_result_write")
    f.write(expected_result_get + "\n")

    return str(f)


@pytest.fixture
def output_filename(tmpdir):
    f = tmpdir.join("start_file_might_have_content")
    f.write("oogabooga\nheyhey")

    return str(f)


@pytest.mark.mergepools
def test_write_index_and_headers_to_file(expected_result_write, infiles,
                                         output_filename):

    infile_index_cols = [0, 1]
    write_infile_headers_to_merged_outfile(infiles, output_filename,
                                           infile_index_cols)
    with open(output_filename) as r, open(expected_result_write) as x:
        actual, expected = r.readline(), x.readline()
        print(actual, "actual")
        print(expected, "expected")
        assert actual == expected
