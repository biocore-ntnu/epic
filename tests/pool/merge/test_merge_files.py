import pytest

import pandas as pd
import numpy as np

from io import StringIO

from epic.pool.merge.merge_files import merge_files


@pytest.fixture
def file18h(tmpdir):
    contents = """Chromosome	Bin	18h_Island	Exp2_18_Polymerase_II	Exp1_18_Polymerase_II	Exp2_18h_Input	Exp1_18h_Input
chr1	10000	True	4	1	2	5
chr1	10050	True	2	1	2	2
chr1	10100	True	1	2	0	2
chr1	10150	True	5	6	2	2
chr1	10200	True	6	1	4	2
chr1	10250	True	2	2	1	3
chr1	10300	True	5	1	4	2
chr1	10350	True	5	4	6	4
chr1	10400	True	1	3	2	1
chrY	59362200	False	1	4	1	0
chrY	59362300	False	1	0	1	0
chrY	59362400	False	1	0	1	1
chrY	59362450	False	1	1	0	1
chrY	59362800	False	0	1	1	0
chrY	59362900	False	1	0	2	1
chrY	59362950	False	1	1	0	0
chrY	59363300	False	4	0	1	0
chrY	59363450	False	3	0	1	0
chrY	59363500	False	0	1	1	0"""

    f = tmpdir.join("file18h")
    f.write(contents)

    return str(f)


@pytest.fixture
def file21h(tmpdir):
    contents = """Chromosome	Bin	21h_Island	Exp2_21_Polymerase_II	Exp1_21_Polymerase_II	Exp2_21h_Input	Exp1_21h_Input
chr1	10000	True	3	3	1	1
chr1	10050	True	3	1	1	2
chr1	10100	True	4	2	0	0
chr1	10150	True	6	6	0	8
chr1	10200	True	4	8	4	4
chr1	10250	True	0	3	2	2
chr1	10300	True	3	5	2	5
chr1	10350	True	6	2	2	2
chr1	10400	True	1	2	1	0
chrY	59362250	False	1	0	0	1
chrY	59362300	False	2	0	1	0
chrY	59362350	False	2	0	0	1
chrY	59362400	False	0	1	0	3
chrY	59362900	False	3	0	0	0
chrY	59362950	False	2	2	0	0
chrY	59363200	False	2	3	0	0
chrY	59363300	False	1	2	1	1
chrY	59363350	False	1	0	1	0
chrY	59363450	False	1	1	0	4"""

    f = tmpdir.join("file21h")
    f.write(contents)

    return str(f)


@pytest.fixture
def file24h(tmpdir):

    contents = """Chromosome	Bin	24h_Island	Exp1_24_Polymerase_II	Exp2_24_Polymerase_II	Exp1_24h_Input	Exp2_24h_Input
chr1	10000	False	5	0	1	1
chr1	10050	True	4	2	4	2
chr1	10100	True	4	2	3	5
chr1	10150	True	1	2	7	3
chr1	10200	True	2	1	2	2
chr1	10250	True	7	2	2	6
chr1	10300	True	1	4	2	1
chr1	10350	True	5	2	4	3
chr1	10400	True	1	1	3	2
chrY	59362450	False	0	2	0	1
chrY	59362750	False	1	1	0	0
chrY	59362800	False	2	0	1	0
chrY	59362900	False	1	0	1	1
chrY	59363000	False	0	1	0	1
chrY	59363100	False	1	0	0	0
chrY	59363300	False	2	2	1	0
chrY	59363350	False	1	0	1	0
chrY	59363450	False	2	1	2	0
chrY	59363500	False	1	0	0	1"""

    f = tmpdir.join("file24h")
    f.write(contents)

    return str(f)


@pytest.fixture
def expected_result():
    pass


@pytest.fixture
def files(file18h, file21h, file24h):
    return [file18h, file21h, file24h]


@pytest.fixture
def output_file(tmpdir):
    f = tmpdir.join("filename")

    return str(f)


def describe_merge_files():
    @pytest.mark.mergepools
    def with_default_arguments(files, output_file, expected_result):

        merge_files(files, output_file)
        result = pd.read_table(output_file, sep="\s+", header=0)
        print(result.to_csv(sep=" "), "actual_result")
        print(result.shape, "actual_result")
        print(expected_result.to_csv(sep=" "), "expected_result")
        print(expected_result.shape, "expected_result")

        assert np.array_equal(result, expected_result)

    @pytest.mark.mergepools
    def with_a_chunksize_of_five(files, output_file, expected_result):

        nb_lines_per_chunk = 5
        merge_files(files, output_file, nb_lines_per_chunk)
        result = pd.read_table(output_file, sep="\s+", header=0)

        assert np.array_equal(result, expected_result)

    @pytest.mark.mergepools
    def with_a_chunksize_of_two(files, output_file, expected_result):

        nb_lines_per_chunk = 2
        merge_files(files, output_file, nb_lines_per_chunk)
        result = pd.read_table(output_file, sep="\s+", header=0)

        assert np.array_equal(result, expected_result)

    @pytest.mark.mergepools
    def with_a_chunksize_of_one(files, output_file, expected_result):

        nb_lines_per_chunk = 1
        merge_files(files, output_file, nb_lines_per_chunk)
        result = pd.read_table(output_file, sep="\s+", header=0)

        assert np.array_equal(result, expected_result)


@pytest.fixture
def expected_result():
    return pd.read_table(
        StringIO(
            """Chromosome	Bin	18h_Island	Exp2_18_Polymerase_II	Exp1_18_Polymerase_II	Exp2_18h_Input	Exp1_18h_Input	21h_Island	Exp2_21_Polymerase_II	Exp1_21_Polymerase_II	Exp2_21h_Input	Exp1_21h_Input	24h_Island	Exp1_24_Polymerase_II	Exp2_24_Polymerase_II	Exp1_24h_Input	Exp2_24h_Input
chr1 10000 1 4 1 2 5 1 3 3 1 1 0 5 0 1 1
chr1 10050 1 2 1 2 2 1 3 1 1 2 1 4 2 4 2
chr1 10100 1 1 2 0 2 1 4 2 0 0 1 4 2 3 5
chr1 10150 1 5 6 2 2 1 6 6 0 8 1 1 2 7 3
chr1 10200 1 6 1 4 2 1 4 8 4 4 1 2 1 2 2
chr1 10250 1 2 2 1 3 1 0 3 2 2 1 7 2 2 6
chr1 10300 1 5 1 4 2 1 3 5 2 5 1 1 4 2 1
chr1 10350 1 5 4 6 4 1 6 2 2 2 1 5 2 4 3
chr1 10400 1 1 3 2 1 1 1 2 1 0 1 1 1 3 2
chrY 59362200 0 1 4 1 0 0 0 0 0 0 0 0 0 0 0
chrY 59362250 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0
chrY 59362300 0 1 0 1 0 0 2 0 1 0 0 0 0 0 0
chrY 59362350 0 0 0 0 0 0 2 0 0 1 0 0 0 0 0
chrY 59362400 0 1 0 1 1 0 0 1 0 3 0 0 0 0 0
chrY 59362450 0 1 1 0 1 0 0 0 0 0 0 0 2 0 1
chrY 59362750 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0
chrY 59362800 0 0 1 1 0 0 0 0 0 0 0 2 0 1 0
chrY 59362900 0 1 0 2 1 0 3 0 0 0 0 1 0 1 1
chrY 59362950 0 1 1 0 0 0 2 2 0 0 0 0 0 0 0
chrY 59363000 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1
chrY 59363100 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0
chrY 59363200 0 0 0 0 0 0 2 3 0 0 0 0 0 0 0
chrY 59363300 0 4 0 1 0 0 1 2 1 1 0 2 2 1 0
chrY 59363350 0 0 0 0 0 0 1 0 1 0 0 1 0 1 0
chrY 59363450 0 3 0 1 0 0 1 1 0 4 0 2 1 2 0
chrY 59363500 0 0 1 1 0 0 0 0 0 0 0 1 0 0 1"""),
        sep="\s+",
        header=0)
