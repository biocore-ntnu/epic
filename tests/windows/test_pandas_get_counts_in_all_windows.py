import pytest

from io import StringIO
import pandas as pd
from numpy import allclose, array_equal

from epic.windows.count.pandas_count_reads_in_windows import pandas_count_reads_in_windows


@pytest.mark.integration
def test_get_count_in_all_windows(expected_result, input_bed_file):

    genome, fragment_size, window_size = "hg19", 150, 200
    keep_duplicates = False

    actual_result = pandas_count_reads_in_windows(
        input_bed_file, genome, fragment_size, window_size, keep_duplicates, 1)

    actual_result = pd.concat(actual_result)

    print("actual_result", actual_result)
    print("expected_result", expected_result)

    assert array_equal(actual_result, expected_result)


@pytest.fixture
def expected_result():

    return pd.read_table(
        StringIO(u"""Count  Chromosome  Bin
2  chr1  39036800
1  chr1  73781000
1  chr1  90059800
1  chr3  55648200
1  chr7  20246600
1  chr7  91135000
2  chr13  100938400
1  chr19  43528800
1  chr19  47108800"""),
        sep=r"\s+")


@pytest.fixture
def input_bed_file(tmpdir):
    bed_file = tmpdir.join("read_file.bed")
    bed_file.write("""
chr1	73781101	73781126	U0	0	+
chr7	91135110	91135135	U0	0	+
chr13	100938475	100938500	U0	0	+
chr3	55648137	55648162	U0	0	+
chr7	20246668	20246693	U0	0	+
chr1	39036822	39036847	U0	0	+
chr19	47109000	47109025	U0	0	-
chr19	43528773	43528798	U0	0	+
chr1	90059861	90059886	U0	0	-
chr1	39036870	39036900	U0	0	-
chr13	100938476	100938500	U0	0	+""")
    return str(bed_file)
