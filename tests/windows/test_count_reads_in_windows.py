import pytest
# from collections import namedtuple

from io import StringIO
import pandas as pd
from numpy import allclose, array_equal, int32

from collections import namedtuple

from epic.windows.count.count_reads_in_windows import (
    count_reads_in_windows, _count_reads_in_windows_paired_end)


@pytest.mark.integration
def test_count_reads_in_windows(expected_result, input_bed_file, args_200_fast):

    actual_result = count_reads_in_windows(input_bed_file, args_200_fast)

    actual_result = pd.concat(actual_result)

    print(actual_result.values, "actual")
    print(expected_result.values, "expected")

    assert array_equal(actual_result, expected_result)


@pytest.fixture
def input_bed_file(tmpdir):
    bed_file = tmpdir.join("read_file.bed")
    bed_file.write("""chr1	39036822	39036847	U0	0	+
chr1	39036870	39036900	U0	0	-
chr1	73781101	73781126	U0	0	+
chr1	90059861	90059886	U0	0	-
chr3	55648137	55648162	U0	0	+
chr7	20246668	20246693	U0	0	+
chr7	91135110	91135135	U0	0	+
chr13	100938475	100938500	U0	0	+
chr13	115816130	115816155	U0	0	+
chr19	43528773	43528798	U0	0	+
chr19	47109000	47109025	U0	0	-""")
    return str(bed_file)


@pytest.fixture
def expected_result(input_bed_file):

    df = pd.read_table(
        StringIO(u"""Count Chromosome Bin
    2 chr1 39036800
    1 chr1 73781000
    1 chr1 90059800
    1 chr3 55648200
    1 chr7 20246600
    1 chr7 91135000
    1 chr13 100938400
    1 chr19 43528800
    1 chr19 47108800"""),
        sep=r"\s+",
        dtype={"Count": int32,
               "Bin": int32})
    df.columns = [input_bed_file, "Chromosome", "Bin"]
    return df


@pytest.fixture
def paired_end(tmpdir):

    bed_file = tmpdir.join("paired_end_file.bed")
    bed_file.write(
        """chr1	33076	33165	chr1	33076	33165	2PJ3LS1:183:C5RR7ACXX:7:1307:11171:75017	0	+	-
chr1	33076	33165	chr1	33076	33165	2PJ3LS1:183:C5RR7ACXX:8:1115:18054:85095	0	+	-
chr1	33660	33744	chr1	33660	33744	2PJ3LS1:183:C5RR7ACXX:2:2201:20777:74838	0	+	-
chr1	33660	33744	chr1	33660	33744	2PJ3LS1:183:C5RR7ACXX:3:1213:9620:83082	0	+	-
chr1	38474	38571	chr1	38474	38571	2PJ3LS1:183:C5RR7ACXX:4:1216:1576:83235	0	+	-
chr1	41985	42043	chr1	42040	42131	2PJ3LS1:183:C5RR7ACXX:3:1214:11778:16074	1	+	-
chr1	23025	23045	chr1	43025	43045	2PJ3LS1:183:C5RR7ACXX:8:1115:18054:85095	0	+	-
chr1	42223	42324	chr1	42285	42386	2PJ3LS1:183:C5RR7ACXX:3:1215:18097:73129	0	+	-""")
    return str(bed_file)


@pytest.fixture
def expected_result_pe(paired_end):
    df = pd.read_table(
        StringIO(u"""Count Chromosome    Bin
2       chr1  33000
1       chr1  33600
1       chr1  38400
1       chr1  42000
1       chr1  42200"""),
        sep="\s+",
        dtype={"Count": int32,
               "Bin": int32})

    df.columns = [paired_end, "Chromosome", "Bin"]
    return df


@pytest.mark.integration
def test_count_reads_in_windows_paired_end(paired_end, expected_result_pe):


    MockNamespace = namedtuple("MockNamespace", ["keep_duplicates", "window_size"])
    args = MockNamespace(False, 200)

    result = _count_reads_in_windows_paired_end(paired_end, args, 249250621,
                                                "chr1")

    print(result)
    print(result.dtypes)
    print(expected_result_pe)
    print(expected_result_pe.dtypes)
    assert array_equal(result, expected_result_pe)
