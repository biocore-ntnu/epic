import pytest

from io import StringIO
import pandas as pd
from numpy import allclose, array_equal, int32

from epic.windows.count.count_reads_in_windows import count_reads_in_windows, _count_reads_in_windows_paired_end

import argparse
parser = argparse.ArgumentParser()

parser.add_argument(
    '--number-cores',
    '-cpu',
    required=False,
    default=1,
    type=int,
    help=
    '''Number of cpus to use. Can use at most one per chromosome. Default: 1.''')

parser.add_argument('--genome',
                    '-gn',
                    required=False,
                    default="hg19",
                    type=str,
                    help='''Which genome to analyze. Default: hg19.''')

parser.add_argument(
    '--keep-duplicates',
    '-k',
    required=False,
    default=False,
    type=bool,
    help=
    '''Keep reads mapping to the same position on the same strand within a library. Default is to remove all but the first duplicate.
                   ''')

parser.add_argument(
    '--window-size',
    '-w',
    required=False,
    default=200,
    type=int,
    help=
    '''Size of the windows to scan the genome. WINDOW_SIZE is the smallest possible island. Default 200.
                   ''')

parser.add_argument(
    '--fragment-size',
    '-fs',
    required=False,
    default=150,
    type=int,
    help=
    '''Size of the sequenced fragment. The center of the the fragment will be taken as half the fragment size. Default 150.
                   ''')

parser.add_argument('--paired-end',
                    '-pe',
                    required=False,
                    default=False,
                    type=bool,
                    help='''(Bam-only) Use paired-end reads.''')

@pytest.mark.integration
def test_count_reads_in_windows(expected_result, input_bed_file):

    args = parser.parse_args(['-fs', '150'])
    args.genome, args.fragment_size, args.window_size = "hg19", 150, 200
    args.keep_duplicates = False
    args.number_cores = 1

    actual_result = count_reads_in_windows(input_bed_file, args)

    actual_result = pd.concat(actual_result)

    print(actual_result)
    print(expected_result)

    assert array_equal(actual_result, expected_result)


@pytest.fixture
def expected_result():

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
    return df


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
def paired_end(tmpdir):

    bed_file = tmpdir.join("paired_end_file.bed")
    bed_file.write(
        """chr1	33076	33165	chr1	33076	33165	2PJ3LS1:183:C5RR7ACXX:7:1307:11171:75017	0	+	-
chr1	33076	33165	chr1	33076	33165	2PJ3LS1:183:C5RR7ACXX:8:1115:18054:85095	0	+	-
chr1	33660	33744	chr1	33660	33744	2PJ3LS1:183:C5RR7ACXX:2:2201:20777:74838	0	+	-
chr1	33660	33744	chr1	33660	33744	2PJ3LS1:183:C5RR7ACXX:3:1213:9620:83082	0	+	-
chr1	38474	38571	chr1	38474	38571	2PJ3LS1:183:C5RR7ACXX:4:1216:1576:83235	0	+	-
chr1	41985	42043	chr1	42040	42131	2PJ3LS1:183:C5RR7ACXX:3:1214:11778:16074	1	+	-
chr1	42223	42324	chr1	42285	42386	2PJ3LS1:183:C5RR7ACXX:3:1215:18097:73129	0	+	-""")
    return str(bed_file)


# @pytest.mark.integration
def test_count_reads_in_windows_paired_end(paired_end):

    genome, fragment_size, window_size = "hg19", 150, 200
    keep_duplicates = False

    # result = _count_reads_in_windows_paired_end(paired_end, genome,
    #                                            fragment_size, window_size,
    #                                            keep_duplicates, "chr1")

    # print(result)
    "Paired end support not implemented yet."
    # assert 0
    # assert result == expected_result
