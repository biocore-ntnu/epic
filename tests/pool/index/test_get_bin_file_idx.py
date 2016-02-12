from collections import OrderedDict
import pytest

from epic.pool.index.get_bin_file_idx import index_bin_files


@pytest.fixture
def input_data(epic_bin_file_path, epic_bin_file_path2):
    return [epic_bin_file_path, epic_bin_file_path2]


@pytest.fixture
def expected_result(epic_bin_file_path, epic_bin_file_path2):
    return OrderedDict([('chr1', [OrderedDict([(epic_bin_file_path, (1, 4)), (
        epic_bin_file_path2, (1, 2))]), OrderedDict([(epic_bin_file_path, (
            5, 7)), (epic_bin_file_path2, (3, 4))]), OrderedDict(
                [(epic_bin_file_path, (8, 10)), (epic_bin_file_path2, (5, 7))
                 ])]), ('chr2', [OrderedDict([(epic_bin_file_path, (11, 14)), (
                     epic_bin_file_path2, (8, 8))]), OrderedDict([(
                         epic_bin_file_path, (15, 17)), (epic_bin_file_path2, (
                             9, 9))]), OrderedDict([(epic_bin_file_path, (
                                 18, 20)), (epic_bin_file_path2, (10, 10))])]),
                        ('chr3', [OrderedDict([(epic_bin_file_path2, (11, 11))
                                               ])])])


@pytest.mark.integration
def test_index_bin_files(input_data, expected_result):

    result = index_bin_files(input_data, 3)
    print(result)

    assert result == expected_result


@pytest.fixture
def epic_bin_file_path(tmpdir):
    epic_bin_file = tmpdir.join("df1_file.txt")
    epic_bin_file.write(
        """Chromosome	Bin	Island	Exp1_3_Polymerase_II	Exp2_3_Polymerase_II	Exp2_3h_Input	Exp1_3h_Input
chr1	10000	False	3	0	4	4
chr1	10050	False	1	1	5	4
chr1	10100	False	3	1	1	3
chr1	10150	False	4	1	3	5
chr1	10200	False	3	0	6	1
chr1	10250	False	1	0	2	1
chr1	10300	False	2	0	3	3
chr1	10350	True	5	1	6	4
chr1	10400	True	3	1	1	1
chr1	10450	True	2	0	1	2
chr2	10100	False	1	0	0	1
chr2	10250	False	1	0	0	0
chr2	10350	False	1	1	0	0
chr2	10450	False	0	2	0	0
chr2	10600	False	0	2	1	0
chr2	11350	False	0	1	1	0
chr2	11600	False	1	0	0	1
chr2	11650	False	0	1	0	0
chr2	11700	False	0	1	1	0
chr2	11750	False	1	0	1	0""")
    return str(epic_bin_file)


@pytest.fixture
def epic_bin_file_path2(tmpdir):
    epic_bin_file = tmpdir.join("df2_file.txt")
    epic_bin_file.write(
        """Chromosome	Bin	Island	Exp1_21_Polymerase_II	Exp2_21_Polymerase_II	Exp2_21h_Input	Exp1_21h_Input
chr1	10000	False	3	0	4	4
chr1	10150	False	4	1	3	5
chr1	10200	False	3	0	6	1
chr1	10300	False	2	0	3	3
chr1	10350	True	5	1	6	4
chr1	10400	True	3	1	1	1
chr1	10450	True	2	0	1	2
chr2	10100	False	1	0	0	1
chr2	11600	False	1	0	0	1
chr2	11700	False	0	1	1	0
chr3	11750	False	1	0	1	0""")
    return str(epic_bin_file)

# @pytest.fixture
# def expected_result():
#     od = OrderedDict()
#     od["df1_file.txt"] = OrderedDict([('chr1', [(0, 3), (4, 6), (7, 9)]), (
#         'chr2', [(10, 13), (14, 16), (17, 19)])])
#     od["df2_file.txt"] = OrderedDict([('chr1', [(0, 1), (2, 3), (4, 6)]), (
#         'chr2', [(7, 7), (7, 8), (8, 9)]), ('chr3', [(10, 10)])])
#     return od
