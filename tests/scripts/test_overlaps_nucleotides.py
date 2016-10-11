# import pytest

# import pandas as pd
# import numpy as np

# import pkg_resources, os
# from collections import defaultdict

# from io import StringIO

# from epic.scripts.overlaps.overlaps_nucleotides import (_create_overlap_matrix_nucleotides,
#                                                         _overlap_matrix_nucleotides,
#                                                         files_to_chromosome_coverage,
#                                                         _nucleotide_overlaps_per_file)

# from epic.config.genomes import (create_genome_size_dict,
#                                  get_effective_genome_length)

# __author__ = "Endre Bakken Stovner https://github.com/endrebak/"
# __license__ = "MIT"

# @pytest.mark.current
# def test_nucleotide_overlaps_per_file(bed_file, extended_rles, expected_result_overlaps_per_file):

#     overlaps = _nucleotide_overlaps_per_file(bed_file, extended_rles, 1)
#     overlaps = overlaps.sort_values(["Main", "Other"])
#     print(overlaps, "ove " * 5)
#     print(expected_result_overlaps_per_file, "exp " * 5)
#     assert np.allclose(expected_result_overlaps_per_file.Overlaps, overlaps.Overlaps)


# @pytest.fixture
# def expected_result_overlaps_per_file():
#     contents = """Chromosome                    Main                    Other  Overlaps
# 0        chr1  more_chromos0h_H3K4me3  more_chromos21h_H3K4me3   2477115
# 1        chr1  more_chromos0h_H3K4me3   more_chromos3h_H3K4me3   2591874
# 2        chr1  more_chromos0h_H3K4me3   more_chromos6h_H3K4me3   2516301
# 3        chr1  more_chromos0h_H3K4me3  more_chromos18h_H3K4me3   2482713
# 4        chr1  more_chromos0h_H3K4me3  more_chromos15h_H3K4me3   2524698
# 5        chr1  more_chromos0h_H3K4me3  more_chromos12h_H3K4me3   2499507
# 6        chr1  more_chromos0h_H3K4me3   more_chromos9h_H3K4me3   2544291
# 7        chr1  more_chromos0h_H3K4me3  more_chromos24h_H3K4me3   2429532
# 8        chr9  more_chromos0h_H3K4me3  more_chromos21h_H3K4me3    220616
# 9        chr9  more_chromos0h_H3K4me3   more_chromos3h_H3K4me3    191840
# 10       chr9  more_chromos0h_H3K4me3   more_chromos6h_H3K4me3    236203
# 11       chr9  more_chromos0h_H3K4me3  more_chromos18h_H3K4me3    227810
# 12       chr9  more_chromos0h_H3K4me3  more_chromos15h_H3K4me3    151074
# 13       chr9  more_chromos0h_H3K4me3  more_chromos12h_H3K4me3    117502
# 14       chr9  more_chromos0h_H3K4me3   more_chromos9h_H3K4me3    161865
# 15       chr9  more_chromos0h_H3K4me3  more_chromos24h_H3K4me3    211024
# 16       chrX  more_chromos0h_H3K4me3  more_chromos21h_H3K4me3    572883
# 17       chrX  more_chromos0h_H3K4me3   more_chromos3h_H3K4me3    573682
# 18       chrX  more_chromos0h_H3K4me3   more_chromos6h_H3K4me3    567290
# 19       chrX  more_chromos0h_H3K4me3  more_chromos18h_H3K4me3    583270
# 20       chrX  more_chromos0h_H3K4me3  more_chromos15h_H3K4me3    588863
# 21       chrX  more_chromos0h_H3K4me3  more_chromos12h_H3K4me3    600049
# 22       chrX  more_chromos0h_H3K4me3   more_chromos9h_H3K4me3    565692
# 23       chrX  more_chromos0h_H3K4me3  more_chromos24h_H3K4me3    557702
# 24       chrY  more_chromos0h_H3K4me3  more_chromos21h_H3K4me3       799
# 25       chrY  more_chromos0h_H3K4me3   more_chromos3h_H3K4me3       799
# 26       chrY  more_chromos0h_H3K4me3   more_chromos6h_H3K4me3       799
# 27       chrY  more_chromos0h_H3K4me3  more_chromos18h_H3K4me3       799
# 28       chrY  more_chromos0h_H3K4me3  more_chromos15h_H3K4me3       799
# 29       chrY  more_chromos0h_H3K4me3  more_chromos12h_H3K4me3       799
# 30       chrY  more_chromos0h_H3K4me3   more_chromos9h_H3K4me3       799
# 31       chrY  more_chromos0h_H3K4me3  more_chromos24h_H3K4me3       799"""
#     df = pd.read_table(StringIO(contents), sep="\s+", header=0, index_col=0).sort_values(["Main", "Other"])
#     return df

# @pytest.fixture
# def expected_result_one_file_overlap():
#     return defaultdict(int, {0: 399544596, 2: 77600, 3: 107400, 4: 135000, 5: 146800, 6: 244600, 7: 250800, 8: 269800, 9: 496200, 10: 3422800})

# @pytest.fixture
# def bed_file(epic_overlap_files_more_chromos):
#     return epic_overlap_files_more_chromos[0]


# @pytest.fixture
# def extended_rles(epic_overlap_files_more_chromos):

#     return files_to_chromosome_coverage(epic_overlap_files_more_chromos, 1)


# @pytest.fixture
# def overlaps(extended_rles, bed_file):

#     return _create_overlap_matrix_nucleotides(bed_file, extended_rles)



# @pytest.mark.unit
# def test__compute_nucleotide_overlap(extended_rles, expected_result_one_file_overlap, bed_file):

#     result = _overlap_matrix_nucleotides(bed_file, extended_rles, 1)
#     assert expected_result_one_file_overlap == result
