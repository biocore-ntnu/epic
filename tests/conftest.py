import pytest

from io import StringIO
from collections import namedtuple

import pandas as pd

from epic.config.genomes import create_genome_size_dict

__author__ = "Endre Bakken Stovner https://github.com/endrebak/"
__license__ = "MIT"

MockNamespace = namedtuple(
    "MockNamespace",
    ["number_cores", "genome", "keep_duplicates", "window_size",
     "fragment_size", "paired_end", "gaps_allowed",
     "false_discovery_rate_cutoff", "effective_genome_fraction", "store_matrix",
     "bigwig", "bed", "chromosome_sizes",
     "treatment", "control", "outfile", "chip_bigwig", "input_bigwig", "log2fc_bigwig", "individual_log2fc_bigwigs"])

# egs = 2290813547.4  # this is the effective genome size used by the original sicer for hg19
egs = 2510169508.0003734
gsd = create_genome_size_dict("hg19")


@pytest.fixture(scope="session")
def args_200_fast():
    return MockNamespace(25, "hg19", False, 200, 150, False, 3, 0.05, egs, False, False, False,
                         gsd, ["examples/test.bed"],
                         ["examples/control.bed"], None, False, False, False, False)


@pytest.fixture(scope="session")
def args_200():
    return MockNamespace(1, "hg19", False, 200, 150, False, 3, 0.05, egs,
                         False, False, False, gsd,
                         ["examples/test.bed"], ["examples/control.bed"], None, False, False, False, False)


# @pytest.fixture(scope="session")
# def args_50():
#     return MockNamespace(1, "hg19", False, 50, 150, False, 3, 0.05, egs, False, False, False,
#                          False, gsd, ["examples/test.bed"],
#                          ["examples/control.bed"], None)


@pytest.fixture(scope="session")
def expected_result_example_input():
    return pd.read_table("examples/expected_results_log2fc.csv", sep=" ", skiprows=1)


@pytest.fixture(scope="session")
def epic_overlap_files():
    return ["examples/epic-overlaps/0h_H3K4me3.regions",
"examples/epic-overlaps/3h_H3K4me3.regions",
"examples/epic-overlaps/6h_H3K4me3.regions",
"examples/epic-overlaps/9h_H3K4me3.regions",
"examples/epic-overlaps/12h_H3K4me3.regions",
"examples/epic-overlaps/15h_H3K4me3.regions",
"examples/epic-overlaps/18h_H3K4me3.regions",
"examples/epic-overlaps/21h_H3K4me3.regions",
"examples/epic-overlaps/24h_H3K4me3.regions"]


@pytest.fixture(scope="session")
def epic_overlap_files_more_chromos():
    return ["examples/epic-overlaps/more_chromos0h_H3K4me3.regions",
"examples/epic-overlaps/more_chromos3h_H3K4me3.regions",
"examples/epic-overlaps/more_chromos6h_H3K4me3.regions",
"examples/epic-overlaps/more_chromos9h_H3K4me3.regions",
"examples/epic-overlaps/more_chromos12h_H3K4me3.regions",
"examples/epic-overlaps/more_chromos15h_H3K4me3.regions",
"examples/epic-overlaps/more_chromos18h_H3K4me3.regions",
"examples/epic-overlaps/more_chromos21h_H3K4me3.regions",
"examples/epic-overlaps/more_chromos24h_H3K4me3.regions"]


@pytest.fixture(scope="session")
def epic_overlap_intermediate_nucleotide_matrixes():
    return ["examples/epic-overlaps/0h_H3K4me3.matrix",
            "examples/epic-overlaps/3h_H3K4me3.matrix",
            "examples/epic-overlaps/6h_H3K4me3.matrix"]


@pytest.fixture(scope="session")
def epic_overlap_intermediate_region_matrixes():
    return ["examples/epic-overlaps/0h_H3K4me3_region.matrix",
            "examples/epic-overlaps/3h_H3K4me3_region.matrix",
            "examples/epic-overlaps/6h_H3K4me3_region.matrix"]



df1_no_enriched = u"""Chromosome Bin Enriched_Gene1_KO.matrix.gz data/align/Gene1_KO_ChIP_1.bed data/align/Gene1_KO_ChIP_2.bed data/align/Gene1_KO_ChIP_3.bed data/align/Gene1_KO_Input_1.bed data/align/Gene1_KO_Input_2.bed data/align/Gene1_KO_Input_3.bed
chr1 10000 0 0 2 0 5 1 3
chr1 10200 0 0 0 0 2 0 0
chr1 10400 0 0 0 0 1 0 1
chr1 10600 0 0 0 0 2 0 3
chr1 10800 0 0 0 0 6 8 5
chr1 11000 0 0 0 0 7 4 8
chr1 11200 0 0 0 0 4 0 1
chr1 11400 0 1 0 0 2 0 0
chr1 12200 0 0 0 0 5 0 1"""

df2_no_enriched = u"""Chromosome Bin Enriched_Gene2.matrix.gz data/align/Gene2_KO_ChIP_1.bed data/align/Gene2_KO_ChIP_2.bed data/align/Gene2_KO_ChIP_3.bed data/align/Gene2_KO_Input_1.bed data/align/Gene2_KO_Input_2.bed data/align/Gene2_KO_Input_3.bed
chr1 10000 0 0 3 0 3 2 1
chr1 10200 0 0 1 1 2 2 1
chr1 10400 0 0 0 0 2 1 1
chr1 10600 0 0 0 0 2 2 1
chr1 10800 0 0 1 0 0 11 4
chr1 11000 0 0 0 0 3 12 12
chr1 11200 0 0 0 0 1 0 3
chr1 11400 0 0 0 0 2 0 0
chr1 12200 0 0 0 0 3 6 8"""

df3_no_enriched = u"""Chromosome Bin Enriched_WT.matrix.gz data/align/WT_ChIP_1.bed data/align/WT_ChIP_2.bed data/align/WT_ChIP_3.bed data/align/WT_Input_1.bed data/align/WT_Input_2.bed data/align/WT_Input_3.bed
chr1 10000 0 0 0 1 2 2 0
chr1 10200 0 0 1 0 2 0 0
chr1 10400 0 0 0 0 2 3 4
chr1 10600 0 0 0 0 0 0 1
chr1 10800 0 0 0 0 2 1 15
chr1 11000 0 0 3 0 1 3 10
chr1 11200 0 0 0 0 0 1 0
chr1 11400 0 0 2 1 0 0 0
chr1 12200 0 0 0 0 3 0 2"""

df3_enriched = u"""Chromosome Bin Enriched_WT.matrix.gz data/align/WT_ChIP_1.bed data/align/WT_ChIP_2.bed data/align/WT_ChIP_3.bed data/align/WT_Input_1.bed data/align/WT_Input_2.bed data/align/WT_Input_3.bed
chr1 10000 0 0 0 1 2 2 0
chr1 10200 1 0 1 0 2 0 0
chr1 10400 1 0 0 0 2 3 4
chr1 10600 1 0 0 0 0 0 1
chr1 10800 0 0 0 0 2 1 15
chr1 11000 0 0 3 0 1 3 10
chr1 11200 0 0 0 0 0 1 0
chr1 11400 0 0 2 1 0 0 0
chr1 12200 0 0 0 0 3 0 2"""


@pytest.fixture(scope="session")
def merge_matrixes_dfs_no_enriched():

    dfs = [pd.read_table(StringIO(c), sep=" ", header=0, index_col=[0, 1]) for c in [df1_no_enriched, df2_no_enriched, df3_no_enriched]]

    names = "Gene1_KO.matrix.gz Gene2_KO.matrix.gz WT.matrix.gz".split()

    dfs = {k: v for k, v in zip(names, dfs)}

    return dfs



@pytest.fixture(scope="session")
def merge_matrixes_dfs_one_enriched():

    dfs = [pd.read_table(StringIO(c), sep=" ", header=0, index_col=[0, 1]) for c in [df1_no_enriched, df2_no_enriched, df3_enriched]]

    names = "Gene1_KO.matrix.gz Gene2_KO.matrix.gz WT.matrix.gz".split()

    dfs = {k: v for k, v in zip(names, dfs)}

    return dfs
