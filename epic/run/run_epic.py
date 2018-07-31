from __future__ import print_function

"""Run whole epic pipeline."""

__author__ = "Endre Bakken Stovner https://github.com/endrebak/"
__license__ = "MIT"

from os.path import dirname, join, basename
from sys import stdout, argv
from itertools import chain
from collections import OrderedDict
from subprocess import call
import logging

from argparse import Namespace
import pandas as pd
from numpy import log2

from typing import Iterable
from natsort import natsorted
from joblib import Parallel, delayed

from epic.windows.count.count_reads_in_windows import (
    count_reads_in_windows, count_reads_in_windows_paired_end)
from epic.statistics.compute_background_probabilites import compute_background_probabilities
from epic.statistics.count_to_pvalue import count_to_pvalue
from epic.statistics.fdr import compute_fdr
from epic.utils.helper_functions import merge_chip_and_input, get_total_number_of_reads, merge_same_files
from epic.windows.cluster.find_islands import find_islands
from epic.matrixes.matrixes import write_matrix_files


def run_epic(args):
    # type: (Namespace) -> pd.DataFrame

    chip_windows = multiple_files_count_reads_in_windows(args.treatment, args)
    input_windows = multiple_files_count_reads_in_windows(args.control, args)

    chip_merged = _merge_files(chip_windows.values(), args.number_cores)

    input_merged = _merge_files(input_windows.values(), args.number_cores)

    chip_merged_sum = sum_columns(chip_merged)
    input_merged_sum = sum_columns(input_merged)

    nb_chip_reads = get_total_number_of_reads(chip_merged_sum)

    nb_input_reads = get_total_number_of_reads(input_merged_sum)

    merged_dfs = merge_chip_and_input(chip_merged_sum, input_merged_sum,
                                      args.number_cores)

    score_threshold, island_enriched_threshold, average_window_readcount = \
        compute_background_probabilities(nb_chip_reads, args)

    dfs = []                    # type: Iterable[pd.DataFrame]
    dfs = count_to_pvalue(merged_dfs, island_enriched_threshold,
                          average_window_readcount, args.number_cores)

    dfs = find_islands(dfs, score_threshold, args)

    logging.info("Done finding islands.")
    logging.info("Concating dfs.")
    df = pd.concat([df for df in dfs if not df.empty])
    logging.info("Labeling island bins.")

    logging.info("Computing FDR.")
    df = compute_fdr(df, nb_chip_reads, nb_input_reads, args)

    # Just in case some ints got promoted to float somewhere
    df[["Start", "End", "ChIP", "Input"]] = df[["Start", "End", "ChIP", "Input"
                                                ]].astype(int)
    # redundancy in below code
    outfile = args.outfile if args.outfile else stdout
    if args.outfile:
        with open(outfile, "w+") as h:
            print("# epic " + " ".join(argv[1:]), file=h)
    else:
        print("# epic " + " ".join(argv[1:]), file=stdout)

    df.to_csv(outfile, index=False, sep=" ", na_rep="NA", mode="a")

    if args.bed:
        df_to_bed(df).to_csv(args.bed, header=False, index=False, sep="\t")

    if (args.store_matrix or args.bigwig or args.chip_bigwig or args.input_bigwig or args.log2fc_bigwig or args.individual_log2fc_bigwigs):
        write_matrix_files(chip_merged, input_merged, df, args)

    return df.reset_index(drop=True) # only returns a value to simplify integration tests


def df_to_bed(df):
    # type: (pd.DataFrame) -> pd.DataFrame

    # '''Chromosome Start End ChIP Input Score Fold_change P FDR
    # chr5 53000 55399 121 13 77.6075622841774 13.655736573980159 6.040968494897508e-92 1.9241805908359603e-91\''

    regions = df["Chromosome Start End FDR".split()].copy()
    regions.insert(4, "Score", log2(df.ChIP/df.Input) * 100)
    regions.loc[regions.Score > 1000, "Score"] = 1000
    regions.loc[regions.Score < 0, "Score"] = 0
    regions.insert(5, "Strand", ".")

    return regions




def sum_columns(dfs):
    # type: (Iterable[pd.DataFrame]) -> List[pd.DataFrame]

    new_dfs = []
    for df in dfs:
        s = df.set_index("Chromosome Bin".split())
        s = s.sum(axis=1)
        s.name = "Count"
        df = s.reset_index()
        new_dfs.append(df)

    return new_dfs


def multiple_files_count_reads_in_windows(bed_files, args):
    # type: (Iterable[str], Namespace) -> OrderedDict[str, List[pd.DataFrame]]
    """Use count_reads on multiple files and store result in dict.

    Untested since does the same thing as count reads."""

    bed_windows = OrderedDict() # type: OrderedDict[str, List[pd.DataFrame]]
    for bed_file in bed_files:
        logging.info("Binning " + bed_file)
        if ".bedpe" in bed_file:
            chromosome_dfs = count_reads_in_windows_paired_end(bed_file, args)
        else:
            chromosome_dfs = count_reads_in_windows(bed_file, args)
        bed_windows[bed_file] = chromosome_dfs

    return bed_windows


def _merge_files(windows, nb_cpu):
    # type: (Iterable[pd.DataFrame], int) -> pd.DataFrame
    """Merge lists of chromosome bin df chromosome-wise.

    windows is an OrderedDict where the keys are files, the values are lists of
    dfs, one per chromosome.

    Returns a list of dataframes, one per chromosome, with the collective count
    per bin for all files.

    TODO: is it faster to merge all in one command?
    """

    # windows is a list of chromosome dfs per file
    windows = iter(windows)  # can iterate over because it is odict_values
    merged = next(windows)

    # if there is only one file, the merging is skipped since the windows is used up
    for chromosome_dfs in windows:
        # merge_same_files merges the chromosome files in parallel
        merged = merge_same_files(merged, chromosome_dfs, nb_cpu)

    return merged
