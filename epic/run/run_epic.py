"""Run whole epic pipeline."""

__author__ = "Endre Bakken Stovner https://github.com/endrebak/"
__license__ = "MIT"

from os.path import dirname
from sys import stdout
from itertools import chain
from collections import OrderedDict
from subprocess import call
import logging

import pandas as pd

from natsort import natsorted

from epic.windows.count.count_reads_in_windows import (
    count_reads_in_windows, count_reads_in_windows_paired_end)
from epic.config.genomes import create_genome_size_dict
from epic.statistics.compute_background_probabilites import compute_background_probabilities
from epic.statistics.count_to_pvalue import count_to_pvalue
from epic.statistics.fdr import compute_fdr
from epic.utils.helper_functions import merge_chip_and_input, get_total_number_of_reads, merge_same_files
from epic.windows.cluster.find_islands import find_islands


def run_epic(args):

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

    dfs = count_to_pvalue(merged_dfs, island_enriched_threshold,
                          average_window_readcount, args.number_cores)

    dfs = find_islands(dfs, score_threshold, args)

    logging.info("Done finding islands.")
    logging.info("Concating dfs.")
    df = pd.concat([df for df in dfs if not df.empty])
    logging.info("Labeling island bins.")

    logging.info("Computing FDR.")
    df = compute_fdr(df, nb_chip_reads, nb_input_reads, args)

    if args.store_matrix:
        print_matrixes(chip_merged, input_merged, df, args)

    # Just in case some ints got promoted to float somewhere
    df[["Start", "End", "ChIP", "Input"]] = df[["Start", "End", "ChIP", "Input"
                                                ]].astype(int)
    df.to_csv(stdout, index=False, sep=" ", na_rep="NA")

    return df.reset_index(
    )  # only returns a value to simplify integration tests


def sum_columns(dfs):

    new_dfs = []
    for df in dfs:
        s = df.set_index("Chromosome Bin".split())
        s = s.sum(axis=1)
        s.name = "Count"
        df = s.reset_index()
        new_dfs.append(df)

    return new_dfs


def enriched_bins(df, args):

    df = df.loc[df.FDR < args.false_discovery_rate_cutoff]

    idx_rowdicts = []
    for _, row in df.iterrows():
        for bin in range(
                int(row.Start), int(row.End) + 2, int(args.window_size)):
            idx_rowdicts.append({"Chromosome": row.Chromosome,
                                 "Bin": bin,
                                 "Island": 1})
    islands = pd.DataFrame.from_dict(idx_rowdicts)
    islands.loc[:, "Chromosome"].astype("category")
    islands.loc[:, "Bin"].astype(int)
    return islands.set_index("Chromosome Bin".split())


def print_matrixes(chip, input, df, args):

    outpath = args.store_matrix
    assert len(chip) == len(input)

    islands = enriched_bins(df, args)

    dir = dirname(outpath)
    if dir:
        call("mkdir -p {}".format(dir), shell=True)

    logging.info("Writing data matrix to file: " + outpath)
    for i, (chip_df, input_df) in enumerate(zip(chip, input)):

        chip_df["Chromosome"] = chip_df["Chromosome"].astype("category")
        chip_df["Bin"] = chip_df["Bin"].astype(int)
        chip_df = chip_df.set_index("Chromosome Bin".split())
        chip_df = islands.join(chip_df, how="right")

        input_df["Chromosome"] = input_df["Chromosome"].astype("category")
        input_df["Bin"] = input_df["Bin"].astype(int)
        input_df = input_df.set_index("Chromosome Bin".split())

        dfm = chip_df.join(input_df, how="outer", sort=False).fillna(0)

        if i == 0:
            header, mode = True, "w+"
        else:
            header, mode = False, "a"

        dfm.astype(int).to_csv(outpath,
                               sep=" ",
                               na_rep="NA",
                               header=header,
                               mode=mode,
                               compression="gzip",
                               chunksize=1e6)


def get_island_bins(df, window_size, genome):
    """Finds the enriched bins in a df."""

    # need these chromos because the df might not have islands in all chromos
    chromosomes = natsorted(list(create_genome_size_dict(genome)))

    chromosome_island_bins = {}
    df_copy = df.reset_index(drop=False)
    for chromosome in chromosomes:
        cdf = df_copy.loc[df_copy.Chromosome == chromosome]
        if cdf.empty:
            chromosome_island_bins[chromosome] = set()
        else:
            island_starts_ends = zip(cdf.Start.values.tolist(),
                                     cdf.End.values.tolist())
            island_bins = chain(*[range(
                int(start), int(end), window_size)
                                  for start, end in island_starts_ends])
            chromosome_island_bins[chromosome] = set(island_bins)

    return chromosome_island_bins


def multiple_files_count_reads_in_windows(bed_files, args):
    """Use count_reads on multiple files and store result in dict.

    Untested since does the same thing as count reads."""

    bed_windows = OrderedDict()
    for bed_file in bed_files:
        logging.info("Binning " + bed_file)
        if args.paired_end:
            chromosome_dfs = count_reads_in_windows_paired_end(bed_file, args)
        else:
            chromosome_dfs = count_reads_in_windows(bed_file, args)
        bed_windows[bed_file] = chromosome_dfs

    return bed_windows


def _merge_files(windows, nb_cpu):
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
