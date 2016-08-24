"""Run whole epic pipeline."""

__author__ = "Endre Bakken Stovner https://github.com/endrebak/"
__license__ = "MIT"

from os.path import dirname, join, basename
from sys import stdout
from itertools import chain
from collections import OrderedDict
from subprocess import call
import logging

import pandas as pd

from natsort import natsorted
from joblib import Parallel, delayed

from epic.windows.count.count_reads_in_windows import (
    count_reads_in_windows, count_reads_in_windows_paired_end)
from epic.config.genomes import create_genome_size_dict
from epic.statistics.compute_background_probabilites import compute_background_probabilities
from epic.statistics.count_to_pvalue import count_to_pvalue
from epic.statistics.fdr import compute_fdr
from epic.utils.helper_functions import merge_chip_and_input, get_total_number_of_reads, merge_same_files
from epic.windows.cluster.find_islands import find_islands
from epic.matrixes.matrixes import write_matrix_files


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

    if args.store_matrix or args.individual_bedgraph or args.bedgraph:
        write_matrix_files(chip_merged, input_merged, df, args)

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


def bedgraph(matrix, args):
    "Create a bedgraph file for ChIP and Input."

    outfolder = args.bedgraph

    call("mkdir -p {}".format(outfolder), shell=True)

    chip_file = join(outfolder, "treatment.bedgraph")
    input_file = join(outfolder, "input.bedgraph")

    c_sum = matrix[args.treatment].sum(1)
    c = c_sum[c_sum > 0]
    c.to_csv(chip_file, sep="\t")

    i_sum = matrix[args.control].sum(1)
    i = i_sum[i_sum > 0]
    i.to_csv(input_file, sep="\t")


def _individual_bedgraphs(matrix, name, outfolder):

    base = basename(name).split(".")

    if len(base) > 2:
        base = "".join(base[:-2])
    else:
        base = "".join(base[:-1])

    outfile = join(outfolder, base + ".bedgraph")
    s = matrix[name]
    nonzeroes_only = s[s != 0]
    nonzeroes_only.to_csv(outfile, sep="\t")


def individual_bedgraphs(matrix, args):
    "Create a bedgraph file for each file used."

    outfolder = args.individual_bedgraph

    call("mkdir -p {}".format(outfolder), shell=True)

    for treatment_file in args.treatment:
        _individual_bedgraphs(matrix, treatment_file, outfolder)

    for control_file in args.control:
        _individual_bedgraphs(matrix, control_file, outfolder)


def _create_matrixes(chromosome, chip, input, islands):

    chip_df = get_chromosome_df(chromosome, chip)
    input_df = get_chromosome_df(chromosome, input)

    chip_df["Chromosome"] = chip_df["Chromosome"].astype("category")
    chip_df["Bin"] = chip_df["Bin"].astype(int)
    chip_df = chip_df.set_index("Chromosome Bin".split())
    chip_df = islands.join(chip_df, how="right")

    input_df["Chromosome"] = input_df["Chromosome"].astype("category")
    input_df["Bin"] = input_df["Bin"].astype(int)
    input_df = input_df.set_index("Chromosome Bin".split())

    dfm = chip_df.join(input_df, how="outer", sort=False).fillna(0)

    return dfm


def create_matrixes(chip, input, df, args):

    "Creates matrixes which can be written to file as is (matrix) or as bedGraph."

    chip = put_dfs_in_chromosome_dict(chip)
    input = put_dfs_in_chromosome_dict(input)
    all_chromosomes = natsorted(set(list(chip.keys()) + list(input.keys())))

    islands = enriched_bins(df, args)

    logging.info("Creating matrixes from count data.")
    dfms = Parallel(n_jobs=args.number_cores)(
        delayed(_create_matrixes)(chromosome, chip, input, islands)
        for chromosome in all_chromosomes)

    return dfms


def print_matrixes(matrixes, args):

    outpath = args.store_matrix

    dir = dirname(outpath)
    if dir:
        call("mkdir -p {}".format(dir), shell=True)

    logging.info("Writing data matrix to file: " + outpath)
    for i, df in enumerate(matrixes):

        if i == 0:
            header, mode = True, "w+"
        else:
            header, mode = False, "a"

        df.astype(int).to_csv(outpath,
                              sep=" ",
                              na_rep="NA",
                              header=header,
                              mode=mode,
                              compression="gzip",
                              chunksize=1e6)


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
