from functools import partial
from logging import info, warning

import pandas as pd
from epic.config.cache_settings import MEMORY
from epic.config.genomes import create_genome_size_dict
from joblib import Parallel, delayed
from numpy import int32

from natsort import natsorted

from epic.windows.count.remove_out_of_bounds_bins import remove_out_of_bounds_bins


# @MEMORY.cache(verbose=0, ignore=["nb_cpus"])
def pandas_count_reads_in_windows(input_file, genome, fragment_size,
                                  window_size, keep_duplicates, nb_cpus):
    """Make table of window counts from bed file."""

    chromosome_size_dict = create_genome_size_dict(genome)

    info("Reading input file {input_file}.".format(**vars()))
    df = pd.read_table(input_file,
                       header=None,
                       usecols=[0, 1, 2, 5],
                       names=["Chromosome", "Start", "End", "Strand"],
                       dtype={"Start": int32,
                              "End": int32},
                       engine="c",
                       converters={"Strand": lambda s: 1 if s == "+" else 0})

    parallel_count_reads = partial(_count_reads_in_windows, fragment_size,
                                   window_size, keep_duplicates)

    # Getting chromosomes from the dict as an optimization;
    # the df might be gigabytes of data
    chromosomes = natsorted(chromosome_size_dict.keys())

    nb_chromosomes = len(chromosomes)
    info("Binning the {nb_chromosomes} chromosomes in {genome}.".format(**vars(
    )))
    # Use joblib parallel to run in parallel to get better error messages
    chromosome_dfs = Parallel(n_jobs=nb_cpus)(delayed(parallel_count_reads)(
        df[df.Chromosome == chromosome], chromosome_size_dict[chromosome])
                                              for chromosome in chromosomes)

    return chromosome_dfs


def _count_reads_in_windows(fragment_size, window_size, keep_duplicates, df,
                            chromosome_size):
    """Bin a chromosome.

    This function is called once per chromosome, usually in parallel.

    (Strand 1 means "+", Strand 0 means "-")
    """

    # Needed so pandas can be sure we are not working on a copy
    # If unsure, it gives an annoying SettingWithCopyWarning
    pd.options.mode.chained_assignment = None

    if not keep_duplicates:
        df = df.drop_duplicates(["Start", "End", "Strand"])

    fragment_size_halved = fragment_size // 2

    # Strand 1 means +, Strand 0 means -
    # Switch end to start for negative strand reads
    df.loc[df["Strand"] == 0, "Start"] = df.End

    # Adjust position of read from start to middle
    df.loc[df.Strand == 1, 'Start'] += fragment_size_halved
    df.loc[df.Strand == 0, 'Start'] -= fragment_size_halved

    # Turn position into a multiple of the window size
    df.Start = df.Start - (df.Start % window_size)

    # For bins we do not consider the end or strand
    df = df.drop(["End", "Strand"], axis=1)

    # Need to sort before counting bins
    df = df.sort_values(by="Start")
    # Count number of reads in each bin
    df["Count"] = df.groupby(["Start"])["Chromosome"].transform("count")
    # Only need one entry per bin
    df = df.drop_duplicates(["Chromosome", "Start"])

    df = df.rename(columns={"Start": "Bin"})

    df = remove_out_of_bounds_bins(df, chromosome_size)
    # Warn user if reads too far out; then something is wrong.
    _check_if_reads_are_too_far_outside_bounds(df, chromosome_size,
                                               window_size)

    df[["Bin", "Count"]] = df[["Bin", "Count"]].astype(int32)

    pd.options.mode.chained_assignment = "warn"

    return df[["Count", "Chromosome", "Bin"]]


def _check_if_reads_are_too_far_outside_bounds(df, chromosome_size,
                                               window_size):
    """If any reads are more than window_size outside of the chromosome endpoint
    this is a sign that something is wrong - most likely they are not using the
    correct genome version."""

    if df.empty:
        return None

    reads_too_far_outside_bounds = (df.Bin - window_size) > chromosome_size
    nb_reads_too_far_outside_bounds = sum(reads_too_far_outside_bounds)
    chromosome = df.iloc[0]["Chromosome"]

    if nb_reads_too_far_outside_bounds:
        warning("There were {nb_reads_too_far_outside_bounds}" \
                " reads too far outside of chromosome bounds " \
                " on chromosome {chromosome}. This probably means" \
                " you are not using the correct genome version.".format(
                    **locals()))
