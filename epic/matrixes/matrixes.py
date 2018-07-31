import sys
import logging
from os.path import dirname, join, basename
from subprocess import call
from itertools import chain
from typing import Iterable, Sequence, Tuple
from argparse import Namespace

import numpy as np
import pandas as pd

from joblib import Parallel, delayed

from natsort import natsorted

from epic.windows.count.remove_out_of_bounds_bins import remove_bins_with_ends_out_of_bounds
from epic.config.genomes import get_genome_size_file


def write_matrix_files(chip_merged, input_merged, df, args):
    # type: (Dict[str, pd.DataFrame], Dict[str, pd.DataFrame], pd.DataFrame, Namespace) -> None

    matrixes = create_matrixes(chip_merged, input_merged, df, args)

    matrix = pd.concat(matrixes, axis=0, sort=False)

    if args.store_matrix:
        print_matrixes(matrix, args)

    # reset and setting index hack to work around pandas bug

    if args.bigwig or args.individual_log2fc_bigwigs or args.chip_bigwig or args.input_bigwig or args.log2fc_bigwig:
        # matrixes = [m.astype(np.float64).reset_index() for m in matrixes if not m.empty]

        matrix = matrix.astype(np.float64)

        matrix = matrix.drop("Enriched", axis=1)
        ends = matrix.index.get_level_values("Bin") + int(args.window_size) - 1
        matrix.insert(0, "End", ends)
        matrix = matrix.set_index("End", append=True)
        matrix = matrix.sort_index(level="Chromosome")

    # TODO: remove out of bounds bins

    if args.bigwig:
        # defer initialization so not run during travis
        from epic.bigwig.create_bigwigs import create_bigwigs
        create_bigwigs(matrix, args.bigwig, args)

    if args.individual_log2fc_bigwigs:
        # defer initialization so not run during travis
        from epic.bigwig.create_bigwigs import create_log2fc_bigwigs
        create_log2fc_bigwigs(matrix, args.individual_log2fc_bigwigs, args)

    if args.chip_bigwig or args.input_bigwig or args.log2fc_bigwig:
        # defer initialization so not run during travis
        from epic.bigwig.create_bigwigs import create_sum_bigwigs
        create_sum_bigwigs(matrix, args)


def _create_matrixes(chromosome, chip, input, islands,
                     chromosome_size, window_size):
    # type: (str, Dict[str, pd.DataFrame], Dict[str, pd.DataFrame], pd.DataFrame, int, int) -> pd.DataFrame

    # print("islands2\n" + islands.head(10).to_csv(sep=" "), file=sys.stderr)

    chip_df = get_chromosome_df(chromosome, chip)
    input_df = get_chromosome_df(chromosome, input)

    try:
        chromo_islands = islands.xs(chromosome, drop_level=False)
    except KeyError:
        return pd.DataFrame(index="Chromosome Bin".split())

    chip_df["Chromosome"] = chip_df["Chromosome"].astype("category")

    # START workaround
    # Should ideally have been just one line: chip_df["Bin"] = chip_df["Bin"].astype(int)
    # Workaround for the following error:
    # ValueError: assignment destination is read-only
    bins = chip_df["Bin"].astype(int)
    chip_df = chip_df.drop("Bin", axis=1)

    chip_df.insert(0, "Bin", bins)

    # print("chip_df1\n", chip_df.head(10).to_csv(sep=" "), file=sys.stderr)

    # END workaround

    chip_df = chip_df.set_index("Chromosome Bin".split())

    # removing duplicates to avoid joining problems
    chip_df = chip_df[~chip_df.index.duplicated(keep='first')]
    chromo_islands = chromo_islands[~chromo_islands.index.duplicated(keep='first')]

    # chromo_islands.to_csv("chromo_islands.csv", sep=" ")
    # chip_df.to_csv("chip_df.csv", sep=" ")
    # print(chromo_islands.head(20).to_csv(sep=" "), file=sys.stderr)

    # print(chromosome)

    # print("chip_df2\n", chip_df.head(10).to_csv(sep=" "), file=sys.stderr)
    # print(chromo_islands.head(10).to_csv(sep=" "), file=sys.stderr)
    chip_df = chromo_islands.join(chip_df, how="outer").fillna(0)
    # print("chip_df3\n", chip_df.head(10).to_csv(sep=" "), file=sys.stderr)

    # print("chip_df", chip_df.tail().to_csv(sep=" "), file=sys.stderr)

    chip_df = chip_df[~chip_df.index.duplicated(keep='first')]
    # print("chip_df4\n", chip_df.head(10).to_csv(sep=" "), file=sys.stderr)

    # print("chip_df", chip_df.tail().to_csv(sep=" "), file=sys.stderr)


    input_df["Chromosome"] = input_df["Chromosome"].astype("category")

    # START workaround
    # Should ideally have been just one line: input_df["Bin"] = input_df["Bin"].astype(int)
    # Workaround for the following error:
    # ValueError: assignment destination is read-only
    bins = input_df["Bin"].astype(int)
    input_df = input_df.drop("Bin", axis=1)

    input_df.insert(0, "Bin", bins)

    input_df = input_df.set_index("Chromosome Bin".split())

    # END workaround

    input_df = input_df[~input_df.index.duplicated(keep='first')]

    dfm = chip_df.join(input_df, how="outer", sort=False).fillna(0)
    # print("dfm1\n", dfm.head(10).to_csv(sep=" "), file=sys.stderr)

    dfm = remove_bins_with_ends_out_of_bounds(dfm, chromosome_size,
                                              window_size)

    dfm = dfm[~dfm.index.duplicated(keep='first')]
    # print("dfm2\n", dfm.head(10).to_csv(sep=" "), file=sys.stderr)

    # print(dfm.tail().to_csv(sep=" "), file=sys.stderr)

    return dfm


def create_matrixes(chip, input, df, args):
    # type: (Iterable[pd.DataFrame], Iterable[pd.DataFrame], pd.DataFrame, Namespace) -> List[pd.DataFrame]
    "Creates matrixes which can be written to file as is (matrix) or as bedGraph."

    genome = args.chromosome_sizes

    chip = put_dfs_in_chromosome_dict(chip)
    input = put_dfs_in_chromosome_dict(input)
    all_chromosomes = natsorted(set(list(chip.keys()) + list(input.keys())))

    # print("df1\n", df, file=sys.stderr)
    islands = enriched_bins(df, args)
    # print("islands1\n", islands, file=sys.stderr)


    logging.info("Creating matrixes from count data.")
    dfms = Parallel(n_jobs=args.number_cores)(delayed(_create_matrixes)(
        chromosome, chip, input, islands, genome[chromosome],
        args.window_size) for chromosome in all_chromosomes)

    return dfms


def print_matrixes(matrix, args):
    # type: (Iterable[pd.DataFrame], Namespace) -> None
    outpath = args.store_matrix

    dir = dirname(outpath)
    if dir:
        call("mkdir -p {}".format(dir), shell=True)

    logging.info("Writing data matrix to file: " + outpath)

    matrix.to_csv(outpath, sep=" ", header=True, compression="gzip")


def get_island_bins(df, window_size, genome, args):
    # type: (pd.DataFrame, int, str, Namespace) -> Dict[str, Set[int]]
    """Finds the enriched bins in a df."""

    # need these chromos because the df might not have islands in all chromos
    chromosomes = natsorted(list(args.chromosome_sizes))

    chromosome_island_bins = {} # type: Dict[str, Set[int]]
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


def put_dfs_in_dict(dfs):
    # type: (Iterable[pd.DataFrame]) -> Dict[str, pd.DataFrame]
    sample_dict = {}
    for df in dfs:

        if df.empty:
            continue

        chromosome = df.head(1).Chromosome.values[0]
        sample_dict[chromosome] = df

    return sample_dict


def put_dfs_in_chromosome_dict(dfs):
    # type: (Iterable[pd.DataFrame]) -> Dict[str, pd.DataFrame]

    chromosome_dict = {}        # type: Dict[str, pd.DataFrame]
    for df in dfs:

        if df.empty:
            continue

        chromosome = df.head(1).Chromosome.values[0]
        chromosome_dict[chromosome] = df

    return chromosome_dict


def get_chromosome_df(chromosome, df_dict):
    # type: (str, Dict[str, pd.DataFrame]) -> pd.DataFrame

    if chromosome in df_dict:
        df = df_dict[chromosome]
    else:
        df = pd.DataFrame(columns="Chromosome Bin".split())


    # print(chromosome, file=sys.stderr)
    # print(df, file=sys.stderr)

    return df


def enriched_bins(df, args):
    # type: (pd.DataFrame, Namespace) -> pd.DataFrame

    df = df.loc[df.FDR < args.false_discovery_rate_cutoff]

    idx_rowdicts = []
    for _, row in df.iterrows():
        for bin in range(
                int(row.Start), int(row.End), int(args.window_size)):
            idx_rowdicts.append({"Chromosome": row.Chromosome,
                                 "Bin": bin,
                                 "Enriched": 1})
    islands = pd.DataFrame.from_dict(idx_rowdicts)
    islands.loc[:, "Chromosome"].astype("category")
    islands.loc[:, "Bin"].astype(int)

    return islands.set_index("Chromosome Bin".split())
