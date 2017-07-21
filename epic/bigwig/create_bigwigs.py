import logging
import numpy as np
from os.path import join, basename, splitext, dirname
from subprocess import call
from argparse import Namespace
import pandas as pd
from typing import Any, Dict, Iterable, List
import pyBigWig

from joblib import Parallel, delayed


def create_bigwigs(matrix, outdir, args):
    # type: (pd.DataFrame, str, Namespace) -> None
    """Create bigwigs from matrix."""
    call("mkdir -p {}".format(outdir), shell=True)
    genome_size_dict = args.chromosome_sizes

    outpaths, data = [], []
    for bed_file in matrix:
        outpath = join(outdir, splitext(basename(bed_file))[0] + ".bw")
        outpaths.append(outpath)

        bed_column = matrix[bed_file]
        data.append(bed_column)

    Parallel(n_jobs=args.number_cores)(delayed(_create_bigwig)(bed_column, outpath, genome_size_dict) for outpath, bed_column in zip(outpaths, data))


def _to_int(l):
    # type: (Iterable[Any]) -> List[int]
    return [int(i) for i in l]


def _create_bigwig(bed_column, outpath, genome_size_dict):
    # type: (pd.Series, str, Dict[str, int]) -> None

    logging.info("Creating biwgwig " + outpath)

    rpkm = 1e6 * bed_column / bed_column.sum()

    rpkm = [float(f) for f in list(rpkm.fillna(0).values.flatten())]

    bed_column = bed_column.reset_index()
    unique_chromosomes = list(bed_column.Chromosome.drop_duplicates())
    chromosomes = list(bed_column.Chromosome)
    starts = _to_int(list(bed_column.Bin))
    ends = _to_int(list(bed_column.End))

    header = [(c, int(genome_size_dict[c])) for c in unique_chromosomes]

    bw = pyBigWig.open(outpath, "w")
    bw.addHeader(header)

    bw.addEntries(chromosomes, starts, ends=ends, values=rpkm)
    bw.close()


def create_sum_bigwigs(matrix, args):

    rpkm_matrix = 1e6 * matrix / matrix.sum()

    chip = rpkm_matrix[args.treatment].sum(axis=1)
    input = rpkm_matrix[args.control].sum(axis=1)

    input_pseudo = input.copy()
    input_pseudo.loc[input_pseudo == 0] = 1
    log2fc = chip / input_pseudo.values

    bigwigs_to_create = []
    if args.chip_bigwig:
        folder = dirname(args.chip_bigwig)
        if folder:
            call("mkdir -p {}".format(folder), shell=True)
        bigwigs_to_create.append([args.chip_bigwig, chip])

    if args.input_bigwig:
        folder = dirname(args.input_bigwig)
        if folder:
            call("mkdir -p {}".format(folder), shell=True)
        bigwigs_to_create.append([args.input_bigwig, input])

    if args.log2fc_bigwig:
        folder = dirname(args.log2fc_bigwig)
        if folder:
            call("mkdir -p {}".format(folder), shell=True)
        bigwigs_to_create.append([args.log2fc_bigwig, log2fc])

    Parallel(n_jobs=args.number_cores)(delayed(_create_bigwig)(bed_column, outpath, args.chromosome_sizes) for outpath, bed_column in bigwigs_to_create)
