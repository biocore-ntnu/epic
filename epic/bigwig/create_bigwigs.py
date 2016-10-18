import logging
import numpy as np
from os.path import join, basename, splitext
from subprocess import call

import pyBigWig

from joblib import Parallel, delayed


def create_bigwigs(matrix, outdir, args):
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
    return [int(i) for i in l]


def _create_bigwig(bed_column, outpath, genome_size_dict):

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


def create_sum_bigwigs(matrix, outdir, args):
    call("mkdir -p {}".format(outdir), shell=True)

    chip = matrix[args.treatment].sum(axis=1)
    input = matrix[args.control].sum(axis=1)

    chip_outpath = join(outdir, "chip_sum" + ".bw")
    input_outpath = join(outdir, "input_sum" + ".bw")

    # _create_bigwig(chip, chip_outpath, args)
    # _create_bigwig(input, input_outpath, args)

    Parallel(n_jobs=args.number_cores)(delayed(_create_bigwig)(bed_column, outpath, args.chromosome_sizes) for outpath, bed_column in zip([chip_outpath, input_outpath], [chip, input]))
