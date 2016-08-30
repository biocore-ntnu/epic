import logging
from os.path import dirname, join, basename
from subprocess import call

import pandas as pd

from joblib import Parallel, delayed

from natsort import natsorted

from epic.windows.count.remove_out_of_bounds_bins import remove_bins_with_ends_out_of_bounds
from epic.config.genomes import get_genome_size_file, create_genome_size_dict


def write_matrix_files(chip_merged, input_merged, df, args):

    matrixes = create_matrixes(chip_merged, input_merged, df, args)

    if args.store_matrix:
        print_matrixes(matrixes, args)

    matrix = pd.concat(matrixes, axis=0)
    matrix = matrix.drop("Island", axis=1)
    ends = matrix.index.get_level_values("Bin") + int(args.window_size) - 1
    matrix.insert(0, "End", ends)
    matrix = matrix.set_index("End", append=True)
    matrix = matrix.sort_index(level="Chromosome")
    # TODO: remove out of bounds bins

    if args.individual_bedgraph or args.individual_bigwig:
        filenames = individual_bedgraphs(matrix, args)
        if args.individual_bigwig:
            individual_bigwigs(filenames, args)

    if args.bedgraph or args.bigwig:
        filenames = bedgraph(matrix, args)
        if args.bigwig:
            bigwig(filenames, args)


def bigwig(filenames, args):

    chromsizes = get_genome_size_file(args.genome)

    outfiles = [
        f.replace(args.bedgraph, args.bigwig + "/", 1).replace(
            ".bedgraph", ".bw").replace("//", "/") for f in filenames
    ]

    call("mkdir -p {}".format(args.bigwig), shell=True)
    for infile, outfile in zip(filenames, outfiles):
        logging.info("Creating bigiwg {}".format(outfile))
        command = "bedGraphToBigWig {} {} {}".format(infile, chromsizes,
                                                     outfile)
        call(command, shell=True)


def individual_bigwigs(filenames, args):

    chromsizes = get_genome_size_file(args.genome)

    outfiles = [f.replace(args.individual_bedgraph,
                          args.individual_bigwig + "/",
                          1).replace(".bedgraph", ".bw").replace("//", "/")
                for f in filenames]

    call("mkdir -p {}".format(args.individual_bigwig), shell=True)
    for infile, outfile in zip(filenames, outfiles):
        command = "bedGraphToBigWig {} {} {}".format(infile, chromsizes,
                                                     outfile)
        logging.info("Creating bigiwg {}".format(outfile))
        call(command, shell=True)


def bedgraph(matrix, args):
    "Create a bedgraph file for ChIP and Input."

    outfolder = args.bedgraph or "temp_bedgraph/"

    call("mkdir -p {}".format(outfolder), shell=True)

    chip_file = join(outfolder, "sum_treatment.bedgraph")
    input_file = join(outfolder, "sum_input.bedgraph")

    c_sum = matrix[args.treatment].sum(1)
    c = c_sum[c_sum > 0]
    logging.info("Writing ChIP bedgraph to {}.".format(chip_file))
    c.astype(int).to_frame().to_csv(chip_file, sep="\t", header=False)

    i_sum = matrix[args.control].sum(1)
    i = i_sum[i_sum > 0]
    logging.info("Writing Input bedgraph to {}.".format(input_file))
    i.astype(int).to_frame().to_csv(input_file, sep="\t", header=False)

    return [chip_file, input_file]


def _individual_bedgraphs(matrix, name, outfolder):

    base = basename(name).split(".")

    if len(base) > 2:
        base = "".join(base[:-2])
    else:
        base = "".join(base[:-1])

    outfile = join(outfolder, base + ".bedgraph")
    s = matrix[name]

    # na only in those where input or chip lacks chromo
    nonzeroes_only = s[s != 0].dropna()
    logging.info("Writing bedgraph for file {} to {}.".format(name, outfile))
    nonzeroes_only.astype(int).to_frame().to_csv(outfile,
                                                 sep="\t",
                                                 na_rep="NA",
                                                 header=False)

    return outfile


def individual_bedgraphs(matrix, args):
    "Create a bedgraph file for each file used."

    outfolder = args.individual_bedgraph or "temp_individual_bedgraph/"

    call("mkdir -p {}".format(outfolder), shell=True)

    outfiles = []
    for treatment_file in args.treatment:
        outfiles.append(_individual_bedgraphs(matrix, treatment_file,
                                              outfolder))

    for control_file in args.control:
        outfiles.append(_individual_bedgraphs(matrix, control_file, outfolder))

    return outfiles


def _create_matrixes(chromosome, chip, input, islands, chromosome_size,
                     window_size):

    chip_df = get_chromosome_df(chromosome, chip)
    input_df = get_chromosome_df(chromosome, input)

    chip_df["Chromosome"] = chip_df["Chromosome"].astype("category")
    chip_df["Bin"] = chip_df["Bin"].astype(int)
    chip_df = chip_df.set_index("Chromosome Bin".split())
    chip_df = islands.join(chip_df, how="right")
    chip_df = chip_df[~chip_df.index.duplicated(keep='first')]

    input_df["Chromosome"] = input_df["Chromosome"].astype("category")
    input_df["Bin"] = input_df["Bin"].astype(int)
    input_df = input_df.set_index("Chromosome Bin".split())
    input_df = input_df[~input_df.index.duplicated(keep='first')]

    dfm = chip_df.join(input_df, how="outer", sort=False).fillna(0)

    dfm = remove_bins_with_ends_out_of_bounds(dfm, chromosome_size,
                                              window_size)

    return dfm


def create_matrixes(chip, input, df, args):

    "Creates matrixes which can be written to file as is (matrix) or as bedGraph."

    genome = create_genome_size_dict(args.genome)

    chip = put_dfs_in_chromosome_dict(chip)
    input = put_dfs_in_chromosome_dict(input)
    all_chromosomes = natsorted(set(list(chip.keys()) + list(input.keys())))

    islands = enriched_bins(df, args)

    logging.info("Creating matrixes from count data.")
    dfms = Parallel(n_jobs=args.number_cores)(delayed(_create_matrixes)(
        chromosome, chip, input, islands, genome[chromosome],
        args.window_size) for chromosome in all_chromosomes)

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


def put_dfs_in_dict(dfs):

    sample_dict = {}
    for df in dfs:

        if df.empty:
            continue

        chromosome = df.head(1).Chromosome.values[0]
        sample_dict[chromosome] = df

    return sample_dict


def put_dfs_in_chromosome_dict(dfs):

    chromosome_dict = {}
    for df in dfs:

        if df.empty:
            continue

        chromosome = df.head(1).Chromosome.values[0]
        chromosome_dict[chromosome] = df

    return chromosome_dict


def get_chromosome_df(chromosome, df_dict):

    if chromosome in df_dict:
        df = df_dict[chromosome]
    else:
        df = pd.DataFrame(columns="Chromosome Bin".split())

    return df


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
