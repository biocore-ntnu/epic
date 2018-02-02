import logging
from collections import defaultdict, OrderedDict
from os.path import basename

import pandas as pd

from docopt import docopt
from natsort import natsorted
from joblib import Parallel, delayed


from epic.merge.compute_bed_bins import merge_bed_bins, compute_bins
from epic.merge.merge_helpers import compute_bin_size, add_new_enriched_bins_matrixes



def main(files, regions, keep_nonenriched, enriched_per_file, nb_cpus):

    dfs = read_dfs(files)
    bin_size = compute_bin_size(dfs)

    if regions:
        dfs = add_new_enriched_bins_matrixes(regions, dfs, bin_size)
        keep_nonenriched = False

    merged_df = merge_matrixes(dfs, keep_nonenriched, regions, enriched_per_file, nb_cpus)

    return merged_df



def enriched_indexes(dfs):

    # find all chromosomes in files
    iter_dfs = iter(dfs.values())
    df = next(iter_dfs)

    current_index = df.filter(like="Enriched_")
    name = current_index.columns[0]

    enriched = current_index.loc[current_index[name] == 1].index

    for df in iter_dfs:

        current_index = df.filter(like="Enriched_")
        name = current_index.columns[0]

        current_index = current_index.loc[current_index[name] == 1].index

        enriched = enriched.union(current_index)

    return enriched


def remove_nonenriched(enriched, dfs):

    new_dfs = OrderedDict()
    for f, df in dfs.items():

        df = df.ix[df.index.intersection(enriched)]

        new_dfs[f] = df

    return new_dfs


def all_chromosomes(dfs):

    # find all chromosomes in files
    chromosomes = set()
    for df in dfs.values():
        chromosomes.update(set(df.index.get_level_values("Chromosome").drop_duplicates()))

    return(chromosomes)


def split_dfs_into_chromosome_dfs(dfs, chromosomes):

    # separate files into chromosomes
    chromosome_dfs = defaultdict(list)
    for f, df in dfs.items():
        for chromosome in chromosomes:
            try:
                chromosome_dfs[chromosome].append(df.xs(chromosome, level="Chromosome", drop_level=False))
            except KeyError: # if file is missing chromosome, insert dummy df to ensure same column ordering
                chromosome_dfs[chromosome].append(pd.DataFrame(columns=dfs[f].columns))

    return chromosome_dfs


def _merge_dfs(dfs):

    return pd.concat(dfs, axis=1).fillna(0)


def merge_dfs(chromosome_dfs, nb_cpu):

    # merge all chromosome dfs horizontally, per df
    merged_chromosome_dfs = Parallel(n_jobs=nb_cpu)(
        delayed(_merge_dfs)(dfs)
        for dfs in chromosome_dfs.values())

    # vertically merge (concat) all chromosome dfs
    merged_df = pd.concat(merged_chromosome_dfs)

    # indexes are shuffled after the intersection, resort
    merged_df = merged_df.reindex(index=natsorted(merged_df.index))

    return merged_df


def merge_matrixes(dfs, keep_nonenriched, regions, enriched_per_file, nb_cpus):

    if not keep_nonenriched and not regions:
        enriched = enriched_indexes(dfs)
        dfs = remove_nonenriched(enriched, dfs)

        if all([df.empty for df in dfs.values()]):
            raise Exception("No enriched bins in your matrixes. Exiting.")

    chromosomes = all_chromosomes(dfs)

    chromosome_dfs = split_dfs_into_chromosome_dfs(dfs, chromosomes)

    merged_df = merge_dfs(chromosome_dfs, nb_cpus)

    enriched_cols = [c for c in merged_df if "Enriched_" in c]

    total_enriched = merged_df[enriched_cols].sum(axis=1)

    merged_df.insert(0, "TotalEnriched", total_enriched)

    merged_df = merged_df.set_index("TotalEnriched", append=True)

    if enriched_per_file:
        merged_df = merged_df.set_index(enriched_cols, append=True)
    else:
        merged_df = merged_df.drop(enriched_cols, axis=1)

    column_order = natsorted(merged_df.columns)

    merged_df = merged_df[merged_df.index.get_level_values("TotalEnriched") > 0]

    return merged_df[column_order]


def read_dfs(files):

    full_path = False
    if not len(files) == len(set([basename(f) for f in files])):
        logging.info("Matrix-files do not have a unique basename. Using full path in header!")
        full_path = True

    dfs = OrderedDict()
    for f in files:
        df = pd.read_table(f, header=0, sep=" ", index_col=[0, 1])

        df = df[~df.index.duplicated(keep='first')]

        columns = list(df.columns)
        file_nick = "Enriched_" + basename(f) if not full_path else "Enriched_" + f
        columns[0] = file_nick
        df.columns = columns

        logging.info("Calling " + f + " " + file_nick + " in matrix file.")
        dfs[f] = df

    return dfs
