import logging
from natsort import natsorted

from joblib import Parallel, delayed

try:
    from functools import lru_cache  # noqa
except ImportError:
    from functools32 import lru_cache  # noqa


def _merge_chip_and_input(chip_df, input_df):

    chip_df_nb_bins = len(chip_df)
    merged_df = chip_df.merge(input_df,
                              how="left",
                              on=["Chromosome", "Bin"],
                              suffixes=[" ChIP", " Input"])
    merged_df = merged_df[["Chromosome", "Bin", "Count ChIP", "Count Input"]]
    merged_df.columns = ["Chromosome", "Bin", "ChIP", "Input"]

    merged_df = merged_df.fillna(0)

    assert len(merged_df) == chip_df_nb_bins

    return merged_df


def merge_chip_and_input(chip_dfs, input_dfs, nb_cpu):

    # should be same length, since missing chromos get empty df
    assert len(chip_dfs) == len(input_dfs)

    logging.info("Merging ChIP and Input data.")
    merged_chromosome_dfs = Parallel(n_jobs=nb_cpu)(
        delayed(_merge_chip_and_input)(chip_df, input_df)
        for chip_df, input_df in zip(chip_dfs, input_dfs))
    return merged_chromosome_dfs


def get_total_number_of_reads(dfs):
    return sum([df.Count.sum() for df in dfs])


def ensure_same_chromosomes_in_list(sample1_dfs, sample2_dfs):

    d1 = create_chromsome_df_map(sample1_dfs)
    d2 = create_chromsome_df_map(sample2_dfs)

    d1, d2 = fill_missing_chromosomes(d1, d2)

    sample1_dfs = [v for (k, v) in natsorted(d1.items())]
    sample2_dfs = [v for (k, v) in natsorted(d2.items())]

    return sample1_dfs, sample2_dfs


def create_chromsome_df_map(dfs):

    sample_dict = {}
    for df in dfs:
        chromosome = df.head(1).Chromosome.values[0]
        sample_dict[chromosome] = df

    return sample_dict


def fill_missing_chromosomes(d1, d2):

    all_chromosomomes = set(d1.keys()).union(d2.keys())

    d1_missing = set(d2.keys()).difference(d1.keys())
    d2_missing = set(d1.keys()).difference(d2.keys())

    for chromosome in d1_missing:
        d1[chromosome] = pd.DataFrame(columns=["Chromosome", "Bin"])

    for chromosome in d2_missing:
        d2[chromosome] = pd.DataFrame(columns=["Chromosome", "Bin"])

    return d1, d2


def merge_same_files(sample1_dfs, sample2_dfs, nb_cpu):

    # if one list is missing a chromosome, we might pair up the wrong dataframes
    # therefore creating dicts beforehand to ensure they are paired up properly
    sample1_dfs, sample2_dfs = ensure_same_chromosomes_in_list(sample1_dfs,
                                                               sample2_dfs)

    logging.info("Merging same class data.")
    merged_chromosome_dfs = Parallel(n_jobs=nb_cpu)(
        delayed(_merge_same_files)(sample1_df, sample2_df)
        for sample1_df, sample2_df in zip(sample1_dfs, sample2_dfs))

    return merged_chromosome_dfs


def _merge_same_files(sample1_df, sample2_df):

    merged_df = sample1_df.merge(sample2_df,
                                 how="outer",
                                 on=["Chromosome", "Bin"])

    return merged_df.fillna(0)
