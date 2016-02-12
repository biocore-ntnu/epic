from functools import partial
from io import BytesIO
from itertools import product
from logging import info
from subprocess import check_output

import pandas as pd
from joblib import Parallel, delayed
from numpy import int32

from natsort import natsorted

from epic.config.genomes import create_genome_size_dict
from epic.windows.count.merge_chromosome_dfs import merge_chromosome_dfs
from epic.windows.count.remove_out_of_bounds_bins import remove_out_of_bounds_bins


def count_reads_in_windows(bed_file,
                           genome,
                           fragment_size,
                           window_size,
                           keep_duplicates,
                           nb_cpus=1):

    chromosome_size_dict = create_genome_size_dict(genome)
    chromosomes = natsorted(list(chromosome_size_dict.keys()))

    parallel_count_reads = partial(_count_reads_in_windows, bed_file,
                                   fragment_size, window_size, keep_duplicates)

    info("Binning chromosomes {}".format(", ".join([c.replace("chr", "")
                                                    for c in chromosomes])))
    chromosome_dfs = Parallel(n_jobs=nb_cpus)(delayed(parallel_count_reads)(
        chromosome_size_dict[chromosome], chromosome,
        strand) for chromosome, strand in product(chromosomes, ["+", "-"]))

    info("Merging the bins on both strands per chromosome.")
    both_chromosome_strand_dfs = [df_pair
                                  for df_pair in _pairwise(chromosome_dfs)]
    merged_chromosome_dfs = Parallel(
        n_jobs=nb_cpus)(delayed(merge_chromosome_dfs)(df_pair)
                        for df_pair in both_chromosome_strand_dfs)

    return merged_chromosome_dfs


def _count_reads_in_windows(bed_file, fragment_size, window_size,
                            keep_duplicates, chromosome_size, chromosome,
                            strand):

    halved_fragment_size = fragment_size // 2
    idx = 1 if strand == "+" else 2  # fragment start indices

    if not keep_duplicates:
        duplicate_handling = " uniq | "
    else:
        duplicate_handling = ""

    if bed_file.endswith(".gz"):
        grep = "zgrep "
    elif bed_file.endswith(".bz2"):
        grep = "bzgrep "
    else:
        grep = "grep "

    command = """{grep} -E '^{chromosome}\\b.*\\{strand}$' {bed_file} |
    cut -f 1-3,6 - | sort -k2,3n | {duplicate_handling}
    LC_ALL=C perl -a -ne '$F[{idx}]{strand}={halved_fragment_size}; # shift fragments
    $F[{idx}]=$F[{idx}]-($F[{idx}] % {window_size}); # turn exact start into bin
    print "@F[0,{idx}]\n"' |
    uniq -c |
    sed -e 's/^[ ]*//'""".format(**locals())

    # debug(command)
    output = check_output(command, shell=True)

    out_table = pd.read_table(
        BytesIO(output),
        header=None,
        sep=" ",
        names=["Count", "Chromosome", "Bin"])

    out_table = remove_out_of_bounds_bins(out_table, chromosome_size)

    out_table[["Bin", "Count"]] = out_table[["Bin", "Count"]].astype(int32)

    return out_table


def _pairwise(iterable):
    col = iter(iterable)
    return zip(col, col)
