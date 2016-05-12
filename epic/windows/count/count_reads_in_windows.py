from functools import partial
from io import BytesIO
from itertools import product
from logging import info
from subprocess import check_output

import pandas as pd
from joblib import Parallel, delayed
from numpy import int32

from natsort import natsorted

from epic.config.genomes import create_genome_size_dict, get_genome_size_file
from epic.windows.count.merge_chromosome_dfs import merge_chromosome_dfs
from epic.windows.count.remove_out_of_bounds_bins import remove_out_of_bounds_bins


def count_reads_in_windows(bed_file, args):

    chromosome_size_dict = create_genome_size_dict(args.genome)
    chromosomes = natsorted(list(chromosome_size_dict.keys()))

    if not args.paired_end:
        parallel_count_reads = partial(_count_reads_in_windows, bed_file, args)
    else:
        parallel_count_reads = partial(_count_reads_in_windows_paired_end,
                                       bed_file, args)

    info("Binning chromosomes {}".format(", ".join([c.replace("chr", "")
                                                    for c in chromosomes])))
    chromosome_dfs = Parallel(n_jobs=args.number_cores)(
        delayed(parallel_count_reads)(chromosome_size_dict[chromosome],
                                      chromosome, strand)
        for chromosome, strand in product(chromosomes, ["+", "-"]))

    info("Merging the bins on both strands per chromosome.")
    both_chromosome_strand_dfs = [df_pair
                                  for df_pair in _pairwise(chromosome_dfs)]
    merged_chromosome_dfs = Parallel(
        n_jobs=args.number_cores)(delayed(merge_chromosome_dfs)(df_pair)
                                  for df_pair in both_chromosome_strand_dfs)

    return merged_chromosome_dfs


def _count_reads_in_windows(bed_file, args, chromosome_size, chromosome,
                            strand):

    halved_fragment_size = args.fragment_size // 2
    idx = 1 if strand == "+" else 2  # fragment start indices

    if not args.keep_duplicates:
        duplicate_handling = " uniq | "
    else:
        duplicate_handling = ""

    if bed_file.endswith(".gz"):
        grep = "zgrep "
    elif bed_file.endswith(".bz2"):
        grep = "bzgrep "
    elif bed_file.endswith(".bam"):
        grep = "bamToBed -i {} | grep".format(bed_file)
        bed_file = ""
    else:
        grep = "grep "

    command = (
        "{grep} -E '^{chromosome}\\b.*\\{strand}$' {bed_file} | "
        "cut -f 1-3,6 - | sort -k2,3n | {duplicate_handling} "
        "LC_ALL=C perl -a -ne '$F[{idx}]{strand}={halved_fragment_size};"  # shift fragments
        "$F[{idx}]=$F[{idx}]-($F[{idx}] % {window_size});"  # turn exact start into bin
        "print \"@F[0,{idx}]\n\"' | "
        "uniq -c | "
        "sed -e 's/^[ ]*//'").format(grep=grep,
                                     chromosome=chromosome,
                                     strand=strand,
                                     bed_file=bed_file,
                                     duplicate_handling=duplicate_handling,
                                     idx=idx,
                                     halved_fragment_size=halved_fragment_size,
                                     window_size=args.window_size)

    output = check_output(command, shell=True)

    out_table = pd.read_table(
        BytesIO(output),
        header=None,
        sep=" ",
        names=["Count", "Chromosome", "Bin"])

    out_table = remove_out_of_bounds_bins(out_table, chromosome_size)

    out_table[["Bin", "Count"]] = out_table[["Bin", "Count"]].astype(int32)

    return out_table


def _count_reads_in_windows_paired_end(bed_file, args, chromosome_size):

    if not args.keep_duplicates:
        duplicate_handling = " uniq | "
    else:
        duplicate_handling = ""

    command = """bamToBed -pe -i {bed_file} |
    grep -E "^{chromosome}\\b.*{chromosome}\\b.*" | # Both chromos must be equal; no chimeras (?)
    cut -f 1-6  | sort -k2,3n -k4,5n | {duplicate_handling} # get chr start end chr start end for PE; sort on location
    LC_ALL=C perl -a -ne 'use List::Util qw[min max]; $start = min($F[1], $F[2]); $end = max($F[4], $F[5]); $middle = $start + int(($end - $start)/2); $bin = $middle - $middle % 200; print "$F[0] $bin \n"' | # Find bin of midpoint between start and ends
    uniq -c |
    sed -e 's/^[ ]*//'
    """.format(**locals())

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
