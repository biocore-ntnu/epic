import pytest

from joblib import Parallel, delayed
from collections import defaultdict
import pandas as pd
import numpy as np

import pkg_resources, os
from natsort import natsorted

from io import StringIO

# from helper.functions

import logging
from rpy2.robjects import r, pandas2ri
pandas2ri.activate()

from rpy2.robjects.packages import importr
importr("S4Vectors")
bioc = importr("GenomicRanges")

length_of_rle = r("function (x) sum(runLength(x))")


__author__ = "Endre Bakken Stovner https://github.com/endrebak/"
__license__ = "MIT"


def nucleotide_overlaps_per_file(all_files, nb_cpu):

    rles = files_to_chromosome_coverage(all_files, nb_cpu)

    nucleotide_overlaps = Parallel(n_jobs=nb_cpu)(delayed(_nucleotide_overlaps_per_file)(
        f, rles) for f in rles)

    return pd.concat(nucleotide_overlaps).sort_values(["Main", "Other"]).reset_index(drop=True)


def _nucleotide_overlaps_per_file(bed_file, extended_rles):

    base_bed = bed_file.split("/")[-1].split(".")[0]
    logging.info("Finding the number of nucleotides in " + base_bed + " overlapping other files.")

    _rle_to_iranges = r("""function(x, c){
    gr = GRanges(c, IRanges(start(x), width=runLength(x)), value=runValue(x))
    gr[elementMetadata(gr)$value==1]
    }""")

    _find_overlaps = r("function(x, o) subjectHits(findOverlaps(subject=x, query=o))")
    _set_metadata = r("function(x, ov, f) {elementMetadata(x) = 0; elementMetadata(x[ov]) = 1; data.frame(Start=start(x), End=end(x), Data=elementMetadata(x))}")

    cvs = extended_rles[bed_file]
    irs = {c: _rle_to_iranges(d, c) for c, d in cvs.items()}

    overlaps = defaultdict(list)
    for f in extended_rles:
        base_bed_other = f.split("/")[-1].split(".")[0]
        cvos = extended_rles[f]
        iros = {c: _rle_to_iranges(d, c) for c, d in cvos.items()}

        overlapping_chromosomes = set(irs.keys()).intersection(iros.keys())
        for c in overlapping_chromosomes:
            ir, iro = irs[c], iros[c]
            fov = _find_overlaps(ir, iro)
            df = _set_metadata(ir, fov, base_bed_other)
            df = pandas2ri.ri2py(df)
            df.insert(2, "Chromosome", c)
            df = df.set_index("Start End Chromosome".split())
            df.columns = [base_bed_other]
            overlaps[c].append(df)

    rowdicts = []
    for c, v in natsorted(overlaps.items()):
        df = pd.concat(v, axis=1).reset_index()
        lengths = df.End - df.Start
        df = df.set_index("Chromosome Start End".split())
        for f in df:
            # print(f)
            idx = df.loc[df[f] != 0, f]
            total_overlaps = lengths.ix[idx].sum()
            rowdict = {"Chromosome": c, "Main": base_bed, "Other": f, "Overlaps": total_overlaps}
            rowdicts.append(rowdict)

    return pd.DataFrame.from_dict(rowdicts)["Chromosome                    Main                    Other  Overlaps".split()]


def overlap_matrix_nucleotides(all_files, nb_cpu):
    rles = files_to_chromosome_coverage(all_files, nb_cpu)

    nucleotide_overlaps = Parallel(n_jobs=nb_cpu)(delayed(_overlap_matrix_nucleotides)(
        f, rles) for f in rles)

    df = pd.concat(nucleotide_overlaps)

    return df.reset_index(drop=True)


def _overlap_matrix_nucleotides(bed_file, extended_rles):

    overlaps = _create_overlap_matrix_nucleotides(bed_file, extended_rles)
    return _counts_runlengths(bed_file, overlaps)


def files_to_chromosome_coverage(all_files, nb_cpu):

    df_to_coverage = r("function(x) coverage(GRanges(x$Chromosome, IRanges(x$Start, x$End)))")

    logging.info("Finding nucleotide coverage of files.")
    coverages = defaultdict(dict)
    for f in all_files:
        df = pd.read_table(f, usecols=[0, 1, 2], header=None, names="Chromosome Start End".split())
        cv = df_to_coverage(df)

        chromosomes = r["names"](cv)
        for chromosome in chromosomes:
            coverages[f][chromosome] = r("function(x, idx) x[[idx]]")(cv, chromosome)

    max_per_chromosome_coverage = r("function(x) sum(runLength(x))")
    remove_duplicate_list_entries = r("function(x) x[unique(names(x))]")

    maxlengths = defaultdict(int)
    for f, data in coverages.items():
        for chromosome, rle in data.items():
            current_len = max_per_chromosome_coverage(rle)[0]
            maxlengths[chromosome] = max(maxlengths[chromosome], current_len)

    extend_rle = r('function(cvg, maxlen) c(cvg,Rle(0,maxlen-length(cvg)))')
    extended_rles = defaultdict(dict)
    for f, d in coverages.items():
        for chromosome, cv in d.items():
            maxlength = maxlengths[chromosome]
            cv = extend_rle(cv, maxlength)
            extended_rles[f][chromosome] = cv

    return extended_rles



def _create_overlap_matrix_nucleotides(bed_file, coverages):

    base_bed_other = bed_file.split("/")[-1]
    logging.info("Processing {} at nucleotide level".format(base_bed_other))

    cvs = coverages[bed_file]

    files = [f for f in coverages if f != bed_file]
    compute_overlap = r("function(x, o) x + (x & o)")
    for f in files:
        # print(f)
        other_cvs = coverages[f]
        overlapping_chromosomes = set(cvs.keys()).intersection(other_cvs.keys())

        for chromosome in overlapping_chromosomes:
            chr_cov = cvs[chromosome]
            other_cov = other_cvs[chromosome]

            overlap = compute_overlap(chr_cov, other_cov)
            cvs[chromosome] = overlap

    return cvs


def _counts_runlengths(bed_file, cvs):

    base_bed = bed_file.split("/")[-1].split(".")[0]

    overlaps = defaultdict(int)
    get_runlength = r("function(x, v) runLength(x[x == v])")
    for chromosome, overlap in cvs.items():
        run_values = set(r["runValue"](overlap))
        for run_value in run_values:
            run_value = int(run_value)
            runlength = get_runlength(overlap, run_value)[0]
            overlaps[run_value] += runlength

    rowdicts = []
    for run_value, run_length in overlaps.items():
        rowdict = {"Main": base_bed, "Other": run_value, "Overlaps": run_length}
        rowdicts.append(rowdict)

    return pd.DataFrame.from_dict(rowdicts)


    # command = "bedtools makewindows -w {} -b {} | bedtools intersect -wo -filenames -a - -b {}".format(window_size, bed_file, all_files_str)
    # output = check_output(command, shell=True).decode()

    # df = pd.read_table(StringIO(output), usecols=[0, 1, 2, 3], header=None, names="Chromosome Start End File".split())

    # df.File = df.File.apply(basename).str.split(".", expand=True).ix[:, 0]

    # df.insert(0, "Main", base_bed_other)

    # return _compute_nucleotide_overlap(df, window_size)



# def make_treatment_control_same_length(treatment, control):

#     treatment_result_dict = dict()
#     control_result_dict = dict()
#     for chromosome in treatment:

#         if chromosome not in control:
#             logging.info("Chromosome " + chromosome +
#                          " missing from input data. Skipping.")
#             continue

#         chr_treatment, chr_control = treatment[chromosome], control[chromosome]

#         treatment_maxlen = r["max"](r["sapply"](list(chr_treatment.values()),
#                                                 length_of_rle))
#         control_maxlen = r["max"](r["sapply"](list(chr_control.values()),
#                                               length_of_rle))

#         maxlen = r["max"](treatment_maxlen, control_maxlen)

#         lapply = r('function(cvg, maxlen) c(cvg,Rle(0,maxlen-length(cvg)))')
