import pytest

from joblib import Parallel, delayed
from collections import defaultdict
import pandas as pd
import numpy as np
from typing import DefaultDict, Dict, Iterable

import pkg_resources, os
from natsort import natsorted

from io import StringIO

# from helper.functions

import logging
from rpy2.robjects import r, pandas2ri
pandas2ri.activate()
from rpy2.robjects.robject import RObject

from rpy2.robjects.packages import importr
importr("S4Vectors")
bioc = importr("GenomicRanges")

from epic.scripts.overlaps.files_to_chromosome_coverage import files_to_chromosome_coverage

__author__ = "Endre Bakken Stovner https://github.com/endrebak/"
__license__ = "MIT"

def overlap_matrix_nucleotides(all_files, nb_cpu):
    # type: (Iterable[str], int) -> pd.DataFrame
    rles = files_to_chromosome_coverage(all_files, nb_cpu)

    nucleotide_overlaps = Parallel(n_jobs=nb_cpu)(delayed(_overlap_matrix_nucleotides)(
        f, rles) for f in rles)

    df = pd.concat(nucleotide_overlaps)

    return df.reset_index(drop=True)


def _overlap_matrix_nucleotides(bed_file, extended_rles):
    # type: (str, Dict[str,Dict[str, RObject]]) -> pd.DataFrame

    overlaps = _create_overlap_matrix_nucleotides(bed_file, extended_rles)
    return _counts_runlengths(bed_file, overlaps)





def _create_overlap_matrix_nucleotides(bed_file, coverages):
    # type: (str, Dict[str,Dict[str, RObject]]) -> Dict[str, RObject]

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
    # type: (str, Dict[str, RObject]) -> pd.DataFrame

    base_bed = bed_file.split("/")[-1].split(".")[0]

    overlaps = defaultdict(int) # type: DefaultDict[str, int]
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
