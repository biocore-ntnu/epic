import pytest

from joblib import Parallel, delayed
from collections import defaultdict
import pandas as pd
import numpy as np
from typing import Dict, Iterable

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

length_of_rle = r("function (x) sum(runLength(x))")


from epic.scripts.overlaps.files_to_chromosome_coverage import (files_to_chromosome_coverage)

__author__ = "Endre Bakken Stovner https://github.com/endrebak/"
__license__ = "MIT"


def nucleotide_overlaps_per_file(all_files, nb_cpu):
    # type: (Iterable[str], int) -> pd.DataFrame

    rles = files_to_chromosome_coverage(all_files, nb_cpu)

    nucleotide_overlaps = Parallel(n_jobs=nb_cpu)(delayed(_nucleotide_overlaps_per_file)(
        f, rles) for f in rles)

    print(nucleotide_overlaps)
    return pd.concat(nucleotide_overlaps).sort_values(["Main", "Other"]).reset_index(drop=True)


def _nucleotide_overlaps_per_file(bed_file, extended_rles):
    # type: (str, Dict[str,Dict[str, RObject]]) -> pd.DataFrame

    base_bed = bed_file.split("/")[-1].split(".")[0]
    logging.info("Finding the number of nucleotides in " + base_bed + " overlapping other files.")

    _find_overlaps = r("""
    function(s, o) {
    runValue(s) = as.logical(runValue(s))
    runValue(o) = as.logical(runValue(o))
    sum(s & o)
    }
    """)

    _find_total = r(" function(s) {runValue(s) = as.logical(runValue(s)); sum(s) }")

    cvs = extended_rles[bed_file]

    rowdicts = []
    for f in extended_rles:

        base_bed_other = f.split("/")[-1].split(".")[0]
        # print("base bed other", base_bed_other)
        cvos = extended_rles[f]

        overlapping_chromosomes = set(cvs.keys()).intersection(cvos.keys())
        overlaps, total = 0, 0
        for c in overlapping_chromosomes:

            cv, cvo = cvs[c], cvos[c]

            fov = _find_overlaps(cv, cvo)[0]
            tot = _find_total(cv)[0]

            overlaps += fov
            total += tot

        ratio = overlaps/total
        rowdict = {"Chromosome": c, "Main": base_bed, "Other": base_bed_other, "Overlaps": ratio}
        rowdicts.append(rowdict)

    return pd.DataFrame.from_dict(rowdicts)["Chromosome Main Other Overlaps".split()]
