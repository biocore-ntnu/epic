


from joblib import Parallel, delayed
from collections import defaultdict
import pandas as pd
import numpy as np

import pkg_resources, os
from natsort import natsorted

from io import StringIO

from typing import DefaultDict, Dict, Iterable
# from helper.functions

import logging
from rpy2.robjects import r, pandas2ri
pandas2ri.activate()
from rpy2.robjects.robject import RObject

from rpy2.robjects.packages import importr

importr("S4Vectors")
bioc = importr("GenomicRanges")



def files_to_chromosome_coverage(all_files, nb_cpu):
    # type: (Iterable[str], int) -> DefaultDict[str,Dict[str, RObject]]

    df_to_coverage = r("function(x) coverage(GRanges(x$Chromosome, IRanges(x$Start, x$End)))")

    logging.info("Finding nucleotide coverage of files.")
    coverages = defaultdict(dict) # type: DefaultDict[str,Dict[str, RObject]]
    for f in all_files:
        df = pd.read_table(f, usecols=[0, 1, 2], header=None, names="Chromosome Start End".split())
        cv = df_to_coverage(df)

        chromosomes = r["names"](cv)
        for chromosome in chromosomes:
            coverages[f][chromosome] = r("function(x, idx) x[[idx]]")(cv, chromosome)

    max_per_chromosome_coverage = r("function(x) sum(runLength(x))")
    remove_duplicate_list_entries = r("function(x) x[unique(names(x))]")

    maxlengths = defaultdict(int) # type: DefaultDict[str, int]
    for f, data in coverages.items():
        for chromosome, rle in data.items():
            current_len = max_per_chromosome_coverage(rle)[0]
            maxlengths[chromosome] = max(maxlengths[chromosome], current_len)

    extend_rle = r('function(cvg, maxlen) c(cvg,Rle(0,maxlen-length(cvg)))')
    extended_rles = defaultdict(dict) # type: DefaultDict[str,Dict[str, RObject]]
    for f, d in coverages.items():
        for chromosome, cv in d.items():
            maxlength = maxlengths[chromosome]
            cv = extend_rle(cv, maxlength)
            extended_rles[f][chromosome] = cv

    return extended_rles
