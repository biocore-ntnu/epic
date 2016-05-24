import logging
from sys import platform
from re import search, IGNORECASE
from io import BytesIO
from subprocess import check_output

import pandas as pd

from epic.config import logging_settings

__author__ = "Endre Bakken Stovner https://github.com/endrebak/"
__license__ = "MIT"


def find_readlength(args):
    """Estimate length of reads based on 1000 first."""

    bed_file = args.treatment[0]

    filereader = "cat "
    if bed_file.endswith(".gz") and search("linux", platform, IGNORECASE):
        filereader = "zcat "
    elif bed_file.endswith(".gz") and search("darwin", platform, IGNORECASE):
        filereader = "gzcat "
    elif bed_file.endswith(".bz2"):
        filereader = "bzgrep "
    elif bed_file.endswith(".bam"):
        filereader = "bamToBed -i "

    command = filereader + "{} | head -1000".format(bed_file)
    output = check_output(command, shell=True)

    df = pd.read_table(
        BytesIO(output),
        header=None,
        usecols=[1, 2],
        sep="\s+",
        names=["Start", "End"])

    avg_readlength = (df.End - df.Start).mean()

    logging.info("Used {} to estimate an average read length of {}".format(
        bed_file, avg_readlength))

    return avg_readlength


def get_closest_readlength(estimated_readlength):
    """Find the predefined readlength closest to the estimated readlength.

    In the case of a tie, choose the shortest readlength."""

    readlengths = [36, 50, 75, 100]
    differences = [abs(r - estimated_readlength) for r in readlengths]
    min_difference = min(differences)
    index_of_min_difference = [i
                               for i, d in enumerate(differences)
                               if d == min_difference][0]

    return readlengths[index_of_min_difference]
