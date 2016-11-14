from natsort import natsorted
from collections import OrderedDict
import pkg_resources
import logging
from typing import Dict

from epic.config import logging_settings
from epic.utils.find_readlength import (find_readlength,
                                        get_closest_readlength)

__author__ = "Endre Bakken Stovner https://github.com/endrebak/"
__license__ = "MIT"


def get_genome_size_file(genome):
    # type: (str) -> str

    genome_names = pkg_resources.resource_listdir("epic", "scripts/chromsizes")
    name_dict = {n.lower().replace(".chromsizes", ""): n for n in genome_names}

    # No try/except here, because get_egs would already have failed if genome
    # did not exist
    genome_exact = name_dict[genome.lower()]

    return pkg_resources.resource_filename(
        "epic", "scripts/chromsizes/{}".format(genome_exact))


def create_genome_size_dict(genome):
    # type: (str) -> Dict[str,int]
    """Creates genome size dict from string containing data."""

    size_file = get_genome_size_file(genome)
    size_lines = open(size_file).readlines()

    size_dict = {}
    for line in size_lines:
        genome, length = line.split()
        size_dict[genome] = int(length)

    return size_dict


def create_genome_size_dict_custom_genome(chromsizes):
    # type: (str) -> OrderedDict[str, int]

    chromosome_lengths = [l.split() for l in open(chromsizes).readlines()]

    od = OrderedDict()          # type: OrderedDict[str, int]

    for c, l in natsorted(chromosome_lengths):
        od[c] = int(l)

    return od


def get_effective_genome_length(genome, read_length):
    # type: (str, int) -> float

    genome_names = pkg_resources.resource_listdir("epic",
                                                  "scripts/effective_sizes")
    name_dict = {n.split("_")[0]: "".join(n.split("_")[:-1])
                 for n in genome_names}

    try:
        genome_exact = name_dict[genome.lower()]
        egf = pkg_resources.resource_string( # type: ignore
            "epic", "scripts/effective_sizes/{}_{}.txt".format(
                genome_exact, read_length)).split()[-1].decode()
    except KeyError:
        genome_list = "\n".join(list(name_dict.keys()))
        logging.error(
            "Genome " + genome +
            " not found.\n These are the available genomes: " + genome_list +
            "\nIf yours is not there, please request it at github.com/endrebak/epic .")

    genome_length = sum(create_genome_size_dict(genome).values())

    logging.info("Using an effective genome fraction of {}.".format(egf))

    assert float(
        egf) < 1, "Something wrong happened, effective genome fraction over 1!"

    egs = float(egf) * genome_length

    return egs
