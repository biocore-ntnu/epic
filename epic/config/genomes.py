import pkg_resources
import logging

from epic.config import logging_settings
from epic.utils.find_readlength import (find_readlength,
                                        get_closest_readlength)

__author__ = "Endre Bakken Stovner https://github.com/endrebak/"
__license__ = "MIT"


def get_genome_size_file(genome):
    return pkg_resources.resource_filename(
        "epic", "scripts/chromsizes/{}.chromsizes".format(genome))


def create_genome_size_dict(genome):
    """Creates genome size dict from string containing data."""

    size_file = get_genome_size_file(genome)
    size_lines = open(size_file).readlines()

    size_dict = {}
    for line in size_lines:
        genome, length = line.split()
        size_dict[genome] = int(length)

    return size_dict


def get_effective_genome_length(genome, read_length):

    egf = pkg_resources.resource_string(
        "epic", "scripts/effective_sizes/{}_{}.txt".format(
            genome, read_length)).split()[-1].decode()

    genome_length = sum(create_genome_size_dict(genome).values())

    logging.info("Using an effective genome fraction of {}.".format(egf))

    return float(egf) * genome_length
