"""This is a module that computes the island threshold.

It is a sanitized and memory efficient version of the SICER version, that
produces results similar to the second or third decimal. I have a less efficient
version that gives even closer results, but I figure that the computations are
very approximate (due to how they are implemented in the original) so I just use
the faster version.
"""

from __future__ import print_function

import logging
import sys

from numpy import log
from itertools import count
from scipy.stats import poisson
from math import log, factorial, pi, exp

from epic.utils.helper_functions import lru_cache
from epic.config.genomes import get_effective_genome_length

from epic.config.cache_settings import MEMORY

import pyximport
pyximport.install()

from epic.statistics.add_to_island_expectations_cython import add_to_island_expectations_cython

WINDOW_P_VALUE = 0.20
BIN_SIZE = 0.001
E_VALUE = 1000
E_VALUE_THRESHOLD = E_VALUE * .0000001
