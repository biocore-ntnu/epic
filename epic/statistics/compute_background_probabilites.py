"""This is a module that computes the island threshold.

It is a sanitized and memory efficient version of the SICER version, that
produces results similar to the second or third decimal. I have a less efficient
version that gives even closer results, but I figure that the computations are
very approximate (due to how they are implemented in the original) so I just use
the faster version.
"""

from __future__ import print_function
from argparse import Namespace
from typing import Tuple
import logging

from epic.statistics.compute_values_needed_for_recurrence import (
    compute_enriched_threshold, compute_gap_factor, compute_boundary)
from epic.statistics.compute_score_threshold import compute_score_threshold


# @MEMORY.cache(verbose=0)
def compute_background_probabilities(total_chip_count, args):
    # type: (int, Namespace) -> Tuple[float, int, float]

    effective_genome_fraction = args.effective_genome_fraction
    logging.debug(str(effective_genome_fraction) + " effective_genome_fraction")
    # move outside of function call

    average_window_readcount = total_chip_count * (
        args.window_size / float(effective_genome_fraction))
    logging.debug(str(args.window_size) + " window size")
    logging.debug(str(total_chip_count) + " total chip count")
    logging.debug(str(average_window_readcount) + " average_window_readcount")

    island_enriched_threshold = compute_enriched_threshold(
        average_window_readcount)
    logging.debug(str(island_enriched_threshold) +
                  " island_enriched_threshold")

    gap_contribution = compute_gap_factor(
        island_enriched_threshold, args.gaps_allowed, average_window_readcount)
    logging.debug(str(gap_contribution) + " gap_contribution")

    boundary_contribution = compute_boundary(
        island_enriched_threshold, args.gaps_allowed, average_window_readcount)
    logging.debug(str(boundary_contribution) + " boundary_contribution")

    genome_length_in_bins = effective_genome_fraction / args.window_size

    score_threshold = compute_score_threshold(
        average_window_readcount, island_enriched_threshold, gap_contribution,
        boundary_contribution, genome_length_in_bins)

    return score_threshold, island_enriched_threshold, average_window_readcount
