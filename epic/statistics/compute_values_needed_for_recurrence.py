from itertools import count
from scipy.stats import poisson

from epic.config.constants import WINDOW_P_VALUE


def compute_enriched_threshold(average_window_readcount):
    """
    Computes the minimum number of tags required in window for an island to be enriched.
    """

    current_threshold, survival_function = 0, 1
    for current_threshold in count(start=0, step=1):
        survival_function -= poisson.pmf(current_threshold,
                                         average_window_readcount)
        if survival_function <= WINDOW_P_VALUE:
            break

    island_enriched_threshold = current_threshold + 1

    return island_enriched_threshold


def compute_gap_factor(island_enriched_threshold, gap_intervals_allowed,
                       poisson_distribution_parameter):

    max_gap_score = 1
    gap_factor = single_gap_factor(island_enriched_threshold,
                                   poisson_distribution_parameter)
    max_gap_score += sum([pow(gap_factor, i)
                          for i in range(1, gap_intervals_allowed + 1)])
    return max_gap_score


def single_gap_factor(island_enriched_threshold,
                      poisson_distribution_parameter):

    poisson_scores = [poisson.pmf(i, poisson_distribution_parameter)
                      for i in range(island_enriched_threshold)]
    return sum(poisson_scores)


def compute_boundary(island_enriched_threshold, gap_intervals_allowed,
                     average):

    single_gap = single_gap_factor(island_enriched_threshold, average)
    single_boundary_score = pow(single_gap, gap_intervals_allowed + 1)
    start_and_end_score = single_boundary_score * single_boundary_score

    return start_and_end_score
