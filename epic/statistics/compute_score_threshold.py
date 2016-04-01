import logging
from scipy.stats import poisson
from numpy import log

from epic.config.constants import BIN_SIZE, E_VALUE_THRESHOLD
from epic.statistics.generate_cumulative_distribution import generate_cumulative_dist
from epic.statistics.add_to_island_expectations import add_to_island_expectations_dict


def compute_score_threshold(average_window_readcount,
                            island_enriched_threshold, gap_contribution,
                            boundary_contribution, genome_length_in_bins):
    """
    What does island_expectations do?
    """

    required_p_value = poisson.pmf(island_enriched_threshold,
                                   average_window_readcount)

    prob = boundary_contribution * required_p_value

    score = -log(required_p_value)

    current_scaled_score = int(round(score / BIN_SIZE))

    island_expectations_d = {}
    island_expectations_d[current_scaled_score] = prob * genome_length_in_bins
    island_expectations_d[
        0] = boundary_contribution * genome_length_in_bins / gap_contribution

    current_max_scaled_score = current_scaled_score

    interval = int(1 / BIN_SIZE)
    partial_cumu = 0
    logging.info("Finding the score required to consider an island enriched.")
    while (partial_cumu > E_VALUE_THRESHOLD or partial_cumu < 1e-100):

        current_scaled_score += interval
        current_max_scaled_score = current_scaled_score - interval
        # logging.debug(island_expectations_d)

        if current_scaled_score > current_max_scaled_score:

            # logging.debug(island_expectations_d)
            island_expectations_d = add_to_island_expectations_dict(
                average_window_readcount, current_max_scaled_score,
                island_enriched_threshold, island_expectations_d,
                gap_contribution)
            partial_cumu = 0.0001
            current_max_scaled_score += 1000

            if max(island_expectations_d) > interval:
                partial_cumu = sum(
                    [val
                     for idx, val in island_expectations_d.items()
                     if idx > current_max_scaled_score - interval])
            else:
                partial_cumu = sum(island_expectations_d.values())

    logging.debug("Computing cumulative distribution.")
    score_threshold = generate_cumulative_dist(island_expectations_d,
                                               current_max_scaled_score + 1)
    logging.info("Enriched score threshold for islands: " + str(
        score_threshold))
    return score_threshold
