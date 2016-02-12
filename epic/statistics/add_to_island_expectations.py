def add_to_island_expectations_dict(
        average_window_readcount, current_max_scaled_score,
        island_eligibility_threshold, island_expectations, gap_contribution):
    """Can probably be heavily optimized.
    Time required to run can be seen from logging info."""

    scaled_score = current_max_scaled_score + E_VALUE
    for index in range(current_max_scaled_score + 1, scaled_score + 1):
        island_expectation = 0.0
        i = island_eligibility_threshold  #i is the number of tags in the added window

        current_island = int(round(index - compute_window_score(
            i, average_window_readcount) / BIN_SIZE))

        while (current_island >= 0):

            if current_island in island_expectations:
                island_expectation += _poisson(
                    i, average_window_readcount) * island_expectations[
                        current_island]
            i += 1
            current_island = int(round(index - compute_window_score(
                i, average_window_readcount) / BIN_SIZE))
        island_expectation *= gap_contribution
        if island_expectation:
            island_expectations[index] = island_expectation

    return island_expectations
