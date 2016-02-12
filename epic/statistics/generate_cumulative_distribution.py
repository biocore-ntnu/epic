from epic.config.constants import E_VALUE, BIN_SIZE


def generate_cumulative_dist(island_expectations_d, total_length):
    """
    Generate cumulative distribution: a list of tuples (bins, hist).
    """

    cumulative = [0] * (total_length + 1)
    partial_sum = 0.0

    island_expectations = []
    for i in range(len(cumulative)):
        if i in island_expectations_d:
            island_expectations.append(island_expectations_d[i])
        else:
            island_expectations.append(0)

    for index in range(1, len(island_expectations) + 1):
        complimentary = len(island_expectations) - index
        partial_sum += island_expectations[complimentary]
        cumulative[complimentary] = partial_sum

    # move to function call
    for index in range(len(cumulative)):
        if cumulative[index] <= E_VALUE:
            score_threshold = index * BIN_SIZE
            break

    return score_threshold
