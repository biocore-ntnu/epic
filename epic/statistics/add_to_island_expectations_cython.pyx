

from numpy import log
from scipy.stats import poisson
from numpy import log, pi, exp
from scipy.misc import factorial
from epic.utils.helper_functions import lru_cache

cdef int E_VALUE
cdef double WINDOW_P_VALUE, E_VALUE_THRESHOLD, BIN_SIZE

WINDOW_P_VALUE = 0.20
BIN_SIZE = 0.001
E_VALUE = 1000
E_VALUE_THRESHOLD = E_VALUE * .0000001

def add_to_island_expectations_cython(double average_window_readcount, int current_max_scaled_score,
                                      int island_eligibility_threshold, island_expectations,
                                      double gap_contribution):

    """Can probably be heavily optimized.
    Time required to run can be seen from logging info."""

    cdef double island_expectation
    cdef int index, i, current_island, scaled_score

    scaled_score = current_max_scaled_score + E_VALUE
    for index in range(current_max_scaled_score + 1, scaled_score+1):
        island_expectation=0.0
        i = island_eligibility_threshold #i is the number of tags in the added window

        current_island = int(round(index - compute_window_score(i, average_window_readcount)/BIN_SIZE))

        while(current_island >= 0):
            if current_island in island_expectations:
                # print("current_island",current_island)
                island_expectation += _poisson(i, average_window_readcount) * island_expectations[current_island]
                # print("island_expectation",island_expectation)
            i += 1
                # print("i",i)
            current_island = int(round(index - compute_window_score(i, average_window_readcount)/BIN_SIZE))
                # print("current_island",current_island)
        island_expectation *= gap_contribution
        if island_expectation:
            island_expectations[index] = island_expectation

    return island_expectations


cdef inline double compute_window_score(i, poisson_parameter):

    # No enrichment; poisson param also average
    if i < poisson_parameter:
        return 0

    p_value = poisson.pmf(i, poisson_parameter)

    if p_value > 0:
        window_score = -log(p_value)
    else:
        # log of zero not defined
        window_score = 1000

    return window_score

cdef inline double _poisson(i, average):
    """
    """
    exponent = -average + i*log(average) - _factln(i)
    return exp(exponent)

cdef inline double _factln(num):
    """
    Computes logfactorial regularly for tractable numbers, uses Ramanujans approximation otherwise.
    """

    if num < 20:
        log_factorial = log(factorial(num))
    else:
        log_factorial = num*log(num) -num + log(num*(1+4*num*(1+2*num)))/6.0 + log(pi)/2

    return log_factorial
