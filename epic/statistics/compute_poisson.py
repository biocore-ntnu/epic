from math import log, factorial, pi, exp
from epic.utils.helper_functions import lru_cache


@lru_cache()
def _factln(num):
    # type: (int) -> float
    """
    Computes logfactorial regularly for tractable numbers, uses Ramanujans approximation otherwise.
    """

    if num < 20:
        log_factorial = log(factorial(num))
    else:
        log_factorial = num * log(num) - num + log(num * (1 + 4 * num * (
            1 + 2 * num))) / 6.0 + log(pi) / 2

    return log_factorial


@lru_cache()
def _poisson(i, average):
    # type: (int, float) -> float
    """
    """
    exponent = -average + i * log(average) - _factln(i)
    return exp(exponent)
