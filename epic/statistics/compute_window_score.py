from epic.utils.helper_functions import lru_cache
from numpy import log
from scipy.stats import poisson


@lru_cache()
def compute_window_score(i, poisson_parameter):

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
