import pytest

import pandas as pd
import numpy as np

from epic.windows.count.pandas_count_reads_in_windows import pandas_count_reads_in_windows
from epic.windows.count.count_reads_in_windows import count_reads_in_windows


# slow so comment out only when needed
def describe_pandas_and_gnu_gives_the_same_results():
    @pytest.mark.slow
    @pytest.mark.parametrize("fragment_size", [150, 200, 50])
    @pytest.mark.parametrize("window_size", [50, 100, 200, 1000])
    @pytest.mark.parametrize("keep_duplicates", [True, False])
    def with_various_arguments(fragment_size, window_size, keep_duplicates):

        bed_file, genome = "examples/test.bed", "hg19"

        gnu_counts = count_reads_in_windows(bed_file, genome, fragment_size,
                                            window_size, keep_duplicates, 1)
        pandas_counts = pandas_count_reads_in_windows(
            bed_file, genome, fragment_size, window_size, keep_duplicates, 1)

        gnu_counts = pd.concat(gnu_counts)
        pandas_counts = pd.concat(pandas_counts)
        print(gnu_counts, gnu_counts.dtypes, "gnu_counts")
        print(pandas_counts, pandas_counts.dtypes, "pandas_counts")

        assert np.array_equal(gnu_counts, pandas_counts)
