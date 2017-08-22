import logging
from epic.config import logging_settings

from scipy.stats.distributions import poisson
from statsmodels.sandbox.stats.multicomp import multipletests

import pandas as pd

def compute_poisson(df, args):

    nb_bins = int(args.effective_genome_fraction/int(args.window_size))

    bad_bins = []
    for fname in df:
        s = df[fname]
        unique_alignments = s.sum()

        average = int(unique_alignments)/nb_bins

        value_counts = s.drop_duplicates().values
        poisson_scores = pd.Series(poisson.sf(value_counts, mu=average))
        poisson_scores = pd.concat([pd.Series(value_counts).to_frame(), poisson_scores], axis=1)
        poisson_scores.columns = ["Value", "Score"]
        poisson_scores = poisson_scores.set_index("Value")

        poisson_scores = pd.Series(index=poisson_scores.index, data=poisson_scores.Score)

        poisson_p_vals = s.replace(poisson_scores.to_dict())

        bonferroni = multipletests(poisson_p_vals, method="bonferroni")[1]
        bonferroni = pd.Series(bonferroni, index=s.index, name="bonferroni")

        bonferroni_df = pd.concat([s, bonferroni], axis=1)

        print(bonferroni_df.head(10).to_csv(sep=" "))

        r = bonferroni_df[bonferroni_df.bonferroni < args.bonferroni]
        logging.info(str(len(r)) + " blacklist-bins found in file " + fname + " out of a total of " + str(len(bonferroni_df)) + " bins (" + str(len(r)/len(bonferroni_df)) + "%)")

        bad_bins.append(r)

    outdf = pd.concat(bad_bins, axis=1).reset_index()
    outdf.insert(1, "End", outdf.Bin + args.window_size - 1)

    return outdf["Chromosome Bin End".split()]
