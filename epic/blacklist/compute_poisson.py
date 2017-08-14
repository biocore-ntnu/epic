from scipy.stats.distributions import poisson
from statsmodels.sandbox.stats.multicomp import multipletests

import pandas as pd

def compute_poisson(df, args):

    bins = int(args.effective_genome_fraction/int(args.window_size))

    bad_bins = []
    for s in df:
        s = df[s]
        unique_alignments = s.sum()

        average = int(unique_alignments)/bins

        value_counts = s.value_counts().fillna(0).index.get_level_values(0)
        poisson_scores = pd.Series(poisson.sf(value_counts, mu=average))
        poisson_scores = pd.concat([pd.Series(value_counts).to_frame(), poisson_scores], axis=1)
        poisson_scores.columns = ["Value", "Score"]
        poisson_scores = poisson_scores.set_index("Value")

        poisson_scores = pd.Series(index=poisson_scores.index, data=poisson_scores.Score)

        poisson_p_vals = s.replace(poisson_scores.to_dict())

        fdr = multipletests(poisson_p_vals, method="fdr_bh")[1]
        fdr = pd.Series(fdr, index=s.index, name="fdr")

        fdr_df = pd.concat([s, fdr], axis=1)

        r = fdr_df[fdr_df.fdr < args.fdr]
        bad_bins.append(r)

    return pd.concat(bad_bins, axis=1).index.to_frame().reset_index(drop=True)
