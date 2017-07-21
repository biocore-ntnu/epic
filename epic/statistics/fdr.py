import pandas as pd
from numpy import log2
from scipy.stats import poisson, rankdata
from argparse import Namespace

def compute_fdr(df, total_chip_reads, total_input_reads, args):
    # type: (pd.DataFrame, int, int, Namespace) -> pd.DataFrame

    total_island_input_reads = df.Input.sum()

    # Hack needed in case we run on test data
    # TODO: why does SICER not need this? Different genome versions?
    # run with FDR=1 on original SICER and get this island:
    # chr7    61606400        61606799        3       2       0.167427550906  1.40719467956   0.1674275 50906
    # does it not show up with epic?
    if total_island_input_reads == 0:
        total_island_input_reads = 2

    scaling_factor = (total_chip_reads * 1.0) / total_input_reads

    zero_controls_multiplier = total_input_reads * 1.0 / args.effective_genome_fraction

    avg_0_denom = (df.End - df.Start + 1) * zero_controls_multiplier
    avg_0_denom[avg_0_denom > 0.25] = 0.25
    avg_0_denom = avg_0_denom * scaling_factor

    avg = df.Input * scaling_factor
    avg[df.Input == 0] = avg_0_denom[df.Input == 0]

    fold_change = df.ChIP / avg
    log2FC = log2(fold_change)
    df.insert(len(df.columns), "Log2FC", log2FC)

    p_vals = pd.Series(poisson.sf(df.ChIP, avg), index=df.index)

    p_vals[df.ChIP <= avg] = 1
    df.insert(len(df.columns), "P", p_vals)

    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * len(df) / ranked_p_values
    fdr[fdr > 1] = 1
    df.insert(len(df.columns), "FDR", fdr)

    df = df[df.FDR < args.false_discovery_rate_cutoff]

    return df
