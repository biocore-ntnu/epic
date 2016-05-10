from epic.config.genomes import get_effective_genome_length
from scipy.stats import poisson, rankdata


def compute_fdr(df, total_chip_reads, total_input_reads, args):

    df.to_csv("for_fdr_test.csv", sep=" ")
    print("total_chip_reads", total_chip_reads)
    print("total_input_reads", total_input_reads)

    total_island_input_reads = df.Input.sum()

    # Hack needed in case we run on test data
    # TODO: why does SICER not need this? Different genome versions?
    # run with FDR=1 on original SICER and get this island:
    # chr7    61606400        61606799        3       2       0.167427550906  1.40719467956   0.1674275 50906
    # does it not show up with epic.
    if total_island_input_reads == 0:
        total_island_input_reads = 2

    scaling_factor = (total_chip_reads * 1.0) / total_input_reads

    effective_genome_size = get_effective_genome_length(args.genome)
    zero_controls_multiplier = total_input_reads * 1.0 / effective_genome_size

    avg_0_denom = (df.End - df.Start + 1) * zero_controls_multiplier
    avg_0_denom[avg_0_denom > 0.25] = 0.25
    avg_0_denom = avg_0_denom * scaling_factor

    avg = df.Input * scaling_factor
    avg[df.Input == 0] = avg_0_denom[df.Input == 0]

    df.P_value = poisson.sf(df.ChIP, avg)
    no_differential_expression = df.ChIP <= avg
    df.loc[no_differential_expression, "P_value"] = 1

    df.Fold_change = df.ChIP / avg

    ranked_p_values = rankdata(df.P_value)
    df.FDR_value = df.P_value * len(df) / ranked_p_values
    fdr_too_high = df.FDR_value > 1
    df.loc[fdr_too_high, "FDR_value"] = 1

    df = df[df.FDR_value < args.false_discovery_rate_cutoff]

    return df
