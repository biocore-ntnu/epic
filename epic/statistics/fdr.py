from epic.config.genomes import get_effective_genome_length
from scipy.stats import poisson, rankdata

def compute_fdr(df, total_chip_reads, total_input_reads, genome, fdr_cutoff):

    total_island_input_reads = df["Input"].sum()

    # Hack needed in case we run on test data
    if total_island_input_reads == 0:
        total_island_input_reads = 2

    scaling_factor = (total_chip_reads * 1.0) / total_input_reads

    effective_genome_size = get_effective_genome_length(genome)
    zero_controls_multiplier = total_input_reads * 1.0 / effective_genome_size

    avg_0_denom = (df["End"] - df["Start"] + 1) * zero_controls_multiplier
    avg_0_denom[avg_0_denom > 0.25] = 0.25
    avg_0_denom = avg_0_denom * scaling_factor

    avg = df["Input"] * scaling_factor
    avg[df["Input"] == 0] = avg_0_denom[df["Input"] == 0]

    df["P value"] = poisson.sf(df["ChIP"], avg)
    no_differential_expression = df["ChIP"] <= avg
    df.loc[no_differential_expression, "P value"] = 1

    df["Fold change"] = df["ChIP"]/avg

    ranked_p_values = rankdata(df["P value"])
    df["FDR value"] = df["P value"] * len(df) / ranked_p_values
    fdr_too_high = df["FDR value"] > 1
    df.loc[fdr_too_high, "FDR value"] = 1

    df = df[df["FDR value"] < fdr_cutoff]

    return df
