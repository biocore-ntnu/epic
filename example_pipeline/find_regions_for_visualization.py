
from sys import argv

import pandas as pd

bed_file = argv[1]

df = pd.read_table(bed_file)

potential_dfs = []
for group, gdf in df.groupby("ClusterID"):

    kinds = gdf.RegionKind

    # want all three types of regions to be represented
    if len(set(kinds)) != 3:
        continue

    if len(gdf) != 5:
        continue

    # max_enriched = gdf.MaxEnrichedCluster.iloc[0]

    # if max_enriched != 3:
    #     continue

    max_median_enriched = gdf.MedianEnrichedRegion.max()

    if max_median_enriched != 3:
        continue

    potential_dfs.append(gdf)

new_df = pd.concat(potential_dfs)

for group, gdf in new_df.groupby("ClusterID"):
    chromosome = gdf.head(1).Chromosome.iloc[0]
    start = gdf.head(1).Start.iloc[0]
    end = gdf.tail(1).End.iloc[0]

    print(chromosome, start, end, group, sep="\t")

# print(new_df.to_csv(sep="\t", index=False))
