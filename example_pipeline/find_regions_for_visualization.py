
from sys import argv

import pandas as pd

bed_file = argv[1]

df = pd.read_table(bed_file)

g1 = df.groupby("ClusterID")
g2 = iter(df.groupby("ClusterID"))

next(g2)

potential_dfs = []
for (group1, gdf1), (group2, gdf2) in zip(g1, g2):

    kinds = gdf1.RegionKind

    if len(set(kinds)) != 3:
        continue

    if len(gdf1) >= 15:
        continue

    if abs(gdf1.tail(1).End.iloc[0] - gdf2.head(1).Start.iloc[0]) > 10000:
        continue

    max_enriched = 0

    for enriched in gdf1.TotalEnriched:
        max_enriched = max(max([int(i) for i in enriched.split(",")]), max_enriched)

    if max_enriched != 3:
        continue

    chromosome = gdf1.head(1).Chromosome.iloc[0]
    start = gdf1.head(1).Start.iloc[0]
    end = gdf2.tail(1).End.iloc[0]

    print(chromosome, start, end, sep="\t")
