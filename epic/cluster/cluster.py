import pandas as pd

from joblib import Parallel, delayed

from natsort import natsorted


def _trunks_flanks_valleys(cdf, trunk_diff, bin_size, distance_allowed):

    dfs = []
    for gid, gdf in cdf.groupby(((cdf.Bin - cdf.Bin.shift()).abs() > distance_allowed).cumsum()):

        enriched_diff = gdf.TotalEnriched >= (gdf.TotalEnriched.max() - trunk_diff)
        enriched_sum = (enriched_diff.shift(1) != enriched_diff).cumsum()
        grpby = gdf.groupby(enriched_sum)
        nb_groups = len(grpby)
        max_value = gdf.TotalEnriched.max()
        chromosome = gdf.Chromosome.head(1).iloc[0]

        dfs2 = []
        # group islands into trunks/flanks/valleys
        for (i, gdf2) in grpby:

            is_trunk = gdf2.head(1).TotalEnriched.iloc[0] >= (max_value - trunk_diff)

            if not is_trunk and (i == 1 or i == (nb_groups)): # first or last, i.e. flank
                status = "flank"
            elif not is_trunk:
                status = "valley"
            else:
                status = "trunk"

            start = str(gdf2.head(1).Bin.iloc[0])
            end = str(gdf2.tail(1).Bin.iloc[0] + bin_size - 1)
            region_id = start + ":" + end

            gdf3 = pd.DataFrame(gdf2.drop("TotalEnriched Chromosome Bin".split(), 1).sum()).T
            gdf3.insert(0, "Kind", status)
            gdf3.insert(0, "RegionID", region_id)
            gdf3.insert(0, "IslandID", gid)
            gdf3.insert(0, "Chromosome", chromosome)

            dfs2.append(gdf3)

        df2 = pd.concat(dfs2, axis=0)
        dfs.append(df2)

    return pd.concat(dfs)


def trunks_flanks_valleys(df, bin_size=200, trunk_diff=1, distance_allowed=200, nb_cpu=1):

        dfs = []

        dfs = Parallel(n_jobs=nb_cpu)(
                delayed(_trunks_flanks_valleys)(
                        df, trunk_diff, bin_size, distance_allowed) for (c, df) in df.groupby("Chromosome"))

        outdf = pd.concat(dfs, axis=0)
        # make all nondata columns one column so r can handle it
        columns_to_join = "IslandID RegionID Kind".split()
        index = outdf.Chromosome.str.cat([outdf[c].astype(str) for c in columns_to_join], sep="_")
        columns_to_drop = columns_to_join + ["Chromosome"]
        outdf = outdf.drop(columns_to_drop, axis=1)

        col_order = ["Index"] + natsorted([c for c in outdf.columns], key=lambda x: x.split("_"))
        outdf.insert(0, "Index", index)

        return outdf.set_index("Index")
