import pandas as pd

from joblib import Parallel, delayed

from epic.bigwig.create_bigwigs import _create_bigwig

def _trunks_flanks_valleys(cdf, trunk_diff, bin_size, distance_allowed):

    dfs = []
    for cid, cdf in cdf.groupby(((cdf.Bin - cdf.Bin.shift()).abs() > distance_allowed).cumsum()):

        enriched_diff = cdf.TotalEnriched >= (cdf.TotalEnriched.max() - trunk_diff)
        enriched_sum = (enriched_diff.shift(1) != enriched_diff).cumsum()
        grpby = cdf.groupby(enriched_sum)
        nb_groups = len(grpby)
        max_value = cdf.TotalEnriched.max()
        chromosome = cdf.Chromosome.head(1).iloc[0]

        cluster_start = str(cdf.head(1).Bin.iloc[0])
        cluster_end = str(cdf.tail(1).Bin.iloc[0] + bin_size - 1)

        dfs2 = []
        # group islands into trunks/flanks/valleys
        for (i, rdf) in grpby:

            is_trunk = rdf.head(1).TotalEnriched.iloc[0] >= (max_value - trunk_diff)

            if not is_trunk and (i == 1 or i == (nb_groups)): # first or last, i.e. flank
                status = "flank"
            elif not is_trunk:
                status = "valley"
            else:
                status = "trunk"

            start = str(rdf.head(1).Bin.iloc[0])
            end = str(rdf.tail(1).Bin.iloc[0] + bin_size - 1)
            region_id = start + ":" + end

            total_enriched = ",".join(rdf.TotalEnriched.astype(str))
            bins = ",".join(rdf.Bin.astype(str))

            # min_enriched = rdf.TotalEnriched.min()
            # median_enriched =  rdf.TotalEnriched.median()

            gdf3 = pd.DataFrame(rdf.drop("TotalEnriched Chromosome Bin".split(), 1).sum()).T
            gdf3.insert(0, "TotalEnriched", total_enriched)
            gdf3.insert(0, "Bins", bins)
            gdf3.insert(0, "MaxEnrichedCluster", max_value)
            # gdf3.insert(0, "MinEnrichedRegion", min_enriched)
            # gdf3.insert(0, "MedianEnrichedRegion", median_enriched)
            gdf3.insert(0, "Start", int(start))
            gdf3.insert(0, "End", int(end))
            gdf3.insert(0, "RegionKind", status)
            gdf3.insert(0, "RegionID", region_id)
            gdf3.insert(0, "ClusterID", "_".join([chromosome, str(cid)]))
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

    return outdf


def create_cluster_bigwig(df, outpath, genome_size_dict, bin_size):

    bins = pd.DataFrame(df.Bins.str.split(',').tolist(), index=df.Chromosome).stack().reset_index().drop("level_1", 1)
    bins.columns = ["Chromosome", "Bin"]
    total_enriched = pd.DataFrame(df.TotalEnriched.str.split(',').tolist(), index=df.Chromosome).stack().reset_index().drop("level_1", 1)
    total_enriched = total_enriched.drop("Chromosome", 1)
    total_enriched.columns = ["values"]

    bw = bins.join(total_enriched)
    bw.loc[:, "Bin"] = bw.Bin.astype(int)
    bw.insert(2, "End", bw.Bin + bin_size - 1)

    bw = bw.set_index(["Chromosome", "Bin", "End"])

    _create_bigwig(bw, outpath, genome_size_dict)





# def create_cluster_bed(bed, bin_size):

#     rowdicts = []
#     for _, row in bed.iterrows():

#         total_enriched = row.TotalEnriched.split(",")
#         bins = row.Bins.split(",")
#         cluster_id = row.Index
#         chromosome = row.Chromosome

#         for _bin, number in zip(bins, total_enriched):
#             _bin, number = int(_bin), int(number)
#             rowdict = {"Chromosome": chromosome, "Start": _bin, "End": _bin +
#                        bin_size - 1, "Name": cluster_id, "Score": 0, "Strand": "."}

#             rowdicts.extend([rowdict] * number)

#     return pd.DataFrame.from_dict(rowdicts)["Chromosome Start End Name Score Strand".split()]
