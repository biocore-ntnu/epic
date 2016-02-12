import pandas as pd


def get_genome_locations_in_file(df_path):

    return pd.read_table(df_path, header=0, usecols=[0, 1])


def get_union_of_all_genome_locations(df_paths):
    """Get the union of the genome locations in all the files."""

    df_paths_iter = iter(df_paths)

    genome_locations_pool = get_genome_locations_in_file(next(df_paths_iter))

    for df_path in df_paths_iter:
        genome_locations = get_genome_locations_in_file(df_path)
        genome_locations_pool = genome_locations_pool.append(
            genome_locations).drop_duplicates()

    return genome_locations_pool.sort_values(["Chromosome", "Bin"
                                              ]).reset_index(drop=True)
