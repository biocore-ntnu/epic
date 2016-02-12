import pandas as pd


def read_chunk(epic_file, index, index_col=[0, 1]):

    start, end = index

    if start or end:
        nb_rows_to_read = end - start + 1
        header = None
        skiprows = start

    # if start, end both zero, just return empty dataframe
    if not (start or end):
        nb_rows_to_read = 0
        header = 0
        skiprows = 1

    return pd.read_table(epic_file,
                         sep="\t",
                         nrows=nb_rows_to_read,
                         skiprows=skiprows,
                         header=None,
                         index_col=index_col)
