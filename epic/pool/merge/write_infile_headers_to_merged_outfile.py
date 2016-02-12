import pandas as pd
from io import StringIO
from os.path import basename


def get_index_and_headers(infiles, index_cols):
    """Get index and header names from files.

    Keyword Arguments:
    infiles -- files to read header from.
    index_cols -- columns containing index
    """

    all_columns = []
    for infile in infiles:
        df = pd.read_table(infile, header=0, index_col=index_cols, sep="\t")

        # replace "Island" header column with "basename(infile)_Island"
        # Important to distinguish the different Island columns later,
        # when merging several files
        basename_infile = basename(infile)
        infile_columns = [c if not c.startswith("Island") else
                          basename_infile + "_Island" for c in df.columns]

        all_columns.extend(infile_columns)
        indexes = df.index.names

    return "\t".join(list(indexes) + all_columns)


def write_infile_headers_to_merged_outfile(infiles, output_filename,
                                           index_cols):

    df_index_and_header = get_index_and_headers(infiles, index_cols)

    with open(output_filename, "w+") as h:
        h.write(df_index_and_header + "\n")
