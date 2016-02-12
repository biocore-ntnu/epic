def append_chunk_to_file(chunk, filename):
    """Append chunk to filename

    Keyword Arguments:
    chunk -- dataframe that should be appended to file
    filename -- path where dataframe chunk should be appended
    """

    # handle, filename = tempfile.mkstemp()
    # close(handle)
    chunk.to_csv(filename, sep="\t", mode="a", header=False)
