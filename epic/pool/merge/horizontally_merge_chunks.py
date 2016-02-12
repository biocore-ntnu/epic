import pandas as pd

from epic.pool.merge.get_chunks import get_chunks


def horizontally_merge_chunks(file_chunk_delimiters):
    """Get chunks from files and merge them.

    Keyword Arguments:

    file_chunk_delimiters -- list of dicts where the keys are files and the values file
    chunk_delimiters. It might look something like this: [
    {"file1.txt": (0, 3), "file2.txt": (0, 1)},
    {"file1.txt": (4, 6), "file2.txt": (5, 7)},
    ...
    ]"""
    horizontally_merged_chunks = []
    for chunk_delimiters in file_chunk_delimiters:
        chunks = get_chunks(chunk_delimiters)
        merged_chunks = pd.concat(chunks, axis=1).fillna(0).astype(int)
        horizontally_merged_chunks.append(merged_chunks)
    return pd.concat(horizontally_merged_chunks)
