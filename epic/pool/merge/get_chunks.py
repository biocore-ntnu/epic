from epic.pool.merge.read_chunk import read_chunk


def get_chunks(file_indexes):
    """Read chunks from each file according to the indexes given.

    Keyword Arguments:

    file_indexes -- dict with file path as the key and
    (start_line, end_line) as the value.
    """

    chunks = []
    for epic_file, (start, end) in file_indexes.items():
        chunk = read_chunk(epic_file, (start, end))
        chunks.append(chunk)
    return chunks
