import logging

from natsort import natsorted

from epic.pool.index.get_bin_file_idx import index_bin_files

from epic.pool.merge.horizontally_merge_chunks import horizontally_merge_chunks
from epic.pool.merge.append_chunk_to_file import append_chunk_to_file
from epic.pool.merge.write_infile_headers_to_merged_outfile import write_infile_headers_to_merged_outfile

from epic.config import logging_settings  # pylint: ignore


def merge_files(infiles,
                output_filename,
                nb_lines_per_chunk=1e5,
                infile_index_cols=[0, 1]):
    """Horizontally concatenate many large files.

    Does not need to use much memory (usage depends on chunksize and
    number of files - smaller chunks/lower number of files equals less RAM
    used).

    Keyword Arguments:
    infiles            -- files to merge
    nb_lines_per_chunk -- int (default 100000)

    """

    logging.info("Files to merge:\n" + "\n".join(infiles))

    bin_file_indexes = index_bin_files(infiles, nb_lines_per_chunk)

    write_infile_headers_to_merged_outfile(infiles, output_filename,
                                           infile_index_cols)

    for chromosome, rowdicts in natsorted(bin_file_indexes.items()):
        logging.info("Merging chunks for chromosome " + chromosome)
        merged_chunk = horizontally_merge_chunks(rowdicts)
        append_chunk_to_file(merged_chunk, output_filename)
