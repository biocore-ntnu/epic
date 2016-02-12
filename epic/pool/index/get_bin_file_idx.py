import logging
from collections import OrderedDict
from os.path import basename

from epic.pool.index.get_genome_locations import get_union_of_all_genome_locations
from epic.pool.index.split_genome_regions_into_chunks import split_genome_regions_into_chunks
from epic.pool.index.find_chunks_start_and_ends import find_chunks_start_and_ends
from epic.pool.index.create_index_for_epic_bin_file import create_index_for_epic_bin_file


def index_bin_files(epic_bin_files, nb_lines_per_chunk=1e5):
    """Create index for each bin file.

    First a main index is created, the union of all genome locations in any
    file.

    This index is split into chunks.

    Each file is given a list of (line_nb_start, line_nb_end) tuples, i.e. file
    regions. The file regions indicated by the line numbers should contain the
    same (but not necessarily all) the genomic regions of a chunk. There should
    be as many index tuples as chunks, and chunk nb. X should correspond to
    index tuple nb. X.

    (This splitting is done to allow for horizontal concatenation of files too
    large for RAM.)

    Keyword Arguments:
    epic_bin_files -- list of filepaths to create index for.
    """
    logging.info("Getting union of all genome bins in all files.")
    main_indexes = get_union_of_all_genome_locations(epic_bin_files)

    logging.info("Creating chunks of size: " + str(nb_lines_per_chunk))
    genome_chunks = split_genome_regions_into_chunks(main_indexes,
                                                     nb_lines_per_chunk)

    chunk_indexes = find_chunks_start_and_ends(genome_chunks)

    bin_file_indexes = OrderedDict()
    for bin_file in epic_bin_files:
        bin_file_index = create_index_for_epic_bin_file(bin_file,
                                                        chunk_indexes)
        bin_file_indexes[bin_file] = bin_file_index

    rowdicts_per_chromosome = _rearrange_dict(epic_bin_files, chunk_indexes,
                                              bin_file_indexes)

    return rowdicts_per_chromosome


def _rearrange_dict(epic_bin_files, chunk_indexes, bin_file_indexes):
    """Rearrange the bin_file_indexes.

    The bin_file_indexes are arranged into dicts of files. We want to
    interleave the files per genomic bin.

    From this (file -> dict chromo to list of chunks per file)

    od = OrderedDict()
    od["df1_file.txt"] = OrderedDict([('chr1', [(0, 3), (4, 6), (7, 9)]), (
        'chr2', [(10, 13), (14, 16), (17, 19)])])
    od["df2_file.txt"] = OrderedDict([('chr1', [(0, 1), (2, 3), (4, 6)]), (
        'chr2', [(7, 7), (7, 8), (8, 9)]), ('chr3', [(10, 10)])])

    To this chr -> list of dicts where files keys, where each dict represents one chunk:

    od = OrderedDict([('chr1', [OrderedDict([(epic_bin_file_path, (1, 4)), (
        epic_bin_file_path2, (1, 2))]), OrderedDict([(epic_bin_file_path, (
            5, 7)), (epic_bin_file_path2, (3, 4))]), OrderedDict(
                [(epic_bin_file_path, (8, 10)), (epic_bin_file_path2, (5, 7))
                 ])]), ('chr2', [OrderedDict([(epic_bin_file_path, (11, 14)), (
                     epic_bin_file_path2, (8, 8))]), OrderedDict([(
                         epic_bin_file_path, (15, 17)), (epic_bin_file_path2, (
                             9, 9))]), OrderedDict([(epic_bin_file_path, (
                                 18, 20)), (epic_bin_file_path2, (10, 10))])]),
                        ('chr3', [OrderedDict([(epic_bin_file_path2, (11, 11))
                                               ])])])


    Keyword Arguments:
    epic_bin_files --
    chunk_indexes  --

    """

    rowdicts_per_chromosome = OrderedDict()
    for chromosome in chunk_indexes:
        logging.info("Rearranging indexes for chromosome: " + chromosome)
        nb_chunks = len(chunk_indexes[chromosome])
        rowdicts = []
        for chunk_nb in range(nb_chunks):
            rowdict = OrderedDict(
            )  # Needs to be ordered for the merger to know the order of the cols
            for bin_file in epic_bin_files:
                if chromosome not in bin_file_indexes[bin_file]:
                    continue

                bin_file_idx = bin_file_indexes[bin_file][chromosome][chunk_nb]

                rowdict[bin_file] = bin_file_idx
            rowdicts.append(rowdict)
        rowdicts_per_chromosome[chromosome] = rowdicts

    return rowdicts_per_chromosome
