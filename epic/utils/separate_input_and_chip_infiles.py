from fnmatch import fnmatch


def separate_input_and_chip_infiles(bed_files, input_pattern="input"):
    """Split list of files depending on whether name contains input_pattern.

    Since docopt cannot take two lists of variadic arguments, splitting it
    based on a pattern is a way to work around this.

    Keyword Arguments:
    input_pattern -- pattern to split on
    bed_files     -- list of files to split
    """

    chip_files, input_files = [], []
    input_pattern = "*" + input_pattern + "*"
    input_pattern = input_pattern.upper()

    for bed_file in bed_files:
        if fnmatch(bed_file.upper(), input_pattern):
            input_files.append(bed_file)
        else:
            chip_files.append(bed_file)

    if len(chip_files) == 0:
        raise ValueError(
            "All files considered input. This means all your files match the"
            " case insensitive pattern: " + input_pattern)
    if len(input_files) == 0:
        raise ValueError(
            "No input files found. This means none of your files match"
            " the case insensitive pattern: " + input_pattern)

    return chip_files, input_files
