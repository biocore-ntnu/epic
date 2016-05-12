import pkg_resources


def get_genome_size_file(genome):
    return pkg_resources.resource_filename(
        "epic", "scripts/chromsizes/{}.chromsizes".format(genome))


def create_genome_size_dict(genome):
    """Creates genome size dict from string containing data."""

    size_file = get_genome_size_file(genome)
    size_lines = open(size_file).readlines()

    size_dict = {}
    for line in size_lines:
        genome, length = line.split()
        size_dict[genome] = int(length)

    return size_dict


EFFECTIVE_GENOME_SIZES = {"hg19": 0.74, "hg18": 0.74}


def get_effective_genome_length(genome):

    genome_size = sum(create_genome_size_dict(genome).values())
    return EFFECTIVE_GENOME_SIZES[genome] * genome_size
