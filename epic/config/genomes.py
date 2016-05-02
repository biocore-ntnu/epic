CHROMOSOME_SIZE_DICT = {"hg19": """chr1	249250621
chr2	243199373
chr3	198022430
chr4	191154276
chr5	180915260
chr6	171115067
chr7	159138663
chr8	146364022
chr9	141213431
chr10	135534747
chr11	135006516
chr12	133851895
chr13	115169878
chr14	107349540
chr15	102531392
chr16	90354753
chr17	81195210
chr18	78077248
chr19	59128983
chr20	63025520
chr21	48129895
chr22	51304566
chrM	16571
chrX	155270560
chrY	59373566""",
                        "hg18": """chr1     247249719
 chr2     242951149
 chr3     199501827
 chr4     191273063
 chr5     180857866
 chr6     170899992
 chr7     158821424
 chr8     146274826
 chr9     140273252
 chr10    135374737
 chr11    134452384
 chr12    132349534
 chr13    114142980
 chr14    106368585
 chr15    100338915
 chr16    88827254
 chr17    78774742
 chr18    76117153
 chr19    63811651
 chr20    62435964
 chr21    46944323
 chr22    49691432
 chrX     154913754
 chrY     57772954
 chrM     16571"""}


def create_genome_size_dict(genome):
    """Creates genome size dict from string containing data."""

    size_lines = CHROMOSOME_SIZE_DICT[genome].splitlines()

    size_dict = {}
    for line in size_lines:
        genome, length = line.split()
        size_dict[genome] = int(length)

    return size_dict


EFFECTIVE_GENOME_SIZES = {"hg19": 0.74, "hg18": 0.74}


def get_effective_genome_length(genome):

    genome_size = sum(create_genome_size_dict(genome).values())
    return EFFECTIVE_GENOME_SIZES[genome] * genome_size
