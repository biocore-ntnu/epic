"""
Snakemake workflow to compute effective genome sizes for various genomes for
various read lengths.

For each genome:
    - downloads .2bit files from UCSC
    - converts to FASTA
    - filters out chromosomes whose names contain "_"
    - computes effective genome size using epic-effective
    - downloads chrom.sizes from UCSC

If running on a cluster, you may want to set $TMPDIR.
"""
__author__ = "Ryan Dale"
__license__ = "MIT"

shell.executable("bash")

from Bio import SeqIO
import os

# On UCSC, some genome FASTAs (e.g., dm3) are provided as tarballs. Others
# (hg38) are gzipped single files. However, all genomes have a .2bit file, so
# we can rely on that for uniformity instead of special-casing a handful of
# genomes.
TWOBIT_PATTERN = 'http://hgdownload.soe.ucsc.edu/goldenPath/{genome}/bigZips/{genome}.2bit'
CHROMSIZES_PATTERN = 'http://hgdownload.soe.ucsc.edu/goldenPath/{genome}/bigZips/{genome}.chrom.sizes'

genomes = ['dm3', 'dm6', 'mm9', 'mm10', 'hg19', 'hg38']
readlengths = [36, 50, 75, 100]

try:
    tmpdir = os.environ['TMPDIR']
except KeyError:
    tmpdir = None


effective_sizes = expand('effective_sizes/{genome}_{readlength}.txt', genome=genomes, readlength=readlengths)
chromsizes = expand('chromsizes/{genome}.chromsizes', genome=genomes)

rule all:
    input: effective_sizes + chromsizes

# Download .2bit file from UCSC
rule download:
    output: 'fasta/{genome}.2bit'
    run:
        url = TWOBIT_PATTERN.format(genome=wildcards.genome)
        shell('wget -O {output}.tmp {url} && mv {output}.tmp {output}')

# Convert to FASTA
rule twoBitToFa:
    input: 'fasta/{genome}.2bit'
    output: 'fasta/{genome}.fa'
    shell:
        'twoBitToFa {input} {output}'

# Assume that most peak-calling will be performed on chromosomes without
# underscores in the name -- so filter out those chromosomes from the FASTA
rule trim_underscore_chroms:
    input: 'fasta/{genome}.fa'
    output: 'trimmed_fasta/{genome}.trimmed.fa'
    log: 'trimmed_fasta/{genome}.trimmed.log'
    run:
        logfile = open(log[0], 'w')
        recs = []
        for rec in SeqIO.parse(input[0], 'fasta'):
            logfile.write('%s' % rec.name)
            if '_' in rec.name:
                logfile.write(': filtered\n')
                continue
            logfile.write('\n')
            recs.append(rec)
        SeqIO.write(recs, output[0], 'fasta')
        logfile.close()


# Compute effective genome size for a particular read length.
rule effective_genome_size:
    input: 'trimmed_fasta/{genome}.trimmed.fa'
    output: 'effective_sizes/{genome}_{readlength}.txt'
    threads: 4
    params:
    shell:
        'epic-effective --read-length={wildcards.readlength} --nb-cpu={threads} --tmpdir {tmpdir} {input} > {output}.tmp && mv {output}.tmp {output}'


# Download chromsizes file
rule chromsizes:
    output: 'chromsizes/{genome}.chromsizes'
    run:
        url = CHROMSIZES_PATTERN.format(genome=wildcards.genome)
        shell('wget -O - {url} | grep -v "_" > {output}.tmp && mv {output}.tmp {output}')

# vim: ft=python
