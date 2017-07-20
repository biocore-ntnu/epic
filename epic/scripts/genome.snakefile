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

# genomes = ["hg38"] # "susScr3 susScr2".split() # ["danRer10"] # ['dm3', 'dm6', 'mm9', 'mm10', 'hg19', 'hg38']
genomes = [s.strip() for s in open("genome_names_no_patches.txt").readlines()]
ucsc_genomes = ['fr2', 'fr1', 'nomLeu3', 'nomLeu2', 'nomLeu1', 'aquChr2', 'rhiRox1', 'gorGor5', 'gorGor4', 'gorGor3', 'chlSab2', 'cavPor3', 'eriEur2', 'eriEur1', 'equCab2', 'equCab1', 'dipOrd1', 'petMar2', 'petMar1', 'braFlo1', 'anoCar2', 'anoCar1', 'galVar1', 'triMan1', 'calJac3', 'calJac1', 'oryLat2', 'geoFor1', 'pteVam1', 'myoLuc2', 'balAcu1', 'micMur2', 'micMur1', 'hetGla2', 'hetGla1', 'oreNil2', 'monDom5', 'monDom4', 'ponAbe2', 'priPac1', 'chrPic1', 'ailMel1', 'susScr3', 'susScr2', 'ochPri3', 'ochPri2', 'ornAna2', 'ornAna1', 'nasLar1', 'oryCun2', 'rn6', 'rn5', 'rn4', 'rn3', 'rheMac8', 'rheMac3', 'rheMac2', 'proCap1', 'sacCer3', 'sacCer2', 'sacCer1', 'strPur2', 'strPur1', 'aplCal1', 'oviAri3', 'oviAri1', 'sorAra2', 'sorAra1', 'choHof1', 'speTri2', 'saiBol1', 'gasAcu1', 'tarSyr2', 'tarSyr1', 'sarHar1', 'echTel2', 'echTel1', 'tetNig2', 'tetNig1', 'nanPar1', 'tupBel1', 'melGal5', 'melGal1', 'macEug2', 'cerSim1', 'xenTro7', 'xenTro3', 'xenTro2', 'xenTro1', 'taeGut2', 'taeGut1', 'danRer10', 'danRer7', 'danRer6', 'danRer5', 'danRer4', 'danRer3', 'ASM18939v1', 'WH8502_v1', 'WH0401_v1', 'WH0402_v1', 'ASM105083v1', 'ASM24754v1', 'RintRC_1']

readlengths = [36, 50, 75, 100]

try:
    tmpdir = os.environ['TMPDIR']
except KeyError:
    tmpdir = "/tmp"

effective_sizes = expand('effective_sizes/{genome}_{readlength}.txt', genome=ucsc_genomes, readlength=readlengths)
chromsizes = expand('chromsizes/{genome}.chromsizes', genome=ucsc_genomes)

rule all:
    input: effective_sizes + chromsizes

# Download .2bit file from UCSC
rule download:
    output: temp('fasta/{genome}.2bit')
    resources: instances = 1
    run:
        url = TWOBIT_PATTERN.format(genome=wildcards.genome)
        shell('wget -O {output}.tmp {url} && mv {output}.tmp {output}')

# Convert to FASTA
rule twoBitToFa:
    input: 'fasta/{genome}.2bit'
    output: temp('fasta/{genome}.fa')
    priority: 2
    shell:
        'twoBitToFa {input} {output}'

# Assume that most peak-calling will be performed on chromosomes without
# underscores in the name -- so filter out those chromosomes from the FASTA
rule trim_underscore_chroms:
    input: 'fasta/{genome}.fa'
    output: temp('trimmed_fasta/{genome}.trimmed.fa')
    log: 'trimmed_fasta/{genome}.trimmed.log'
    priority: 3
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
    threads: 25
    resources: instances = 1
    priority: 4
    params:
    shell:
        'epic-effective --read-length={wildcards.readlength} --nb-cpu={threads} --tmpdir {tmpdir} {input} > {output}.tmp && mv {output}.tmp {output}'


# Download chromsizes file
rule chromsizes:
    output: 'chromsizes/{genome}.chromsizes'
    resources: instances = 1
    run:
        url = CHROMSIZES_PATTERN.format(genome=wildcards.genome)
        shell('wget -O - {url} | grep -v "_" > {output}.tmp && mv {output}.tmp {output}')

# vim: ft=python
