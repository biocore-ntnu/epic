To support additional genomes that are available on UCSC, add the genome to the
list in `genome.snakefile`. Prepare the environment:

```
conda create -n effective-genome-size -c bioconda --file requirements.txt python=3
```

Activate the environment:

```
source activate effective-genome-size
```

Dry-run the workflow:

```
snakemake -n -s genome.snakefile
```

Set the $TMPDIR env variable if you want to use something other than /tmp.

Run the workflow with 8 cores:

```
snakemake -s genome.snakefile -j 8
```

NOTE: jellyfish, which is used for counting k-mers, uses a lot of RAM. For
example, k=100 in hg19 consumes 90GB of RAM and takes about 40 mins running on
4 cores.
