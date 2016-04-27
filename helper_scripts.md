# Helper scripts

Various helpful utilities will be added to epic. Currently it includes:

#### epic-effective

A perennial question on bioinformatics sites is how to compute the effective
genome size for a genome. epic includes a script called `epic-effective` to do
just this. It can use multiple cores. Please share your results on the issue
tracker.

And the script to compute the effective genome size has the following CLI:

```
epic-effective
Compute the effective genome size from a fasta file.

(Visit github.com/endrebak/epic for examples and help.)

Usage:
    epic-effective --read-length=LEN [--nb-cpu=CPU] FILE
    epic-effective --help

Arguments:
    FILE                      Fasta genome
    -r LEN --read-length LEN  length of reads

Options:
    -h --help                 show this help message
    -n CPU --nb-cpu CPU       number of cores to use [default: 1]
```
