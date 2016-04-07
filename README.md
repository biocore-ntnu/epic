# epic

epic is a software package for finding medium to diffusely enriched domains in
chip-seq data. It is a fast, parallel and memory-efficient implementation of the
SICER algorithm. By running epic on a set of data ("ChIP") files and a set of
control ("Input") files, epic is able to quickly identify differentially
enriched regions.

epic is an improvement over the original SICER by being faster, more memory
efficient, multicore, and significantly much easier to install and use.

The MIT-licensed code is available at https://github.com/endrebak/epic

## Install

epic is available for python2.7 and above. It can be installed from the Python
Package Index with `pip install bioepic` or by cloning the repo at
https://github.com/endrebak/epic

## Usage

(Might be slightly out of date.)

```
epic

Diffuse domain ChIP-Seq caller based on SICER.
(Visit github.com/endrebak/epic for examples and help.)

Usage:
    epic [--input-string=STR] [--fragment-size=FRG] [--window-size=WIN] [--gaps-allowed=GAP]
         [--nb-cpu=CPU] [--fdr-cutoff=FDR] [--keep-duplicates] [--pandas-only] [--genome=GEN] FILE...
    epic --help

Arguments:
    FILE                        a list of chip _and_ input files

Options:
    -h --help                   show this help message
    -i STR --input-string=STR   case insensitive string used to distinguish input/control-files [default: input]
    -f FRG --fragment-size=FRG  estimated length of dna fragments [default: 150]
    -w WIN --window-size=WIN    size of bins [default: 200]
    -g GAP --gaps-allowed=GAP   number of gaps allowed [default: 3]
    -v GEN --genome=GEN         genome-version to use  [default: hg19]
    -q FDR --fdr-cutoff         false-discovery rate cutoff [default: 1.0]
    -n CPU --nb-cpu             number of cpus to use [default: 1]
    -k --keep-duplicates        do not delete duplicate reads with equal chromosome, start and end coordinates
    -p --pandas-only            use pandas for all taks (does not rely on GNU coretools, but is slower)

Note:
    The suggested settings for the different types of modifications are as following:
        H3K27me3: --window-size=200 --gaps=3
        H3K4me3:  --window-size=200 --gaps=1
```

## Brag brag brag

#### Actively developed

Will be maintained and further updated.

We hope to make further refinements to the actual algorithm and make it even better.

#### Functionality

epic accepts several input and ChIP files at the same time and accepts
block-zipped and gzipped input files.

Works on files of any size.

Works on all Python versions 2.7/3+.

#### Speed

epic can use one core per chromosome, which should give a speedup of < ~22-25
(differs by species) by itself. In addition, epic uses the Python science stack,
including Pandas, for almost all tasks, which means each core runs heavily
optimized C, Fortran and Cython code for further speed gains.

#### Memory

epic streams the data instead of loading it all into memory, which should result
in a much smaller memory footprint.

#### Usage

Instead of needing eleven command line arguments to run, epic contains sensible
defaults and only needs the files it is to analyze as parameters.

epic can be run from whichever location with files found anywhere on the disk.

## Version

This is a pre-alpha release. Please do aggressively report issues, quirks, complaints and anything that just feels slightly off to the issue tracker. Also please ask questions and make docrequests - there are loads of neat stuff I have not documented.

I have been using epic with great success.

## License

MIT

## Requirements

Python data science stack and a fairly recent version of Pandas (0.17 >=).
Python 2.7 or 3+

## Quickstart

```
pip install bioepic
# you only need git clone to get the test data
git clone https://github.com/endrebak/epic.git
epic -i control epic/examples/test.bed epic/examples/control.bed
```

<!-- ## Helper scripts -->

<!-- Various helpful utilities will be added to epic. Currently it includes: -->

<!-- #### epic-effective -->

<!-- A perennial question on bioinformatics sites is how to compute the effective -->
<!-- genome size for a genome. epic includes a script called `epic-effective` to do -->
<!-- just this. It can use multiple cores, but is not optimized for memory -->
<!-- consumption so run it on a bioinformatics cluster/server. Please share -->
<!-- your results on the issue tracker. Currently only finds effective genome size -->

## TODO

* Change to argparse so that can take separate ChIP/Input lists.
* Fix "bug" that prints chromosome start and end in output as floats.
* Print command line args exactly as recieved.
* Change output format to use tab as delimiter and underscore as within name delimiter ("P_value", not "P value")
* Add bam and paired end support

## NAQ/Various

#### Difference from the original SICER

Note that the island enriched threshold computation produces results that are < ~10% more conservative than in the original SICER.
This is due to numerics (summing many very small numbers is done in both implementations, albeit slightly differently).

This gives a different cutoff than the original, but produces virtually identical results (since epic and SICER produces the same candidate island list, with the same order, but epic selects slightly fewer islands from this list).

#### Why another ChIP-Seq domain caller?

MACS2 is great for narrow peaks, but epic performs better on diffuse domains. For medium size domains, such as PolII, our tests indicate that both perform about equally well, but epic uses only a fraction of the time.

#### Why not SICER?

SICER is a wonderful piece of software, but advances in the Python data science libraries has made it possible to implement it much more efficiently. Furthermore, SICER was not made to handle the mountains of data we have now; it simply cannot run on very large datasets due to (sensible) restrictions in the original implementation.

#### When is your paper coming out?

Dunno. We do not want to write a methods paper, but rather just include a section about epic in an appropriate biology paper sometime.

#### Credit

Chongzhi Zang, Dustin E. Schones, Chen Zeng, Kairong Cui, Keji Zhao and Weiqun Peng for the original SICER. Please consider citing their paper (*in addition* to our eventual paper) if you use epic. And if you use any (helper) scripts in SICER that are not included in epic you should of course cite the SICER paper!

Most of the improvements in epic were possible due to Python Science libraries that were not available when SICER was originally written. Thanks to the Pandas developers!

Endre Bakken Stovner for the implementation of epic.

#### Why the name epic?

It stands for electronic pic [sic] caller or epigenome cartographer, whichever you prefer. Or perhaps it isn't just another bogus bioinformatics acronym. Hope you find the name fitting.

#### Other pieces of software you might prefer

* [SICER](http://home.gwu.edu/~wpeng/Software.htm) - great diffuse domain ChIP-Seq caller (which epic is based on.)
* [SICERpy](https://github.com/dariober/SICERpy) - a wrapper around SICER for convenience/parallelism.
* [csaw](https://github.com/LTLA/csaw) - R package. Uses an approach to island finding that complements epic very well. Requires more statistical sophistication and programming skill to use.
* [MACS2](https://github.com/taoliu/MACS) - my preferred peak caller.

#### Bug in the original SICER?

I would appreciate it if anyone can send me the output islands with FDR from running the original SICER on the example data they provide. I seem to remember finding a bug in the original, but cannot be bothered to get SICER running again. (It wasn't critical or anything just giving slightly inaccurate counts for islands.)
