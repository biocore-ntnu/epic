
# epic: diffuse domain ChIP-Seq caller based on SICER

[![Build Status](https://travis-ci.org/endrebak/epic.svg?branch=master)](https://travis-ci.org/endrebak/epic) [![Coverage Status](https://coveralls.io/repos/github/endrebak/epic/badge.svg?branch=master)](https://coveralls.io/github/endrebak/epic?branch=master) [![Code Health](https://landscape.io/github/endrebak/epic/master/landscape.svg?style=flat)](https://landscape.io/github/endrebak/epic/master)  [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/bioepic/README.html)

epic is a software package for finding medium to diffusely enriched domains in
chip-seq data. It is a fast, parallel and memory-efficient implementation of the
incredibly popular SICER algorithm. By running epic on a set of data ("ChIP")
files and control ("Input") files, epic is able to quickly identify
differentially enriched regions.

epic is an improvement over the original SICER by being faster, more memory
efficient, multicore, and significantly much easier to install and use.

The MIT-licensed code is available at https://github.com/endrebak/epic
Please post issues, questions, docrequests or whatever is on your mind on the [issue tracker](https://github.com/endrebak/epic/issues).

<!-- markdown-toc start - Don't edit this section. Run M-x markdown-toc-generate-toc again -->
**Table of Contents**
- [epic: diffuse domain ChIP-Seq caller based on SICER](#epic-diffuse-domain-chip-seq-caller-based-on-sicer)
    - [Install](#install)
    - [Changelog](#changelog)
    - [Improvements](#improvements)
    - [Version](#version)
    - [License](#license)
    - [Requirements](#requirements)
    - [Helper scripts](#helper-scripts)
    - [TODO](#todo)
    - [Usage](#usage)
    - [Credit](#credit)
    - [NAQ/Various](#naqvarious)

<!-- markdown-toc end -->

## Install

epic is available for python2.7 and above. It can be installed from the Python
Package Index with `pip install bioepic`, with bioconda using `conda install bioepic` ([bioconda setup instructions](http://bioconda.github.io/index.html#setup)) or by cloning the repo at
https://github.com/endrebak/epic

## Changelog

```
# 0.0.7 (06.05.16)
- fix default value of keep duplicates (was True, now is False)
- fix max-gap bug (max gap used was actually max-gap - 1)
- (internal: no more functions with > 5 arguments)
- (the difference between SICER and epic I thought was a due to a bug, was caused
  by them having different genome versions as default)

# 0.0.6 (10.04.16)
- add script to compute the effective genome size (thanks to Heng Li for alerting me to the existence of jellyfish2)

# 0.0.5 (07.04.16)
- remove Pandas only option. Now unix is required.
- change to argparse instead of docopt (thanks Dario Beraldi)
- separate args for control and input files (thanks Dario Beraldi)
- print cl-args exactly as recieved (thanks to Istvan Albert for pointing out glitch)
- add version info
- add support for bam files (requires bedtools)
- fix "bug" that prints chromosome start and end in output as floats.
- change output format to use tab as delimiter and underscore as within name delimiter ("P_value", not "P value")
```

## Quickstart

```
$ pip install bioepic
$ # you only need git clone to get the test data
$ git clone https://github.com/endrebak/epic.git
$ # -t is treatment files, -i input files
$ epic -t epic/examples/test.bed -i epic/examples/control.bed
```

## Improvements

#### Actively developed

Will be maintained and further updated.

We hope to make further refinements to the actual algorithm and make it even better.

#### Functionality

epic accepts several input and ChIP files at the same time and accepts bed files
(both block-zipped and gzipped) and bam files.

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

## License

MIT

## Requirements

Python data science stack and a fairly recent version of Pandas (0.17 >=).
Python 2.7 or 3+.
Various unix tools found on all major distributions.

For the effective genome size script [jellyfish2](https://github.com/gmarcais/Jellyfish) is required.
To read bam files [bedtools](https://github.com/arq5x/bedtools2) is required.

## Helper scripts

See [this page](helper_scripts.md) for the various helper scripts that are a part of epic.

## TODO

* Clean up tests
* Add paired end support
* Add more examples of usage
* Add bigwig output for islands
* Improve logging messages
* Write proper docs
* Explain the effective genome size and test how much it matters
* Test island threshold difference between SICER and epic on different datasets

## Usage

(Might be slightly out of date.)

```
usage: epic [-h] --treatment TREATMENT [TREATMENT ...] --control CONTROL
            [CONTROL ...] [--number-cores NUMBER_CORES] [--genome GENOME]
            [--keep-duplicates KEEP_DUPLICATES] [--window-size WINDOW_SIZE]
            [--gaps-allowed GAPS_ALLOWED] [--fragment-size FRAGMENT_SIZE]
            [--false-discovery-rate-cutoff FALSE_DISCOVERY_RATE_CUTOFF]
            [--version]

Diffuse domain ChIP-Seq caller based on SICER. (Visit github.com/endrebak/epic
for examples and help.)

optional arguments:
  -h, --help            show this help message and exit
  --treatment TREATMENT [TREATMENT ...], -t TREATMENT [TREATMENT ...]
                        Treatment (pull-down) file(s) in
                        bam/bed/bed.gz/bed.bz2 format.
  --control CONTROL [CONTROL ...], -c CONTROL [CONTROL ...]
                        Control (input) file(s) in bam/bed/bed.gz/bed.bz2
                        format.
  --number-cores NUMBER_CORES, -cpu NUMBER_CORES
                        Number of cpus to use. Can use at most one per
                        chromosome. Default: 1.
  --genome GENOME, -gn GENOME
                        Which genome to analyze. Default: hg19.
  --keep-duplicates KEEP_DUPLICATES, -k KEEP_DUPLICATES
                        Keep reads mapping to the same position on the same
                        strand within a library. Default is to remove all but
                        the first duplicate.
  --window-size WINDOW_SIZE, -w WINDOW_SIZE
                        Size of the windows to scan the genome. WINDOW_SIZE is
                        the smallest possible island. Default 200.
  --gaps-allowed GAPS_ALLOWED, -g GAPS_ALLOWED
                        Multiple of window size used to determine the gap
                        size. Must be an integer. Default: 3.
  --fragment-size FRAGMENT_SIZE, -fs FRAGMENT_SIZE
                        Size of the sequenced fragment. The center of the the
                        fragment will be taken as half the fragment size.
                        Default 150.
  --false-discovery-rate-cutoff FALSE_DISCOVERY_RATE_CUTOFF, -fdr FALSE_DISCOVERY_RATE_CUTOFF
                        Remove all islands with an FDR below cutoff. Default
                        1.0, that is, all islands included.
  --version, -v         show program's version number and exit
```

## Credit

Chongzhi Zang, Dustin E. Schones, Chen Zeng, Kairong Cui, Keji Zhao and Weiqun Peng for the original SICER. Please consider citing their paper (*in addition* to our eventual paper) if you use epic. And if you use any (helper) scripts in SICER that are not included in epic you should of course cite the SICER paper!

Most of the improvements in epic were possible due to Python Science libraries that were not available when SICER was originally written. Thanks to the Pandas developers!

#### Author

Endre Bakken Stovner

#### Contributors

* Pål Sætrom (algorithmic/theoretical discussions, endless patience)
* Dario Beraldi (argparsing, bam support)
* Ryan Dale (bioconda, ideas, genome info script)

## NAQ/Various

Answers to some questions no-one has ever asked me.

#### Why is the SICER algorithm so great?

The wonderful thing about the SICER algorithm is that is very careful about dropping windows with few reads in them. All ChIP-seq callers I know have some preprocessing step where this is done liberally. SICER pools these windows together and gives them a composite score, allowing very long stretches of very diffuse signal to be detected.

#### Why another ChIP-Seq domain caller?

MACS2 is great for narrow peaks, but epic performs better on diffuse domains. For medium size domains, such as PolII, our tests indicate that both perform about equally well, but epic uses only a fraction of the time.

#### Why not SICER?

SICER contains a great algorithm and is a wonderful piece of software, but advances in the Python data science libraries has made it possible to implement it much more efficiently. Furthermore, SICER was not made to handle the mountains of data we have now; it simply cannot run on very large datasets due to (sensible) restrictions in the original implementation.

#### When is your paper coming out?

Dunno. We do not want to write a methods paper, but rather just include a section about epic in an appropriate biology paper sometime.

#### Why the name epic?

It stands for electronic pic [sic] caller or epigenome cartographer, whichever you prefer. Or perhaps it isn't just another bogus bioinformatics acronym. Hope you find the name fitting.

But suggestions for better names accepted. On paper I liked the epi/epic/epigenetics link but now when I hear it it sounds so boastful I cringe. exorcised sounds like a slight on the original software... Mad MACS?

#### Which other ChIP-Seq callers do you use?

* [SICER](http://home.gwu.edu/~wpeng/Software.htm) - great diffuse domain ChIP-Seq caller (which epic is based on.)
* [SICERpy](https://github.com/dariober/SICERpy) - a wrapper around SICER for convenience/parallelism. Stole some good ideas from there.
* [csaw](https://github.com/LTLA/csaw) - R package. Uses an approach to island finding that complements epic very well. Requires more statistical sophistication and programming skill to use.
* [MACS2](https://github.com/taoliu/MACS) - my preferred peak caller.
