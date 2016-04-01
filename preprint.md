Quickly find enriched diffuse domains in ChIP-seq data using epic
------------------------------------------------------------------

epic is a software package for finding medium to diffusely enriched domains in
chip-seq data. It is a fast, parallel and memory-efficient implementation of the
SICER algorithm with some new features (1). By running epic on a set of data
("ChIP") files and a set of control ("Input") files, epic is able to quickly
identify differentially enriched regions.

epic is an improvement over the original SICER by being faster, more memory
efficient, multicore, and being significantly much easier to install and use.

The MIT-licensed code is available at https://github.com/endrebak/epic

1. Introduction

ChIP-Seq is a method to analyze how proteins and DNA interact. The final result
of a ChIP-seq analysis is the DNA that bound to the protein under study (called
the "ChIP" library). Modern ChIP-studies often use of a control library of
unspecific DNA background ("Input"). ChIP-seq callers work by finding the
regions with more ChIP than Input signal and then perform statistical tests to
see whether the region should be considered enriched for ChIP signal.

1.1 Quickstart

$ pip install bioepic
$ # you only need git clone to get the test data
$ git clone https://github.com/endrebak/epic.git
$ epic -i control epic/examples/test.bed epic/examples/control.bed

2. Methods

2.1 Installation

epic is available for python2.7 and above. It can be installed from the Python
Package Index with `pip install bioepic` or by cloning the repo at
https://github.com/endrebak/epic

For the mathematical details of how the algorithm works, please see (1).

3. Improvements

In this section we highlight the areas of SICER where we feel epic improves
on the original software.

3.1 Speed

The original SICER could only use one core, while epic can use one core per
chromosome, which should give a speedup of ~22-25 (differs by species) by
itself. In addition, epic uses the Python science stack, including Pandas, for
almost all tasks, which means each core runs heavily optimized C, Fortran and
Cython code instead of the (largely) interpreted Python code which SICER uses.

One suggested way to find optimal parameters for window-based broad domain
ChIP-seq callers is to repeat the analyses and see which parameters result in
the most promising islands (2, 3) (without looking at p-values, of course.)
This is a task made feasible by the speedup offered by epic.

3.2 Memory

epic streams the data instead of loading it all into memory, which should result
in a much smaller memory footprint.

3.3 Usage

SICER needs eleven command line arguments to run, while epic contains sensible
defaults and only needs the files it is to analyze as input.

SICER is finicky about the location and naming of the input files, while epic is
very easy to use.

epic can be run from anywhere, while SICER was hard to run from other places
than the install directory.

3.4 Functionality

epic can be easily installed using the python package manager using `pip install
bioepic`.

epic works on both Python 2 and 3 and should also work on Windows with a special
--pandas-only flag (experimental), while SICER only works on Python2/Unix
machines.

The regular SICER script can only compare one ChIP file to one Input file, while
epic can compare an arbitrary number of ChIP files to Input files.

The original SICER contained hard-coded restrictions on the size of the files it
could analyze (for reasons of computational efficiency), that no longer make
sense in light of the explosive growth in the size of biological datasets.

epic implements a test-suite to 1) ensure correct behavior, 2) greatly lessen
the time and effort required to accept community contributions, and 3) serve as
documentation for the implementation.

3.5 New helper scripts

(Not implemented yet): epic contains a helper script to compute the "Effective genome size" for any
genome. The effective genome size differs by genome, and is a required parameter
for SICER, which makes it hard to use on non-standard genomes (anything but
hg19/hg38/mm9/mm10).

4. Features missing from epic

SICER contains a wealth of functionality, much of which has not been implemented
in epic since the authors did not use them. From such a perspective epic might
be a downgrade in certain ways. The authors are open to well-considered
suggestions for new functionality.

5. Gammel intro om ChIP-Seq (for dårlig/unødvendig)

ChIP-seq is performed by 1) binding proteins to DNA, 2) sonically shearing the
DNA into pieces and 3) capture the protein of interest, to finally 4) sequence
the DNA attached to the captured proteins. After the DNA has been sequenced, it
is 5) mapped to the genome, which shows where the protein bound. Since ChIP-seq
is performed on a population of cells, there will be some stochasticity
involved, and the proteins might bind differently in different cells, but the
more sequences align to a particular place in the genome, the more likely the
protein is to perform a biological function in that particular genomic location
in the cell or tissue of interest. In modern ChIP-seq experiments a control
library called "Input" is used; this is just the DNA collected from the same
type of cell using the same crosslinking and shearing, but without the
immunoprecipitation. This gives a background signal of unspecific DNA. Regions
where the chromatin-immunoprecipitated ("ChIP") signal is much stronger than the
background ("Input") signal are considered enriched, which means that it might
be a region of significance for the protein under study.


1. “A clustering approach for identification of enriched domains from histone
   modification ChIP-Seq data” Chongzhi Zang, Dustin E. Schones, Chen Zeng,
   Kairong Cui, Keji Zhao, and Weiqun Peng, Bioinformatics 25, 1952 - 1958
   (2009)
2. "csaw: a Bioconductor package for differential binding analysis of ChIP-seq
   data using sliding windows" Aaron T.L. Lun and Gordon K. Smyth
3. "Spatial clustering for identification of ChIP-enriched regions (SICER) to
   map regions of histone methylation patterns in embryonic stem cells." Xu S,
   Grullon S, Ge K, Peng W.
