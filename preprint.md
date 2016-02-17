Quickly find enriched diffused domains in ChIP-seq data using epic
------------------------------------------------------------------

epic is a software package for finding medium to diffusely enriched domains in
chip-seq data. It is a fast, parallel and memory-efficient implementation of the
SICER algorithm (1). By running epic on a set of data ("ChIP") files and a set
of control ("Input") files, epic is able to quickly identify differentially
enriched regions.

epic is an improvement over the original SICER by being faster, more memory
efficient, having the ability to use multiple cores, and being significantly
much easier to install and use.

One suggested way to find optimal parameters for window-based broad domain
ChIP-seq callers is to repeat the analyses and see which parameters result in
the most promising islands (2, 3) (without looking at p-values, of course.)

Furthermore, the original SICER contained hard-coded restrictions on the size of
the files it could analyze (for reasons of computational efficiency), that no
longer make sense in light of the explosive growth in the size of biological
datasets. epic also corrects a few (non-critical) bugs in the original
implementation of SICER and implements a test-suite to 1) ensure correct
behavior, 2) greatly lessen the time and effort required to accept community
contributions, and 3) serve as documentation for the implementation.

1. Introduction

ChIP-Seq is a method to analyze how proteins and DNA interact. The final result
of a ChIP-seq analysis is the DNA that bound to the protein under study (called
the "ChIP" library). Modern ChIP-studies often use of a control library of
unspecific DNA background ("Input"). ChIP-seq callers work by finding the
regions with more ChIP than Input signal and then perform statistical tests to
see whether the region should be considered enriched for ChIP signal.

2. Methods

2.1 Installation

epic is available for python2.7 and above. It can be installed from the Python
Package Index with `pip install bio-epic` or by cloning the repo at
https://github.com/endrebak/epic





For the mathematical details of how the algorithm works, ple

1. “A clustering approach for identification of enriched domains from histone
   modification ChIP-Seq data” Chongzhi Zang, Dustin E. Schones, Chen Zeng,
   Kairong Cui, Keji Zhao, and Weiqun Peng, Bioinformatics 25, 1952 - 1958
   (2009)
2. "csaw: a Bioconductor package for differential binding analysis of ChIP-seq
   data using sliding windows" Aaron T.L. Lun and Gordon K. Smyth
3. "Spatial clustering for identification of ChIP-enriched regions (SICER) to
   map regions of histone methylation patterns in embryonic stem cells." Xu S,
   Grullon S, Ge K, Peng W.

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
