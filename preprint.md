Quickly find enriched diffused domains in ChIP-seq data using epic
------------------------------------------------------------------

epic is a software package for finding diffusely enriched domains in chip-seq
data. It is a fast, parallel and memory-efficient implementation of the SICER
algorithm (1). By running epic on a set of data ("ChIP") files and a set of
control ("Input") files, epic is able to quickly identify differentially
enriched regions.

epic is an improvement over the original SICER by being faster, more memory
efficient, having the ability to use multiple cores, and being significantly
much easier to install and use. Furthermore, the original SICER contained
hard-coded restrictions on the size of the files it could analyze (for reasons
of computational efficiency), that no longer make sense in light of the
explosive growth in the size of biological datasets. epic also corrects a few
(non-critical) bugs in the original implementation of SICER and implements a
test-suite to 1) ensure correct behavior, 2) greatly lessen the time and effort
required to accept community contributions.

Introduction

ChIP-Seq is a method to analyze how proteins and DNA interact. ChIP-seq is
performed by 1) binding proteins to DNA, 2) sonically shearing the DNA into
pieces and 3) picking out the protein of interest, to finally 4) sequence the
DNA attached to the captured proteins. After the DNA has been sequenced, it is
5) mapped to the genome, which leads to a map of where the protein bound. Since
ChIP-seq is performed on a population of cells, there will be some stochasticity
involved in where the protein is bound in any particular cell, but the more
sequences align to a particular place in the genome (i.e. the region is
enriched), the more likely the protein is to perform a biological function in
that particular genomic location in the cell or tissue of interest. In modern
ChIP-seq experiments a control library called "Input" is used; this is just the
DNA collected from the same type of cell using the same crosslinking and
shearing, but without the immunoprecipitation. This gives a background signal of
unspecific DNA. Regions where the chromatin-immunoprecipitated ("ChIP") signal
is much stronger than the background ("Input") signal are considered enriched.

This is how ChIP-seq peak/domain callers work: finding the regions with more
ChIP than Input signal and then performing statistical tests according to some
model, to see whether the region should be considered enriched for ChIP signal
(meaning that the enriched location is a place where the protein of interest is
likely to interact.)

For the mathematical details of how the algorithm works, ple

1. “A clustering approach for identification of enriched domains from histone
   modification ChIP-Seq data” Chongzhi Zang, Dustin E. Schones, Chen Zeng,
   Kairong Cui, Keji Zhao, and Weiqun Peng, Bioinformatics 25, 1952 - 1958
   (2009)
