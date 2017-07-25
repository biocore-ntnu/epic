epic-merge
==========

epic-merge is used to merge two or more output-matrixes from epic. The result is
a count-matrix containing several experiments, so that producing an input-file
for differential analysis with state-of-the-art tools such as limma is
simplified.

By default, epic-merge only keeps the bins within a region epic considered
enriched, but you can tell it to keep all bins - although this obviously is much
more memory and time-consuming. epic-merge can also take one or more bed-files
with regions to consider enriched instead of those found by epic. This allows
you to use the epic toolsuite for differential analysis of experiments, even if
you found the enriched regions with other region callers such as macs2.

The problem epic-merge solves (relatively) quickly and efficiently is that of
horizontal concatenation.

Options
~~~~~~~

epic-merge has the following flags:

* **-m, --matrixes**

   The epic count-matrixes to be merged.

* **-r, --regions**

   One or more bed-files with regions to consider enriched. For each bin,
   epic-merge can keep a column with how many times (i.e. in how many files) it
   was considered enriched, so you should keep the peaks/regions in separate
   files if they indeed were from different experiments.

* **-e, --enriched-per-file**

   Whether to keep a column with info about in how many experiments a bin was
   considered enriched or not. The downstream tool epic-cluster can use this
   info.

* **-o, --output**

   Path to write the (gzipped) merged matrix to.

* **-cpu, --number-cores**

   The number of CPUs to use. Can at most use one per chromosome.


Example1: merging two output-matrixes from epic.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here we show how to merge two output-matrixes from epic; one from timepoint 0
and one from timepoint 12 of a time-series experiment.

.. code-block:: bash

   zcat H3K27_me3_12h.gz | head -5
   # Chromosome Bin Enriched H3K27me3/Exp1_12h_H3K27me3.bed H3K27me3/Exp2_12h_H3K27me3.bed Input/Exp1_12h_Input.bed Input/Exp2_12h_Input.bed
   # chr1 10000 0 0 0 2 5
   # chr1 10200 0 2 0 0 1
   # chr1 10400 0 0 0 1 0
   # chr1 13200 0 0 0 1 0
   zcat H3K27_me3_0h.gz | head -5
   # Chromosome Bin Enriched H3K27me3/Exp1_0h_H3K27me3.bed H3K27me3/Exp2_0h_H3K27me3.bed Input/Exp1_0h_Input.bed Input/Exp2_0h_Input.bed
   # chr1 10000 0 0 1 2 5
   # chr1 10200 0 0 2 1 1
   # chr1 10400 0 1 0 0 0
   # chr1 13200 0 0 0 0 1

   epic-merge -cpu 25 -m H3K27_me3_0h.gz H3K27_me3_12h.gz -o 0_12.gz

   Chromosome Bin TotalEnriched H3K27me3/Exp1_0h_H3K27me3.bed H3K27me3/Exp1_12h_H3K27me3.bed H3K27me3/Exp2_0h_H3K27me3.bed H3K27me3/Exp2_12h_H3K27me3.bed Input/Exp1_0h_Input.bed Input/Exp1_12h_Input.bed Input/Exp2_0h_Input.bed Input/Exp2_12h_Input.bed
   chr1 264400 1.0 1.0 1.0 5.0 3.0 0.0 0.0 0.0 0.0
   chr1 264600 1.0 3.0 3.0 3.0 5.0 1.0 2.0 0.0 1.0
   chr1 265200 1.0 3.0 6.0 2.0 2.0 4.0 3.0 2.0 0.0
   chr1 265400 1.0 2.0 1.0 2.0 0.0 1.0 2.0 1.0 3.0

Example2: Using custom region files to extract bins.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   head test1.bed
   # chr1	10050	10100
