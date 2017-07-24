epic-merge
==========

epic-merge is used to merge two or more output-matrixes from epic. The result is
a count-matrix containing several experiments, so that producing an input-file
for differential analysis with state-of-the-art tools such as limma is
simplified.

By default, epic-merge only keeps the bins within a region epic considered
enriched, but you can tell it to keep all bins - although this obviously is much
more memory and time-consuming. epic-merge can also take one or more bed-files
with regions to consider enriched, instead of those found by epic. This allows
you to do differential analysis of experiments, even if you found the enriched
regions with other tools such as macs2.

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

   zcat

Example2: Using custom region files to extract bins.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
