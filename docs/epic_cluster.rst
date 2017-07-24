epic-cluster
============

When doing differential analysis of multiple ChIP-Seq experiments, one problem
is to choose the regions that should be tested for differential enrichment.
This is the problem epic-cluster tries to solve.

epic-cluster takes one merged-matrix created by the epic-merge tool and clusters
the genomic bins within it to produce regions that should be tested for
differential enrichment.

Several algorithms might be added to the tool; for now epic-merge contains one,
which we have called trunks, flanks and valleys.

Algorithm: trunks, flanks and valleys
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This algorithm first looks for genomic bins that are at most
`--distance-allowed` apart. These bins are grouped together. Then these bins are
further divdied into subregions based on the number of experiments in which the
bins within the subregion were considered enriched.

Example
~~~~~~~

We start with an epic-merge matrix that only includes two experiments:

.. code-block:: bash

  zcat
  Chromosome Bin TotalEnriched Exp1_0h_H3K27me3.bed Exp1_12h_H3K27me3.bed Exp2_0h_H3K27me3.bed Exp2_12h_H3K27me3.bed Exp1_0h_Input.bed Exp1_12h_Input.bed Exp2_0h_Input.bed Exp2_12h_Input.bed
  chr1 264400 1.0 1.0 1.0 5.0 3.0 0.0 0.0 0.0 0.0
  chr1 264600 1.0 3.0 3.0 3.0 5.0 1.0 2.0 0.0 1.0
  chr1 265200 1.0 3.0 6.0 2.0 2.0 4.0 3.0 2.0 0.0
  chr1 265400 1.0 2.0 1.0 2.0 0.0 1.0 2.0 1.0 3.0
  chr1 265800 1.0 7.0 12.0 6.0 7.0 0.0 2.0 4.0 3.0
  chr1 268400 1.0 5.0 5.0 6.0 5.0 6.0 3.0 3.0 3.0
  chr1 268600 1.0 2.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
  chr1 269000 1.0 4.0 1.0 0.0 2.0 1.0 3.0 1.0 0.0
  chr1 269200 1.0 0.0 1.0 3.0 2.0 0.0 0.0 2.0 1.0

Since it only includes two experiments, the TotalEnriched column can only
include 1 or 2 (which means the bin was enriched in 1 or 2 experiments) [#]_.

.. [#] If epic-merge was used with the option `--keep-nonenriched`, the
       TotalEnriched column can also include 0.

.. code-block:: bash



Options
~~~~~~~

* **-m, --matrix**

   A single epic-merge output matrix.

* **-t, --trunk-diff**

   The difference in number of enriched regions that decides when to consider a
   subregion a flank or valley, instead of a trunk.

* **-d, --distance-allowed**

   The number of nucleotides before two bins are considered to be in separate
   clusters.

* **-b, --bin-size**

   The bin-size used in the input matrix file. Auto-inferred by default.

* **-cpu, --number-cores**

   The number of cores to use. Can at most use one per chromosome
