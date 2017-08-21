.. epic documentation master file, created by
   sphinx-quickstart on Wed Jul 19 11:24:44 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

epic: find diffusely enriched domains in ChIP-Seq data
================================

epic is a modern reimplementation of the extremely popular SICER algorithm for
finding broad, diffusely enriched regions in ChIP-Seq data. epic focuses on
speed, innovation and ease of use.

In addtion to the classic SICER algorithm, epic contains an array of tools to
help you do differential analysis of ChIP-Seq data from multiple conditions.

Novel features
--------------

* epic creates output files that allow you to leverage state-of-the-art
  statistical software to do rigorous differential analyses of experimental
  conditions
* epic creates bigwigs and bed files that allow you to visualize and explore
  the results in genome browsers

Improvements
------------

epic contains a slew of improvements

* Python 2.7 and 3+ compatible
* much faster and multicore
* easy to install and use; exists both in Bioconda and PyPI
* accepts both single- and paired-end reads
* can analyse multiple ChIP and input files at the same time (and even mix
  paired and single end files in the same analysis)
* sensible defaults - only needs the files to analyse as command line args
* effective genome size autoselected based on read length (autodetected) and genome
* epic can compute the effective genome fraction for you, allowing you to
  analyze any genome for which you have a fasta file
* (soon) metadata for all genomes in UCSC is available so you only need to give
  the genome name for the species you wish to analyse
* accepts custom genomes and assemblies
* automatically tested, which safely and easily allows contributions
* continually used, updated and maintained by the author


-----


.. toctree::
   :caption: Getting Started
   :hidden:
   :maxdepth: 2

   installation
   quick_start

.. toctree::
   :caption: epic
   :hidden:
   :maxdepth: 2

   basic_intro
   options
   output_files

.. toctree::
   :caption: epic-tools
   :hidden:
   :maxdepth: 2

   epic_merge
   epic_cluster
   epic_count
   epic_blacklist


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
