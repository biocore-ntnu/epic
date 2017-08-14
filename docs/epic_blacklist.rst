epic-blacklist
==========

epic-blacklist takes one or more ChIP-files that contain data unrelated to the
experiment (ChIP from another species for example) and finds bins where a lot of
reads still bind. Bins with a statistically significant number of reads is
computed according to a Poisson model. These bins are written as a bed-file to
stdout.

* **-i, --infiles**

  One or more bed/bedpe files to count reads in.

* **-o, --outfile**

   File to write results to. By default sent to stdout.

* **-cpu, --number-cores**

   The number of cores epic should use. Can at most take advantage of 1 core per
   strand per chromosome (i.e. 46 for humans). Default: 1

* **-gn, --genome**

   Which genome to analyze. By default hg19.

* **-k, --keep-duplicates**

   Keep reads mapping to the same position on the same strand within a library.
   The default is to remove all but the first duplicate (this is done once per
   file, not for all files collectively.)

* **-fs, --fragment-size**

   (Only used for single-end files) Size of the sequenced fragment. The center of
   the fragment will be used to calculate which window a read ended up in. So
   reads are shifted by fragment-size/2. Default 150.

* **-cs, --chromsizes**

   Set the chromosome lengths yourself in a file with two
   columns: chromosome names and sizes. Useful to analyze
   custom genomes, assemblies or simulated data. Only
   chromosomes included in the file will be analyzed.

* **-f, --fdr**

   Cut-off to consider a bin enriched (Default: 0.05)

* **-egf, --effective-genome-fraction**

   Use a different effective genome fraction than the one included in epic. Or
   include an egf for custom genomes that are not a part of epic. Should be a
   number between 0 and 1. Autoinferred by sampled read-length and genome by
   default.
