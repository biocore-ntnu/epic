Options
=======

epic allows for many flags to denote (optional) output files or to change the
execution of epic.

* **-t, --treatment**

   One or more ChIP files (bed or bedpe format)

* **-c, --control**

   One or more input files (bed or bedpe format)

* **-o, --outfile**

   File to write results to. By default sent to stdout.

* **-l, --log**

   File to write log messages to. Also written to stderr by default.

* **-cpu, --number-cores**

   The number of cores epic should use. Can at most take advantage of 1 core per
   strand per chromosome (i.e. 46 for humans). Default: 1

* **-gn, --genome**

   Which genome to analyze. By default hg19.

* **-k, --keep-duplicates**

   Keep reads mapping to the same position on the same strand within a library.
   The default is to remove all but the first duplicate (this is done once per
   file, not for all files collectively.)

* **-w, --window-size**

   Size of the windows (bins) used to scan the genome. This is also the smallest
   possible enriched region you can get. Default 200.

* **-g, --gaps-allowed**

   How many non-enriched windows in a row can be part of the same enriched
   region. If the number of gaps between two enriched windows is higher than this
   number, they are considered separate regions. Default: 3

* **-fs, --fragment-size**

   (Only used for single-end files) Size of the sequenced fragment. The center of
   the fragment will be used to calculate which window a read ended up in. So
   reads are shifted by fragment-size/2. Default 150.

* **-fdr, --false-discovery-rate-cutoff**

   Remove all regions with an FDR below cutoff. Note: this also affects which
   windows are considered enriched in the optional matrix output and which
   regions are included in the optional bed output.

* **-egs, --effective-genome-size**

   Use a different effective genome fraction than the one included in epic. Or
   include an egs for custom genomes that are not a part of epic. Should be a
   number between 0 and 1. Autoinferred by sampled read-length and genome by
   default.

* **-cs, --chromsizes**

   Set the chromosome lengths yourself in a file with two
   columns: chromosome names and sizes. Useful to analyze
   custom genomes, assemblies or simulated data. Only
   chromosomes included in the file will be analyzed.

* **-sm, --store-matrix**

   The path in which to store a gzipped matrix of read counts per window. One
   column for each ChIP and Input file.

* **-bw, --bigwig**

   The folder in which to store a RPKM-normalized bigwig for each file in the
   dataset. The bigwig shows how epic sees the data. It shows all bins, not just
   those in enriched regions.

* **-b, --bed**

   A summary bed file of all regions, for display in the UCSC genome browser or
   for use in downstream analyses with e.g. bedtools. The score field is
   log2(#ChIP/#Input) * 100 capped at a 1000.

* **-v, --version**

   Show version.
