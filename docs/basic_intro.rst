Basic intro
===========

epic is used to find differentially enriched regions in ChIP-Seq data. It is
especially well suited for proteins that bind diffusely over longer regions (such
as many histone modifications).

The epic domain caller has three features.

* finding enriched regions in ChIP-Seq data
* creating files that allow you to visualize the ChIP-Seq signal in genome browsers
* creating a matrix of read counts for differential analysis between
  experimental conditions

Let's look at these in turn.

Finding enriched regions in ChIP-Seq data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An enriched region is a region where the amount of ChIP signal is significantly
higher than the amount of input (background) signal.

The list of enriched regions is the main output-file from epic and is always
included (matrix and visualization files are optional).

In the `quick-start <quick_start.html>`_ example we found enriched regions and stored them in.

Now let us see what such a file looks like. Here the first five lines are shown.

.. code-block:: bash

    # epic -c control.bedpe -t examples/test.bed -o results.csv
    Chromosome Start End ChIP Input Score Fold_change P FDR
    chr1 23568400 23568599 2 0 13.54464274405703 2534.5007148630584 8.184731894908448e-11 6.319374393278151e-10
    chr1 26401200 26401399 2 0 13.54464274405703 2534.5007148630584 8.184731894908448e-11 6.319374393278151e-10
    chr1 33054800 33055399 2 0 14.288252844518933 844.8335716210196 2.207263610333191e-09 3.85690272963484e-09

The first line merely contains the command used to produce the output. The
second line contains the column names. From the third line out, the data we are
actually interested in is shown.

The three first columns give the location of the enriched region. The ChIP and
Input columns contain the number of ChIP and Input reads in those regions,
respectively.

The last four columns contain statistical information about how enriched the region is.

Visualization
~~~~~~~~~~~~~

To allow for visualization of the data in genome browsers, epic can output both
bed and bigwig files that allow you to explore and analyze the data.

The -bw flag takes a folder to store the bigwigs in (if the folder does not
exist, it is created for you).

.. code-block:: bash

   epic -t test.bed -c control.bedpe -o enriched_regions.csv -bw bigwigs

Now one bigwig for each analysed file is stored in the folder bigwigs

.. code-block:: bash

   ls bigwigs/
   # control.bw test.bw

Count-matrix
~~~~~~~~~~~~

For downstream statistical analyses, epic can output a matrix of read-counts.

You can ask for one with the -sm flag. (Remember to give it a .gz extension as
the results are stored as a gzipped file!)

.. code-block:: bash

   epic -t test.bed -c control.bedpe -o enriched_regions.csv -sm matrix.gz
