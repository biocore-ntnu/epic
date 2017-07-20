Quick-start
================================

To use epic, you need at least one ChIP and one input bed file.
Let's download some test data from github:

.. code-block:: bash

   wget https://raw.githubusercontent.com/biocore-ntnu/epic/master/examples/control.bedpe
   wget https://raw.githubusercontent.com/biocore-ntnu/epic/master/examples/test.bed

Now you can run epic on these files with

.. code-block:: bash

   epic -t test.bed -c control.bedpe -o enriched_regions.csv

The -t means treatment (ChIP), while the -c means control (input). The -o is the
path where you want to store the results. If no output path is given, the
results are written to stdout.
