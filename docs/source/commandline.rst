.. _commandline:

Command Line Interface
======================

Usage
******

Once installed, some of Espalier's methods can be run from the command line using the ``espalier`` executable.

To get an overview of the command line features, at the command line type ``espalier --help``.
::

	$ espalier --help
	Usage: espalier [OPTIONS] COMMAND [ARGS]...

	Options:
	  -h, --help  Show this message and exit.

	Commands:
	  maf        Compute maximum agreement forest (MAF) between two...
	  reconcile  Reconcile two discordant trees through their maximum...
	  spr        Compute subtree-prune-and-regraft distance between two trees...

Computing Maximum Agreement Forests
***********************************

Maximum agreement forests (MAFs) can be computed using the ``maf`` command. Here, we will use the trees *maf_example_tree1.tre* and *maf_example_tree2.tre* found in the *examples* folder. The basic usage requires that we pass the Newick files containing the two trees as arguments to the ``maf`` command:
::

	$ espalier maf maf_example_tree1.tre maf_example_tree2.tre
	Computing MAF...
	MAF contains 3 subtrees:
	(1:2.03543736451416,2:2.03543736451416):721.1792750257996;

	5:723.2147123903137;

	((3:21.0348587795342,(10:3.56730088556742,(4:1.6955616440807,8:1.6955616440807):1.87173924148672):17.46755789396678):702.1798536107794,((6:54.32114231628793,7:54.32114231628793):501.67017755551933,9:555.9913198718073):167.22339251850644);

The subtrees in the MAF are returned as Newick tree strings, one on each line.

By default, subtrees are written to *stdout* but subtrees can be written to a file instead by supplying a third argument specifying an output file:
::

	$ espalier maf maf_example_tree1.tre maf_example_tree2.tre maf_subtrees.txt

Finally, the subtrees in the MAF can be plotted using the ``--plot`` option.
::

	$ espalier maf maf_example_tree1.tre maf_example_tree2.tre --plot
	Computing MAF...
	MAF contains 3 subtrees:

	-------------------------
	Maximum agreement forest:
	-------------------------

	/---------------------------------------------------------------------------- 1
	+                                                                              
	\---------------------------------------------------------------------------- 2
	                                                                               
	                                                                               

	---5

	                   /-------------------------------------------------------- 3 
	/------------------+                                                           
	|                  |                  /------------------------------------- 10
	|                  \------------------+                                        
	|                                     |                  /------------------ 4 
	+                                     \------------------+                     
	|                                                        \------------------ 8 
	|                                                                              
	|                                                        /------------------ 6 
	|                                     /------------------+                     
	\-------------------------------------+                  \------------------ 7 
	                                      |                                        
	                                      \------------------------------------- 9

Computing SPR distances
***********************

Subtree-prune-and-regraft (SPR) distances between trees can also be computed from the command line. We simply use the ``spr`` command, passing the two tree files as arguments:
::

	$ espalier spr maf_example_tree1.tre maf_example_tree2.tre
	Computing SPR distance...
	SPR distance = 2

Tree reconciliation
*******************

Two discordant trees can be reconciled through their MAF using the ``reconcile`` command. More information about how the trees are reconciled can be found in the *Tree Reconciliation* section of the :ref:`primer`. The ``reconcile`` command requires three positional arguments: the first two arguments provide the pair of trees to be reconciled and the third argument a Fasta sequence file. Topological discordances between the two trees strongly supported by the sequence data will be retained while those not supported will be eliminated.

Here we will reconcile a pair of ML trees for two different genomic regions based on the sequence data from the second genomic region. We can think of the first tree as a reference or consensus tree and the second tree as an alternative local tree topology reconstructed for the second genomic region. Our goal then is to reconcile the conflicts between the reference tree and the local ML tree that are not strongly supported by the sequence data.
::

	$ espalier reconcile reconciler_MLTree_r1.tre reconciler_MLTree_r2.tre reconciler_r2.fasta
	Likelihood: -680.88
	(((4:0.105918,((8:0.102793,3:0.102793):0.0031250000000000028,5:0.105918):5.299999999985872e-07):2.86982,10:2.97574):1.00409,((7:0.309185,2:0.309185):0.04671170000000002,(1:0.309185,(6:0.1545925,9:0.1545925):0.1545925):0.04671170000000002):3.6239383);

Usually the ``reconcile`` command returns a set of reconciled trees. The reconciled trees are sorted in descending order based on the likelihood of the sequence data given each tree. However, in this example, only one reconciled tree with significant support was found.

Reconciled trees are written to *stdout* by default, but trees can be written to a file by supplying a fourth argument specifying an output file:
::

	$ espalier reconcile reconciler_MLTree_r1.tre reconciler_MLTree_r2.tre reconciler_r2.fasta reconciled_trees.txt

One important optional parameter is ``--lbr``, which sets the lower-bound-ratio used to constrain the search path explored by the reconciliation algorithm to only include reconciled trees strongly supported by the sequence data. We can however use a smaller lower-bound-ratio to broaden the search path to include reconciled trees that are less strongly supported by the sequence data. Here setting the ``--lbr`` parameter to 0.0001 instead of 0.1 returns two plausible reconciled trees instead of one: 
::

	espalier reconcile reconciler_MLTree_r1.tre reconciler_MLTree_r2.tre reconciler_r2.fasta --lbr 0.0001

	Tree log likelihood: -680.88
	(((4:0.105918,((8:0.102793,3:0.102793):0.0031250000000000028,5:0.105918):5.299999999985872e-07):2.86982,10:2.97574):1.00409,((7:0.309185,2:0.309185):0.04671170000000002,(1:0.309185,(6:0.1545925,9:0.1545925):0.1545925):0.04671170000000002):3.6239383);

	Tree log likelihood: -689.09
	((7:0.309185,(2:0.309185,(1:0.309185,(6:0.1545925,9:0.1545925):0.1545925):0.0):0.0):3.67065,((4:0.105918,((8:0.102793,3:0.102793):0.0031250000000000028,5:0.105918):5.299999999985872e-07):2.86982,10:2.97574):1.00409);

ARG reconstruction
******************

ARG reconstruction is not currently supported through the command line interface. This feature should be coming soon in Espalier version 0.2.


