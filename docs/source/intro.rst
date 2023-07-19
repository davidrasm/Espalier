Introduction
============
 
Espalier is a Python package for working with discordant phylogenetic trees using maximum agreement forests. The package can be used to compute maximum agreement forests, compute subtree-prune-and-regraft (SPR) distances, reconcile discordant trees and reconstruct approximate ancestral recombination graphs (ARGs).

For more information on the algorithms behind Espalier please see:

`Rasmussen, D.A. and Guo, F. Espalier: Efficient tree reconciliation and ARG reconstruction using maximum agreement forests, Systematic Biology, 2023; syad040 <https://doi.org/10.1093/sysbio/syad040>`_


Motivation
**********

**Espalier** *(noun): the ancient agricultural art of controlling woody plant growth by pruning and attaching branches to a common frame such as a trellis.*

Phylogenetic trees reconstructed from sequence data for different genomic regions often have conflicting or *discordant* tree topologies. Discordance may arise biologically from recombination or horizontal transfers of genetic material where individuals inherit genomic material from different parents. In this case, phylogenetic relationships among individuals/taxa may vary across the genome. 

Espalier uses maximum agreement forests (MAFs) to identify and reconcile topological discordance between trees. Given two discordant trees, a MAF identifies a set of topologically concordant subtrees present in both trees. MAFs are therefore very useful for working with discordant trees as they allow us to identify what relationships agree and which conflict between trees.

MAFs can also help identify which lineages recombined. To see this, consider the topological effects of a recombination event where an individual inherits genetic material from different parents to the left and right of a recombination breakpoint. In a phylogenetic tree, a lineage that has undergone recombination may *attach* to different parent lineages on either side of a breakpoint. The topological effect of recombination is therefore analogous to a *subtree-prune-and-regraft* (SPR) move where a subtree is cut from one branch and regrafted to another branch elsewhere in the tree. 

The subtrees in a MAF therefore identify lineages that may have moved due to recombination. MAFs are especially useful in this regard because the number of subtrees cut to obtain the MAF (i.e the SPR distance) tells us the minimum number of recombination events necessary to explain the discordance between trees. SPR distances can therefore be more useful than other distance metrics like Robinson-Fould distances in the presence of recombination because distance is quantified in terms of topological rearrangements due to subtree transfers.

However, discordance may also arise from errors or uncertainty in reconstructing trees from limited sequence data.Espalier therefore allows discordances between trees arising from phylogenetic noise to be reconciled through MAFs. Starting with a MAF, we can regraft subtrees that were cut to obtain the MAF to their respective branches in the original trees and then select the resulting tree that maximizes the likelihood of the sequence data. The reconciliation approach employed by Espalier repeats this process iteratively until all subtrees have been regrafted in a process we call *iterative regrafting*. The resulting reconciled trees may therefore be a topological hybrid of the starting trees. Reconciling trees in this way therefore eliminates conflicts most likely attributable to phylogenetic noise while retaining conflicts that are strongly supported by the sequence data.

Of course, in the presence of true topology-altering recombination events, the ancestral relationships among sampled individuals can no longer be represented by any single tree topology but rather requires a mosaic of phylogenetic histories. Ancestral recombination graphs (ARGs) provide an ideal way to capture this mosaic as a network of linked coalescent and recombination events. An ARG can also be thought of as a sequence of *local* trees representing the phylogenetic history over different regions of the genome with additional recombination nodes linking lineages that have undergone recombination between neighboring local trees.

Espalier efficiently reconstructs approximate ARGs using (surprise!!) MAFs. Starting with an initial candidate local tree for each region of the genome, the iterative-regrafting algorithm is used to first reconcile discordances between neighboring local trees that are not strongly supported by the sequence data. Espalier then uses a dynamic programming approach based on the Viterbi algorithm to select a path of trees along the genome that maximizes the overall likelihood of the sequence data while minimizing discordance. After the tree path is chosen, MAFs are used yet again to identify the recombination events necessary to explain any remaining discordance between neighboring local trees and these events are inserted to obtain a fully connected ARG. While this approach is heuristic, Espalier performs very well at reconstructing ARGs even when compared against more exact but computationally expensive methods [RasmussenGuo2023]_.   

Installation/Usage
******************

The easiest way to install Espalier is through `pip <https://pypi.org/project/pip/>`_:
::

	$ pip install Espalier

pip will install the other Python packages required by Espalier. However, Espalier also relies on `RAxML NG <https://github.com/amkozlov/raxml-ng>`_. See further requirements below.

If Espalier was installed correctly, the command
::

	$ espalier --help

should return the options for Espalier's :ref:`commandline`. 

Requirements
************

Espalier requires Python3 as well as several Python packages to be installed:

* numpy

* pandas

* dendropy

* tskit

* biopython

* scipy

* click

These dependencies will be automatically installed by pip.

Espalier also requires but is not packaged with RAxML-NG. Instructions for installing RAxML-NG can be found `here <https://github.com/amkozlov/raxml-ng>`_.

.. [RasmussenGuo2023] Rasmussen, D.A. and Guo, F. Espalier: Efficient tree reconciliation and ARG reconstruction using maximum agreement forests, Systematic Biology, 2023; syad040.
