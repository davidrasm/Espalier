"""
Created on Wed Dec 22 15:23:09 2021

Example of how to use ARGBuilder to reconstruct a simulated ARG

@author: david
"""

from Espalier.Reconciler import Reconciler
from Espalier.RAxML import RAxMLRunner
from Espalier.ARGBuilder import ARGBuilder
from Espalier.viz import PlotTanglegrams

# Set up temp dir for temp tree output
import os
import datetime
import shutil
import numpy as np
import dendropy
from dendropy.calculate import treecompare

temp_dir = 'temp-{date:%Y-%m-%d_%H:%M:%S}/'.format(date=datetime.datetime.now())
if not os.path.isdir(temp_dir):
    os.mkdir(temp_dir)

sim_path = './sim1/'    
if os.path.isdir(sim_path):
    shutil.rmtree(sim_path)
    os.mkdir(sim_path)
else:
    os.mkdir(sim_path)

"""
    Run simulation
"""

# import ARGSimulator
# import Utils

# sim = 1
# sample_size = 10
# Ne = 1.0 / 2
# genome_length = 1000
# rec_rate = 0.2 / genome_length # rate per site
# mut_rate = 0.1

# # Sim ARG
# ts = ARGSimulator.sim_ARG(sample_size=sample_size,Ne=Ne,length=genome_length,recombination_rate=rec_rate,min_breakpoints=2)
# breaks = ts.breakpoints(as_array=True)
# segments = len(breaks) - 1 # number of non-recombinant segments between breakpoints
# print("Segments: " + str(segments))

# # Write local tree and sim seq files
# tree_files = [sim_path + "example_tree" + str(i) + ".tre" for i in range(segments)]
# seq_files = [sim_path + "example_tree" + str(i) + ".fasta" for i in range(segments)]
# for tr_num, tree in enumerate(ts.trees()):
#     seq_length = round(tree.interval[1] - tree.interval[0])
#     print("Segment " + str(tr_num) + " length: " + str(seq_length))
#     with open(tree_files[tr_num], "w") as text_file:
#         print(tree.newick(), file=text_file)
#     ARGSimulator.sim_seqs(tree_files[tr_num],seq_files[tr_num],mut_rate=mut_rate,seq_length=seq_length)

# # Reconstruct local ML trees in RAxML
# raxml = RAxMLRunner(raxml_path='raxml-ng',lsd_path='lsd',temp_dir=temp_dir)
# ML_tree_files = [sim_path + "example_MLTree" + str(i) + ".tre" for i in range(segments)]
# tip_date_file = temp_dir + 'dates-lsd.txt'
# Utils.write_tip_dates(tree_files[0], tip_date_file)
# rate_file = temp_dir + 'rate.txt'
# Utils.write_rate_file(mut_rate, rate_file)
# for i in range(segments):
#     raxml.get_dated_raxml_tree(seq_files[i], ML_tree_files[i], tip_date_file, rate_file)


# # Use "best" true tree as ref
# ref = Utils.get_consensus_tree(ML_tree_files)
# taxa = ref.taxon_namespace
# true_tree_dists = np.zeros((segments,segments))
# for seg_i in range(segments):   
#     tree_i = dendropy.Tree.get(file=open(tree_files[seg_i], 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
#     for seg_j in range(segments):
#         if seg_i != seg_j:
#             tree_j = dendropy.Tree.get(file=open(tree_files[seg_j], 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)    
#             true_tree_dists[seg_i][seg_j] = treecompare.symmetric_difference(tree_i, tree_j)
# best_tree = np.argmin(np.sum(true_tree_dists,axis=0))


"""
    Example code for primer - some duplicated from above for clarity
"""

# Plot tanglegram for local trees in true ARG
segments = 4 # i.e. number of genomic regions
tree_files = ["ARG_example_tree" + str(i) + ".tre" for i in range(segments)]
tree0 = dendropy.Tree.get(file=open(tree_files[0], 'r'), schema="newick", rooting="default-rooted")
taxa = tree0.taxon_namespace # get taxon namespace from first tree
#tree_files = [sim_path + "example_tree" + str(i) + ".tre" for i in range(segments)]
tanglegram_fig_name = 'trueARG-tanglegram.png' 
PlotTanglegrams.plot(tree_files, tanglegram_fig_name, numerical_taxa_names=True)

# Get ML tree files and corresponding seq files
ML_tree_files = ["ARG_example_MLTree" + str(i) + ".tre" for i in range(segments)]
#ML_tree_files = [sim_path + "example_MLTree" + str(i) + ".tre" for i in range(segments)]
seq_files = ["ARG_example_tree" + str(i) + ".fasta" for i in range(segments)]
#seq_files = [sim_path + "example_tree" + str(i) + ".fasta" for i in range(segments)]

# Plot tanglegram for local ML trees
tanglegram_fig_name = 'localMLTree-tanglegram.png' 
PlotTanglegrams.plot(ML_tree_files, tanglegram_fig_name, numerical_taxa_names=True)

#Initialize callable instances of Espalier objects
lower_bound_ratio = 0.1 # lower bound ratio for reconciliation algorithm
rec_rate = 0.0002 # recombination rate per site
prior_gamma = rec_rate # decay rate for weighting site likelihoods at sites external to a genomic region
raxml = RAxMLRunner(raxml_path='raxml-ng',lsd_path='lsd',temp_dir=temp_dir)
reconciler = Reconciler(raxml, lower_bound_ratio=lower_bound_ratio,prior_gamma=prior_gamma,temp_dir=temp_dir)
argb = ARGBuilder(reconciler,raxml)

# Get a reference tree
#ref = Utils.get_consensus_tree(ML_tree_files)
#ref = dendropy.Tree.get(file=open(tree_files[best_tree], 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
ref = dendropy.Tree.get(file=open('ARG_example_consensus_ref.tre', 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
#ref.suppress_unifurcations(update_bipartitions=False)
#ref.write(path='ARG_example_consensus_ref.tre',schema='newick',suppress_annotations=True,suppress_rooting=True)

# Run reconstruction
tree_path = argb.reconstruct_ARG(ML_tree_files,seq_files,ref,rec_rate,temp_dir=temp_dir)

# Write local trees in reconstructed ARG to files
ARG_tree_files =  ["ARG_example_ARGLocalTree" + str(i) + ".tre" for i in range(segments)]
#ARG_tree_files =  [sim_path + "example_ARGLocalTree" + str(i) + ".tre" for i in range(segments)]
for idx,tr in enumerate(tree_path):
    tr.write(path=ARG_tree_files[idx],schema='newick',suppress_annotations=True,suppress_rooting=True)
tanglegram_fig_name = 'ARGLocalTree-tanglegram.png' 
PlotTanglegrams.plot(ARG_tree_files, tanglegram_fig_name, numerical_taxa_names=True)

