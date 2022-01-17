#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 23 12:54:07 2021

@author: david
"""
import dendropy
from Espalier.Reconciler import Reconciler
from Espalier.RAxML import RAxMLRunner
from Espalier import MAF

#from Espalier import Utils
#import ARGSimulator
#import msprime
import os
import datetime



# Set up temp dir for temp tree output
temp_dir = 'temp-{date:%Y-%m-%d_%H:%M:%S}/'.format(date=datetime.datetime.now())
if not os.path.isdir(temp_dir):
    os.mkdir(temp_dir)


# Initialize Espalier objects
raxml = RAxMLRunner(raxml_path='raxml-ng',lsd_path='lsd',temp_dir=temp_dir)
reconciler = Reconciler(raxml,lower_bound_ratio=0.1,prior_gamma=0.0,temp_dir=temp_dir)

# Set up files
tree_file = 'reconciler_true.tre'
seq_file_r1 = 'reconciler_r1.fasta'
seq_file_r2 = 'reconciler_r2.fasta'
ML_tree_file_r1 = 'reconciler_MLTree_r1.tre'
ML_tree_file_r2 = 'reconciler_MLTree_r2.tre'    

"""
    Run simulation
"""
run_sim = False
if run_sim:
    seq_length = 1000 
    ts = msprime.simulate(sample_size=10, Ne=1, length=seq_length, recombination_rate=0, record_full_arg=True)
    tree = ts.first()
    with open(tree_file, "w") as text_file:
        print(tree.newick(), file=text_file)
    
    taxa = dendropy.TaxonNamespace()
    ref = dendropy.Tree.get(file=open(tree_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
    
    # Simulate sequence data for alternate trees
    mut_rate = 0.1
    ARGSimulator.sim_seqs(tree_file,seq_file_r1,mut_rate=mut_rate,seq_length=900)
    ARGSimulator.sim_seqs(tree_file,seq_file_r2,mut_rate=mut_rate,seq_length=100)
    
    # Reconstuct alt tree from seq data"
    tip_date_file = temp_dir + 'dates-lsd.txt'
    Utils.write_tip_dates(tree_file, tip_date_file)
    rate_file = temp_dir + 'rate.txt'
    Utils.write_rate_file(mut_rate, rate_file)
    raxml.get_dated_raxml_tree(seq_file_r1, ML_tree_file_r1, tip_date_file, rate_file)
    raxml.get_dated_raxml_tree(seq_file_r2, ML_tree_file_r2, tip_date_file, rate_file)

"""
    Example code for primer
"""
taxa = dendropy.TaxonNamespace()
tree_r1 = dendropy.Tree.get(file=open(ML_tree_file_r1, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
tree_r2 = dendropy.Tree.get(file=open(ML_tree_file_r2, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)

print("Tree 1:")
tree_r1.print_plot()

print("Tree 2:")
tree_r2.print_plot()

print('SPR distance:', str(MAF.get_spr_dist(tree_r1,tree_r2)))


# Reconcile ML trees through their MAF
maf = MAF.get_maf_4cut(tree_r1,tree_r2)
sampled_trees = reconciler(tree_r1,tree_r2,maf,seq_file_r2)
sampled_trees.sort(key=lambda x: x.like) # sort by ascending likelihood
sampled_trees.reverse() # now sorted in descending order
rec_tree = sampled_trees[0].tree # sampled trees are sorted such that first tree will have highest likelihood
#rec_tree.write(
#    path=reconciled_tree_file,
#    schema='newick',
#    suppress_annotations=True)

print("Reconciled tree:")
rec_tree.print_plot()

print('SPR distance:', str(MAF.get_spr_dist(tree_r1,rec_tree)))