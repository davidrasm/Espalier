"""
Created on Wed May 20 09:36:28 2020

Test Reconciler to see if we can correctly reconcile discordant trees

@author: david
"""

from Espalier.sim import ARGSimulator
from Espalier.sim import TreeRandomizer
from Espalier import MAF
from Espalier import Utils
from Espalier.RAxML import RAxMLRunner
from Espalier.Reconciler import Reconciler
from Espalier import Utils
import dendropy
from dendropy.calculate import treecompare
import msprime
import numpy as np
import pandas as pd
import os
import datetime


"Simulate tree without recombination"
path = './reconciler-randomSPRTest-16xPriorGamma-100reps/'
ref_tree_file = path + 'randomSPRs_refTree.tre'
ref_seq_file = path + 'randomSPRs_refSeqs.tre'
alt_tree_file = path + 'randomSPRs_altTree.tre'
ML_tree_file = path + 'randomSPRs_MLTree.tre'
alt_seq_file = path + 'randomSPRs_altSeqs.fasta'
concat_seq_file = path + 'randomSPRs_concatSeqs.fasta'
reconciled_tree_file = path + 'randomSPRs_reconciled.tre'

"Sim params"
breakpoints = np.arange(100,1000,100)
n_sims = len(breakpoints)
reps = 100
sample_size = 20
Ne= 1/2
mut_rate = 0.1
seq_length=1000
plot = False

if not os.path.isdir(path):
    os.mkdir(path)
results_file = path + "reconciler-test_results.csv"

"Set up temp dir for temp tree output"
temp_dir = 'temp-{date:%Y-%m-%d_%H:%M:%S}/'.format(date=datetime.datetime.now())
if not os.path.isdir(temp_dir):
    os.mkdir(temp_dir)

"Performance metrics"
results_cols = ['Simulation','Breakpoint Location','Method', 'RF Distance','SPR Distance']
start_num = 1
if start_num == 1:
    results_df = pd.DataFrame(columns=results_cols)
else:
    results_df = pd.read_csv(results_file)

"Compute appropriate prior_gamma"
expected_tree_length = 0.0
for k in range(2,sample_size+1):
    expected_tree_length += k * 2*Ne / (k*(k-1)/2)
prob_recomb_per_site = expected_tree_length * 0.001
prior_gamma = 16*prob_recomb_per_site

"Initialize callable instances of Espalier objects"
raxml = RAxMLRunner(raxml_path='raxml-ng',lsd_path='lsd',temp_dir=temp_dir)
reconciler = Reconciler(raxml,lower_bound_ratio=0.1,prior_gamma=prior_gamma,temp_dir=temp_dir)

"Run sims"
#for sim in range(0,n_sims):
sim = 0
for bk in breakpoints:
    for rp in range(reps):
    
        sim += 1
        print("-"*20)
        print("Simulation = " + str(sim))
        print("-"*20)
        
        "Simulate ref tree"
        ts = msprime.simulate(sample_size=sample_size, Ne=Ne, length=seq_length, recombination_rate=0, record_full_arg=True)
        tree = ts.first()
        with open(ref_tree_file, "w") as text_file:
            print(tree.newick(), file=text_file)
        taxa = dendropy.TaxonNamespace()
        ref = dendropy.Tree.get(file=open(ref_tree_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
        
        """
        Perform random SPR moves on tree
        """
        moves = 1 # number of SPR moves
        alt = ref
        for i in range(moves):
            alt = TreeRandomizer.random_SPR(alt,exclude_root_children=True,plot=False)
        alt.taxon_namespace = taxa
        alt.reconstruct_taxon_namespace() #reindex_subcomponent_taxa()
        alt.write(path=alt_tree_file,schema='newick',suppress_annotations=True,suppress_rooting=True)
    
        "Need to actually simulate sequence data for alternate tree"
        ref_length = bk
        alt_legnth = seq_length - ref_length
        ARGSimulator.sim_seqs(ref_tree_file,ref_seq_file,mut_rate=mut_rate,seq_length=ref_length)
        ARGSimulator.sim_seqs(alt_tree_file,alt_seq_file,mut_rate=mut_rate,seq_length=alt_legnth)
        
        Utils.concate_aligns([ref_seq_file,alt_seq_file],concat_seq_file)
        
        "Reconstuct alt tree from seq data"
        tip_date_file = temp_dir + 'dates-lsd.txt'
        Utils.write_tip_dates(alt_tree_file, tip_date_file)
        rate_file = temp_dir + 'rate.txt'
        Utils.write_rate_file(mut_rate, rate_file)
        raxml.get_dated_raxml_tree(alt_seq_file, ML_tree_file, tip_date_file, rate_file)
        ML_tree = dendropy.Tree.get(file=open(ML_tree_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
        
        if plot:
            print("Reference tree:")
            ref.print_plot()
        
        if plot:
            print("Alternate tree:")
            alt.print_plot()
        
        "Get maximum agreement forest"
        maf = MAF.get_maf_4cut(ref,ML_tree,plot=False)
        if plot:
            print("MAF:")
            MAF.plot_maf(maf)
        
        "Reconcile ML tree for alt and ref"
        sampled_trees = reconciler(ref,ML_tree,maf,concat_seq_file,start_pos=bk,end_pos=seq_length)
        sampled_trees.sort(key=lambda x: x.like) # sort by ascending likelihood
        sampled_trees.reverse() # now sorted in descending order
        rec_tree = sampled_trees[0].tree
        rec_tree.taxon_namespace = taxa
        rec_tree.reconstruct_taxon_namespace() #reindex_subcomponent_taxa()
        #rec_tree.write(
        #    path=reconciled_tree_file,
        #    schema='newick',
        #    suppress_annotations=True)
        
        if plot:
            print("Reconciled tree:")
            rec_tree.print_plot()
            
        rec_RF_dist = treecompare.symmetric_difference(rec_tree,alt) # RF distance between reconciled tree and alt tree
        ML_RF_dist = treecompare.symmetric_difference(ML_tree,alt) # RF distance between reconciled tree and alt tree
        rec_SPR_dist = MAF.get_spr_dist(rec_tree,ref) # SPR distance between reconciled tree and reference
        ML_SPR_dist = MAF.get_spr_dist(ML_tree,ref) # SPR distance betweeen ML tree and reference
       
        "Add sim results to df for ML tree"
        results = {'Simulation':sim,
                   'Breakpoint Location':bk,
                   'Method':'ML',
                   'RF Distance':ML_RF_dist,
                   'SPR Distance':ML_SPR_dist}
        results_df = results_df.append(results, ignore_index=True)
        
        results = {'Simulation':sim,
                   'Breakpoint Location':bk,
                   'Method':'Reconciled',
                   'RF Distance':rec_RF_dist,
                   'SPR Distance':rec_SPR_dist}
        results_df = results_df.append(results, ignore_index=True)
                
    results_df.to_csv(results_file,index=False)

