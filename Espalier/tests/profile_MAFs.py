"""
Created on Tue May 12 10:39:44 2020

Test maximum agreement forest algorithm with random trees created by SPR moves in reference

@author: david
"""

from Espalier.sim import TreeRandomizer
from Espalier import MAF
import dendropy
import msprime
import random
import time
import numpy as np

def sim_trees(sims,plot=False):
    
    sample_size = 100
    max_cuts = 5
    spr_moves = []
    ref_trees = []
    alt_trees = []
    for i in range(sims):
        
        print("Sim: " + str(i))
        
        ts = msprime.simulate(sample_size=sample_size, Ne=100, length=1e4, recombination_rate=0, record_full_arg=True)
        tree = ts.first()
        with open("maf_randomSPRs_refTree_temp.tre", "w") as text_file:
            print(tree.newick(), file=text_file)
        ref_file = "maf_randomSPRs_refTree_temp.tre"
    
        taxa = dendropy.TaxonNamespace()
        ref = dendropy.Tree.get(file=open(ref_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
    
        #Perform random SPR moves on tree#
        moves = random.randint(1,max_cuts) # number of SPR moves
        alt = ref
        for i in range(moves):
            #print("SPR move = ", str(i))
            alt = TreeRandomizer.random_SPR(alt,exclude_root_children=True,plot=False)
        
        spr_moves.append(moves)
        ref_trees.append(ref)
        alt_trees.append(alt)
    
    return ref_trees, alt_trees

# Simulate pairs of discordant trees
sims = 100
ref_trees, alt_trees = sim_trees(sims)
  
# Run maximum agreement forest algorithm"
run_times = []
for ref,alt in zip(ref_trees,alt_trees):
    tic = time.perf_counter()
    maf = MAF.get_maf_4cut(ref,alt,plot=False)
    toc = time.perf_counter()
    elapsed = toc - tic
    run_times.append(elapsed)
    print(f"Elapsed time: {elapsed:0.4f} seconds")
    
avg_time = np.mean(run_times)
print(f"Average elapsed time: {avg_time:0.4f} seconds")
    
    
    