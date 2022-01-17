"""
Created on Tue Jan 19 10:04:21 2021

    Test if ARGs simulated by msprime meet expected theoretical properties in terms of heights and numbers of breakpoints

@author: david
"""
import msprime
import tskit
import numpy as np
#import matplotlib.pyplot as plt
#import seaborn as sns

def basic_sim(Ne=100,plot=False):

    "Simulate diploid tree"
    tree_sequence = msprime.simulate(sample_size=2, Ne=Ne)
    tree = tree_sequence.first()
    if plot: print(tree.draw(format="unicode"))
    
    "Travere tree back to root"
    u = 2
    while u != tskit.NULL:
        if plot: print("node {}: time = {}".format(u, tree.time(u)))
        if plot: print("node {}: branch length = {}".format(u, tree.branch_length(u)))
        height = tree.time(u)
        u = tree.parent(u)
        
    return height

def test_tree_heights():
    
    """
        msprime assumes diploid pop sizes so for haploid pops need to set Ne = Ne / 2
    """
    sims = 1000
    Ne = 1.0 / 2 # divide by two because msprime assumes individuals are diploid
    tree_heights = np.zeros([sims,1])
    for s in range(sims):
        tree_heights[s] = basic_sim(Ne=Ne)
    print(np.mean(tree_heights))

def sim_ARG(sample_size=10,Ne=100,length=1e3,recombination_rate=5e-6):
    
    "Simulate local trees in ARG using msprime"
    ts = msprime.simulate(sample_size=sample_size, Ne=Ne, length=length, recombination_rate=recombination_rate, record_full_arg=True)
    breaks = len(ts.breakpoints(as_array=True)) - 2 # -2 because tskit counts ends as breakpoints
    
    return breaks

def expected_tree_length(sample_size,Ne):
    
    """
        sample_size = sample size
        Ne = haploid pop size
    """
    
    "Expected total tree length"
    tree_length = 0.0
    for k in range(2,sample_size+1):
        tree_length += k * Ne / (k*(k-1)/2)
    return tree_length
    

def test_breakpoint_count():
    
    sample_size = 10
    Ne = 1.0 / 2 # divide by two because msprime assumes individuals are diploid
    rho = 1.0
    genome_length = 1000
    rho_per_site = rho / genome_length
    
    T_length = expected_tree_length(sample_size,2*Ne)
    
    expected_breaks = T_length * rho
    
    mu = 0.05 * genome_length
    expected_muts = T_length * mu
    expected_muts_per_block = expected_muts / (expected_breaks+1)
    print()
    
    
    sims = 1000
    breaks = np.zeros([sims,1])
    for s in range(sims):
        breaks[s] = sim_ARG(sample_size=sample_size,Ne=Ne,length=genome_length,recombination_rate=rho_per_site)
    
    print ('Expected breakpoints: %s' % expected_breaks)
    print ('Mean sim breakpoints: %s' % np.mean(breaks))


#test_tree_heights()

test_breakpoint_count()
