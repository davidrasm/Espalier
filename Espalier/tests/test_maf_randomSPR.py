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
import pandas as pd

def batch_test(plot=False):
    
    sample_size = 100
    sims = 100
    max_cuts = 1
    sample_sizes = [sample_size]*sims
    spr_moves = []
    spr_dists_4cut = []
    spr_dists_3approx = []
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
            alt = TreeRandomizer.random_SPR(alt,plot=False)
    
        if plot:    
            print("Reference tree:")
            ref.print_plot()
            print()
            print("Alternate tree:")
            alt.print_plot()
            #alt.write(path='maf_randomSPRs_altTree_test1.tre',schema='newick',suppress_annotations=True)
    
        "Get maximum agreement forest"
        maf = MAF.get_maf_4cut(ref,alt,plot=False)
        #maf = MAF.get_maf_4cut(ref,alt,plot=False)
        
        "Get 3-approx maximum agreement forest"
        maf3 = MAF.get_maf_3approx(ref,alt,plot=False)
        
        if plot:
            MAF.plot_maf(maf)
        
        spr_moves.append(moves)
        spr_dists_4cut.append(len(maf)-1)
        spr_dists_3approx.append(len(maf3)-1)
    
    "Put results for both methods in one row"
    #data = {'TrueSPRMoves':spr_moves,'SPRDists4Cut':spr_dists_4cut,'SPRDists3Approx': spr_dists_3approx}   
    #df = pd.DataFrame(data)
    
    "Split methods into different results"
    data_4cut = {'TrueSPRMoves':spr_moves,'SPRDist':spr_dists_4cut,'SampleSizes':sample_sizes,'Method':['4Cut']*sims}
    df = pd.DataFrame(data_4cut)
    
    data_3approx = {'TrueSPRMoves':spr_moves,'SPRDist':spr_dists_3approx,'SampleSizes':sample_sizes,'Method':['3Approx']*sims}
    df_3approx = pd.DataFrame(data_3approx)
    df = df.append(df_3approx)
    
    df.to_csv("test_MAF_randomSPRS_4cut_vs_3approx_maxSPR1_s100_results.csv",index=False)

"Batch test"
batch = False
if batch:
    batch_test()
else:

    "Tree files"
    ref_file = "maf_example_tree1.tre"
    alt_file = "maf_example_tree2.tre"

    "Simulate tree without recombination"
    sim = True
    if sim:
        ts = msprime.simulate(sample_size=10, Ne=100, length=1e4, recombination_rate=0, record_full_arg=True)
        tree = ts.first()
        print("-" * 20)
        print(tree.draw(format="unicode"))
        with open(ref_file, "w") as text_file:
            print(tree.newick(), file=text_file)
    
    taxa = dendropy.TaxonNamespace()
    ref = dendropy.Tree.get(file=open(ref_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
    
    
    """
    Perform random SPR moves on tree
    """
    moves = 2 # number of SPR moves
    alt = ref
    for i in range(moves):
        print("SPR move = ", str(i))
        alt = TreeRandomizer.random_SPR(alt,plot=False)
    
    print("Reference tree:")
    ref.print_plot()
    
    print("Alternate tree:")
    alt.print_plot()
    alt.write(path=alt_file,schema='newick',suppress_annotations=True)
    
    "Get maximum agreement forest"
    maf = MAF.get_maf_4cut(ref,alt,plot=False)
    MAF.plot_maf(maf)
    
    
    