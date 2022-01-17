"""
Created on Fri Jul 31 13:45:35 2020

@author: david
"""
from Espalier.Reconciler import Reconciler
from Espalier.RAxML import RAxMLRunner
from Espalier.ARGBuilder import ARGBuilder
from Espalier.sim import ARGSimulator
from Espalier.sim import ARGweaver
from Espalier import MAF
from Espalier import Utils
import dendropy
from dendropy.calculate import treecompare
import pandas as pd
import glob
import math
import os
import shutil
import argparse
import datetime
    
new_sims = True # simulate new trees/seqs?
argweaver = new_sims # compare with ARGWeaver?
plot = False

"Set defaults if not passing arguments"
n_sims = 1
mut_rate = 0.1
sample_size = 20
prior_gamma = 0.01 # exp decay rate on site likelihood prior
lower_bound_ratio = 0.05 # lower bound used in branch and bound algorithm (was 0.01)
path = "./arg-sampler-test_k20_mu0.01_prevTrees_gammaNatural/"
restart_sim = 1

parser = argparse.ArgumentParser(description='Run simulations to check performance of ARG sampler')
parser.add_argument('--sims', type=int, default=n_sims, help='number of sims to run')
parser.add_argument('--sample_size', type=int, default=sample_size, help='the sample size if a .trees file is not given')
parser.add_argument('--mut_rate', '-mu', type=float, default=mut_rate, help='the mutation rate')
parser.add_argument('--prior_gamma', '-gamma', type=float, default=prior_gamma, help='exp decay rate on site likelihood prior')
parser.add_argument('--lower_bound_ratio', '-lbr', type=float, default=lower_bound_ratio, help='lower bound used in branch and bound algorithm')
parser.add_argument('--out_dir', '-d',type=str,default=path,help='the path to the output directory')
parser.add_argument('--restart_sim', type=int, default=restart_sim, help='sim number to start with if restarting')
args = parser.parse_args()

n_sims = args.sims
path = args.out_dir
if not os.path.isdir(path):
    os.mkdir(path)
results_file = path + "arg-sampler-test_k" + str(args.sample_size) + "_mu" + str(args.mut_rate) + "_bestTrueReference_results.csv"
#trees_file = path + "arg-sampler-test_k" + str(args.sample_size) + "_mu" + str(args.mut_rate) + "_trueTrees_noRecCost_trees.csv"

"Set up temp dir for temp tree output"
temp_dir = 'temp-{date:%Y-%m-%d_%H:%M:%S}/'.format(date=datetime.datetime.now())
if not os.path.isdir(temp_dir):
    os.mkdir(temp_dir)

"""
Control parameters in ARG sampler that need to be tuned:
    -rho recomb prob --> set this to true prob for now
    -prior ratio on alt/ref in reconciler --> now using locally adaptive prior
    -number of local trees to sample
    
Testing regimes:
    -Low,Med,High mutation rates (and recomb rates)
"""

"Sim params"
sample_size = args.sample_size
Ne =  1.0 / 2 # divide by two because msprime assumes individuals are diploid
genome_length = 1e3
recomb_rate = 1.0 # rate per genome
rho = recomb_rate / genome_length # rate per site
mut_rate = args.mut_rate

"prior"
prior_gamma = args.prior_gamma
lower_bound_ratio = args.lower_bound_ratio

"Or set natural gamma decay rate based on prob of no recombination"
expected_tree_length = 0.0
for k in range(2,sample_size+1):
    expected_tree_length += k * 2*Ne / (k*(k-1)/2)
prob_recomb_per_site = expected_tree_length * rho
prob_mut_per_site = expected_tree_length * mut_rate
prior_gamma = prob_recomb_per_site

"Initialize callable instances of Espalier objects"
raxml = RAxMLRunner(raxml_path='raxml-ng',lsd_path='lsd',temp_dir=temp_dir)

reconciler = Reconciler(raxml, lower_bound_ratio=lower_bound_ratio,prior_gamma=prior_gamma,temp_dir=temp_dir)

argb = ARGBuilder(reconciler,raxml,max_sample_size=10)

"Just to check exp decay approximation"
#import numpy as np
#prob_no_recomb_per_site = 1 - prob_recomb_per_site
#mult_prob = prob_no_recomb_per_site**1
#exp_prob = np.exp(-prob_recomb_per_site*1)

"Performance metrics"
start_num = args.restart_sim
if start_num == 1:
    sims = []
    intervals = []
    methods = []
    RF_dists = []
    weighted_RF_dists = []
    SPR_dists = []
    true_SPR_dists = []
    best_trees = []
    selected_trees = []
    #sim_trees = []
    prop_zero_branches = []
    informative_sites = []
else:
    results_df = pd.read_csv(results_file)
    sims = results_df['Simulation'].tolist()
    intervals = results_df['LocalTree'].tolist()
    methods = results_df['Method'].tolist()
    RF_dists = results_df['RF Dist'].tolist()
    weighted_RF_dists = results_df['Weighted RF Dist'].tolist()
    SPR_dists = results_df['SPR Dist'].tolist()
    true_SPR_dists = results_df['True SPR Dist'].tolist()
    
    #trees_df = pd.read_csv(trees_file)
    best_trees = results_df['Best Tree'].tolist()
    selected_trees = results_df['Selected Tree'].tolist()
    #sim_trees = trees_df['Simulation'].tolist()
    prop_zero_branches = results_df['Proportion Zero Branches'].tolist()
    informative_sites = results_df['Informative Sites'].tolist()

"Simulate ARG in msprime"    
min_seg_length = 20
for sim in range(start_num,n_sims+1): 
    
    print("-"*20)
    print("Simulation = " + str(sim))
    print("-"*20)
    
    sim_path = path + 'sim' + str(sim).rjust(3, '0') + '/'
    
    sim_flag = True
    while sim_flag:
        
        if new_sims:
        
            "Create sim dir or delete old one if it exists"
            if os.path.isdir(sim_path):
                shutil.rmtree(sim_path)
                os.mkdir(sim_path)
            else:
                os.mkdir(sim_path)
            
            sim_flag = False
            ts = ARGSimulator.sim_ARG(sample_size=sample_size,Ne=Ne,length=genome_length,recombination_rate=rho,min_breakpoints=2)
            breaks = ts.breakpoints(as_array=True)
            segments = len(breaks) - 1 # number of non-recombinant segments between breakpoints
            for tr_num, tree in enumerate(ts.trees()):
                seq_length = round(tree.interval[1] - tree.interval[0])
                if seq_length < min_seg_length:
                    sim_flag = True
            if sim_flag:
                continue
            
            print("Non-recombinant segments = " + str(segments))
        
            "Write local tree and seq files"
            tree_files = [sim_path + "sim" + str(sim) + "_tree" + str(i) + ".tre" for i in range(segments)]
            seq_files = [sim_path + "sim" + str(sim) + "_tree" + str(i) + ".fasta" for i in range(segments)]
            #seg_lengths = []
            informative_site_counts = []
            for tr_num, tree in enumerate(ts.trees()):
                seq_length = round(tree.interval[1] - tree.interval[0])
                #seg_lengths.append(seq_length)
                with open(tree_files[tr_num], "w") as text_file:
                    print(tree.newick(), file=text_file)
                ARGSimulator.sim_seqs(tree_files[tr_num],seq_files[tr_num],mut_rate=mut_rate,seq_length=seq_length)
                informative_site_counts.append(Utils.informative_sites(seq_files[tr_num]))
            
            "Reconstruct local ML trees in RaXML"
            reconstruct = True
            ML_tree_files = [sim_path + "sim" + str(sim) + "_MLTree" + str(i) + ".tre" for i in range(segments)]
            tip_date_file = temp_dir + 'dates-lsd.txt'
            Utils.write_tip_dates(tree_files[0], tip_date_file)
            rate_file = temp_dir + 'rate.txt'
            Utils.write_rate_file(mut_rate, rate_file)
            if reconstruct:
                for i in range(segments):
                    if informative_site_counts[i] > 0:
                        #ARGSampler.get_raxml_tree(seq_files[i], ML_tree_files[i])
                        #RAxML.get_dated_raxml_tree_ng(seq_files[i], ML_tree_files[i], tip_date_file, rate_file, temp_dir)
                        raxml.get_dated_raxml_tree(seq_files[i], ML_tree_files[i], tip_date_file, rate_file)
        
        else:
            sim_flag = False
            fasta_file_list = glob.glob(sim_path + '*.fasta') # get all smc files in dir as list
            informative_site_counts = []
            for seq_file in fasta_file_list:
                informative_site_counts.append(Utils.informative_sites(seq_file))
            segments = len(fasta_file_list) - 1
            tree_files = [sim_path + "sim" + str(sim) + "_tree" + str(i) + ".tre" for i in range(segments)]
            seq_files = [sim_path + "sim" + str(sim) + "_tree" + str(i) + ".fasta" for i in range(segments)]
            ML_tree_files = [sim_path + "sim" + str(sim) + "_MLTree" + str(i) + ".tre" for i in range(segments)]
    
        "Get reference tree from consensus of reconstructed trees - set non-reconstructed trees to consensus"
        reconstructed_trees = [ML_tree_files[x] for x in range(segments) if informative_site_counts[x] > 0]
        "If reconstructred trees is empty sim flag = True"
        if not reconstructed_trees:
            sim_flag = True # repeat sim -- no trees were reconstructed
        nonreconstructed_trees = [ML_tree_files[x] for x in range(segments) if informative_site_counts[x] == 0]
        ref = Utils.get_consensus_tree(reconstructed_trees)
        consensus_tree_file = sim_path + "sim" + str(sim) + "_consensus.tre"
        ref.write(path=consensus_tree_file,schema='newick',suppress_annotations=True,suppress_rooting=True) 
        for tree_file in nonreconstructed_trees:
            ref.write(path=tree_file,schema='newick',suppress_annotations=True,suppress_rooting=True) 
    
        "Get priors based on segment/total length ratio"
        #prior_tree_ratio = []
        #max_seg_length = max(seg_lengths)
        #for i in range(segments):
        #    prior_tree_ratio.append(seg_lengths[i] / max_seg_length)
        

        """
            Use "best" true tree as ref
        """
        import numpy as np
        taxa = ref.taxon_namespace
        true_tree_dists = np.zeros((segments,segments))
        for seg_i in range(segments):   
            tree_i = dendropy.Tree.get(file=open(tree_files[seg_i], 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
            for seg_j in range(segments):
                if seg_i != seg_j:
                    tree_j = dendropy.Tree.get(file=open(tree_files[seg_j], 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)    
                    true_tree_dists[seg_i][seg_j] = treecompare.symmetric_difference(tree_i, tree_j)
        best_tree = np.argmin(np.sum(true_tree_dists,axis=0))
        ref = dendropy.Tree.get(file=open(tree_files[best_tree], 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
        ref.suppress_unifurcations(update_bipartitions=False)
        
        "Run the ARG sampler to get a tree path"
        out_report = sim_path + "arg-sampler-sim" + str(sim) + "-report.txt"
        tree_path = argb.reconstruct_ARG(ML_tree_files,seq_files,ref,rho,temp_dir=temp_dir)
        tree_source_path = argb.tree_source_path
    
        "Compare tree path with ARGweaver"
        ARGweaver_tree_files = [sim_path + "sim" + str(sim) + "_ARGweaverTree" + str(i) + ".tre" for i in range(segments)]
        if argweaver:
            concat_fasta = sim_path + "sim" + str(sim) + "_concatenated.fasta"
            ARGweaver.concat_alignments(seq_files,concat_fasta)
            smc_path = sim_path #'./argweaver/sim' + str(sim) + '/'
            smc_files = smc_path + 'arg-out'
            sim_flag = ARGweaver.run_arg_sampler(concat_fasta,smc_files,2*Ne,rho,mut_rate,verbose=False)
            if sim_flag:
                continue
            smc_list = glob.glob(smc_path + '*.smc') # get all smc files in dir as list
            tree_df,taxa_map = ARGweaver.parse_smc(smc_list)
            for seg in range(segments):
                position = math.floor(breaks[seg]) + 1 # a position within this segment
                consensus = ARGweaver.local_consensus_tree(tree_df,position)
                for lf in consensus.leaf_node_iter():
                    lf.taxon.label = taxa_map[lf.taxon.label] # Reset taxa names back to original labels
                consensus.write(path=ARGweaver_tree_files[seg],schema='newick',suppress_annotations=True)
                
    "Test performance"
    for seg, tree in enumerate(tree_path):
        

        """
            For ARG sampler method
        """
        sims.append(sim)
        intervals.append(seg)
        methods.append('ARG Sampler')
        
        "Get RF (NNI) dist to true tree"
        true_tree = dendropy.Tree.get(file=open(tree_files[seg], 'r'), schema="newick", rooting="default-rooted", taxon_namespace=tree.taxon_namespace)
        RF_dist = treecompare.symmetric_difference(tree, true_tree)
        RF_dists.append(RF_dist)
        
        "Get weighted RF dist to true tree"
        wRF_dist = treecompare.weighted_robinson_foulds_distance(tree, true_tree)
        weighted_RF_dists.append(wRF_dist)
        
        "Get SPR dist from tree to previous tree"
        if seg > 0:
            SPR_dist = MAF.get_spr_dist(tree,tree_path[seg-1])
        else:
            SPR_dist = 0
        SPR_dists.append(SPR_dist)
        
        "Get SPR between true trees"
        if seg > 0:
            true_prev_tree = dendropy.Tree.get(file=open(tree_files[seg-1], 'r'), schema="newick", rooting="default-rooted", taxon_namespace=tree.taxon_namespace)
            true_SPR_dist = MAF.get_spr_dist(true_prev_tree,true_tree)
        else:
            true_SPR_dist = 0
        true_SPR_dists.append(true_SPR_dist) 
        
        """
            For ML trees
        """
        sims.append(sim)
        intervals.append(seg)
        methods.append('ML Tree')
        
        "Get RF (NNI) dist from ML tree to true tree"
        ML_tree = dendropy.Tree.get(file=open(ML_tree_files[seg], 'r'), schema="newick", rooting="default-rooted", taxon_namespace=tree.taxon_namespace)
        ML_dist = treecompare.symmetric_difference(ML_tree, true_tree)
        RF_dists.append(ML_dist)
        
        "Get weighted RF dist from ML tree to true tree"
        wRF_dist = treecompare.weighted_robinson_foulds_distance(ML_tree, true_tree)
        weighted_RF_dists.append(wRF_dist)
        
        "Get SPR dist from ML tree to previous ML tree"
        prev_ML_tree = dendropy.Tree.get(file=open(ML_tree_files[seg-1], 'r'), schema="newick", rooting="default-rooted", taxon_namespace=tree.taxon_namespace)
        if seg > 0:
            ML_SPR_dist = MAF.get_spr_dist(ML_tree,prev_ML_tree)
        else:
            ML_SPR_dist = 0
        SPR_dists.append(ML_SPR_dist)
        true_SPR_dists.append(true_SPR_dist) 
        
        """
            For ARGWeaver trees
        """
        sims.append(sim)
        intervals.append(seg)
        methods.append('ARGWeaver')
        
        "Get RF (NNI) dist from local ARGWeaver tree to true tree"
        ARGweaver_tree = dendropy.Tree.get(file=open(ARGweaver_tree_files[seg], 'r'), schema="newick", rooting="default-rooted", taxon_namespace=tree.taxon_namespace)
        ARGweaver_dist = treecompare.symmetric_difference(ARGweaver_tree, true_tree)
        RF_dists.append(ARGweaver_dist)
        
        "Get weighted RF dist from local ARGWeaver tree to true tree"
        wRF_dist = treecompare.weighted_robinson_foulds_distance(ARGweaver_tree, true_tree)
        weighted_RF_dists.append(wRF_dist)
        
        "Get SPR dist from ARGweaver tree to previous ARGweaver tree"
        prev_ARGweaver_tree = dendropy.Tree.get(file=open(ARGweaver_tree_files[seg-1], 'r'), schema="newick", rooting="default-rooted", taxon_namespace=tree.taxon_namespace)
        if seg > 0:
            ARGweaver_SPR_dist = MAF.get_spr_dist(ARGweaver_tree,prev_ARGweaver_tree)
        else:
            ARGweaver_SPR_dist = 0
        SPR_dists.append(ARGweaver_SPR_dist) 
        true_SPR_dists.append(true_SPR_dist)
        
        """
            For consensus tree
        """
        sims.append(sim)
        intervals.append(seg)
        methods.append('Consensus')
        
        "Get RF (NNI) dist from consensus tree to true local tree"
        consensus_tree = dendropy.Tree.get(file=open(consensus_tree_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=tree.taxon_namespace)
        consensus_dist = treecompare.symmetric_difference(consensus_tree, true_tree)
        RF_dists.append(consensus_dist)
        
        "Get weighted RF dist from consensus tree to true tree"
        wRF_dist = treecompare.weighted_robinson_foulds_distance(consensus_tree, true_tree)
        weighted_RF_dists.append(wRF_dist)
        
        SPR_dists.append(0) # dist to prev tree is zero because consensus tree is constant 
        true_SPR_dists.append(true_SPR_dist)
        
        """
            Record best vs. picked trees
        """
        tree_options = ['Selected','ML','Consensus']
        seg_RF_dists = [RF_dist,ML_dist,consensus_dist]
        best_index = seg_RF_dists.index(min(seg_RF_dists))
        best_tree = tree_options[best_index]
        selected_tree = tree_source_path[seg]
        #sim_trees.append(sim)
        best_trees.append(best_tree)
        best_trees.append(best_tree)
        best_trees.append(best_tree)
        best_trees.append(best_tree)
        selected_trees.append(selected_tree)
        selected_trees.append(selected_tree)
        selected_trees.append(selected_tree)
        selected_trees.append(selected_tree)
        
        "Compute proprotion of zero-length branches"
        total_branches = len(ML_tree.edges())
        zero_branches = 0
        for edge in ML_tree.postorder_edge_iter():
            if not edge.length is None:
                if edge.length == 0:
                    zero_branches += 1
        prop_zero = round(zero_branches / total_branches, 3)
        prop_zero_branches.append(prop_zero)
        prop_zero_branches.append(prop_zero)
        prop_zero_branches.append(prop_zero)
        prop_zero_branches.append(prop_zero)
        
        informative_sites.append(informative_site_counts[seg])
        informative_sites.append(informative_site_counts[seg])
        informative_sites.append(informative_site_counts[seg])
        informative_sites.append(informative_site_counts[seg])

    data = {'Simulation':sims,
            'LocalTree':intervals,
            'Method':methods, 
            'RF Dist': RF_dists,
            'Weighted RF Dist': weighted_RF_dists,
            'SPR Dist': SPR_dists,
            'True SPR Dist': true_SPR_dists,
            'Best Tree': best_trees,
            'Selected Tree': selected_trees,
            'Proportion Zero Branches': prop_zero_branches,
            'Informative Sites': informative_sites}
    df = pd.DataFrame(data)
    df.to_csv(results_file,index=False)
    
    #tree_data = {'Simulation':sim_trees,
    #        'Best Tree': best_trees,
    #        'Selected Tree': selected_trees,
    #        'Proportion Zero Branches': prop_zero_branches} 
    
    #tree_df = pd.DataFrame(tree_data)
    #tree_df.to_csv(trees_file,index=False)

