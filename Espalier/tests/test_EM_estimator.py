"""
Created on Fri Jul 31 13:45:35 2020

Remeber: Tree path is set to true trees for testing right now

@author: david
"""
from Espalier.Reconciler import Reconciler
from Espalier.RAxML import RAxMLRunner
from Espalier.ARGBuilder import ARGBuilder
from Espalier.SCARLikelihood import SCAR
from Espalier.sim import ARGSimulator
import Utils
import dendropy
from dendropy.calculate import treecompare
import pandas as pd
import glob
import os
import shutil
import argparse
import numpy as np
import datetime

new_sims = True # simulate new trees/seqs?
rerun_sims = True
argweaver = new_sims # compare with ARGWeaver?
plot = False

redo = False

"Set defaults if not passing arguments"
n_sims = 1
mut_rate = 0.1
sample_size = 10
prior_gamma = 0.01 # exp decay rate on site likelihood prior
lower_bound_ratio = 0.05 # lower bound used in branch and bound algorithm (was 0.01)
path = "./EM-test_k20_mu0.1_L10000_1sim/"
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
results_file = path + "EM-test_k" + str(args.sample_size) + "_mu" + str(args.mut_rate) + "_bestTrueReference_results.csv"
#trees_file = path + "arg-sampler-test_k" + str(args.sample_size) + "_mu" + str(args.mut_rate) + "_trueTrees_noRecCost_trees.csv"

"Set up temp dir for temp tree output"
temp_dir = 'temp-{date:%Y-%m-%d_%H:%M:%S}/'.format(date=datetime.datetime.now())
if not os.path.isdir(temp_dir):
    os.mkdir(temp_dir)

"Sim params"
sample_size = args.sample_size
Ne =  1.0 / 2 # divide by two because msprime assumes individuals are diploid
genome_length = 1e4 # 1e3
recomb_rate_per_genome = 0.1 # rate per genome
rec_rate = recomb_rate_per_genome / genome_length # rate per site
mut_rate = args.mut_rate
M = [[0]]

"If varying rec rate by simulation"
#sim_rec_rates = np.linspace(0.05,2.0,num=n_sims) # used these for 1000bp genome
#sim_rec_rates = np.logspace(-2.0, 1.0, base = 10, num=n_sims)
#sim_rec_rates = np.linspace(0.1,5.0,num=n_sims) # rates used for sims in paper with 10000bp genome
sim_rec_rates = np.linspace(0.5,0.5,num=n_sims) # rates used for sims in paper with 10000bp genome

"Trellis builder params"
prior_gamma = args.prior_gamma
lower_bound_ratio = args.lower_bound_ratio

"Or set natural gamma decay rate based on prob of no recombination"
expected_tree_length = 0.0
for k in range(2,sample_size+1):
    expected_tree_length += k * 2*Ne / (k*(k-1)/2)
prob_recomb_per_site = expected_tree_length * rec_rate
prob_mut_per_site = expected_tree_length * mut_rate
prior_gamma = prob_recomb_per_site

"Just to check exp decay approximation"
#prob_no_recomb_per_site = 1 - prob_recomb_per_site
#mult_prob = prob_no_recomb_per_site**1
#exp_prob = np.exp(-prob_recomb_per_site*1)

"Initialize callable instances of Espalier objects"
raxml = RAxMLRunner(raxml_path='raxml-ng',lsd_path='lsd',temp_dir=temp_dir)

reconciler = Reconciler(raxml,lower_bound_ratio=lower_bound_ratio,prior_gamma=prior_gamma,temp_dir=temp_dir)

argb = ARGBuilder(reconciler,raxml,max_sample_size=10)

"Initialize SCAR model class"
bounds = (0.0,0.001)
scar_model = SCAR(rec_rate,M,Ne,genome_length,bounds=bounds)

"Initialize results dataframe"
results_cols = ['Simulation','True Rec Rate','Rec Rate Estimate','True Recomb Events','Discordant Recomb Events','Inferred Recomb Events']

"Performance metrics"
start_num = args.restart_sim
if start_num == 1 and not redo:
    results_df = pd.DataFrame(columns=results_cols)
    #sims = []
else:
    results_df = pd.read_csv(results_file)
    #sims = results_df['Simulation'].tolist()

"Simulate ARG in msprime"    
min_seg_length = 20
 
"If rerunning sims"
#rerun_sims = [9, 13]
#rerun_sims = [x+1 for x in rerun_sims]
#for sim in rerun_sims:

"If running new sims"
for sim in range(start_num,n_sims+1): 

    print("-"*20)
    print("Simulation = " + str(sim))
    print("-"*20)
    
    "Put SCAR params in dict"
    rec_rate = sim_rec_rates[sim-1] / genome_length # rate per site
    scar_model.rec_rate = rec_rate
    #scar_params = {'Ne': 2*Ne, 'rec_rate': rec_rate, 'M': M, 'genome_length': genome_length}
    
    "Update prior gamma"
    prob_recomb_per_site = expected_tree_length * rec_rate
    prior_gamma = prob_recomb_per_site
    
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
            ts = ARGSimulator.sim_ARG(sample_size=sample_size,Ne=Ne,length=genome_length,recombination_rate=rec_rate,min_breakpoints=2)
            ts.dump(sim_path + 'recomb_placement_tables')
            true_recomb_events = len(ts.tables.nodes.flags[ts.tables.nodes.flags==131072]) / 2 # number of recomb nodes divided by two since each is present twice in nodes tables
            discordant_recomb_events = ARGSimulator.count_topo_changes(ts)
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
                        raxml.get_dated_raxml_tree(seq_files[i], ML_tree_files[i], tip_date_file, rate_file)
        else:
            sim_flag = False
            fasta_file_list = glob.glob(sim_path + '*.fasta') # get all smc files in dir as list
            informative_site_counts = []
            for seq_file in fasta_file_list:
                informative_site_counts.append(Utils.informative_sites(seq_file))
            segments = len(fasta_file_list)
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
        
        "Run EM algorithm"
        out_report = sim_path + "arg-sampler-sim" + str(sim) + "-report.txt"
        
        try:
        
            true_trees = None # tree_files (should be None or true trees will be included in state space)
            #ts_est, rec_rate_est, inferred_recomb_events = ARGSampler.run_EM(ML_tree_files,seq_files,ref,scar_params,lower_bound_ratio=lower_bound_ratio,prior_gamma=prior_gamma,midpoint=False,use_viterbi=True,true_trees=true_trees,report=out_report,temp_dir=temp_dir)
            #tree_path_est, rho_est, inferred_recomb_events  = ARGSampler.run_EM_noconvert(ML_tree_files,seq_files,ref,scar_params,lower_bound_ratio=lower_bound_ratio,prior_gamma=prior_gamma,midpoint=False,use_viterbi=True,true_trees=true_trees,report=out_report,temp_dir=temp_dir)
            ts_est, rec_rate_est, inferred_recomb_events = argb.run_EM(ML_tree_files,seq_files,ref,scar_model)
            
        except Exception as e: 
                
            print(e)
            #sim_flag = True # triggers sim to run again
            sim_flag = rerun_sims # set False if we just want to exit
            
            print("-"*20)
            print("Rerunning simulation = " + str(sim))
            print("-"*20)
            
        else:
            
            #inferred_recomb_events = len(ts_est.tables.nodes.flags[ts_est.tables.nodes.flags==131072]) / 2 # number of recomb nodes divided by two since each is present twice in nodes tables
            results = {'Simulation':sim,
                       'True Rec Rate':rec_rate,
                       'Rec Rate Estimate':rec_rate_est,
                       'True Recomb Events':true_recomb_events,
                       'Discordant Recomb Events':discordant_recomb_events,
                       'Inferred Recomb Events':inferred_recomb_events}
            if redo:
                results_df.loc[sim-1,'Rec Rate Estimate'] = rec_rate_est
                results_df.loc[sim-1,'True Recomb Events'] = true_recomb_events
                results_df.loc[sim-1,'Discordant Recomb Events'] = discordant_recomb_events
                results_df.loc[sim-1,'Inferred Recomb Events'] = inferred_recomb_events  
            else:
                results_df = results_df.append(results, ignore_index=True)

    results_df.to_csv(results_file,index=False)
            
            
            