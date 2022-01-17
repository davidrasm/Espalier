##########################
##
## Espalier: A Python package for tree reconciliation and reconstructing ARGs using maximum agreement forests.
##
## Copyright 2021-2022 David A. Rasmussen (drasmus@ncsu.edu)
##
## If using Espalier or this code, please cite:
##
##      Rasmussen, D.A. and Guo, F. Espalier: Efficient tree reconciliation and ARG reconstruction using maximum agreement forests. 2022.
##
############################

from Espalier import MAF
import msprime
import dendropy
import pyvolve

def sim_ARG(sample_size=10,Ne=100,length=1e3,recombination_rate=5e-6,min_breakpoints=1,max_breakpoints=1000,plot=False):
    
    """
        Simulate ARGs with local trees in a tree sequence using msprime.
        
        Parameters:     
           sample_size (int): Number of tips/samples in ARG
           Ne (float): Haploid effective population size
           length (int): Genome length
           recombination_rate (float): recombination rate per site per lineage
           min_breakpoints (int): minumum number of allowed breakpoints in simulated ARG
           max_breakpoints (int): maximum number of allowed breakpoints in simulated ARG
           plot (boolean): display ARG as local trees and ts tables
           
        Returns:     
           ts (tskit.TreeSequence): TreeSequence representation of simulated ARG
    """   
    
    breaks = 0
    while breaks < min_breakpoints or breaks > max_breakpoints: 
        ts = msprime.simulate(sample_size=sample_size, Ne=Ne, length=length, recombination_rate=recombination_rate, record_full_arg=True)
        breaks = len(ts.breakpoints(as_array=True)) - 2 # -2 because tskit counts ends as breakpoints
    
    for tr_num, tree in enumerate(ts.trees()):
        if plot:
            print("-" * 20)
            print("tree {}: interval = {}".format(tree.index, tree.interval))
            print(tree.draw(format="unicode"))
            
    # Count topo changes in tree sequences
    topo_changes = count_topo_changes(ts)
    
    if plot:
        print(ts.tables.nodes)
        print()
        print(ts.tables.edges)
        print()
        print("Recombination breakpoints: " + str(breaks))
        print()
        print("Topological changes in tree sequence: " + str(topo_changes))
    
    return ts

def sim_seqs(tree_file,seq_file,mut_rate=0.001,seq_length=1000,freqs=[0.25,0.25,0.25,0.25],kappa=2.75):
    
    """
        Simulate sequences for local trees in ARG under a HKY model in pyvolve
        
        Parameters:     
           tree_file (str): Newick file with local tree for which sequences will be simulated
           seq_file (str): Output fasta file for simulated sequences
           
        Optional keyword arguements:   
           mut_rate (float): Mutation rate per site
           seq_length (int): Genome length
           freqs (list[float]): Equilibrium nucleotide frequencies
           kappa (float): transition/transversion ratio for HKY model  
        
    """
    
    # Read in phylogeny along which Pyvolve should simulate seqs
    my_tree = pyvolve.read_tree(file = tree_file, scale_tree = mut_rate) # scale_tree sets absolute mutation rate
    #pyvolve.print_tree(my_tree) # Print the parsed phylogeny
    
    # Parameterize HKY substitution model for sim
    nuc_model = pyvolve.Model( "nucleotide", {"kappa":kappa, "state_freqs":freqs})
    
    # Define a Partition object which evolves set # of positions
    my_partition = pyvolve.Partition(models = nuc_model, size = seq_length)
    
    # Define an Evolver instance to evolve a single partition
    my_evolver = pyvolve.Evolver(partitions = my_partition, tree = my_tree) 
    
    # Evolve sequences with custom file names 
    # my_evolver(ratefile = "AMR_ratefile.txt", infofile = "AMR_infofile.txt", seqfile = "AMR-seqsim.fasta" )
    # **ratefile** is a custom name for the "site_rates.txt" file. Provide None or False to suppress file creation.
    # **infofile** is a custom name for the "site_rates_info.txt" file. Provide None or False to suppress file creation.
    my_evolver(seqfile = seq_file, ratefile = None, infofile = None)

def count_topo_changes(ts):
    
    """
        Counts number of topological changes between local trees in a TreeSequence.
        
        Parameters:     
           ts (tskit.TreeSequence): Tree sequence
           
        Returns:     
           changes (int): number of topological changes
    """   
    
    changes = 0
    prev_tree = None
    for tr_num, tree in enumerate(ts.trees()):
        if tr_num > 0:
            taxa = dendropy.TaxonNamespace()
            ref = dendropy.Tree.get(data=prev_tree, schema="newick", rooting="default-rooted", taxon_namespace=taxa)
            alt = dendropy.Tree.get(data=tree.newick(), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
            if MAF.test_discordance(ref,alt):
                changes += 1
        prev_tree = tree.newick()
        
    return changes


    


