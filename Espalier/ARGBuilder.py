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
from Espalier import Reconciler
from Espalier import TreeOps
from Espalier.Viterbi import viterbi
from Espalier.Reconciler import ReconciledTree
from Espalier import Dendro2TSConverter
from Espalier import Utils
import dendropy
import numpy as np
from scipy.stats import poisson
from Bio import SeqIO
from operator import attrgetter
from collections import namedtuple
import logging

class ARGBuilder(object):
    
    """
        ARGBuilder reconstructs ARGs from reconciled trees
    """
    
    def __init__(self,reconciler,raxml,max_sample_size=10,**kwargs):
        
        '''             
            Parameters:
                reconciler (Reconciler): instance of Reconciler class used to reconcile local trees
                raxml (RAxML): instance of RAxML runner used to compute seq likelihoods
                max_sample_size (int): maximum number of reconciled trees to sample per genomic region
                
            Optional keyword arguments:
                temp_dir (str): directory to store temp output

        '''
        
        self.reconciler = reconciler
        self.raxml = raxml
        self.max_sample_size = max_sample_size
        self.temp_dir = kwargs.get('temp_dir', './')

    def reconstruct_ARG(self,local_tree_files,seq_files,ref,rec_rate,**kwargs):
        
        """
            
            Reconstruct an approximate ARG from a set of local tree files.
            Local trees are first reconciled against a reference using iterative-regrafting,
            creating a trellis of candidate trees for each genomic region.
            A Viterbi algorithm is then run to select a path of local trees in the ARG.
            
            TODO: Remove true_trees option?
        
            Parameters:     
               local_tree_files (list): list of newick tree files containing initial local trees for each genomic region
               seq_files (list): list of fasta sequence files containing alignment for each genomic region
               ref (dendropy.Tree): reference tree to reconcile local trees against
               rec_rate (float): assumed recombination rate per site
               
            Optional keyword arguments:
                true_trees: true local trees -- only used for debugging
                report (boolean): prints report logging performance if True
               
            Returns:     
               tree_path (list): Tree path containing local trees in ARG
        """
        
        # Or initialize in __init__?
        report = kwargs.get('report', None)
        true_trees = kwargs.get('true_trees', None)
        
        # Get genomic segment info for each genomic region
        self._get_genome_segments(seq_files)
        n_segments = len(self.segments)
    
        # Build trellis of candidate trees by reconciling local trees        
        tree_trellis, sources_array = self._build_trellis(local_tree_files,seq_files,ref)
                
        # Append ref to the sampled trees for each region
        for loc in range(n_segments):
            tree_trellis[loc].append(ref) # add ref to set of sampled trees
            sources_array[loc].append('Consensus')
            
            # Append true trees if provided (only used for debugging purposes)
            # TODO: to be removed
            # if true_trees:
            #     true_tree_file = true_trees[loc] # local tree
            #     true_tree = dendropy.Tree.get(file=open(true_tree_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
            #     true_tree.suppress_unifurcations(update_bipartitions=False)
            #     tree_trellis[loc].append(true_tree) # add ref to set of sampled trees
            #     sources_array[loc].append('True')
            
        # Get seq likelihood for all trees in tree_trellis
        like_array = self._get_like_array(seq_files,tree_trellis)
        
        # Get recombination transition probs between trees based on SPR distances
        trans_probs = self._get_recomb_trans_probs(tree_trellis,rec_rate)
        
        # Select a path of local trees to form the basis of the ARG
        tree_path, opt = viterbi(tree_trellis, trans_probs, like_array)
            
        # Get source of trees in tree_path
        self.tree_source_path = [sources_array[loc][opt[loc]] for loc in range(n_segments)]
        logging.debug("Tree source path: " + ' '.join(self.tree_source_path))
            
        # Write report to file -- no longer used
        # if report:
        #     txt=open(report,"w")
        #     txt.write("Total number of genomic segments: " + str(n_segments) + "\n")
        #     for loc in range(n_segments):
        #         txt.write("-"*20 + "\n")    
        #         txt.write("Segment = " + str(loc) + "\n")
        #         txt.write("-"*20 + "\n")
        #         txt.write("Segment length = " + str(seg_lengths[loc]) + "\n")
        #         txt.write("Subtrees in MAF = " + str(MAF_sizes[loc]) + "\n") 
        #         #txt.write("Prior alt/ref tree ratio = " + str(f"{prior_tree_ratio[loc]:0.2f}") + "\n")
        #         #txt.write("Reconciled trees proposed = " + str(num_trees_proposed[loc]) + "\n")
        #         txt.write("Reconciled trees accepted = " + str(num_trees_accepted[loc]) + "\n")
        #         txt.write("\n")
        #     txt.close()
        
        return tree_path
    
    
    class EMError(Exception):
        pass
    
    def run_EM(self,local_tree_files,seq_files,ref,coal_model,iters=10,max_attempts=5,min_rec_rate=1e-06,**kwargs):
        
        """
            
            EM (Expectation-Maximization) algorithm returns maximum likelihood estimates 
            of demographic params while jointly reconstructing an approximate ARG.
            In the Expectation step, an ARG is sampled using the Viterbi tree path sampler.
            In the Maximization step the demographic params are optimized,
            to find their MLEs under a demographic (e.g. coalescent) model.
            
            TODO: Remove true_trees option and tree_source_path from return?
        
            Parameters:     
               local_tree_files (list): list of newick tree files containing local trees
               seq_files (list): list of fasta sequence file containing alignment for each genome region
               ref (dendropy.Tree): reference tree to reconcile local trees against
               coal_model (SCAR): Coalescent model used to compute likelihood of ARG
               
            Optional keyword arguments:
                iters (int): number of EM iterations to run
                max_attempts (int): maximum number of attempts to sample valid tree_path for ARG in each EM iteration
                true_trees: true local trees -- only used for debugging
                report (boolean): prints report logging performance if True
                min_rec_rate (float): threshold for minimum allowable recombination rate estimated in maximization step
               
            Returns:     
               tree_path (list): Tree path containing local trees in ARG
               rec_rate (float): MLE of recombination rate
               inferred_recomb_events (int): number of recombination events in reconstructed ARG

        """
        
        "Sample an ARG given starting local trees and a reference/consensus" 
        "Smooth local tree sequencing by reconciling against reference"
        
        # Or initialize in __init__?
        report = kwargs.get('report', None)
        true_trees = kwargs.get('true_trees', None)
        
        # Get initial rec_rate from coal_model
        rec_rate = coal_model.rec_rate # recombination rate per site
        
        # Get genomic segment info for each genomic region
        self._get_genome_segments(seq_files)
        tree_intervals = [(s.start,s.end) for s in self.segments]
        n_segments = len(self.segments)
    
        # Build trellis of candidate trees by reconciling local trees        
        tree_trellis, sources_array = self._build_trellis(local_tree_files,seq_files,ref)
        
        # Append ref to the sampled trees for each region
        for loc in range(n_segments):
            tree_trellis[loc].append(ref.clone(depth=2)) # add ref to set of sampled trees
            sources_array[loc].append('Consensus')
            
            # Append true trees if provided (only used for debugging purposes)
            # if true_trees:
            #     true_tree_file = true_trees[loc] # local tree
            #     true_tree = dendropy.Tree.get(file=open(true_tree_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
            #     true_tree.suppress_unifurcations(update_bipartitions=False) # maybe this should be True?
            #     tree_trellis[loc].append(true_tree.clone(depth=2)) # add ref to set of sampled trees
            #     sources_array[loc].append('True')
        
        # Get seq likelihood for all trees in tree_trellis
        like_array = self._get_like_array(seq_files,tree_trellis)
        
        # Get rooted SPR distance between adjacent local trees in trellis
        rSPR_array = _get_rSPR_array(tree_trellis)
        
        # Run EM algorithm iteratively sampling tree paths and opitimizing recomb rate or other demographic params
        rec_rate_samples = []
        for i in range(iters):        
        
            logging.info("-" * 20)
            logging.info("EM step: " + str(i+1))
            logging.info("-" * 20)
            
            step_completed = False # tracks whether current EM iter has sucessfully completed or not
            failed_attempts = 0 # counter for failed attempts to convert tree_path to ARG
        
            while not step_completed:
        
                try:
            
                    ###
                        # Expectation step: sample new tree path and convert to tree sequence
                    ###
                
                    logging.info("Expectation step: Sampling new tree path")
                    trans_probs = self._get_recomb_trans_probs(tree_trellis,rec_rate,rSPR_array=rSPR_array)
                    tree_path, opt = viterbi(tree_trellis, trans_probs, like_array)
                    
                    # Get source of trees in tree_path"
                    tree_source_path = [sources_array[loc][opt[loc]] for loc in range(n_segments)]
                    logging.info("Tree source path: " + ' '.join(tree_source_path))
                
                    # Reconcile tree node heights for equivalent nodes in neighboring trees
                    logging.info("Reconciling linked heights")
                    tree_path = Reconciler.reconcile_linked_heights(tree_path)
                    
                    # Jitter internal node times so no two coalescent events occur at the same time in the same tree
                    tree_path = _jitter_coal_times(tree_path,displace_dt=0.0001)
                    
                    # Add required rec nodes to trees in path
                    logging.info("Adding recombination events to trees in path")
                    tree_path, total_recs_added = add_path_rec_nodes(tree_path)
                    
                    # Convert tree path with recombination nodes to tskit TreeSequence
                    ts = Dendro2TSConverter.convert(tree_path,tree_intervals)
                    inferred_recomb_events = len(ts.tables.nodes.flags[ts.tables.nodes.flags==131072]) / 2 # number of recomb nodes divided by two since each is present twice in nodes tables
                    
                    ###
                        # Maximization step: find params that optimize likelihood
                    ###
                    
                    logging.info("Maximization step: Optimizing model params")
                    rec_rate = coal_model.opt_MLE(ts)
                    rec_rate_samples.append(rec_rate)
                    logging.info("Recombination rate estimate: %s", f'{rec_rate:.6f}')
                                
                except Exception as e: 
                    
                    print(e)
                    failed_attempts += 1
                
                else:
                    
                    step_completed = True
                    
                finally:
                                   
                    if failed_attempts == max_attempts:
                        
                        # Can either raise error
                        #raise EMError("Too many failed attempts to sample valid tree path")
                        
                        # Or crudely estimate rec rate
                        step_completed = True
                        logging.warning("Could not find valid tree path. Crudely estimating recombination rate as # events over tree length")
                        
                        inferred_recomb_events = total_recs_added / 2 # number of recomb events divided by two since each event is added twice
                        tree_lengths = [tr.length() for tr in tree_path]
                        mean_tree_length = np.mean(tree_lengths)
                        rec_rate_per_genome = inferred_recomb_events / mean_tree_length
                        rec_rate = rec_rate_per_genome / coal_model.genome_length
                        rec_rate = max(min_rec_rate, rec_rate)
                        logging.info("Recombination rate estimate: %s", f'{rec_rate:.6f}')
                                
        return ts, rec_rate, inferred_recomb_events
    
    def _build_trellis(self,local_tree_files,seq_files,ref):
    
        """
            Build trellis of candidate trees by reconciling local trees againt reference 
        """
        
        logging.info("Building tree trellis")
        
        # Get taxon_namespace from ref
        taxa = ref.taxon_namespace
        
        # Get concatenated sequence 
        temp_concat_file = self.temp_dir + 'temp_concat_seqs.fasta'
        Utils.concate_aligns(seq_files,temp_concat_file)
        
        # Old stuff logged for ARG reconstruction report
        #MAF_sizes = []
        #num_trees_accepted = []
        
        tree_trellis = []
        previous_trees = dendropy.TreeList()
        sampled_trees = dendropy.TreeList()
        sources_array = []
        
        # Get total length of genome
        total_length = self.segments[-1].end
        
        # Reconcile local trees against reference
        n_segments = len(self.segments)
        for loc in range(n_segments):
            
            #print("-"*20)
            logging.info("Segment = " + str(loc))
            #print("-"*20)
            
            # Get genomic segment info
            seg = self.segments[loc]
            
            # Get local tree for this segment"
            alt_file = local_tree_files[loc] # local tree
            local_tree = dendropy.Tree.get(file=open(alt_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
            #seq_file = seq_files[loc] # path + "disentangler_test1_tree" + str(loc) + ".fasta"
    
            # Get locally adapted consensus for previous segement
            #if loc > 0:
            #   ref = get_consensus_tree(sampled_trees)
            #else:
            #   ref = global_ref
            
            # Get maximum agreement forest
            maf = MAF.get_maf_4cut(ref,local_tree) # switched ref/local_tree order so subtrees in MAF are pruned from ref
            
            # Get set of candidate trees using the iterative branch and bound regrafting algorithm to search for reconciled trees with high likelihoods
            sampled_rec_trees = self.reconciler(ref,local_tree,maf.clone(),temp_concat_file,seg.start,seg.end)
            
            # Compute genomic distances for discounting likelihood contribution from external sites
            left_dists = np.abs(np.arange(0,seg.start) - seg.start) # distances to left of start_pos
            right_dists = np.abs(np.arange(seg.end,total_length) - seg.end) # distances to right of start_pos
            
            # Add sampled trees from previous interval
            new_sample_count = len(sampled_rec_trees)
            for tr in previous_trees:
                tree_file = self.temp_dir + 'temp_ref.tre'
                tr.write(path=tree_file,schema='newick',suppress_annotations=True,suppress_rooting=True)
                site_likes = self.raxml.get_site_likelihoods(tree_file,temp_concat_file)
                local_like = np.sum(site_likes[seg.start:seg.end])  
                left_site_likes = site_likes[:seg.start] * np.exp(-self.reconciler.prior_gamma*left_dists)
                right_site_likes = site_likes[seg.end:] * np.exp(-self.reconciler.prior_gamma*right_dists)
                global_like = np.sum(left_site_likes) + np.sum(right_site_likes)
                tr_like = local_like + global_like
                sampled_rec_trees.append(ReconciledTree(parent = None, tree = tr, index = 'Previous', like = tr_like, source = 'Previous')) # index 'Previous' is assigned just to let us know this was a previously sampled tree
            
            # Sample trees with highest overall likelihood
            if len(sampled_rec_trees) > self.max_sample_size:
                sampled_rec_trees.sort(key=lambda x: x.like)
                sampled_rec_trees = sampled_rec_trees[-self.max_sample_size:] # take last x nodes sorted in ascending order
            sampled_sources = [s.source for s in sampled_rec_trees]
        
            # Count how many sampled trees are from previous interval"
            previous_count =  0
            for n in sampled_rec_trees:
                if n.index == 'Previous':
                    previous_count += 1
            
            # Put sampled trees in TreeList - we deep clone so downstream ops on individual local trees do not impact other trees
            sampled_trees = dendropy.TreeList()
            sampled_trees.taxon_namespace = ref.taxon_namespace
            for nd in sampled_rec_trees:
                sampled_trees.append(nd.tree.clone(depth=2)) # clone tree before putting in tree_array

            # Append ML tree
            sampled_trees.append(local_tree) # append ML tree for segment
            sampled_sources.append('ML')           
            
            logging.info(f"Newly sampled trees : {new_sample_count}" + "; Previously sampled: " + f"{previous_count}")

            #################################
            
            tree_trellis.append(sampled_trees)
            previous_trees = sampled_trees
            sources_array.append(sampled_sources)
            
            # Stats for report
            #MAF_sizes.append(len(maf))
            #num_trees_proposed.append(counter)
            #num_trees_accepted.append(len(sampled_trees))
            
        return tree_trellis, sources_array
    
    def _get_genome_segments(self,seq_files):
    
        """
            Get genomic coordinates and length of each genome segment (i.e. region/interval)
        """

        Segment = namedtuple('Segment', 'start end length')
        self.segments = []
        seg_lengths = []
        n_segments = len(seq_files)
        for loc in range(n_segments):
            seq_file = seq_files[loc]
            seg_length = len(list(SeqIO.parse(seq_file, "fasta"))[0])
            seg_start = int(np.sum(seg_lengths))
            seg_end = seg_start + seg_length # this is actually seg_end + 1 but this makes indexing sites easier
            seg_lengths.append(seg_length)
            self.segments.append(Segment(start=seg_start, end=seg_end, length=seg_length))
    
    
    def _get_like_array(self,seq_files,tree_array):
    
        """
            Compute likelihood for all trees in tree_array
        """
        
        like_array = []
        segments = len(seq_files)
        for loc in range(segments):
            tree_likes = []
            seq_file = seq_files[loc]
            for tr in tree_array[loc]:
                temp_tree_file = self.temp_dir + "temp_ref.tre"
                tr.write(path=temp_tree_file,schema='newick',suppress_annotations=True,suppress_rooting=True)
                logL = self.raxml.get_tree_likelihood(temp_tree_file,seq_file)
                tree_likes.append(logL)
            like_array.append(tree_likes)
        
        return like_array
    
    def _get_recomb_trans_probs(self,tree_trellis,rec_rate,poisson_density=False,rSPR_array=None):
    
        """
            Compute transition probabilities between all pairs of adjacent local trees in tree_trellis
            based on minimum number of recombination events required to explain discordance.
            The transition probs can be computed either assuming that recombination is a Poisson process over an entire genomic region,
            or recombination rates are modeled as instantaneuous events occuring at exact breakpoints.  
            Recombination rate here is the rate per site.
            
            Can provide pre-computed rSPR array so SPR dists are not recomputed each time trans probs are computed.
            
            Parameters:     
               tree_trellis (list): list of TreeLists containting candidate tree states for each genomic region
               rec_rate (float): assumed recombination rate per site
               
            Optional keyword arguments:
                poisson_density: Models recombination as a Poisson process if True
                rSPR_array (3D array-like): pre-computed rSPR distances between pairs of local trees
    
            Returns:     
                trans_probs (3D array-like): transition probabilities where entry n,i,j contains prob of transitioning from tree i in region n-1 to tree j in region n
            
        """
        
        trans_probs = []
        trans_probs.append(None)
        n_segments = len(tree_trellis)
        seg_lengths = [s.length for s in self.segments]
        
        for loc in range(1,n_segments):
            
            # Compute recomb transition probs between trees in previous loc and current loc
            trees = tree_trellis[loc]
            prev_trees = tree_trellis[loc-1]
            probs = np.zeros((len(prev_trees),len(trees))) # matrix with rows corresponding to trees in T-1 and columns in T
            for i, tree_i in enumerate(prev_trees):
                for j, tree_j in enumerate(trees):
                    
                    # Get min number of rec events necessary to reconcile trees
                    if rSPR_array:
                        rSPR = rSPR_array[loc][i][j]
                    else:
                        rSPR = MAF.get_spr_dist(tree_i,tree_j)
                                    
                    # Modeling recombination events as a Poisson process
                    if poisson_density:
                        seg_length = seg_lengths[loc-1]
                        tree_length = tree_i.length()
                        lambda_rec = rec_rate * seg_length * tree_length # recombintation itensity over genomic segment/block
                        probs[i][j] = poisson.pmf(rSPR, lambda_rec)
                    else: # Modeling as instantaneous events
                        seg_length = 1
                        tree_length = tree_i.length()
                        lambda_rec = rec_rate * seg_length * tree_length # recombintation itensity over genomic segment/block
                        probs[i][j] = lambda_rec**rSPR           
                    
            trans_probs.append(probs)
        
        return trans_probs

def _get_rSPR_array(tree_trellis):
    
    """
        Compute rooted subtree-prune-regraft distance between all pairs of adjacent local trees in tree_trellis
    """
    
    rSPR_array = []
    rSPR_array.append(None)
    segments = len(tree_trellis)
    for loc in range(1,segments):
        trees = tree_trellis[loc]
        prev_trees = tree_trellis[loc-1]
        rSPRs = np.zeros((len(prev_trees),len(trees)))
        for i, tree_i in enumerate(prev_trees):
            for j, tree_j in enumerate(trees):
                rSPRs[i][j] = MAF.get_spr_dist(tree_i,tree_j)
        rSPR_array.append(rSPRs)
    
    return rSPR_array

def _jitter_coal_times(tree_path,displace_dt=0.0001):

    """
        Jitter internal node heights (coalescent times) in the same tree that occur at the exact same time 
        by adding a small displacement dt to the length of their children.
        This is only necessary if ARG tree_path is being converted to tskit.TreeSequence. 
    """    

    for tree in tree_path:
        tree.calc_node_ages(ultrametricity_precision=False)
        coal_times = tree.internal_node_ages(ultrametricity_precision=False)
        if len(set(coal_times)) < len(coal_times): # some most not be unique
            existing_times = [] 
            for node in tree.postorder_internal_node_iter():
                while node.age in existing_times:
                    for child in node.child_node_iter():
                        child.edge.length = child.edge.length + displace_dt
                    tree.calc_node_ages(ultrametricity_precision=False)
                existing_times.append(node.age)
                
    return tree_path
    


def add_rec_node(tree,attachment_edge,recomb_time,midpoint=False):
    
    """
        Add a recombination node to tree
        Only used by add_recombination_events()
        TODO: Move to Reconciler or TreeOps?
    """
    
    # Create a new 'graft' node to attach sub_tree
    attachment_parent_node = attachment_edge.tail_node
    attachment_child_node = attachment_edge.head_node
    deep_rooted = False
    if not attachment_parent_node:
        
        # Grafting above root of tree such that graft node will be new root
        logging.debug('Grafting recombinant node above root - deep rooting tree')
        
        tree = TreeOps.add_deep_root(tree)
        deep_rooted = True 
        attachment_parent_node = attachment_edge.tail_node
        attachment_child_node = attachment_edge.head_node
        
    # else: I'm not sure why this else statement was there
        
    attachment_edge_length = attachment_child_node.edge_length
    graft_node = attachment_parent_node.new_child() # this will represent the recombination node

    # Extract child subtree descending from attachment edge and prune it from tree    
    child_subtree = TreeOps.extract_subtree(tree,attachment_child_node)
    tree.prune_subtree(attachment_child_node, update_bipartitions=False, suppress_unifurcations=False) # supress_unifurcations was True, which was why we were losing recombinant nodes we already added

    # Get and set length of graft edge
    if midpoint: # attach to midpoint of attachment edge
        graft_node.edge.length = attachment_edge_length / 2.0
    else: # figure out length of edge from attachment height
        tree.calc_node_ages(ultrametricity_precision=False)
        graft_node_edge_length = attachment_parent_node.age - recomb_time
        if graft_node_edge_length >= 0.0:
            graft_node.edge.length = graft_node_edge_length
        else:
            logging.warning('Graft edge length is less than zero!!')
            graft_node.edge.length = 0.0

    # Reattach child subtree to rec/graft node
    if midpoint:
        child_subtree_root = child_subtree.seed_node
        child_subtree_root.edge.length = attachment_edge_length / 2.0
        graft_node.add_child(child_subtree_root)
    else:
        TreeOps.regraft_subtree(child_subtree,graft_node,recomb_time)
    
    if deep_rooted:
         tree = TreeOps.remove_deep_root(tree)
    
    return tree


def find_recombinant_edge(tree,subtree_taxa):
    
    """
        Find edge where subtree attaches in tree
        Only used by add_recombination_events()
        TODO: Move to Reconciler or TreeOps?
    """
    
    # Find edge where subtree attaches in alt tree 
    # Note: this will be the parent of the first edge that includes all subtree_taxa in its leaf set"
    for edge in tree.postorder_edge_iter():
        edge_taxa = set([lf.taxon.label for lf in edge.head_node.leaf_iter()])
        if subtree_taxa.issubset(edge_taxa): # subtree_taxa are subset of all edge_taxa
            break
    
    attachment_edge = edge
    parent_node = edge.tail_node # parent of recombinant/attachment edge
    child_node = edge.head_node # child of recombinant edge
    child_nodes = parent_node.child_nodes()
    #if len(child_nodes) < 2:
        #print('WTF!') # I think this can only happen if we've already added a recombination node as a unifurcation
    
    # TODO: can remove this since sibling nodes no longer place constraints on recombination event time
    #if child_nodes[0] is child_node:
    #    sibling_node = child_nodes[1]
    #else:
    #    sibling_node = child_nodes[0]
    
    # Find constraints on timing of recomb event
    tree.calc_node_ages(ultrametricity_precision=False)
    parent_time = parent_node.age # max recombination time is height/age of parent node
    child_time = child_node.age

    return attachment_edge, parent_time, child_time,


def find_recombination_events(ref,alt,maf,ref_rec_nodes,alt_rec_nodes):
    
    """
        Find necessary recombination events to reconcile two discordant trees based on their MAF
        Recombination events that need to be added are stored as a named tuple
        Each rec event is added to both local trees as required by tskit TreeSequence 
        Note: this has been corrected with updated node height constraints
        
        TODO: Should we just throw an exception here if min_recomb_time < max_recomb_time since we know the local trees are incompatible?
    """
    
    RecNode = namedtuple('RecNode', 'edge_bitmask recomb_time')
    
    maf.pop() # remove last/largest connected component
    for sub_tree in reversed(maf):
    
        # Get taxa set in subtree
        subtree_taxa = set([lf.taxon.label for lf in sub_tree.leaf_node_iter()])
            
        # Find edge where subtree attaches in alt (and then find equivalent edge in ref)
        attachment_edge_alt, parent_time_alt, child_time_alt = find_recombinant_edge(alt,subtree_taxa)
        
        # Find edge where subtree attaches in ref (and then find equivalent edge in alt)
        attachment_edge_ref, parent_time_ref, child_time_ref  = find_recombinant_edge(ref,subtree_taxa)
        
        # Pick a recombination event time that satisfies time constraints in both alt and ref trees
        max_recomb_time = min(parent_time_alt,parent_time_ref) # event must be below parent nodes in either tree
        min_recomb_time = max(child_time_alt,child_time_ref) # event must be above child time in either tree
    
        #assert max_recomb_time >= min_recomb_time
        if max_recomb_time >= min_recomb_time:
            # Add random rec event time based on min and max allowed rec times
            recomb_time = np.random.uniform(low=min_recomb_time,high=max_recomb_time) # draw a uniform time between these constraints
            ref_rec_nodes.append(RecNode(edge_bitmask=attachment_edge_ref.split_bitmask, recomb_time=recomb_time))
            alt_rec_nodes.append(RecNode(edge_bitmask=attachment_edge_alt.split_bitmask, recomb_time=recomb_time))
        else:
            logging.warning("WARNING: Could not add recombination nodes due to time constraint violations!")
    
    return ref_rec_nodes,alt_rec_nodes


def add_path_rec_nodes(tree_path):
    
    """
        Add recombination events required to reconcile discordances between local trees in tree_path
    """
    
    segments = len(tree_path)
    
    # Create list to hold rec_nodes for each tree in path
    rec_nodes = [[] for x in range(segments)]

    # Find recombination events between local trees in path
    for loc in range(segments-1):
        ref = tree_path[loc]
        alt = tree_path[loc+1]
        ref_rec_nodes = rec_nodes[loc]
        alt_rec_nodes = rec_nodes[loc+1]
        maf = MAF.get_maf_4cut(ref,alt)
        logging.debug("Adding " + str(len(maf)-1) + " recombination events between segments " + str(loc) + " and " + str(loc+1))
        ref_rec_nodes,alt_rec_nodes = find_recombination_events(ref,alt,maf,ref_rec_nodes,alt_rec_nodes)
    
    # Add recombination events to local trees in path
    total_recs_added = 0
    for tree_idx, tree in enumerate(tree_path):
        
        bipartition_dict = tree.split_bitmask_edge_map
        
        # Sort rec_nodes in descending order by time so we add older events first
        tree_rec_nodes = rec_nodes[tree_idx]
        tree_rec_nodes = sorted(tree_rec_nodes, key=attrgetter('recomb_time'), reverse=True)
        total_recs_added += len(tree_rec_nodes)
        for rn in tree_rec_nodes:
                        
            # Add rec event
            rec_edge = bipartition_dict.get(rn.edge_bitmask,None)
            tree = add_rec_node(tree,rec_edge,rn.recomb_time)
            
            # Re-encode bipartitions in case bipartitions were lost after adding rec_node
            tree.encode_bipartitions(suppress_unifurcations=False)
            bipartition_dict = tree.split_bitmask_edge_map
            
    logging.info("Added " + str(total_recs_added) + " recombination nodes to tree path")
            
    return tree_path, total_recs_added
    
    