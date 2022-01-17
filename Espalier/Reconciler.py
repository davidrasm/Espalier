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

from Espalier import TreeOps
import numpy as np
import math
from Bio import SeqIO
import logging

# For profiling only
# import time

class Reconciler(object):
    
    """
        A callable class used to reconcile topological differeces between pairs of trees
    """
    
    def __init__(self,raxml,lower_bound_ratio=0.1,prior_gamma=0.0,**kwargs):
        
        '''             
            Parameters: 
                raxml (RAxML): instance of RAxML runner used to compute seq likelihoods
                lower_bound_ratio (float): lower bound ratio used by branch-and-bound algorithm
                prior_gamma (float): rate at which likelihood contribution from external sites decays exponentially
               
            Optional keyword arguments:
                temp_dir (str): directory to store temp output
        '''
        
        self.raxml = raxml
        self.lower_bound_ratio = lower_bound_ratio        
        self.prior_gamma = prior_gamma
        self.temp_dir = kwargs.get('temp_dir', './')
    
    def __call__(self,ref,alt,maf,seq_file,start_pos=None,end_pos=None):
        
        """
        
            Reconcile a pair of trees (ref and alt) through their MAF using iterative branch and bound regrafting
            Produces a set of reconciled trees with high likelihood given the sequence data
            If reconciling a local tree for a given genome region, can specify start and end position of region
        
            Parameters:     
               ref (dendropy.Tree): input ref tree (T1)
               alt (dendropy.Tree): input alt tree (T2)
               maf (dendropy.TreeList): maximum agreement forest between ref and alt
               seq_file (str): fasta sequence file containing alignment
               
            Optional keyword arguments:
                start_pos (int): start position in seq alignment
                end_pos (int): end position in seq alignment
               
            Returns:     
               sampled_trees (list(ReconciledTrees)): List of reconciled trees
        """
                        
        # Get dictionary of sequence ids : seq records
        # TODO: Could just pass seq_dict from init so we don't need to repeatedly import alignment
        my_records = list(SeqIO.parse(seq_file, "fasta"))
        
        # Remove underscores from taxon names so compatible with dendropy taxon names
        #seq_dict = {record.id:record for record in my_records}
        seq_dict = {record.id.replace('_',' '):record for record in my_records}
        total_length = len(seq_dict[next(iter(seq_dict))].seq)
        
        # Assume entire align in seq_file is genomic region if start/end positions not given
        if not start_pos:
            start_pos = 0
        if not end_pos:
            end_pos = total_length
        
        # Compute genomic distances for discounting likelihood contribution from external sites"
        left_dists = np.abs(np.arange(0,start_pos) - start_pos) # distances to left of start_pos
        right_dists = np.abs(np.arange(end_pos,total_length) - end_pos) # distances to right of start_pos
        
        # Reverse order MAF subtrees so we start with largest connected component
        maf = list(reversed(maf))
        
        logging.info("Starting regrafting")
        logging.info("Subtrees in MAF: %s", str(len(maf)))
        
        # Create root node in search tree and add to queue
        root_tree = maf[0].clone(depth=2)
        root = ReconciledTree(parent = None, tree = root_tree, index = 0)
        queue = [] # queue of "open" nodes we still need to process
        queue.append(root)
        
        sampled_trees = [] # sampled trees represent regrafted (full) trees
        
        nodes_visited = 0
        ref_count = 0
        alt_count = 0
        while queue:
            
            parent = queue.pop()
            
            maf_index = parent.index+1
            
            nodes_visited += 1
            if nodes_visited % 10 == 0:
                logging.info("Processing next search node : %s", str(nodes_visited))
            
            if maf_index == len(maf): # this is a full tree
            
                if parent.tree.__len__() != ref.__len__():
                    logging.warning("Trees of unequal size found!!")
                parent.source = 'Reconciled'
                sampled_trees.append(parent)
                
                logging.debug("Reached tip node in search tree. Output tree #: %s", str(len(sampled_trees)))
            
            else: # node contains partially reconciled tree
            
                sub_tree = maf[maf_index]
                rec_tree = parent.tree # should we clone rec_tree here?
                    
                # Find attachment edge (in rec_tree) for both alt and ref trees
                logging.debug("Finding attachment edges")
                attachment_edge_alt = find_attachment_edge(rec_tree,alt,sub_tree)
                rec_tree_clone = rec_tree.clone(depth=2)
                attachment_edge_ref = find_attachment_edge(rec_tree_clone,ref,sub_tree)
                

                # Attach sub_tree to rec_tree by adding a new "graft" node along attachment edge        
                logging.debug("Attaching subtrees")
                sub_tree_clone = sub_tree.clone(depth=2) #Should we also deep clone sub_tree?
                sub_tree_2 = sub_tree.clone(depth=2)
                rec_tree_ref = attach_subtree(rec_tree_clone,sub_tree_clone,ref,attachment_edge_ref)
                rec_tree_alt = attach_subtree(rec_tree,sub_tree_2,alt,attachment_edge_alt) # rec_tree_ref.clone(depth=2)
                
                if rec_tree_alt.__len__() > 3: # min tree size for RAxML is 4
                    
                    logging.debug("Computing site likelihoods")
                
                    # Write alt reconciled tree file
                    alt_tree_file = self.temp_dir + 'temp_alt.tre'
                    rec_tree_alt.write(path=alt_tree_file,schema='newick',suppress_annotations=True,suppress_rooting=True)
                    
                    # Write ref reconciled tree file
                    ref_tree_file = self.temp_dir + 'temp_ref.tre'
                    rec_tree_ref.write(path=ref_tree_file,schema='newick',suppress_annotations=True,suppress_rooting=True)
                    
                    # Write fasta seq file for taxa currently in rec_tree
                    temp_seq_file = self.temp_dir + 'temp_seqs.fasta'
                    temp_records = [seq_dict[nd.taxon.label] for nd in rec_tree_alt.leaf_node_iter()]
                    SeqIO.write(temp_records, temp_seq_file, "fasta")
                    
                    # Get site likelihoods (for all sites) given both alt and ref reconciled trees 
                    # tic = time.perf_counter()
                    alt_site_likes = self.raxml.get_site_likelihoods(alt_tree_file,temp_seq_file)
                    ref_site_likes = self.raxml.get_site_likelihoods(ref_tree_file,temp_seq_file)
                    # toc = time.perf_counter()
                    # elapsed = toc - tic
                    # print(f"Elapsed time for site likelihoods: {elapsed:0.4f} seconds")
                    
                    # Compute local (internal) site likelihoods for genomic region                    
                    alt_local_like = np.sum(alt_site_likes[start_pos:end_pos])  
                    ref_local_like = np.sum(ref_site_likes[start_pos:end_pos])
        
                    # Compute likelihoods for sites to left and right of genomic region and discount by distance
                    left_alt_site_likes = alt_site_likes[:start_pos] * np.exp(-self.prior_gamma*left_dists)
                    right_alt_site_likes = alt_site_likes[end_pos:] * np.exp(-self.prior_gamma*right_dists)
                    left_ref_site_likes = ref_site_likes[:start_pos] * np.exp(-self.prior_gamma*left_dists)
                    right_ref_site_likes = ref_site_likes[end_pos:] * np.exp(-self.prior_gamma*right_dists)
                    
                    # Sum (log) likelihoods for all sites outside of genomic region
                    alt_global_like = np.sum(left_alt_site_likes) + np.sum(right_alt_site_likes)  
                    ref_global_like = np.sum(left_ref_site_likes) + np.sum(right_ref_site_likes) 
                
                else: # if tree size is 3 or smaller cannot compute likelihoods
                    
                    alt_local_like = 0 
                    ref_local_like = 0
                    alt_global_like = 0 
                    ref_global_like = 0
        
                # Compute total likelihood of each tree for all sites
                alt_like = alt_local_like + alt_global_like
                ref_like = ref_local_like + ref_global_like
                logging.debug('Alt / ref likelihood = ' + f'{alt_like:.2f}' + ' / ' + f'{ref_like:.2f}')
                
                logging.debug("Updating search tree")
                
                # Add reconciled trees to queue
                if ref_like > alt_like:
                    child = ReconciledTree(parent = parent, tree = rec_tree_ref, index = maf_index, like = ref_like)
                    queue.append(child)
                    ref_count += 1
                    ratio = math.exp(alt_like - ref_like)
                    if ratio > self.lower_bound_ratio:
                        child = ReconciledTree(parent = parent, tree = rec_tree_alt, index = maf_index, like = alt_like)
                        queue.append(child)
                        alt_count += 1
                        logging.debug("Added lower like node to search tree. Like ratio = " + f'{ratio:.2f}')
                else:
                    child = ReconciledTree(parent = parent, tree = rec_tree_alt, index = maf_index, like = alt_like)
                    queue.append(child)
                    alt_count += 1
                    ratio = math.exp(ref_like - alt_like)
                    if ratio > self.lower_bound_ratio:
                        child = ReconciledTree(parent = parent, tree = rec_tree_ref, index = maf_index, like = ref_like)
                        queue.append(child)
                        ref_count += 1
                        logging.debug("Added lower like node to search tree. Like ratio = " + f'{ratio:.2f}')
        
        logging.info("Alt / ref reattachments = " + f'{alt_count}' + ' / ' + f'{ref_count}')
        
        return sampled_trees


class ReconciledTree:
    """
        Simple class for storing reconciled trees along with their metadata used for post-processing
    """
    def __init__(self, parent=None, tree=None, index=-1, like=0.0, source=None):
        self.parent = parent
        self.tree = tree
        self.index = index
        self.like = like
        self.source = source

def find_attachment_edge(rec_tree,guide_tree,subtree):
    
    """
    
        Find edge in (reconciled) rec_tree where subtree attaches.
        Equivalent to the subtree's position in guide_tree.
            
        Proceeds in two steps:
            
        Step 1: Find attachment edge in guide_tree
        Attachment edge will be the first edge that includes all subtree_taxa in its leaf set
        
        Step 2: Find edge equivalent to the attachment edge in rec_tree
        Note: this edge will be the parent of the edge leading to the subtree after regrafting
        Hence we find this edge by the parent node's taxa set using edge.tail
        Equivalent edge will have all taxa in attachment_edge_taxa_set that are also taxa in the rec_tree
        But there may be no equivalent edge to attach to if none of the taxa in attachment_edge_taxa_set have been added to rec_tree yet
        So if taxa_set_intersect is null set move to parent of attachment edge until union is not null set
    
        Parameters:     
           rec_tree (dendropy.Tree): partially reconciled tree to which subtree is attached
           guide_tree (dendropy.Tree): full tree used to determine where to attach subtree
           subtree (dendropy.Tree): subtree to attach
           
        Returns:     
            attachment_edge (dendropy.Edge): edge to attach subtree to
    """
    
    # Get taxa set in subtree and rectree we are attaching"
    subtree_taxa = set([lf.taxon.label for lf in subtree.leaf_node_iter()])
    rectree_taxa = set([lf.taxon.label for lf in rec_tree.leaf_node_iter()])
    
    # Step 1: Find attachment edge in guide_tree
    for edge in guide_tree.postorder_edge_iter():
        edge_taxa = set([lf.taxon.label for lf in edge.head_node.leaf_iter()])
        if subtree_taxa.issubset(edge_taxa): # subtree_taxa are subset of all edge_taxa
            break
    if not edge.tail_node:
        print("Error finding attachment edge: edge has no parent/tail")
    
    # Step 2: Find equivalent edge in rec_tree to attach subtree 
    attachment_edge_taxa_set = set([lf.taxon.label for lf in edge.tail_node.leaf_iter()]) # tail_node because we need parent edge of node we will regraft to 
    taxa_set_intersection = attachment_edge_taxa_set & rectree_taxa
    while not taxa_set_intersection: # while intersect is null set
        edge = edge.tail_node.edge # get parent edge
        attachment_edge_taxa_set = set([lf.taxon.label for lf in edge.tail_node.leaf_iter()])
        taxa_set_intersection = attachment_edge_taxa_set & rectree_taxa
    for edge in rec_tree.postorder_edge_iter():
        edge_taxa = set([lf.taxon.label for lf in edge.head_node.leaf_iter()])
        if not taxa_set_intersection.difference(edge_taxa): 
            break # if there is no difference
    attachment_edge = edge # edge subtree will be attached to in rec_tree
    
    return attachment_edge

def attach_subtree(rec_tree,sub_tree,attachment_tree,attachment_edge,midpoint=False):
    
    """
    
        Attach sub_tree to rec_tree along attachment_edge by regrafting.
        Attachment tree (should be called guide tree) is only used to determine height of node subtree will be regrafted to.
        
        Parameters:     
           rec_tree (dendropy.Tree): partially reconciled tree to which subtree is attached
           sub_tree (dendropy.Tree): subtree to attach
           guide_tree (dendropy.Tree): full tree used to determine regraft node height
           attachment_edge (dendropy.Edge): edge to attach subtree to
           
        Optional keyword arguments:
            midpoint = If true, subtrees will be regrafted to node at midpoint of attachment_edge 
            
        Returns:     
            rec_tree (dendropy.Tree): Reconciled tree with subtree regrafted to attachment edge
    """
    
    min_edge_length = 0.0 # minimum branch length
    
    # Create a new 'graft' node to attach sub_tree and sister of sub_tree
    attachment_parent_node = attachment_edge.tail_node
    attachment_sister_node = attachment_edge.head_node
    deep_rooted = False
    if not attachment_parent_node:
        
        # Grafting above root of rec tree such that graft node will be new root
        # This can happen if cut one of the root's child edges
        
        logging.debug("Deep rooting tree")
        
        rec_tree = TreeOps.add_deep_root(rec_tree) # reroot == True by default now
        deep_rooted = True
        attachment_parent_node = attachment_edge.tail_node
        attachment_sister_node = attachment_edge.head_node
    
    # Get temporal constraints on graft node height:
    rec_tree.calc_node_ages(ultrametricity_precision=False)
    attachment_parent_height = attachment_parent_node.age
    max_height =  attachment_parent_height # Graft node must be below attrachment parent node
    sister_subtree_height = attachment_sister_node.age # But above the root node of the grafted sub_tree and sister_subtree
    sub_tree.calc_node_ages(ultrametricity_precision=False)
    subtree_root = sub_tree.seed_node # TreeOps.get_root_node(sub_tree)
    subtree_height = subtree_root.age
    min_height = max(sister_subtree_height,subtree_height)
    
    if attachment_parent_height < sister_subtree_height:
        logging.warning("Parent height is below sister subtree height")
    if attachment_parent_height < subtree_height:
        logging.warning("Parent height is below subtree height")
    
    #if attachment_parent_height < min_height:
        #label_fn = lambda nd: nd.taxon.label if nd.is_leaf() else str(f"{nd.age:0.2f}")
        #rec_tree.print_plot(show_internal_node_labels=True,node_label_compose_fn=label_fn)
        #sub_tree.print_plot(show_internal_node_labels=True,node_label_compose_fn=label_fn)
        #attachment_tree.calc_node_ages(ultrametricity_precision=False)
        #attachment_tree.print_plot(show_internal_node_labels=True,node_label_compose_fn=label_fn)
        #logging.warning("Parent height is below sister or subtree height")

    # Add new graft node to attach sub_tree to
    attachment_edge_length = attachment_sister_node.edge_length
    graft_node = attachment_parent_node.new_child()

    # Alternative way to exact - does not remove dummy 'new_child' node we've grafted above
    #sister_subtree = rec_tree.extract_tree_with_taxa(taxa=set([lf.taxon for lf in attachment_sister_node.leaf_iter()]))
    
    # Extract sister subtree descending from attachment edge and prune it from rec_tree
    sister_subtree = TreeOps.extract_subtree(rec_tree,attachment_sister_node)
    rec_tree.prune_subtree(attachment_sister_node, update_bipartitions=False, suppress_unifurcations=True)

    # Get and set length of graft edge
    if midpoint:
        graft_node.edge.length = attachment_edge_length / 2.0
    else:
        subtree_taxon = sub_tree.leaf_nodes()[0].taxon.label # get one leaf in subtree
        sister_subtree_taxon = sister_subtree.leaf_nodes()[0].taxon.label # get one leaf in sister subtree
        attachment_tree.calc_node_ages(ultrametricity_precision=False)
        mrca = attachment_tree.mrca(taxon_labels=[subtree_taxon,sister_subtree_taxon]) #Should this be tree specific (rec/alt?)
        graft_node_age = mrca.age
        
        #Check constraints on graft node height/age"
        if graft_node_age < min_height or graft_node_age > max_height:
            graft_node_age = min_height + (max_height - min_height)/2.0 # center between max and min
        graft_node.edge.length = attachment_parent_height - graft_node_age #mrca.edge_length

    # Reattach sister subtree to graft node
    if midpoint:
        sister_subtree_root = sister_subtree.seed_node
        sister_subtree_root.edge.length = attachment_edge_length / 2.0
        graft_node.add_child(sister_subtree_root)
    else:
        TreeOps.regraft_subtree(sister_subtree,graft_node,graft_node_age,min_edge_length=min_edge_length)
    
    
    # Attach sub_tree to graft node"
    if midpoint:
        subtree_root = sub_tree.seed_node
        subtree_root.edge.length = attachment_edge_length / 2.0
        graft_node.add_child(subtree_root)
    else:
        TreeOps.regraft_subtree(sub_tree,graft_node,graft_node_age,min_edge_length=min_edge_length)
    
    
    # Remove deep root if deep rooted
    if deep_rooted:
        rec_tree = TreeOps.remove_deep_root(rec_tree)
        
    # Double check that there are no negative branch lengths
    for edge in rec_tree.postorder_edge_iter():
        if edge.length:
            if edge.length < 0.0:
                edge.length = min_edge_length
                logging.debug("Negative branch length in reconciled tree")
    
    # Recheck rooting
    if not rec_tree.is_rooted:
        logging.warning("Reconciled tree has lost rooting")   

    return rec_tree

def reconcile_linked_heights(tree_path, min_branch_length=0.01):
    
    """
        Reconcile node heights for equivalent nodes across a tree path using the "lock and chain" method.
        Equivalent nodes in neighboring local trees will be assigned the same height.
        Here, equivalent means two nodes have same descendent subtree as defined by their edge bipartitions
        Optional: minimum branch length threshold can be set so no branches have zero length
        
        TODO: Find heights that maximize the likelihood of the sequence data across the linked blocks rather than just setting median heights.
        
        Parameters:     
            tree_path (TreeList or list(dendropy.Tree)): tree_path
               
        Optional keyword arguments:
            min_branch_length (float): min branch length threshold
        
    """
    
    path_length = len(tree_path)
    bipartition_maps = []
    for tree in tree_path:
        tree.calc_node_ages(ultrametricity_precision=False)
        tree.encode_bipartitions() # need to encode bipartitions before we can access them
        bipartition_maps.append(tree.split_bitmask_edge_map) # create dict the maps edges to edge (split) bitmasks
    
    # Set state of all nodes to unlocked. Once locked a nodes age cannot be reset, which prevents repeatedly updating the same node's age
    for tree_idx, tree in enumerate(tree_path):
        for node in tree.postorder_internal_node_iter():
            node.locked = False
    
    for tree_idx, tree in enumerate(tree_path):
        
        logging.debug("Starting node height reconciliation on tree # " + str(tree_idx))
        
        for node in tree.postorder_internal_node_iter():
            
            if not node.locked:
            
                # Find equivalent edges in adjacent trees until there is no such equivalent edge
                edge = node.edge
                linked_trees = [tree] # list of linked trees with equivalent edge
                linked_heights = [node.age]
                linked_edges = [edge]
                min_constraint = max([nd.age for nd in node.child_node_iter()])
                linked_min_constraints = [min_constraint]
                for idx in range(tree_idx+1,path_length):
                    edge = bipartition_maps[idx].get(edge.split_bitmask,None)
                    if edge:
                        linked_trees.append(tree_path[idx]) #(idx)
                        linked_heights.append(edge.head_node.age)
                        linked_edges.append(edge)
                        min_constraint = max([nd.age for nd in edge.head_node.child_node_iter()])
                        linked_min_constraints.append(min_constraint)
                    else:
                        break
        
                # Set new node height/age to median of linked hights
                new_height = np.median(linked_heights)
                min_constraint = max(linked_min_constraints) + min_branch_length # added so node has to be above closest child by min_branch_length
                if new_height < min_constraint:
                    new_height = min_constraint
                
                # Reset edge lengths in all linked trees
                for edge_idx, edge in enumerate(linked_edges):
                    parent = edge.head_node
                    parent.locked = True
                    for child in parent.child_node_iter():
                        new_edge_length = new_height - child.age
                        if new_edge_length < 0:
                            logging.warning("Child edge length negative after height reconciliation!!")
                        child.edge.length = new_edge_length
                    
                    # Recompute node ages here to reflect new edge lengths
                    linked_trees[edge_idx].calc_node_ages(ultrametricity_precision=False)

    return tree_path
