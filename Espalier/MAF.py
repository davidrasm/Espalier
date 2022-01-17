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
import dendropy

def get_maf_4cut(ref_in,alt_in,plot=False):

    """
    
        Computes a MAF between two trees using the 4-cut algorithm. 
        A discordant edge is identified at each iteration and
        four possible edges surrounding that edge are cut.
        If a given edge in ALT is discordant with children A,C
        we find the sibling of A in REF and call this B,
        and then find the sibling of C in REF and call this D.
        B or D may not exist in ALT
        The edge (A,B,C or D) that minimizes discordance in the pruned trees is then cut.
        IMPORTANT: Subtrees in forest will be cut from ALT such that branch lengths in the MAF will reflect those in ALT
    
        Parameters:     
           ref_in (dendropy.Tree): input ref tree (T1)
           alt_in (dendropy.Tree): input alt tree (T2) 
           
        Optional keyword arguments:
            plot (boolean): plots subtrees in MAF if True
           
        Returns:     
           forest (dendropy.TreeList): the MAF as a list of trees
    """
    
    #Deep copy original input trees so they are not modified
    ref = ref_in.clone(depth=2)
    alt = alt_in.clone(depth=2)
    
    # Deep root trees by adding dummy node as outgroup
    #if deep_root:
        #ref,alt = deep_root_pair(ref,alt)
    
    # Encode bipartitions of leaf taxa set in each tree
    ref.encode_bipartitions()    
    alt.encode_bipartitions()
    
    # Create dict mapping edge/split bitmasks (integers) to their corresponding edge 
    ref_bipartition_dict = ref.split_bitmask_edge_map 
    #for key, value in ref_bipartition_dict.items(): print(key, value)
    
    # Check for discordance between trees
    discordant = _update_concordance(alt,ref_bipartition_dict)
        
    # Run MAF algorithm
    forest = dendropy.TreeList() # forest to hold subtrees in MAF
    counter = 0
    while discordant:
                    
        # Find first discordant edge
        for edge in alt.postorder_edge_iter():
            if not edge.concordant: break
            
        parent = edge.head_node #first discordant sibling pair (edge) in t1
        edge_cut_list = parent.child_edges()
        
        # Update edge bipartition dicts
        ref_bipartition_dict = ref.split_bitmask_edge_map
        alt_bipartition_dict = alt.split_bitmask_edge_map
        
        # Get sibling of A and C in ALT
        edge_A = edge_cut_list[0]
        edge_C = edge_cut_list[1]
        
        # Get sibling edge B 
        ref_edge_A = ref_bipartition_dict.get(edge_A.leafset_bitmask,None)
        parent_edge_A = ref_edge_A.tail_node
        child_edges = parent_edge_A.child_edges()
        if child_edges[0] is ref_edge_A:
            sibling_edge = child_edges[1]
        else:
            sibling_edge = child_edges[0]
        alt_edge_B = alt_bipartition_dict.get(sibling_edge.leafset_bitmask,None)
        if alt_edge_B:
            edge_cut_list.append(alt_edge_B)
        
        # Get sibling edge D
        ref_edge_C = ref_bipartition_dict.get(edge_C.leafset_bitmask,None)
        parent_edge_C = ref_edge_C.tail_node
        child_edges = parent_edge_C.child_edges()
        if child_edges[0] is ref_edge_C:
            sibling_edge = child_edges[1]
        else:
            sibling_edge = child_edges[0]
        alt_edge_D = alt_bipartition_dict.get(sibling_edge.leafset_bitmask,None)
        if alt_edge_D:
            edge_cut_list.append(alt_edge_D)
        
        # Cut each edge in edge_cut_list and determine which minimizes discordance"
        D = [] # array to hold discordance metric
        alt_cut_trees = []
        ref_cut_trees = []
        components = [] # list to hold cut subtrees
        for edge in edge_cut_list:
            
            alt_cut = alt.clone(depth=2) #Deep copy tree
            ref_cut = ref.clone(depth=2) #Deep copy tree
            
            alt_cut_bipartition_dict = alt_cut.split_bitmask_edge_map
            ref_cut_bipartition_dict = ref_cut.split_bitmask_edge_map
            bitmask = edge.leafset_bitmask
            alt_cut_edge = alt_cut_bipartition_dict.get(bitmask,None)
            ref_cut_edge = ref_cut_bipartition_dict.get(bitmask,None)
            
            # Make sure we're getting the correct edges
            #edge_taxa = set([lf.taxon.label for lf in edge.head_node.leaf_iter()])
            #alt_cut_edge_taxa = set([lf.taxon.label for lf in alt_cut_edge.head_node.leaf_iter()])
            #ref_cut_edge_taxa = set([lf.taxon.label for lf in ref_cut_edge.head_node.leaf_iter()])

            # Extract subtree from alt and place in subtree components"
            child_nd = alt_cut_edge.head_node
            component = TreeOps.extract_subtree(alt_cut,child_nd)
            components.append(component)
        
            # Prune subtree from alt_cut"
            alt_cut.prune_subtree(child_nd, update_bipartitions=True, suppress_unifurcations=True)
        
            # Find and prune corresponding edge in ref_cut tree
            ref_cut_child = ref_cut_edge.head_node
            ref_cut.prune_subtree(ref_cut_child, update_bipartitions=True, suppress_unifurcations=True)
            
            # Update concordance
            ref_cut_bipartition_dict = ref_cut.split_bitmask_edge_map
            discordant = _update_concordance(alt_cut,ref_cut_bipartition_dict)
            D.append(_count_discordant_edges(alt_cut,ref_cut_bipartition_dict))
            alt_cut_trees.append(alt_cut)
            ref_cut_trees.append(ref_cut)
            
        # Decide which cut minimizes D and add cut subtree to forest
        cut = D.index(min(D))
        alt = alt_cut_trees[cut] # Set alt for next iter
        ref = ref_cut_trees[cut] # Set ref for next iter
        cut_tree = components[cut]
        forest.append(components[cut])
        
        counter += 1
        if plot:
            print("-" * 20)
            print("MAF iteration: " + str(counter))
            print("Reference tree:")
            ref.print_plot()
            print()
            print("Alternate tree:")
            alt.print_plot()
            print()
            print("Cut tree:")
            if len(cut_tree.nodes()) > 1:
                cut_tree.print_plot()
            else:
                print('---' + cut_tree.nodes()[0].taxon.label)
            print()
            print("-" * 20)
           
        # Make sure discordant is updated 
        ref_bipartition_dict = ref.split_bitmask_edge_map 
        discordant = _update_concordance(alt,ref_bipartition_dict)
    
    # Remove dummy node at deep root
    #if deep_root:
        #alt = remove_deep_root(alt)
    
    forest.append(alt) # append what's left of alt (should at least be a rooted tree with two edges)
    return forest

def get_maf_3approx(ref_in,alt_in,deep_root=True,plot=False):
    
    
    """
    
        Computes a MAF between two trees using the 3-approximation of Whidden and Zeh (2009)
        If a given edge in ALT is discordant with children A,C
        We find sibling of A in REF and call this B
        Then cut edge A,B and C placing resulting sub
        Note that default is to deep_root tree otherwise root edge is often cut and algorithm fails
    
        Parameters:     
           ref_in (dendropy.Tree): input ref tree (T1)
           alt_in (dendropy.Tree): input alt tree (T2) 
           
        Optional keyword arguments:
            deep_root (boolean): adds deep root both ALT and REF if true
            plot (boolean): plots subtrees in MAF if True
           
        Returns:     
           forest (dendropy.TreeList): the MAF as a list of trees
    """
    
    # Deep copy original trees so we don't lose them
    ref = ref_in.clone(depth=2)
    alt = alt_in.clone(depth=2)
    
    # Deep root trees by adding dummy node as outgroup
    if deep_root:
        ref,alt = TreeOps.deep_root_pair(ref,alt)
    
    # Encode bipartitions of leaf taxa set in each tree
    ref.encode_bipartitions()
    alt.encode_bipartitions()
    
    # Create dict mapping edge/split bitmasks (integers) to their corresponding edge 
    ref_bipartition_dict = ref.split_bitmask_edge_map #dictionary mapping split bitmasks (integers) to their corresponding edges through the Tree.split_bitmask_edge_map attribute
    #for key, value in ref_bipartition_dict.items(): print(key, value)
    
    # Check for discordance between trees
    discordant = _update_concordance(alt,ref_bipartition_dict)
        
    # Run MAF algorithm
    forest = dendropy.TreeList() # forest to hold subtrees in MAF
    counter = 0
    while discordant:
                    
        # Find first discordant edge
        for edge in alt.postorder_edge_iter():
            if not edge.concordant: break
            
        parent = edge.head_node #first discordant sibling pair (edge) in t1
        edge_cut_list = parent.child_edges()
        
        # Update edge bipartition dicts
        ref_bipartition_dict = ref.split_bitmask_edge_map
        alt_bipartition_dict = alt.split_bitmask_edge_map
        
        # Get sibling of A and C in ALT
        edge_A = edge_cut_list[0]
        edge_C = edge_cut_list[1]
        
        # Get sibling edge B
        ref_edge_A = ref_bipartition_dict.get(edge_A.leafset_bitmask,None)
        parent_edge_A = ref_edge_A.tail_node
        child_edges = parent_edge_A.child_edges()
        if child_edges[0] is ref_edge_A:
            sibling_edge = child_edges[1]
        else:
            sibling_edge = child_edges[0]
        alt_edge_B = alt_bipartition_dict.get(sibling_edge.leafset_bitmask,None)
        if alt_edge_B:
            edge_cut_list.append(alt_edge_B)
            
        if not alt_edge_B:
            # Get sibling edge D in place of B
            ref_edge_C = ref_bipartition_dict.get(edge_C.leafset_bitmask,None)
            parent_edge_C = ref_edge_C.tail_node
            child_edges = parent_edge_C.child_edges()
            if child_edges[0] is ref_edge_C:
                sibling_edge = child_edges[1]
            else:
                sibling_edge = child_edges[0]
            alt_edge_D = alt_bipartition_dict.get(sibling_edge.leafset_bitmask,None)
            if alt_edge_D:
                edge_cut_list.append(alt_edge_D)
        
        # Cut all three edges in edge_cut_list
        for edge in edge_cut_list:
            
            alt_bipartition_dict = alt.split_bitmask_edge_map
            ref_bipartition_dict = ref.split_bitmask_edge_map
    
            bitmask = edge.leafset_bitmask
            alt_cut_edge = alt_bipartition_dict.get(bitmask,None)
            ref_cut_edge = ref_bipartition_dict.get(bitmask,None)
            
            # Make sure we're getting the correct edges
            #edge_taxa = set([lf.taxon.label for lf in edge.head_node.leaf_iter()])
            #alt_cut_edge_taxa = set([lf.taxon.label for lf in alt_cut_edge.head_node.leaf_iter()])
            #ref_cut_edge_taxa = set([lf.taxon.label for lf in ref_cut_edge.head_node.leaf_iter()])

            # Extract subtree from alt and place in subtree forest
            child_nd = alt_cut_edge.head_node
            component = TreeOps.extract_subtree(alt,child_nd)
            forest.append(component)
        
            # Prune subtree from alt_cut
            alt.prune_subtree(child_nd, update_bipartitions=True, suppress_unifurcations=True)
        
            # Find and prune corresponding edge in ref_cut tree
            ref_cut_child = ref_cut_edge.head_node
            ref.prune_subtree(ref_cut_child, update_bipartitions=True, suppress_unifurcations=True)
        
        counter += 1
        if plot:
            print("-" * 20)
            print("MAF iteration: " + str(counter))
            print("Reference tree:")
            ref.print_plot()
            print()
            print("Alternate tree:")
            alt.print_plot()
            print()
            print("Cut tree:")
            if len(component.nodes()) > 1:
                component.print_plot()
            else:
                print('---' + component.nodes()[0].taxon.label)
            print()
            print("-" * 20)
            
        # Make sure discordant is updated
        ref_bipartition_dict = ref.split_bitmask_edge_map 
        discordant = _update_concordance(alt,ref_bipartition_dict)
    
    # Remove dummy node at deep root
    #if deep_root:
        #alt = remove_deep_root(alt)
    
    forest.append(alt) # append what's left of alt (should at least be a rooted tree with two edges)
    return forest

def _update_concordance(alt,ref_bipartition_dict):
    
    """
        
        Determine if each edge in ALT is (dis)concordant based on  edge bipartitions present in REF
        Returns discordant = True if any edge in ALT is discordant
    
        Parameters:     
           alt (dendropy.Tree): input ALT tree (T1)
           ref_bipartition_dict (dict): dict mapping edge/split bitmasks (integers) to their corresponding edge in REF (T2) 
           
        Returns:     
           discordant (boolean): True if any edge in ALT/REF is discordant
    """
    
    # Determine if each edge and its child node is concordant/discordant
    discordant = False
    for edge in alt.postorder_edge_iter():
        bitmask = edge.leafset_bitmask
        ref_edge = ref_bipartition_dict.get(bitmask,None)
        if ref_edge: # if edge bipartition exists in reference
            edge.concordant = True
            edge.head_node.concordant = True # child node of edge
        else:
            edge.concordant = False
            edge.head_node.concordant = False
            discordant = True
            
    return discordant

def _count_discordant_edges(alt,ref_bipartition_dict):
    
    """
        
        Count discordant edges between ALT and REF trees
        Assumes concordant status is known from _update_concordance()
    
        Parameters:     
           alt (dendropy.Tree): input ALT tree (T1)
           ref_bipartition_dict (dict): dict mapping edge/split bitmasks (integers) to their corresponding edge in REF (T2) 
           
        Returns:     
           count (int): Number of discordant edges between trees
    """
    
    count = 0
    for edge in alt.postorder_edge_iter():
        bitmask = edge.leafset_bitmask
        ref_edge = ref_bipartition_dict.get(bitmask,None)
        if not ref_edge: # if edge bipartition does not exist in reference 
            count += 1
            
    return count

def get_maf(ref,alt):

    """
    
        Short-cut alias for get_maf_4cut()
        Computes a MAF between two trees using the 4-cut algorithm.
    
        Parameters:     
           ref (dendropy.Tree): input ref tree (T1)
           alt (dendropy.Tree): input alt tree (T2) 
           
        Returns:     
           forest (dendropy.TreeList): the MAF as a list of trees
    """
    
    forest = get_maf_4cut(ref,alt)
    
    return forest

def get_spr_dist(ref,alt):
    
    """
    
        Compute subtree prune and regraft (SPR) distance between ref and alt trees using the 4-cut MAF algorithm.
    
        Parameters:     
           ref (dendropy.Tree): input ref tree (T1)
           alt (dendropy.Tree): input alt tree (T2) 
           
        Returns:     
           spr_dist (int): SPR distance
    """
    
    if ref.taxon_namespace is not alt.taxon_namespace:
        alt.taxon_namespace = ref.taxon_namespace
        alt.reindex_subcomponent_taxa()
    maf = get_maf_4cut(ref,alt)
    spr_dist = len(maf) - 1
    
    return spr_dist


def test_discordance(ref,alt):
    
    """
    
        Wrapper function to test discordance between a pair of trees REF and ALT
    
        Parameters:     
           ref (dendropy.Tree): input ref tree (T1)
           alt (dendropy.Tree): input alt tree (T2) 
           
        Returns:     
           discordant (boolean): True if input trees are discordant
    """
    
    ref.encode_bipartitions() # need to encode bipartitions before can access them    
    alt.encode_bipartitions()
    ref_bipartition_dict = ref.split_bitmask_edge_map #dictionary mapping split bitmasks (integers) to their corresponding edges through the Tree.split_bitmask_edge_map attribute
    discordant = _update_concordance(alt,ref_bipartition_dict)
    return discordant

def plot_maf(maf):
        
    print("-" * 25)
    print("Maximum agreement forest:")
    print("-" * 25)
    for tr in maf:
        print()
        if len(tr.nodes()) > 1:
            tr.print_plot()
        else:
            print('---' + tr.nodes()[0].taxon.label)
            
