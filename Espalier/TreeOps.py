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

import dendropy
import logging

"""
    Functions for operating on tree topologies including SPR moves
"""

def extract_subtree(tree,node):
    
    """
    
        Extract subtree starting at 'node' from tree.
        Resulting subtree will be rooted on 'node'.
    
        Parameters:     
           tree (dendropy.Tree): input tree
           node (dendropy.Node): node to extract at
           
        Returns:     
           subtree (dendropy.TreeList): extracted subtree
    """
    
    for nd in tree.preorder_node_iter():
        nd.remove = False
    for nd in node.preorder_iter():
        nd.remove = True
    node_filter_fn = lambda nd: nd.remove
    subtree = tree.extract_tree(node_filter_fn=node_filter_fn)
    
    return subtree


def regraft_subtree(subtree,graft_node,graft_node_age,min_edge_length=0.0):
    
    """
        
        Regraft subtree to tree at graft_node.
        Length of edge connecting subtree to graft_node is determined by graft_node_ag.e
    
        Parameters:     
           subtree (dendropy.Tree): subtree to regraft
           graft_node (dendropy.Node): node to regraft to
           graft_node_age (float): age (height) of graft node
        
        Optional keyword arguements:
            min_edge_lenth (float): minimum allowable edge/branch length
        
        Returns:     
           None
    """
    
    subtree.calc_node_ages(ultrametricity_precision=False)
    subtree_root = subtree.seed_node
    subtree_height = subtree_root.age
        
    # Edge length is difference between mrca age and height of subtree
    subtree_root_length = graft_node_age - subtree_height #subtree_root.age
    if subtree_root_length >= 0.0:
        subtree_root.edge.length = subtree_root_length
    else:
        logging.debug('Subtree edge length is less than zero!!')
        subtree_root.edge.length = min_edge_length
    graft_node.add_child(subtree_root)

def add_deep_root(tree,reroot=True):
    
    """
    
        Add deep root above current root by adding dummy node as new root.
    
        Parameters:     
           tree (dendropy.Tree): input tree
           
        Returns:     
            new_tree (dendropy.Node): deep rooted tree
    """
    
    new_tree = dendropy.Tree(taxon_namespace=tree.taxon_namespace)
    dummy = new_tree.seed_node.new_child() # dummy node above new root
    dummy.edge.length = 1
    curr_root = tree.seed_node # get_root_node(tree)
    new_tree.seed_node.add_child(curr_root)
    
    if reroot:
        new_root = new_tree.seed_node #get_root_node(new_tree)
        new_tree.reroot_at_node(new_root, update_bipartitions=False)    
    
    # Add taxon label for dummy node
    new_tree.taxon_namespace.add_taxon(dendropy.Taxon("dummy"))
    dummy.taxon = new_tree.taxon_namespace.get_taxon("dummy")
    
    return new_tree

def remove_deep_root(tree):
    
    """
    
        Remove deep root above current root by adding dummy node as new root
    
        Parameters:     
           tree (dendropy.Tree): input tree
           
        Returns:     
            tree (dendropy.Node): tree with deep root removed
    """
    
    tree.prune_taxa_with_labels(["dummy"])
    
    return tree


def deep_root_pair(ref,alt):
    
    """
    
        Deep root a pair of trees preserving taxon namespace by writing and reading in new newick trees
    
        Parameters:     
           ref (dendropy.Tree): first input tree T1
           alt (dendropy.Tree): second input tree T2
           
        Returns:     
            ref_deep_rooted (dendropy.Tree): ref tree with deep root
            alt_deep_rooted (dendropy.Tree): ref tree with deep root
    """
    
    ref_deep_rooted = add_deep_root(ref)
    alt_deep_rooted = add_deep_root(alt)
    
    # Write out and read in newick trees to update taxa namespaces
    temp_ref_file = 'temp_ref.tre'
    temp_alt_file = 'temp_alt.tre'
    ref_deep_rooted.write(path=temp_ref_file,schema='newick',suppress_annotations=True)
    alt_deep_rooted.write(path=temp_alt_file,schema='newick',suppress_annotations=True)
    taxa = dendropy.TaxonNamespace()
    ref_deep_rooted = dendropy.Tree.get(file=open(temp_ref_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
    alt_deep_rooted = dendropy.Tree.get(file=open(temp_alt_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)

    return ref_deep_rooted, alt_deep_rooted

#################################################
    # Deprecated functions slated for removal
#################################################

# def add_deep_root_with_ages(tree,reroot=True):

#     "Get age of curr root in tree"
#     curr_root = get_root_node(tree)
#     tree.calc_node_ages(ultrametricity_precision=False)
#     curr_root_age = curr_root.age
    
#     "Add deep root above current root by adding dummy node at new root"
#     new_tree = dendropy.Tree(taxon_namespace=tree.taxon_namespace)
#     dummy = new_tree.seed_node.new_child() # dummy node above new root
#     dummy.edge.length = curr_root_age
    
#     new_tree.seed_node.add_child(curr_root)
#     curr_root.edge.length = 0.0
    
#     if reroot:
#         new_root = get_root_node(new_tree)
#         new_tree.reroot_at_node(new_root, update_bipartitions=False)    
    
#     "Add taxon label for dummy node"
#     new_tree.taxon_namespace.add_taxon(dendropy.Taxon("dummy"))
#     dummy.taxon = new_tree.taxon_namespace.get_taxon("dummy")
    
#     return new_tree

# def get_root_node(tree):
    
#     """
#         TODO: replace all instances with tree.seed_node
#         Get root node of tree
    
#         Parameters:     
#            tree (dendropy.Tree): input tree
           
#         Returns:     
#             root (dendropy.Node): root node
#     """
    
#     for nd in tree.preorder_node_iter():
#         root = nd
#         break
    
#     seed_node = tree.seed_node
    
#     assert seed_node == root
    
#     return root

# def get_root_edge(tree):
    
#     """
    
#         Get root edge of tree
    
#         Parameters:     
#            tree (dendropy.Tree): input tree
           
#         Returns:     
#             root (dendropy.Node): root edge
#     """
    
#     for ed in tree.preorder_edge_iter():
#         root = ed
#         break
#     return root


