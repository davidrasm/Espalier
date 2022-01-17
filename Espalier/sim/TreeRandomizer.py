#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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

from Espalier.TreeOps import extract_subtree
import dendropy
import msprime
import random

def random_SPR(tree,exclude_root_children=False,plot=False):
    
    """
        Randomize topology of tree by random subtree prune regraft moves
    """
    
    tr = tree.clone(depth=2) #Deep copy tree
    tr.calc_node_ages(ultrametricity_precision=False)
    edge_list = tr.edges()
    root_edge = tr.seed_node.edge
        
    # Check if child edge of root are terminal branches - exclude them if they are
    root_child_edges = root_edge.head_node.child_edges()
    forbidden = [root_edge]
    if exclude_root_children:
        forbidden.append(root_child_edges[0])
        forbidden.append(root_child_edges[1])
    else:
        if root_child_edges[0].head_node.is_leaf():
            forbidden.append(root_child_edges[1])
        if root_child_edges[1].head_node.is_leaf():
            forbidden.append(root_child_edges[0])
    
    # Pick random edge for subtree to prune
    cut_edge = random.choice(edge_list)
    while cut_edge in forbidden:
        print("Cut edge is root")
        cut_edge = random.choice(edge_list)
    subtree_node = cut_edge.head_node
    parent = cut_edge.tail_node
    child_edges = parent.child_edges()
    if child_edges[0] is cut_edge:
        sister_edge = child_edges[1]
    else:
        sister_edge = child_edges[0]
    
    # Extract and prune subtree
    subtree = extract_subtree(tr,subtree_node)
    tr.prune_subtree(subtree_node, update_bipartitions=False, suppress_unifurcations=True)
    if plot:
        if len(subtree.nodes()) > 1:
            subtree.print_plot()
        else:
            print('---' + subtree.nodes()[0].taxon.label)
        tr.print_plot()
    
    # Pick random edge to regraft subtree
    edge_list = tr.edges()
    root_edge = tr.seed_node.edge
    graft_edge = random.choice(edge_list)
    forbidden = [root_edge,sister_edge] # Graft edge cannot be root or sister to cut edge
    constraint_flag = True
    while constraint_flag:
        
        constraint_flag = False
        
        #print("Attachment edge is the same or sister to cut edge")
        graft_edge = random.choice(edge_list)
        if graft_edge in forbidden:
            constraint_flag = True
        
        # Check height constraints
        if not graft_edge.tail_node:
            attachment_parent_height = tr.seed_node.age
        else:
            attachment_parent_height = graft_edge.tail_node.age
        subtree_height = subtree_node.age
        if attachment_parent_height < subtree_height:
            constraint_flag = True
    
    # Regraft subtree by creating a new 'graft' node to attach sub_tree and sister of sub_tree
    attachment_parent_node = graft_edge.tail_node
    #if not attachment_parent_node: # this should never happen
    #    parent_leaf_set = set([lf.taxon.label for lf in graft_edge.head_node.leaf_iter()])
    #    print(parent_leaf_set)
    attachment_sister_node = graft_edge.head_node
    attachment_edge_length = attachment_sister_node.edge_length
    graft_node = attachment_parent_node.new_child()
    
    # Make sure graft node height is above grafted subtree and sister height!!!
    attachment_parent_height = attachment_parent_node.age
    subtree_height = subtree_node.age
    sister_subtree_height = attachment_sister_node.age
    max_height = attachment_parent_height
    min_height = max(sister_subtree_height,subtree_height)
    
    if attachment_parent_height < subtree_height:
        print("WARNING: Parent height is below subtree height")
    
    # Prune sister subtree descending from attachment edge
    tr.prune_subtree(attachment_sister_node, update_bipartitions=False, suppress_unifurcations=True)
    
    # Get and set length of graft edge - currently using midpoint because subtree could attach anywhere
    graft_node_age = min_height + (max_height - min_height)/2.0
    graft_node.edge.length = attachment_parent_height - graft_node_age
    
    # Reattach sister subtree to graft node
    attachment_sister_node.edge.length = graft_node_age - sister_subtree_height
    graft_node.add_child(attachment_sister_node)
    
    # Edge length is difference between graft node age and height of subtree
    subtree_node.edge.length = graft_node_age - subtree_height
    graft_node.add_child(subtree_node)
        
    if plot:
        print("New SPR tree")
        tr.print_plot()
    
    return tr

if  __name__ == '__main__':

    # Simulate tree without recombination"
    sim = True
    if sim:
        ts = msprime.simulate(sample_size=6, Ne=100, length=1e4, recombination_rate=0, record_full_arg=True)
        tree = ts.first()
        print("-" * 20)
        print(tree.draw(format="unicode"))
        with open("maf_randomSPRs_refTree_test1.tre", "w") as text_file:
            print(tree.newick(), file=text_file)
    
    ref_file = "maf_randomSPRs_refTree_test1.tre"    
    taxa = dendropy.TaxonNamespace()
    ref = dendropy.Tree.get(file=open(ref_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
    
    
    # Perform random SPR moves on tree
    moves = 2 # number of SPR moves
    alt = ref
    for i in range(moves):
        print("SPR move = ", str(i))
        alt = random_SPR(alt,plot=True)
    
    print("Reference tree:")
    ref.print_plot()
    
    print("Alternate tree:")
    alt.print_plot()
    #alt.write(path='maf_randomSPRs_altTree_test1.tre',schema='newick',suppress_annotations=True)
    
    
