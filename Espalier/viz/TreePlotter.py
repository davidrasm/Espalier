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

import dendropy

def plot_tree_path(tree_path,label=None):
    
    """
        
        Plot trees in tree_path as ASCII text plots
        
        Parameters:     
            tree_path (list[dendropy.Tree]): list of trees in path
            label (str): label type for nodes in tree
              
    """
    
    if label:
        if label == 'age':
            label_fn = lambda nd: nd.taxon.label if nd.is_leaf() else str(f"{nd.age:0.2f}")
        elif label == 'bitmask':
            label_fn = lambda nd: str(nd.edge.split_bitmask)
        else:
            raise Exception("%s is not an accepted node labeling function", label)
        
    for loc, tree in enumerate(tree_path):
        print("-" * 20)
        print("Tree #" + str(loc))
        print("-" * 20)
        tree.print_plot(show_internal_node_labels=True,node_label_compose_fn=label_fn)