#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 14:52:10 2020

11-16-2021: add_recombination_events was modified to test whether we need to
include sibling child node constraints when adding recombination times 

@author: david
"""
import tskit
import dendropy
import numpy as np
import TreeReconciler
import ARGSimulator
from operator import attrgetter
from collections import namedtuple
import os
import shutil

def add_rec_node(tree,attachment_edge,recomb_time,midpoint=False):
    
    "Create a new 'graft' node to attach sub_tree and sister of sub_tree"
    attachment_parent_node = attachment_edge.tail_node
    attachment_child_node = attachment_edge.head_node
    deep_rooted = False
    if not attachment_parent_node:
        
        print('Grafting recombinant node above root')
        
        "Grafting above root of rec tree such that graft node will be new root"
        "This can happen if cut one of the root's child edges"
        "This should be working now but not very well tested!!"
        tree = TreeReconciler.add_deep_root(tree) # reroot == True by default now
        deep_rooted = True 
        attachment_parent_node = attachment_edge.tail_node
        attachment_child_node = attachment_edge.head_node
        
    else:
        
        attachment_edge_length = attachment_child_node.edge_length
        graft_node = attachment_parent_node.new_child() # this will represent the recombination node
    
        "Extract child subtree descending from attachment edge and prune it from tree"
        for node in tree.preorder_node_iter():
            node.remove = False
        for node in attachment_child_node.preorder_iter():
            node.remove = True
        node_filter_fn = lambda nd: nd.remove
        child_subtree = tree.extract_tree(node_filter_fn=node_filter_fn)
        tree.prune_subtree(attachment_child_node, update_bipartitions=False, suppress_unifurcations=False) # supress_unifurcations was True, which was why we were losing recombinant nodes we already added

        "Get and set length of graft edge"
        if midpoint: # attach to midpoint of attachment edge
            graft_node.edge.length = attachment_edge_length / 2.0
        else: # figure out length of edge from attachment height
            tree.calc_node_ages(ultrametricity_precision=False)
            "This sets the height of the recombination node/event below the parent node"
            graft_node_edge_length = attachment_parent_node.age - recomb_time
            if graft_node_edge_length >= 0.0:
                graft_node.edge.length = graft_node_edge_length
            else:
                print('WARNING: Graft edge length is less than zero!!')
                graft_node.edge.length = 0.0
    
        "Reattach child subtree to graft node"
        child_subtree_root = TreeReconciler.get_root_node(child_subtree)
        if midpoint:
            child_subtree_root.edge.length = attachment_edge_length / 2.0
        else:
            "Edge length is difference between graft node age and height/age of child_subtree"
            child_subtree.calc_node_ages(ultrametricity_precision=False)
            child_subtree_edge_length = recomb_time - child_subtree_root.age
            if child_subtree_edge_length >= 0.0:
                child_subtree_root.edge.length = child_subtree_edge_length
            else:
                print('WARNING: Child subtree edge length is less than zero!!')
                child_subtree_root.edge.length = 0.0
        graft_node.add_child(child_subtree_root)
        
        if deep_rooted:
             tree = TreeReconciler.remove_deep_root(tree)
    
    return tree

def add_path_rec_nodes(tree_path,rec_nodes):
    
    for tree_idx, tree in enumerate(tree_path):
        bipartition_dict = tree.split_bitmask_edge_map
        
        "Sort rec_nodes in descending order by time so we add older events first"
        tree_rec_nodes = rec_nodes[tree_idx]
        tree_rec_nodes = sorted(tree_rec_nodes, key=attrgetter('recomb_time'), reverse=True)

        for rn in tree_rec_nodes:
            rec_edge = bipartition_dict.get(rn.edge_bitmask,None)
            tree = add_rec_node(tree,rec_edge,rn.recomb_time)
            
    return tree_path

"Find edge where subtree attaches in alt and then find equivalent edge in ref"
def find_recombinant_edge(tree,subtree_taxa):
    
    "Find edge where subtree attaches in alt tree" 
    "Note: this will be the parent of the first edge that includes all subtree_taxa in its leaf set"
    for edge in tree.postorder_edge_iter():
        edge_taxa = set([lf.taxon.label for lf in edge.head_node.leaf_iter()])
        if subtree_taxa.issubset(edge_taxa): # subtree_taxa are subset of all edge_taxa
            "All subtree_taxa are in edge_taxa so subtree must attach here!"
            break
    
    attachment_edge = edge
    parent_node = edge.tail_node # parent of recombinant/attachment edge
    child_node = edge.head_node # child of recombinant edge
    child_nodes = parent_node.child_nodes()
    if len(child_nodes) < 2:
        print('WTF!') # I think this can only happen if we've already added a recombination node as a unifurcation
    if child_nodes[0] is child_node:
        sibling_node = child_nodes[1]
    else:
        sibling_node = child_nodes[0]
    
    "Find constraints on timing of recomb event"
    tree.calc_node_ages(ultrametricity_precision=False)
    parent_time = parent_node.age # max recombination time is height/age of parent node
    child_time = child_node.age
    sibling_time = sibling_node.age

    return attachment_edge, parent_time, child_time, sibling_time

"""
Add necessary recombination events to reconcile two trees based on MAF
"""
def add_recombination_events(ref,alt,maf,ts,loc,path_recomb_edges,ref_rec_nodes,alt_rec_nodes,plot=False):
    
    "Deep copy original trees so we don't lose them"
    #ref = ref_in.clone(depth=2)
    #alt = alt_in.clone(depth=2)
    
    RecNode = namedtuple('RecNode', 'edge_bitmask recomb_time')
    
    maf.pop() # remove last/largest connected component
    for sub_tree in reversed(maf):
    
        "Get taxa set in subtree"
        subtree_taxa = set([lf.taxon.label for lf in sub_tree.leaf_node_iter()])
            
        "Find edge where subtree attaches in alt and then find equivalent edge in ref"
        attachment_edge_alt, parent_time_alt, child_time_alt, sibling_time_alt = find_recombinant_edge(alt,subtree_taxa)
        
        "Find edge where subtree attaches in ref and then find equivalent edge in alt"
        attachment_edge_ref, parent_time_ref, child_time_ref, sibling_time_ref  = find_recombinant_edge(ref,subtree_taxa)
        
        "Check if we've picked the right recombinant edge"
        chosen_recomb_edge = attachment_edge_alt.split_bitmask
        #ts_tree = ts.at_index(loc) # oddly can't access added instance variables of tree objects this way
        true_recomb_edge = path_recomb_edges[loc][0]
        if chosen_recomb_edge != true_recomb_edge:
            print("WARNING: chosen recombinant edge is not true recombinant edge!!")
        #print("Chosen recombinant edge: " + str(chosen_recomb_edge))
        #print("True recombinant edge: " + str(true_recomb_edge))
        
        
        "Pick a recombination event time that satisfies time constraints in both alt and ref trees"
        max_recomb_time = min(parent_time_alt,parent_time_ref) # event must be below parent nodes in either tree
        
        min_recomb_times = [child_time_alt,child_time_ref]
        
        """
            Commenting out lines below removes sibling time constraints
        """
        #if sibling_time_alt < max_recomb_time:
        #    min_recomb_times.append(sibling_time_alt)
        #if sibling_time_ref < max_recomb_time:
        #    min_recomb_times.append(sibling_time_ref)
        min_recomb_time = max(min_recomb_times)        
        
        #min_recomb_time = max(min_recomb_time_ref,min_recomb_time_alt) # event must be above parent's children in either tree
        #assert max_recomb_time >= min_recomb_time
        if max_recomb_time <= min_recomb_time:
            for tr_num, tree in enumerate(ts.trees()):
                if tr_num == loc or tr_num == loc+1:
                    print("-" * 20)
                    print("tree {}: ".format(tr_num))
                    print(tree.draw(format="unicode"))
            print("Max time constraint: " + str(max_recomb_time))
            print("Min time constraint: " + str(min_recomb_time))
            assert max_recomb_time >= min_recomb_time
        recomb_time = np.random.uniform(low=min_recomb_time,high=max_recomb_time) # draw a uniform time between these constraints 
        
        """
            Add single node as a unification representing recombination to both alt and ref tree        
            We will now add rec_nodes last so bifurcations and edge_bitmasks
            We will use a named tuple so we can easilty sort these by times 
        """
        #ref = add_rec_node(ref,attachment_edge_ref,recomb_time)
        #alt = add_rec_node(alt,attachment_edge_alt,recomb_time)
        ref_rec_nodes.append(RecNode(edge_bitmask=attachment_edge_ref.split_bitmask, recomb_time=recomb_time))
        alt_rec_nodes.append(RecNode(edge_bitmask=attachment_edge_alt.split_bitmask, recomb_time=recomb_time))
        
        if plot:
            print("Reference tree:")
            ref.print_plot()
            print("Alternate tree:")
            alt.print_plot()
        
    #return ref, alt
    return ref_rec_nodes,alt_rec_nodes


def relabel_taxa(tree_file):
    
    taxa = dendropy.TaxonNamespace()
    tree = dendropy.Tree.get(file=open(tree_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
    for lf in tree.leaf_node_iter():
        lf.taxon.label = str(int(lf.taxon.label) - 1) # to make consistent with taxon lables from msprime
    tree.write(path=tree_file,schema='newick',suppress_annotations=True)
    
def add_edge_bitmasks(dendro_tree, ts_tree, ts):
        
    """
        For each edge in ts table find equivelent edge in dendropy tree using the edge's taxon leaf set
    """
    ts_edge_bitmasks = []
    dendro_tree.encode_bipartitions(suppress_unifurcations=False)
    for edge in ts.tables.edges:
        u = edge.child
        leaf_set = set([str(v) for v in ts_tree.leaves(u)])
        equiv_edge = None
        for edge in dendro_tree.postorder_edge_iter():
            edge_taxa = set([lf.taxon.label for lf in edge.head_node.leaf_iter()])
            if leaf_set==edge_taxa: # subtree_taxa are subset of all edge_taxa
                equiv_edge = edge
                break
        
        if equiv_edge:
            ts_edge_bitmasks.append(edge.split_bitmask)
        else:
            ts_edge_bitmasks.append(None)
    
    ts_tree.edge_bitmasks = ts_edge_bitmasks # maybe attach to individual trees in ts so we can iterate through ts without indexing a seperate list

"""
    Each tree can have multiple recombination nodes/events
    Here, we try to find the recombinant edge associated with a particular breakpoint
    This breakpoint is assumed to be the right pos of each tree's inteval
"""
def add_recombinant_ts_edges(tree,ts):
    
    edge_parents = np.array(ts.tables.edges.parent,ndmin=1)
    edge_ids = np.arange(len(edge_parents))  
    left_pos = tree.interval[0]
    right_pos = tree.interval[1]

    recomb_edges = []
    for node_id, node in enumerate(ts.tables.nodes):
        if node.flags == 131072: # if recombination flag
            "Find row(s) in edge table"
            for edge_idx in edge_ids[edge_parents == node_id]:
                edge_left = ts.tables.edges[edge_idx].left
                edge_right = ts.tables.edges[edge_idx].right
                if left_pos >= edge_left and right_pos == edge_right: 
                    recomb_edges.append(tree.edge_bitmasks[edge_idx])

    #tree.recomb_edges = recomb_edges
    return recomb_edges

"Simulate a tree series in ms prime"
sim = True
plot = False

"Sim params"
sample_size = 10 
Ne = 1
genome_length = 1e3
rho = 2e-4

path = 'recomb_placement_path_noSiblingConstraints'
if not os.path.isdir(path):
    os.mkdir(path)
    
sim_test_num = 1
sim_path = path + '/sim' + str(sim_test_num).rjust(3, '0') + '/'

if sim:
    "Create sim dir or delete old one if it exists"
    if os.path.isdir(sim_path):
        shutil.rmtree(sim_path)
        os.mkdir(sim_path)
    else:
        os.mkdir(sim_path)
    ts = ARGSimulator.sim_ARG(sample_size=sample_size,Ne=Ne,length=genome_length,recombination_rate=rho,min_breakpoints=2)
    ts.dump(sim_path + 'recomb_placement_tables')
else:
    ts = tskit.load(sim_path + 'recomb_placement_tables')

breaks = ts.breakpoints(as_array=True)
segments = len(breaks) - 1 # number of non-recombinant segments between breakpoints

"Print the simulated ts tables"
#if plot:
print("the simulated ts tables")
print(ts.tables.nodes)
print(ts.tables.edges)
for tree in ts.trees():
    print("-" * 20)
    print("tree {}: interval = {}".format(tree.index, tree.interval))
    print(tree.draw(format="unicode"))

sim = 1
tree_files = [sim_path + "tree" + str(i) + ".tre" for i in range(segments)]
tree_path = []
taxa = dendropy.TaxonNamespace()
for tr_num, tree in enumerate(ts.trees()):
    tree_file = tree_files[tr_num]
    with open(tree_file, "w") as text_file:
        print(tree.newick(), file=text_file)
    relabel_taxa(tree_file)
    dendro_tree = dendropy.Tree.get(file=open(tree_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
    dendro_tree.suppress_unifurcations(update_bipartitions=True) # If we don't suppress unifurfactions can get recomb positions from dendropy nodes
    tree_path.append(dendro_tree)

"Record true recombinant edges in ts trees to compare with the ones we pick"
path_recomb_edges = []
for tr_num, tree in enumerate(ts.trees()):
    #print("Tree number: " + str(tr_num))
    "Get edge mapping between tskit and dendro tree"
    add_edge_bitmasks(tree_path[tr_num], tree, ts)
    recomb_edges = add_recombinant_ts_edges(tree, ts)
    path_recomb_edges.append(recomb_edges)

"Create list to hold rec_nodes for each tree in path"
rec_nodes = [[] for x in range(segments)]

"Add recombination events to each tree in path"
print("Adding recombination events to trees in path")
for loc in range(segments-1):
    ref = tree_path[loc]
    alt = tree_path[loc+1]
    ref_rec_nodes = rec_nodes[loc]
    alt_rec_nodes = rec_nodes[loc+1]
    maf = TreeReconciler.get_maf_4cut(ref,alt,plot=False)
    print("Adding " + str(len(maf)-1) + " recombination events between segments " + str(loc) + " and " + str(loc+1))
    ref_rec_nodes,alt_rec_nodes = add_recombination_events(ref,alt,maf,ts,loc,path_recomb_edges,ref_rec_nodes,alt_rec_nodes,plot=False)
    
"Add rec nodes to trees in path"
tree_path = add_path_rec_nodes(tree_path,rec_nodes)

"Write tree in path to newick files"
write_tree_path = True
if write_tree_path:
    for loc, tree in enumerate(tree_path):
        tree_file = tree_files[loc].replace('tree','PathTree')
        tree.write(path=tree_file,schema='newick',suppress_annotations=True,suppress_rooting=True)

