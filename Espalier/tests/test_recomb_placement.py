#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 14:52:10 2020

@author: david
"""
import msprime
import tskit
import dendropy
import numpy as np
import TreeReconciler

def sim_ts(desired_breakpoints):
    
    "Simulate tree series with exact num of breakpoints and discordant trees"
    discordant = False
    breaks = 0
    while breaks != desired_breakpoints + 2 or not discordant:
        newicks = []
        ts = msprime.simulate(sample_size=10, Ne=100, length=1e3, recombination_rate=5e-6, record_full_arg=True)
        breaks = len(ts.breakpoints(as_array=True))
        prev_tree = None
        for tr_num, tree in enumerate(ts.trees()):
            if tr_num > 0:
                "Test for discordance"
                taxa = dendropy.TaxonNamespace()
                ref = dendropy.Tree.get(data=prev_tree, schema="newick", rooting="default-rooted", taxon_namespace=taxa)
                ref.suppress_unifurcations(update_bipartitions=True) # get rid of all unifications
                alt = dendropy.Tree.get(data=tree.newick(), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
                alt.suppress_unifurcations(update_bipartitions=True) # get rid of all unifications
                discordant =  TreeReconciler.test_discordance(ref,alt)
            prev_tree = tree.newick()
            newicks.append(tree.newick())  
        print("Breaks: ", breaks-2, "Discordant: ", discordant)
          
    "Write tree files"
    tree_file_A = "recomb_placement_test1_treeA.tre"
    with open(tree_file_A, "w") as text_file:
        print(newicks[0], file=text_file)
    tree_file_B = "recomb_placement_test1_treeB.tre"
    with open(tree_file_B, "w") as text_file:
            print(newicks[1], file=text_file)
            
    return ts

"Find edge where subtree attaches in alt and then find equivalent edge in ref"
def find_recombinant_edge_defunct(alt_tree,ref_tree,subtree_taxa):
    
    "Get bipart_index_leaf_dict - maps bipartition indexes (keys) to taxon lables (values)"
    #bipart_index_leaf_dict = TreeReconciler.get_bipart_indexes(alt_tree)
    #for key, value in bipart_index_leaf_dict.items(): print(key, value)
    
    "Find edge where subtree attaches in alt tree" 
    "Note: this will be the parent of the first edge that includes all subtree_taxa in its leaf set"
    for edge in alt_tree.postorder_edge_iter():
        
        "Old way"
        #edge_leafset_bitstring = edge.bipartition.leafset_as_bitstring()
        #edge_taxa = set([bipart_index_leaf_dict[pos] for pos, char in enumerate(edge_leafset_bitstring) if char == '1'])
        
        "New way"
        edge_taxa = set([lf.taxon.label for lf in edge.head_node.leaf_iter()])
        
        "Could also try this:"
        #edge_taxa.is_leafset_nested_within(subree_taxa.leafset)
        
        if subtree_taxa.issubset(edge_taxa): # subtree_taxa are subset of all edge_taxa
            "All subtree_taxa are in edge_taxa so subtree must attach here!"
            break
    
    parent_node = edge.tail_node # parent of attachment edge    
    #attachment_edge_leafset = parent_node.leafset_as_bitstring()
    #attachment_edge_taxa_set = set([bipart_index_leaf_dict[pos] for pos, char in enumerate(attachment_edge_leafset) if char == '1'])
    attachment_edge_taxa_set = set([lf.taxon.label for lf in parent_node.leaf_iter()])
    
    #"Find height of attachment"
    #alt_tree.calc_node_ages(ultrametricity_precision=False)
    #attachment_height = edge.tail_node.age # should this be the age of the head or tail node?
    
    "Find constraints on timing of recomb event"
    alt_tree.calc_node_ages(ultrametricity_precision=False)
    max_recomb_time = parent_node.age # max recombination time is height/age of parent node
    child_nodes = parent_node.child_nodes()
    child0_age = child_nodes[0].age
    child1_age = child_nodes[1].age
    min_recomb_time = max(child0_age,child1_age) # min recombination time is max of parent nodes children
    #time_constraints = (max_recomb_time, min_recomb_time)
    
    "Find edge equivalent to the attachment edge in ref_tree"
    "Note: this edge will have all taxa in attachment_edge_taxa_set that ARE NOT in subtree"
    taxa_set_intersection = attachment_edge_taxa_set - subtree_taxa
    for edge in ref_tree.postorder_edge_iter():
        edge_taxa = set([lf.taxon.label for lf in edge.head_node.leaf_iter()])
        if not taxa_set_intersection.difference(edge_taxa): 
            "If there is no difference:"
            break
    attachment_edge = edge # edge subtree will be attached to in rec_tree
    
    return attachment_edge, max_recomb_time, min_recomb_time

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
        tree.prune_subtree(attachment_child_node, update_bipartitions=False, suppress_unifurcations=True)

        "Get and set length of graft edge"
        if midpoint: # attach to midpoint of attachment edge
            graft_node.edge.length = attachment_edge_length / 2.0
        else: # figure out length of edge from attachment height
            tree.calc_node_ages(ultrametricity_precision=False)
            "This sets the height of the recombination node/event below the parent node"
            graft_node.edge.length = attachment_parent_node.age - recomb_time
    
        "Reattach child subtree to graft node"
        child_subtree_root = TreeReconciler.get_root_node(child_subtree)
        if midpoint:
            child_subtree_root.edge.length = attachment_edge_length / 2.0
        else:
            "Edge length is difference between graft node age and height/age of child_subtree"
            child_subtree.calc_node_ages(ultrametricity_precision=False)
            child_subtree_root.edge.length = recomb_time - child_subtree_root.age
        graft_node.add_child(child_subtree_root)
        
        if deep_rooted:
             tree = TreeReconciler.remove_deep_root(tree)
    
    return tree

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
    parent_node = edge.tail_node # parent of attachment edge
    
    "Find constraints on timing of recomb event"
    tree.calc_node_ages(ultrametricity_precision=False)
    max_recomb_time = parent_node.age # max recombination time is height/age of parent node
    child_nodes = parent_node.child_nodes()
    min_recomb_time = max(child_nodes[0].age,child_nodes[1].age) # min recombination time is max of parent nodes children

    return attachment_edge, max_recomb_time, min_recomb_time

"""
Add necessary recombination events to reconcile two trees based on MAF
"""
def add_recombination_events(ref,alt,maf,plot=False):
    
    "Deep copy original trees so we don't lose them"
    #ref = ref_in.clone(depth=2)
    #alt = alt_in.clone(depth=2)
    
    maf.pop() # remove last/largest connected component
    for sub_tree in reversed(maf):
    
        "Get taxa set in subtree"
        subtree_taxa = set([lf.taxon.label for lf in sub_tree.leaf_node_iter()])
            
        "Find edge where subtree attaches in alt and then find equivalent edge in ref"
        #attachment_edge_ref, max_recomb_time_ref, min_recomb_time_ref = find_recombinant_edge(alt,ref,subtree_taxa)
        attachment_edge_alt, max_recomb_time_alt, min_recomb_time_alt = find_recombinant_edge(alt,subtree_taxa)
        #print("Ref attachment edge leaf taxa")
        #print([lf.taxon.label for lf in attachment_edge_ref.head_node.leaf_iter()])
        
        "Find edge where subtree attaches in ref and then find equivalent edge in alt"
        #attachment_edge_alt, max_recomb_time_alt, min_recomb_time_alt  = find_recombinant_edge(ref,alt,subtree_taxa)
        attachment_edge_ref, max_recomb_time_ref, min_recomb_time_ref  = find_recombinant_edge(ref,subtree_taxa)
        #print("Alt attachment edge leaf taxa")
        #print([lf.taxon.label for lf in attachment_edge_alt.head_node.leaf_iter()])
         
        "Pick a recombination event time that satisfies time constraints in both alt and ref trees"
        max_recomb_time = min(max_recomb_time_ref,max_recomb_time_alt) # event must be below parent nodes in either tree
        min_recomb_time = max(min_recomb_time_ref,min_recomb_time_alt) # event must be above parent's children in either tree
        assert max_recomb_time >= min_recomb_time
        recomb_time = np.random.uniform(low=min_recomb_time,high=max_recomb_time) # draw a uniform time between these constraints 
        
        """
        Add single node as a unification representing recombination to both alt and ref tree        
        """
        ref = add_rec_node(ref,attachment_edge_ref,recomb_time)
        alt = add_rec_node(alt,attachment_edge_alt,recomb_time)
        
        if plot:
            print("Reference tree:")
            ref.print_plot()
            print("Alternate tree:")
            alt.print_plot()
        
    return ref, alt

def relabel_taxa(tree_file):
    
    taxa = dendropy.TaxonNamespace()
    tree = dendropy.Tree.get(file=open(tree_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
    for lf in tree.leaf_node_iter():
        lf.taxon.label = str(int(lf.taxon.label) - 1) # to make consistent with taxon lables from msprime
    tree.write(path=tree_file,schema='newick',suppress_annotations=True)

"Simulate a tree series in ms prime"
sim = True
plot = False
desired_breakpoints = 1 # of desired breakpoints in tree series
if sim:
    ts = sim_ts(desired_breakpoints) 
    ts.dump('recomb_placement_test1_tables')
else:
    ts = tskit.load('recomb_placement_test1_tables')

"Load in trees"
tree_file_A = "recomb_placement_test1_treeA.tre"
tree_file_B = "recomb_placement_test1_treeB.tre"

"Relabel -- can move inside sim condition"
relabel_taxa(tree_file_A)
relabel_taxa(tree_file_B)

"Print the simulated ts tables"
#if plot:
    #print("the simulated ts tables")
print(ts.tables.nodes)
print(ts.tables.edges)
for tree in ts.trees():
    print("-" * 20)
    print("tree {}: interval = {}".format(tree.index, tree.interval))
    print(tree.draw(format="unicode"))

"""
    If we don't suppress unifurfactions can get recomb positions from dendropy nodes'
"""
       
"Get trees as dendropy objects"
taxa = dendropy.TaxonNamespace()
ref = dendropy.Tree.get(file=open(tree_file_A, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
print("Reference tree:")
ref.print_plot()
ref.suppress_unifurcations(update_bipartitions=True) 

alt = dendropy.Tree.get(file=open(tree_file_B, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
print("Alternate tree:")
alt.print_plot()
alt.suppress_unifurcations(update_bipartitions=True)
    
"Get maximum agreement forest"
maf = TreeReconciler.get_maf_4cut(ref,alt,plot=False)
print("Trees in MAF: " + str(len(maf)))
#TreeReconciler.plot_maf(maf)

"Add recombination events as unifurcations to trees"
ref, alt = add_recombination_events(ref,alt,maf,plot=True)

tree_file_A = "recomb_events_added_test1_treeA.tre"
tree_file_B = "recomb_events_added_test1_treeB.tre"
ref.write(path=tree_file_A,schema='newick',suppress_annotations=True)
alt.write(path=tree_file_B,schema='newick',suppress_annotations=True)

