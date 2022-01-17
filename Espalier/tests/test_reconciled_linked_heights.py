"""
Created on Thu Sep  2 10:43:04 2021

Test node height reconcilation between trees in path using "lock and chain" method

@author: david
"""
from Espalier import Reconciler
import dendropy

sim_path = './recomb_placement_path/sim001/'
segments = 4
tree_files = [sim_path + "PathTree" + str(i) + ".tre" for i in range(segments)]
taxa = dendropy.TaxonNamespace()
tree_path = []
for seg in range(segments):   
    tree = dendropy.Tree.get(file=open(tree_files[seg], 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
    tree_path.append(tree)


Reconciler.reconcile_linked_heights(tree_path)
        
"Plot original tree with node root dists"
plot = True
if plot:
    age_label_fn = lambda nd: nd.taxon.label if nd.is_leaf() else str(f"{nd.age:0.2f}")
    for tree_idx, tree in enumerate(tree_path):
        print("-" * 20)
        print("Tree #" + str(tree_idx))
        print("-" * 20)
        orig_tree = dendropy.Tree.get(file=open(tree_files[tree_idx], 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
        orig_tree.calc_node_ages(ultrametricity_precision=False)
        print("Original tree:")
        orig_tree.print_plot(show_internal_node_labels=True,node_label_compose_fn=age_label_fn)
        print("New tree:")
        tree.print_plot(show_internal_node_labels=True,node_label_compose_fn=age_label_fn)
        #tree.print_plot(plot_metric='age',show_internal_node_labels=True,node_label_compose_fn=age_label_fn)
        
        