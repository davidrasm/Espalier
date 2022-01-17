"""
Created on Wed Dec 22 15:23:09 2021

Compute a MAF between two trees using a simple example

@author: david
"""
from Espalier import MAF
import dendropy

tree1_file = './maf_example_tree1.tre'
taxa = dendropy.TaxonNamespace()
tree1 = dendropy.Tree.get(file=open(tree1_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)

tree2_file = './maf_example_tree2.tre'
tree2 = dendropy.Tree.get(file=open(tree2_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)

print("Tree 1:")
tree1.print_plot()

print("Tree 2:")
tree2.print_plot()

# Get maximum agreement forest
maf = MAF.get_maf(tree1,tree2)

# Plot MAF
from Espalier.MAF import plot_maf
plot_maf(maf)

