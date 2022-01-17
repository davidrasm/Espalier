"""
Created on Wed Dec 22 15:23:09 2021

Compute SPR distance between two trees from the simple_maf.py example

@author: david
"""
from Espalier import MAF
import dendropy

tree1_file = './maf_example_tree1.tre'
taxa = dendropy.TaxonNamespace()
tree1 = dendropy.Tree.get(file=open(tree1_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)

tree2_file = './maf_example_tree2.tre'
tree2 = dendropy.Tree.get(file=open(tree2_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)

# Get SPR distance
spr_dist = MAF.get_spr_dist(tree1,tree2)
print("SPR distance: " + str(spr_dist))

# Get maximum agreement forest
maf = MAF.get_maf(tree1,tree2)
maf_spr_dist = len(maf) - 1
print("SPR distance from MAF: " + str(maf_spr_dist))


