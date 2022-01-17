"""
Created on Tue Jul 14 09:00:27 2020

@author: david
"""

from Espalier import TreeOps
import dendropy

ref_file = "./maf_randomSPRs_refTree_test1.tre"
#alt_file = "maf_2TopoChanges_test1_treeB.tre"
taxa = dendropy.TaxonNamespace()
ref = dendropy.Tree.get(file=open(ref_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa) 
alt = ref.clone(depth=2)

print("Reference tree:")
ref.print_plot()
print()

new_ref_tree = TreeOps.add_deep_root(ref,reroot=True)
print("Deep rooted tree:")
new_ref_tree.print_plot()
print()

if not new_ref_tree.is_rooted:
    print("Tree has lost rooting")
else:
    print("Tree is rooted")

"Does adding dummy preserve concordance?"
#from Espalier import MAF
#new_ref_tree,new_alt_tree = TreeOps.deep_root_pair(ref,alt)
#discordant = MAF.test_discordance(new_ref_tree,new_alt_tree)
#print(f"Trees discordant: {discordant}")

"Test removal"
derooted_tree = TreeOps.remove_deep_root(new_ref_tree)
print("De-rooted tree:")
derooted_tree.print_plot()
print()

if not derooted_tree.is_rooted:
    print("Tree has lost rooting")
else:
    print("Tree is rooted")

"Get maximum agreement forest"
#maf = TreeReconciler.get_maf_4cut(ref,alt,plot=False)