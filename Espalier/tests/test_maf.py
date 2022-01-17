"""
Created on Tue May 12 10:39:44 2020

Test maximum agreement forest algorithm with random trees created by SPR moves in reference

@author: david
"""
from Espalier import MAF
import dendropy

ref_file = "./randomSPRs_refTree_test1.tre"
alt_file = "./randomSPRs_altTree_test1.tre"

taxa = dendropy.TaxonNamespace()
ref = dendropy.Tree.get(file=open(ref_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
alt = dendropy.Tree.get(file=open(alt_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)

print("Reference tree:")
ref.print_plot()
print()

print("Alternate tree:")
alt.print_plot()
print()

"Get maximum agreement forest"
maf = MAF.get_maf_4cut(ref,alt)

"Or using 3-cut algorithm"
#maf = MAF.get_maf_3approx(ref,alt)


print()
print("Maximum agreement forest:")
for tr in maf:
    print()
    if len(tr.nodes()) > 1:
        tr.print_plot()
    else:
        print('---' + tr.nodes()[0].taxon.label)