"""
Created on Tue May 12 10:39:44 2020

Test core functions and imports in Espalier

@author: david
"""
from Espalier import MAF
from Espalier import TreeOps
from Espalier.Reconciler import Reconciler
from Espalier.RAxML import RAxMLRunner
from Espalier.ARGBuilder import ARGBuilder
from Espalier.Dendro2TSConverter import convert
import dendropy
from pathlib import Path
import os
import datetime

"""
    Test MAF
"""
examples_dir = Path(__file__).parent.parent.parent / 'examples' # the parent dir of wherever we are
ref_file = examples_dir / "maf_example_tree1.tre"
alt_file = examples_dir / "maf_example_tree2.tre"

taxa = dendropy.TaxonNamespace()
ref = dendropy.Tree.get(file=open(ref_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
alt = dendropy.Tree.get(file=open(alt_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
maf = MAF.get_maf_4cut(ref,alt)

print()
print("Maximum agreement forest:")
for tr in maf:
    print()
    if len(tr.nodes()) > 1:
        tr.print_plot()
    else:
        print('---' + tr.nodes()[0].taxon.label)
        
"""
    Test Reconciler   
"""

# Set up temp dir for temp tree output
temp_dir = 'temp-{date:%Y-%m-%d_%H:%M:%S}/'.format(date=datetime.datetime.now())
if not os.path.isdir(temp_dir):
    os.mkdir(temp_dir)

# Initialize Espalier objects
raxml = RAxMLRunner(raxml_path='raxml-ng',lsd_path='lsd',temp_dir=temp_dir)
reconciler = Reconciler(raxml,lower_bound_ratio=0.1,prior_gamma=0.0,temp_dir=temp_dir)

# Set up files
seq_file_r2 = examples_dir / 'reconciler_r2.fasta'
ML_tree_file_r1 = examples_dir / 'reconciler_MLTree_r1.tre'
ML_tree_file_r2 = examples_dir / 'reconciler_MLTree_r2.tre'   

taxa = dendropy.TaxonNamespace()
tree_r1 = dendropy.Tree.get(file=open(ML_tree_file_r1, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
tree_r2 = dendropy.Tree.get(file=open(ML_tree_file_r2, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)

# Reconcile ML trees through their MAF
maf = MAF.get_maf_4cut(tree_r1,tree_r2)
sampled_trees = reconciler(tree_r1,tree_r2,maf,seq_file_r2)