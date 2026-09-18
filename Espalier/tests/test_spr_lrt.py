#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 09:11:33 2026

Test likelihood ratio test (LRT) for statistical support of SPR/reassortment events

@author: davidrasmussen
"""

#from Espalier.sim import ARGSimulator
#from Espalier.sim import TreeRandomizer
from Espalier import MAF
#from Espalier import Utils
from Espalier.RAxML import RAxMLRunner
from Espalier.SPRTester import SPRTester
#from Espalier import Utils
import dendropy
#from dendropy.calculate import treecompare
#import msprime
#import numpy as np
import pandas as pd
import os
import datetime

"Plot during dubugging"
plot = False

"File names for input trees and alignments"
path = '../tswv-test/'
ref_tree_file = path + 'trees/TSWV-NC-segS-n30_ML.tre'
ref_seq_file = path + 'alignments/TSWV-NC-segS-n30.fasta'
alt_tree_file = path + 'trees/TSWV-NC-segM-n30_ML.tre'
alt_seq_file = path + 'alignments/TSWV-NC-segM-n30.fasta'

"Load in trees"
taxa = dendropy.TaxonNamespace()
ref_tree = dendropy.Tree.get(file=open(ref_tree_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
alt_tree = dendropy.Tree.get(file=open(alt_tree_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)

# if not os.path.isdir(path):
#     os.mkdir(path)
# results_file = path + "reassortment_ltr_test_results.csv"

"Set up temp dir for temp tree output"
temp_dir = 'temp-{date:%Y-%m-%d_%H:%M:%S}/'.format(date=datetime.datetime.now())
if not os.path.isdir(temp_dir):
    os.mkdir(temp_dir)
    
"Initialize callable instances of Espalier objects"
raxml = RAxMLRunner(raxml_path='raxml-ng',lsd_path='lsd',temp_dir=temp_dir)
tester = SPRTester(raxml,temp_dir=temp_dir)

"Get maximum agreement forest"
#maf = MAF.get_maf_4cut(alt_tree,ref_tree,plot=False)
maf = MAF.get_maf_4cut(ref_tree,alt_tree,plot=False)
if plot:
    print("MAF:")
    MAF.plot_maf(maf)

"""
    Run likelihood ratio test to identify statistically well supported recombination/reassortment events
    from the SPR events identified in a MAF
"""
results_df = tester(ref_tree,ref_seq_file,alt_tree,alt_seq_file,maf)
results_df.to_csv("test_spr_lrt_results.csv",index=False)


