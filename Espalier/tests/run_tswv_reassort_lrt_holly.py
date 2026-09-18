#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 5 08:11:33 2026

Run likelihood ratio test (LRT) for statistical support of reassortment events
in Holly's TSWV isolates

@author: davidrasmussen
"""

from Espalier import MAF
from Espalier.RAxML import RAxMLRunner
from Espalier.SPRTester import SPRTester
import dendropy
import pandas as pd
import os
import datetime

"Plot during dubugging"
plot = True

"File names for input trees and alignments"
path = '../tswv-test/holly-rb/'
S_tree_file = path + '6_18_26_S.tre'
S_seq_file = path + '6_18_26_S.fas'
M_tree_file = path + '6_18_26_M.tre'
M_seq_file = path + '6_18_26_M.fas'
L_tree_file = path + '6_18_26_L.tre'
L_seq_file = path + '6_18_26_L.fas'

"Load in trees"
taxa = dendropy.TaxonNamespace()
S_tree = dendropy.Tree.get(file=open(S_tree_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
M_tree = dendropy.Tree.get(file=open(M_tree_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
L_tree = dendropy.Tree.get(file=open(L_tree_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)

"Set up temp dir for temp tree output"
temp_dir = 'temp-{date:%Y-%m-%d_%H:%M:%S}/'.format(date=datetime.datetime.now())
if not os.path.isdir(temp_dir):
    os.mkdir(temp_dir)
    
"Initialize callable instances of Espalier objects"
raxml = RAxMLRunner(raxml_path='raxml-ng',lsd_path='lsd',temp_dir=temp_dir)
tester = SPRTester(raxml,temp_dir=temp_dir)

"""
    Run LRT for S and M segments
"""
maf = MAF.get_maf_4cut(S_tree, M_tree, plot=False)
if plot:
    print("MAF:")
    MAF.plot_maf(maf)
results_df = tester(S_tree,S_seq_file,M_tree,M_seq_file,maf)
results_df.to_csv("seg_S-M_lrt_results.csv",index=False)

"""
    Run LRT for M and L segments
"""
maf = MAF.get_maf_4cut(M_tree, L_tree, plot=False)
if plot:
    print("MAF:")
    MAF.plot_maf(maf)
results_df = tester(M_tree,M_seq_file,L_tree,L_seq_file,maf)
results_df.to_csv("seg_M-L_lrt_results.csv",index=False)

"""
    Run LRT for S and L segments
"""
maf = MAF.get_maf_4cut(S_tree, L_tree, plot=False)
if plot:
    print("MAF:")
    MAF.plot_maf(maf)
results_df = tester(S_tree,S_seq_file,L_tree,L_seq_file,maf)
results_df.to_csv("seg_S-L_lrt_results.csv",index=False)

