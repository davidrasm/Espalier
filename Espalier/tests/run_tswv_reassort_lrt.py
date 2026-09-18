#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 5 08:11:33 2026

Run likelihood ratio test (LRT) for statistical support of reassortment events

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
path = '../tswv-test/original-data/'
#S_tree_file = path + 'tswv-S-mcc-NC109.tre'
S_seq_file = path + 'TSWV-Tr-NCGLOBAL-SegS-NC109.fasta'
#M_tree_file = path + 'tswv-M-mcc-NC109.tre'
M_seq_file = path + 'TSWV-Tr-NCGLOBAL-SegM-NC109.fasta'
#L_tree_file = path + 'tswv-L-mcc-NC109.tre'
L_seq_file = path + 'TSWV-Tr-NCGLOBAL-SegL-NC109.fasta'


"Set up temp dir for temp tree output"
temp_dir = 'temp-{date:%Y-%m-%d_%H:%M:%S}/'.format(date=datetime.datetime.now())
if not os.path.isdir(temp_dir):
    os.mkdir(temp_dir)
    
"Initialize callable instances of Espalier objects"
raxml = RAxMLRunner(raxml_path='raxml-ng',lsd_path='lsd',temp_dir=temp_dir)
tester = SPRTester(raxml,temp_dir=temp_dir)

"Set tree files"
path = '../tswv-test/'
#S_tree_file = path + 'tswv-S-mcc-NC109.tre'
#M_tree_file = path + 'tswv-M-mcc-NC109.tre'
#L_tree_file = path + 'tswv-L-mcc-NC109.tre'
S_tree_file = path + 'trees/TSWV-NC-segS-n109_ML.tre'
M_tree_file = path + 'trees/TSWV-NC-segM-n109_ML.tre'
L_tree_file = path + 'trees/TSWV-NC-segL-n109_ML.tre'

"Or infer ML trees"
# raxml.get_raxml_tree(S_seq_file, S_tree_file)
# raxml.get_raxml_tree(M_seq_file, M_tree_file)
# raxml.get_raxml_tree(L_seq_file, L_tree_file)

"Load in trees"
taxa = dendropy.TaxonNamespace()
S_tree = dendropy.Tree.get(file=open(S_tree_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
M_tree = dendropy.Tree.get(file=open(M_tree_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
L_tree = dendropy.Tree.get(file=open(L_tree_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)

"""
    Run LRT for S and M segments
"""
# maf = MAF.get_maf_4cut(S_tree, M_tree, plot=False)
# if plot:
#     print("MAF:")
#     MAF.plot_maf(maf)
# results_df = tester(S_tree,S_seq_file,M_tree,M_seq_file,maf)
# results_df.to_csv("tswv_NC109_seg_S-M_lrt_results.csv",index=False)

"""
    Run LRT for M and L segments
"""
maf = MAF.get_maf_4cut(M_tree, L_tree, plot=False)
if plot:
    print("MAF:")
    MAF.plot_maf(maf)
results_df = tester(M_tree,M_seq_file,L_tree,L_seq_file,maf)
results_df.to_csv("tswv_NC109_seg_M-L_lrt_results.csv",index=False)

"""
    Run LRT for S and L segments
"""
# maf = MAF.get_maf_4cut(S_tree, L_tree, plot=False)
# if plot:
#     print("MAF:")
#     MAF.plot_maf(maf)
# results_df = tester(S_tree,S_seq_file,L_tree,L_seq_file,maf)
# results_df.to_csv("tswv_NC109_seg_S-L_lrt_results.csv",index=False)

