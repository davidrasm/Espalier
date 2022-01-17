##########################
##
## Espalier: A Python package for tree reconciliation and reconstructing ARGs using maximum agreement forests.
##
## Copyright 2021-2022 David A. Rasmussen (drasmus@ncsu.edu)
##
## If using Espalier or this code, please cite:
##
##      Rasmussen, D.A. and Guo, F. Espalier: Efficient tree reconciliation and ARG reconstruction using maximum agreement forests. 2022.
##
############################

import dendropy
from Bio import AlignIO
from Bio import SeqIO
from dendropy.calculate import treecompare

"""
    Helper functions for running Espalier
"""

def write_tip_dates(tree_file,date_file):
    
    """
        Write tip labels and dates to file for running LSD
    
        Parameters:     
            tree_file (str): Newick tree file for tree
            date_file (str): Date file to write dates
    
    """
    
    taxa = dendropy.TaxonNamespace()
    tree = dendropy.Tree.get(file=open(tree_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
    tip_date_dict = {tx:0.0 for tx in taxa.labels()}
    txt=open(date_file,"w")
    txt.write(str(len(taxa.labels())) + '\n')
    for k,v in tip_date_dict.items():
        txt.write(k + '\t' + str(v) + '\n')
    txt.close()
    
def write_rate_file(rate,rate_file):
    
    """
        Write molecular clock rate to file for running LSD
        
        Parameters:     
            rate (float): Molecular clock rate
            rate_file (str): File to write rate to
    """
    
    txt=open(rate_file,"w")
    txt.write(str(rate) + '\n')
    txt.close()
    
def informative_sites(seq_file):
    
    """
        Count number of phylogenetically informative sites (with SNPs) in alignment.
    
        Parameters:     
            seq_file (str): Fasta sequence file containing alignment
            
        Returns:
            i_sites (int): number of informative sites
    
    """
    
    align = AlignIO.read(seq_file, "fasta")
    L = len(align[0])
    i_sites = 0
    for i in range(L):
        if len(set(align[:,i])) > 1:
            i_sites += 1
    return i_sites

def concate_aligns(seq_files,file_out):
    
    """
        Concatenate alignments horizontally such that the 5' end of each sequence
        is concatenated to the 3' end of the sequences in the previous alignments.
        
        TODO: rename concat_aligns?
    
        Parameters:     
            seq_files (list): List of Fasta sequence files containing alignments
            file_out (str): Fasta file to write concatenated alignment
    
    """

    for idx, file in enumerate(seq_files):
        seq_dict = SeqIO.to_dict(SeqIO.parse(file, "fasta")) # one line alternative
        if idx == 0:
            concat_seqs = seq_dict
        else:
            for key in concat_seqs:
                concat_seqs[key] += seq_dict[key]            
    
    # Convert to list and write seq records to fasta file
    concat_records = [concat_seqs[key] for key in concat_seqs]
    SeqIO.write(concat_records, file_out, "fasta")
    
def get_consensus_tree(tree_list,root=True):
    
    """
        Get consensus maximum clade credibility tree from a list of trees using DendroPy
        
        Parameters:     
            tree_list (list[dendropy.Tree] or dendropy.TreeList): List of trees from which consensus will be computed

        Optional keyword aguements:
            root (boolean): Consensus tree will be rooted if True.
            
    """
    
    if isinstance(tree_list, dendropy.TreeList):
        trees = tree_list
    else: 
        trees = dendropy.TreeList()
        for tree_file in tree_list: trees.read(path=tree_file,schema='newick',rooting="default-rooted")
    #consensus = trees.consensus(min_freq=0.95) # does not preserve branch lengths
    consensus = trees.maximum_product_of_split_support_tree()
    #print("\nTree {} maximizes the product of split support (log product = {}): {}".format(trees.index(mcct), mcct.log_product_of_split_support, mcct))
    if not consensus.is_rooted and root:
        consensus.reroot_at_midpoint()
    
    return consensus 
    
def remove_seqs(file_in,file_out,drop_list):
    
    """
        Remove particular sequences in drop_list from file_in and write new seqs to file_out.
    
        Parameters:     
            file_in (str): Input fasta sequence file containing original seqs
            file_out (str): Output fasta sequence file with seqs removed
            drop_list (list): List of sequence/taxa names to remove.
    
    """
    
    records = SeqIO.parse(file_in, "fasta")
    new_records = []
    for rec in records:
        if rec.name not in drop_list:
            new_records.append(rec)
    SeqIO.write(new_records,file_out, "fasta")
    
def tree_in_set(test_tree,trees):
    
    """
        Test if test_tree is in list of trees
        
        Parameters:
            test_tree (dendropy.Tree): Tree to test if in tree list.
            trees (dendropy.TreeList): List of trees to query
            
        Returns:
            True if test_tree in trees
    """
    
    for tr in trees:
        dist = treecompare.symmetric_difference(test_tree, tr) # # Unweighted Robinson-Foulds distance
        if dist == 0:
            return True
    
    
if  __name__ == '__main__':
    
    # Test write_tip_dates
    #tree_file = './sim_trees/original_test_trees/disentangler_test1_tree0.tre'
    #date_file = 'dates-lsd2.txt'
    #write_tip_dates(tree_file,date_file)
    
    # Test write_rate_file
    #rate = 0.001
    #rate_file = 'rate2.txt'
    #write_rate_file(rate,rate_file)
    
    # Test informative sites
    seq_file = './sim_trees/original_test_trees/disentangler_test1_tree0.fasta'
    i_sites = informative_sites(seq_file)
    print(i_sites)
    