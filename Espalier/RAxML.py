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

import numpy as np
import subprocess
import sys
import os
import re
import dendropy
from io import StringIO
import logging

class RAxMLRunner(object):
    
    """
        Wrappers for reconstructing trees and computing likelhoods using RAxML
    """
    
    def __init__(self, **kwargs):
        
        '''             
            Parameters:
               
            Optional keyword arguments:
                temp_dir (str): directory to store temp output
                threads (int): number of threads RAxML can use
                model (str): substitution model with specified rates e.g. GTR{1.0/2.0/1.0/1.0/2.0/1.0}+FU{0.25/0.25/0.25/0.25}+G4m{4.0}
                seed (int): random number seed
                raxml_path (str): path to RAxML application
                lsd_path (str): path to LSD application
        '''
        
        self.temp_dir = kwargs.get('temp_dir', './')
        self.threads = kwargs.get('threads', 1)
        self.model = kwargs.get('model', 'GTR{1.0/1.0/1.0/1.0/1.0/1.0}')
        self.seed = kwargs.get('seed', 12345)
        self.raxml_path = kwargs.get('raxml_path', 'raxml-ng')
        self.lsd_path = kwargs.get('lsd_path', 'lsd')

    def get_tree_likelihood(self,tree_file,seq_file):
        
        '''
            Compute log likelihood of sequnece alignment in seq_file given tree in tree_file
         
            Parameters:     
               tree_file (str): input newick tree file
               seq_file (str): input fasta file containing seq alignment
               
            Returns:
                logL (float): log likelihood of seq data given tree
        '''
        
        logL = 0
    
        temp_prefix = self.temp_dir + 'temp'
        
        # Evaluatae likelihood of seq data given tree in RAxML
        cmd_args = [self.raxml_path,
                    '--evaluate --msa', seq_file,
                    '--threads', str(self.threads),
                    '--model', self.model,
                    '--tree', tree_file,
                    '--prefix', temp_prefix,
                    '--opt-model', 'off',
                    '--opt-branches', 'off', 
                    '--redo']
        cmd = ' '.join(cmd_args)

        try:
            output = subprocess.check_output(cmd, shell=True,stderr=subprocess.STDOUT)
            #sys.stdout.write(output)
        except subprocess.CalledProcessError as exc:
            logging.error(exc.output)
            logging.error('Execution of "%s" failed!\n' % cmd)
            sys.exit(1)
        
        # Parse log likelihood from raxml log file
        f = open(self.temp_dir + 'temp.raxml.log')
        line = f.readline()
        while line:
            if "Final LogLikelihood:" in line:
                logL = float(line.split()[2])
            line = f.readline()
        f.close()
        
        # Clean up
        try:
            os.remove(temp_prefix + '.raxml.bestModel')
            os.remove(temp_prefix + '.raxml.bestTree')
            os.remove(temp_prefix + '.raxml.bestTreeCollapsed')
            os.remove(temp_prefix + '.raxml.log')
            os.remove(temp_prefix + '.raxml.rba')
            os.remove(temp_prefix + '.raxml.siteLH')
            os.remove(temp_prefix + '.raxml.startTree')
        except OSError:
            pass
        
        return logL
    
    def get_site_likelihoods(self,tree_file,seq_file):   
        
        """
            Compute per site log likelihoods of seq data given tree using raxml
            This allows site likes to be computed without optimization of model params or branch lengths
            
            Parameters:     
               tree_file (str): input newick tree file
               seq_file (str): input fasta file containing seq alignment
               
            Returns:
                site_likes (numpy.array): site-specific log likelihood of seq data given tree
            
        """
        
        temp_prefix = self.temp_dir + 'temp'

        # Evaluatae likelihood of seq data given tree in RAxML
        cmd_args = [self.raxml_path,
                    '--sitelh --msa', seq_file,
                    '--threads', str(self.threads),
                    '--model', self.model,
                    '--tree', tree_file,
                    '--prefix', temp_prefix,
                    '--opt-model', 'off',
                    '--opt-branches', 'off',
                    '--force', 'msa_allgaps', # compte likelihoods even for sites with all gaps
                    '--redo']
        cmd = ' '.join(cmd_args)

        try:
            output = subprocess.check_output(cmd, shell=True,stderr=subprocess.STDOUT)
            #sys.stdout.write(output)
        except subprocess.CalledProcessError as exc:
            print(exc.output)
            print('Execution of "%s" failed!\n' % cmd)
            sys.exit(1)
            
        # Read in and parse site-specific likelihoods    
        f = open(temp_prefix + '.raxml.siteLH')
        f.readline() # first line is meta data
        txt = f.readline()
        f.close()
        #txt = txt.split('\t')[1]
        txt = txt.replace('tree1    ','') # very annoyingly raxml-ng inserts spaces here
        site_likes = np.loadtxt(StringIO(txt))
        
        # Clean up
        try:
            os.remove(temp_prefix + '.raxml.bestModel')
            os.remove(temp_prefix + '.raxml.bestTree')
            os.remove(temp_prefix + '.raxml.bestTreeCollapsed')
            os.remove(temp_prefix + '.raxml.log')
            os.remove(temp_prefix + '.raxml.rba')
            os.remove(temp_prefix + '.raxml.siteLH')
            os.remove(temp_prefix + '.raxml.startTree')
        except OSError:
            pass
        
        return site_likes
       
    
    def _parse_best_model(self,best_model_file):
    
        """
            Parse model from *.raxml.bestModel output
        """    
    
        f = open(best_model_file)
        line = f.readline()
        f.close()
        model = line.split(',')[0]
        
        return model
    
    
    def get_dated_raxml_tree(self,seq_file,tree_file,tip_date_file,rate_file,dirty=True,parse_model=True,root=True,outgroup_taxon=None):
        
        """
            Reconstruct ML phylogeny from input alignment file and date tree using LSD
            If root is True, tree is rooted on midpoint unless outgroup_taxon is specified
            Substition model is assumed to be GTR+G for tree reconstruction
            
            Parameters:     
               seq_file (str): Input alignment file
               tree_file (str): Output newick tree file
               tip_date_file (str): Input file with tip sampling dates/times
               rate_file (str): Input file with molecular clock rate (subst per site per unit time)             
               
            Optional keyword arguements:  
               dirty (boolean): Perform quick and dirty tree search from a single starting tree?
               parse_model (boolean): Parse best model from tree reconstruction output and set as current model
               root (boolean): Root unrooted ML tree before writing to output file?
               outgroup_taxon: taxon name label for taxa to root tree on
            
            Returns:     
               None (writes dated ML tree to tree_file)
        """  
        
        temp_prefix = self.temp_dir + 'temp'
        
        # Set tree search mode
        search_mode = ''
        if dirty:
            search_mode = '--search1' 
        
        # Reconstruct ML tree in RAxML
        cmd_args = [self.raxml_path,
            search_mode,
            '--msa', seq_file,
            '--model', 'GTR+G',
            '--prefix', temp_prefix,
            '--threads', str(self.threads),
            '--seed', str(self.seed),
            '--redo']
        cmd = ' '.join(cmd_args)
        
        try:
            output = subprocess.check_output(cmd, shell=True,stderr=subprocess.STDOUT)
            #sys.stdout.write(output)
        except subprocess.CalledProcessError as exc:
            logging.error(exc.output)
            logging.error('Execution of "%s" failed!\n' % cmd)
            sys.exit(1)
        
        
        # Parse best model from RAxML output
        if parse_model:
            best_model_file = temp_prefix + '.raxml.bestModel'
            self.model = self._parse_best_model(best_model_file)
        
        # Rename best ML tree as output tree_file
        try:
            os.rename(temp_prefix + '.raxml.bestTree', tree_file)
        except OSError:
            pass
        
        # Root ML tree using midpoint or, if specified, an outgroup
        if root:
            taxa = dendropy.TaxonNamespace()
            tree = dendropy.Tree.get(file=open(tree_file, 'r'), schema="newick", taxon_namespace=taxa)
            if outgroup_taxon:
                outgroup_node = tree.find_node_with_taxon_label(outgroup_taxon)
                new_root_node = outgroup_node.parent_node
                tree.reroot_at_node(new_root_node,collapse_unrooted_basal_bifurcation=False)
                tree.prune_taxa_with_labels([outgroup_taxon])
            else:
                tree.reroot_at_midpoint()
            tree.write(path=tree_file,schema='newick',suppress_annotations=True,suppress_rooting=True)  
        
        # Run least-squares dating in LSD to get a time-calibrated ML tree
        temp_lsd_prefix = self.temp_dir + 'temp.lsd'
        #cmd = lsd_path + ' -i ' + tree_file + ' -d ' + tip_date_file + ' -c -w ' + rate_file + ' -o ' + temp_lsd_prefix
        cmd_args = [self.lsd_path,
            '-i', tree_file,
            '-d', tip_date_file,
            '-c', # constrain internal node times by tip times
            '-w', rate_file,
            '-o', temp_lsd_prefix]
        cmd = ' '.join(cmd_args)
        try:
            output = subprocess.check_output(cmd, shell=True,stderr=subprocess.STDOUT)
            #sys.stdout.write(output)
        except subprocess.CalledProcessError as exc:
            logging.error(exc.output)
            logging.error('Execution of "%s" failed!\n' % cmd)
            sys.exit(1)
        
        # Rename dated ML tree as output tree_file
        # This may be in a newick or nexus file depending on the LSD version
        dated_newick_file = temp_lsd_prefix + '.date.newick'
        if os.path.isfile(dated_newick_file):
            os.rename(dated_newick_file, tree_file)
        else:
            dated_nexus_file = temp_lsd_prefix + '.date.nexus'
            if os.path.isfile(dated_nexus_file): # extract from nexus
                f = open(dated_nexus_file)
                line = f.readline()
                while line:
                    if "tree 1 =" in line:
                        tree = line.split()[3]
                    line = f.readline()
                tree = re.sub("[\[].*?[\]]", "", tree) # remove metadata in brackets
                f.close()
                nwk=open(tree_file,"w")
                nwk.write(tree + '\n')
                nwk.close()
            else:
                logging.error("ERROR: Cannot find dated newick or nexus file from LSD!!")
        
        # Clean up
        try:
            #os.remove(temp_dir + 'temp.raxml.reduced.phy')
            os.remove(temp_prefix + '.raxml.rba')
            os.remove(temp_prefix + '.raxml.startTree')
            os.remove(temp_prefix + '.raxml.bestModel')
            os.remove(temp_prefix + '.raxml.log')
            os.remove(temp_lsd_prefix  + '.newick')
            os.remove(temp_lsd_prefix  + '.nexus')
            os.remove(temp_lsd_prefix)
            os.remove(temp_prefix + '.raxml.bestTreeCollapsed')
        except OSError as e:  ## if failed, report it back to the user ##
            print ("Error: %s - %s." % (e.filename, e.strerror))
    
if __name__ == '__main__':
    
    raxml = RAxMLRunner(raxml_path='raxml-ng',temp_dir='./')
    
    # Test tree reconstruction using raxml-ng
    seq_file = 'temp_seqs.fasta'
    tree_file = 'temp_ML.tre'
    tip_date_file = 'dates-lsd.txt'
    rate_file = 'rate.txt'
    raxml.get_dated_raxml_tree(seq_file,tree_file,tip_date_file,rate_file)
    
    # Test tree likelihood calculation
    tree_file = 'temp_ML.tre'
    seq_file = 'temp_seqs.fasta'
    L = raxml.get_tree_likelihood(tree_file,seq_file)
    print("Log likelihood: " + f'{L:.2f}')
    
    # Test site-specific likelihoods"
    site_L = raxml.get_site_likelihoods(tree_file,seq_file)
    print("Site-specific log likelihoods: ") 
    #print(site_L)
    
    