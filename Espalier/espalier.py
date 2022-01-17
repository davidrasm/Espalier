#!/usr/bin/env python

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

from Espalier.MAF import get_maf_4cut
from Espalier.MAF import plot_maf
from Espalier.Reconciler import Reconciler
from Espalier.RAxML import RAxMLRunner
import click
import dendropy
import os
import datetime

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(context_settings=CONTEXT_SETTINGS)
def main():
    pass


@main.command()
@click.argument('tree1',type=click.Path(exists=True))
@click.argument('tree2',type=click.Path(exists=True))
@click.argument('file_out', type=click.File('w'),default='-',required=False)
@click.option('--plot', is_flag=True, help='Plots subtrees in MAF')
def maf(tree1,tree2,file_out,plot):
    
    """
        Compute maximum agreement forest (MAF) between two discordant trees. Tree files should be in Newick format and rooted.
        
        Arguements:     
               tree1 (str): Newick tree file for tree1
               tree2 (str): Newick tree file for tree2
               file_out (str): Output file to write subtrees in MAF. Default '-' writes to stdout.
    """
    
    # Read in trees
    taxa = dendropy.TaxonNamespace()
    tree1 = dendropy.Tree.get(file=open(tree1, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
    tree2 = dendropy.Tree.get(file=open(tree2, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)

    # Compute MAF
    click.echo('Computing MAF...') 
    forest = get_maf_4cut(tree1,tree2)
    click.echo('MAF contains ' + str(len(forest)) + ' subtrees:')
        
    # Write subtrees to file_out
    for tr in forest:
        click.echo(tr.as_string(schema="newick",suppress_annotations=True,suppress_rooting=True), file=file_out)
    
    if plot:
        plot_maf(forest)


@main.command()
@click.argument('tree1',type=click.Path(exists=True))
@click.argument('tree2',type=click.Path(exists=True))
def spr(tree1,tree2):
    
    """
    
        Compute subtree-prune-and-regraft distance between two trees through their maximum agreement forest. 
        Tree files should be in Newick format and rooted.
        
        Arguements:     
               tree1 (str): Newick tree file for tree1
               tree2 (str): Newick tree file for tree2
    """
    
    # Read in trees
    taxa = dendropy.TaxonNamespace()
    tree1 = dendropy.Tree.get(file=open(tree1, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
    tree2 = dendropy.Tree.get(file=open(tree2, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
    
    # Compute SPR distance from MAF
    click.echo('Computing SPR distance...') 
    forest = get_maf_4cut(tree1,tree2)
    d_spr = len(forest) - 1
    click.echo("SPR distance = " + str(d_spr))

@main.command()
@click.argument('tree1',type=click.Path(exists=True))
@click.argument('tree2',type=click.Path(exists=True))
@click.argument('seq_file',type=click.Path(exists=True))
@click.argument('file_out', type=click.File('w'),default='-',required=False)
@click.option('--lbr', default=0.1, help='Lower bound ratio for iterative regrafting')
@click.option('--best', is_flag=True, help='Only return reconciled tree with highest likelihood')
@click.option('--temp_dir', default='temp-{date:%Y-%m-%d_%H:%M:%S}/'.format(date=datetime.datetime.now()), help='Directory to store temp output')
@click.option('--raxml_path', default='raxml-ng', help='Path to RAxML-NG executable')
@click.option('--lsd_path', default='lsd', help='Path to LSD executable')
#@click.option('--prior_gamma', default=0.0, help='Exponential decay rate for likelihood contribution from external sites')
def reconcile(tree1,tree2,seq_file,file_out,lbr,best,temp_dir,raxml_path,lsd_path):
    
    """
        Reconcile two discordant trees through their maximum agreement forest (MAF). Tree files should be in Newick format and rooted.
        Discordances will be reconciled if not strongly supported by sequence data in seq_file
        Returns set of reconciled trees unless sorted by the log likelihood of the sequence data given each reconciled tree.
        Does not currently allow for infomration from external sites to be considered when reconciling local trees.
        
        Arguements:     
               tree1 (str): Newick tree file for tree1
               tree2 (str): Newick tree file for tree2
               seq_file (str): Fasta sequence file used to evaluate support for discordances.
               file_out (str): Output file to write reconciled trees. Default '-' writes to stdout.

    """
    
    # Set up temp dir for temp tree output
    if not os.path.isdir(temp_dir):
        os.mkdir(temp_dir)

    # Initialize Espalier objects
    raxml = RAxMLRunner(raxml_path='raxml-ng',lsd_path='lsd',temp_dir=temp_dir)
    reconciler = Reconciler(raxml,lower_bound_ratio=lbr,prior_gamma=0.0,temp_dir=temp_dir)

    # Read in trees
    taxa = dendropy.TaxonNamespace()
    tree1 = dendropy.Tree.get(file=open(tree1, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
    tree2 = dendropy.Tree.get(file=open(tree2, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
    
    # Reconcile trees through their MAF
    forest = get_maf_4cut(tree1,tree2)
    reconciled = reconciler(tree1,tree2,forest,seq_file)
    reconciled.sort(key=lambda x: x.like) # sort by ascending likelihood
    reconciled.reverse() # now sorted in descending order

    # Write reconciled trees to file_out
    if best:
        rec_tree = reconciled[0].tree # trees are sorted such that first tree will have highest likelihood
        click.echo(rec_tree.as_string(schema="newick",suppress_annotations=True,suppress_rooting=True), file=file_out)    
    else:
        for rec_tree in reconciled:
            tr = rec_tree.tree
            like = rec_tree.like
            click.echo('Tree log likelihood: ' + f'{like:.2f}', file=file_out)
            click.echo(tr.as_string(schema="newick",suppress_annotations=True,suppress_rooting=True), file=file_out)

if __name__ == '__main__':
    main()