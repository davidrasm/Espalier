#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 10:39:45 2022

@author: david
"""

from Espalier.Dendro2TSConverter import convert
import dendropy
import tskit

sim_path = './recomb_placement_path/sim004/'
segments = 7
genome_length = 1e3
tree_files = [sim_path + "PathTree" + str(i) + ".tre" for i in range(segments)]
taxa = dendropy.TaxonNamespace()
tree_path = []
for seg in range(segments):   
    tree = dendropy.Tree.get(file=open(tree_files[seg], 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
    tree_path.append(tree)
path_length = len(tree_path)

"Get interval info from TreeSequence as we would always know this anyways"
ts = tskit.load(sim_path + 'recomb_placement_tables')

"Get tree intervals from original ts"
tree_intervals = []
for seg in range(segments):
    left_pos = ts.at_index(seg).interval.left
    right_pos = ts.at_index(seg).interval.right
    tree_intervals.append((left_pos,right_pos))

"Covert tree_path to TreeSequence"
#tree_interals = [(i,j) for i,j in zip(seg_starts,seg_ends)]
merged_ts = convert(tree_path,tree_intervals)
  
"""
    Print temp TreeSequence as SVG for visualization
"""
#svg_size = (800, 250)
#ts.draw_svg(path='original_tree_sequence.svg',size=svg_size,time_scale="log_time", x_scale="treewise")
#merged_ts.draw_svg(path='converted_tree_sequence.svg',size=svg_size,time_scale="log_time", x_scale="treewise")

"""
    Get original tskit tables and trees to compare against merged
"""

print("Original nodes table")
print(ts.tables.nodes)
print()
print("Converted nodes table")
print(merged_ts.tables.nodes)
print("-" * 20)
print()

print("Original edges table")
print(ts.tables.edges)
print()
print("Converted edges table")
print(merged_ts.tables.edges)
print("-" * 20)
print()