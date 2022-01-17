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

import tskit
import numpy as np
import pandas as pd
import copy
import logging

def convert(tree_path,tree_intervals):
    
    """
        Main conversion routine for converting tree_path with recombination events into tskit TreeSequence.
        Local trees are first converted into node and edge tables, which are sequentially merged across genome regions (segments)
        After merging, duplicate nodes/edges are removed from these tables.
        Two checks are performed on tables to ensure a valid ARG can be constructed:
        1) A "post-order" check to ensure all nodes have one and only one parent
        2) A "pre-order" check to ensure no node has more than two children
        Additional "hidden" nodes may be added during checks to satisfy parent-child requirements
        While these two checks are necessary to ensure a valid ARG, they are not sufficient and conversion may still fail.
        
        Important: all unifurcations in local trees will be interpreted as recombination events
        
        Parameters:     
            tree_path (list(dendropy.Tree)): tree path
            tree_intervals (list(tuple)): list of genomic intervals for each local tree specified as a tuple (start,end) 
                 
        Returns:     
           merged_ts (tskit.TreeSequence): tree sequence representing a connected ARG

    """
    
    # Get local tree and ts for first interval
    segments = len(tree_path)
    tree = tree_path[0]
    tip_count = tree.__len__()
    left_pos = tree_intervals[0][0]
    right_pos = tree_intervals[0][1]
    total_length = tree_intervals[-1][-1]
    
    # Convert tree to node and edge dataframes
    merged_node_df, tree = tree2nodesdf(tree)
    id_dict = dict(zip(merged_node_df['unique_id'].astype('int'), merged_node_df.index))
    merged_edge_df = tree2edgedf(tree,left_pos,right_pos,id_dict)
    
    logging.debug("Merging node and edge tables")
    
    # Iteratively add nodes and edges from each local tree into merged dataframes 
    for loc in range(1,segments):
        
        tree = tree_path[loc]
        left_pos = tree_intervals[loc][0] #ts.first().interval.left
        right_pos = tree_intervals[loc][1] #ts.first().interval.right
        
        # Get node df for next tree and merge node tables
        next_node_df, tree = tree2nodesdf(tree, prev_node_df=merged_node_df)
        merged_node_df = merged_node_df.append(next_node_df,ignore_index=True)
        merged_node_df.drop_duplicates(ignore_index=True,inplace=True)
        merged_node_df.sort_values('metadata', ignore_index=True, inplace=True) # first sort on metadata (sample/taxon labels)
        merged_node_df.sort_values('time', ignore_index=True, inplace=True) # then sort on times
        
        # Create dictionary to map unique ids to new tskit node ids
        id_dict = dict(zip(merged_node_df['unique_id'].astype('int'), merged_node_df.index))
        
        # Reindex edge df to update parent/child ids based on mapping in id_dict
        # Explaination: Once we've merged the node tables we may have added new coalescent/recomb nodes
        # So we need to get a new id_dict to map unique_ids to new node indexes
        # BUT merged_edge_tables still has old node ids so we need to reindex parent/child ids in edge table before we merge 
        merged_edge_df = reindex_edgedf(merged_edge_df,id_dict)
        
        # Get edge df for next trees and merge edge tables
        next_edge_df = tree2edgedf(tree,left_pos,right_pos,id_dict)
        merged_edge_df = merged_edge_df.append(next_edge_df)
        
        # Check if we can assemble TableCollection and TreeSequence after each merger
        #merged_tables = df2TreeTables(merged_edge_df,merged_node_df,total_length)
        #merged_tables.edges.squash() # squash sorts and combines "equivalent" adjacent edges
        #merged_ts = merged_tables.tree_sequence()
    
    # Check to make sure all child node times are below parent times 
    #check_time_constraints(merged_edge_df,merged_node_df)
    
    # Check to make sure each child has only one parent over each genomic interval
    check_parent_child_intervals(merged_edge_df,merged_node_df)

    merged_tables = df2TreeTables(merged_edge_df,merged_node_df,total_length)
    #merged_tables.sort() # squash will sort edges below
    merged_tables.edges.squash() # squash sorts and combines "equivalent" adjacent edges
    merged_ts = merged_tables.tree_sequence()
    
    # Work with TableCollection for checking instead of TreeSequence b/c TS are not mutable
    tables = merged_ts.dump_tables() # return copy of tables we'll work with
    edges_df, nodes_df = treeTables2df(tables)
    
    # Run post-order traversal checking child-parent relationships
    nodes_df, edges_df = postorder_check(edges_df,nodes_df,tip_count)
    
    # Run pre-order traversal checking parent-child relationships
    nodes_df, edges_df = preorder_check(edges_df,nodes_df,tip_count)
    
    # Convert back to tskit tables
    edges_df.drop_duplicates(inplace=True) # need to drop dups or squash will throw conflicting child/parent error
    merged_tables = df2TreeTables(edges_df,nodes_df,total_length)
    merged_tables.edges.squash() # squash sorts and combines "equivalent" adjacent edges
    
    # Attemp to simplify ARG?
    # drop_list = []
    # for nd_idx, nd in enumerate(tables.nodes):
    #     if tables.nodes[nd_idx].flags == 0: # if coal node
    #         nd_children = list(set(edges_df.loc[edges_df['parent'] == nd_idx, 'child'].tolist()))
    #         if len(nd_children) < 2:
    #             drop_list.append(nd_idx)
    # node_set = set(list(range(len(tables.nodes))))
    # nodes_to_retain = list(node_set.difference(set(drop_list)))
    # #nodes_to_retain = list(range(10)) + [nd_idx for nd_idx, nd in enumerate(tables.nodes) if tables.nodes[nd_idx].flags == 131072]
    # merged_tables.simplify(samples=nodes_to_retain) # simplifies but retains samples and nodes_to_retain
    
    merged_ts = merged_tables.tree_sequence()  
        
    return merged_ts


def postorder_check(edges_df,nodes_df,tip_count):

    """
        Run postorder ARG check making sure all nodes have one and only one parent 
        Adds additional "hidden" nodes required to have proper TreeSequence
    """
    
    logging.info("Post-order ARG traversal checking child-parent relationships")   
     
    # Post-order traversal: length of nodes_df can change while iterating due to adding nodes 
    nd_idx = 0
    while nd_idx < len(nodes_df.index):
        
        #print("Current node: " + str(nd_idx))
        
        # Get parents of current node
        nd_parents = list(set(edges_df.loc[edges_df['child'] == nd_idx, 'parent'].tolist()))
        
        #print("Parents: ", nd_parents)
        #if nd_idx == 28:
            #print()
        
        # Add necessary nodes if node has more than one parent in ARG
        if len(nd_parents) > 1:
            
            # Determine if parents are coalescent or recombination nodes
            parent_types = [nodes_df.at[x,'flags'] for x in nd_parents]
            #coal_parents = [y for x,y in enumerate(nd_parents) if parent_types[x] == 0] # not used
            recomb_parents = [y for x,y in enumerate(nd_parents) if parent_types[x] == 131072]
            
            # Determine order of parents by their node times and get the most recent parent
            parent_times = [nodes_df.at[x,'time'] for x in nd_parents]
            mrp = nd_parents[np.argmin(parent_times)] # most recent parent
            mrp_type = nodes_df.at[mrp,'flags']
            
            if mrp_type == 131072: # recomb event
                
                if len(recomb_parents) > 2:
                    logging.error('Warning: node has more than two recombination node parents')
            
                # Find recombinant parents
                recomb_parent_edges_1 = edges_df.loc[(edges_df['parent'] == recomb_parents[0]) & (edges_df['child'] == nd_idx)]
                recomb_parent_edges_2 = edges_df.loc[(edges_df['parent'] == recomb_parents[1]) & (edges_df['child'] == nd_idx)]
                if recomb_parent_edges_1.left.min() < recomb_parent_edges_2.left.min():
                    left_recomb_parent = recomb_parents[0]
                    left_recomb_edge = recomb_parent_edges_1
                    right_recomb_parent = recomb_parents[1]
                elif recomb_parent_edges_1.left.min() > recomb_parent_edges_2.left.min(): # replace with genome length param
                    left_recomb_parent = recomb_parents[1]
                    left_recomb_edge = recomb_parent_edges_2
                    right_recomb_parent = recomb_parents[0]
                rec_breakpoint = left_recomb_edge['right'].values[0]
                
                # Split edges missing recombination event by adding recombination node between parent and child
                deeper_parents = copy.deepcopy(nd_parents)
                deeper_parents.remove(left_recomb_parent)
                deeper_parents.remove(right_recomb_parent)
                
                # Check if adding more recent rec node between node and deeper parents createss a conflicting parent-child edge
                # Current code cannot handle this condition
                conflict = check_conflicting_edges(edges_df,deeper_parents,left_recomb_parent,nd_idx)
                if conflict:
                    logging.error("Adding left recombinant parent conflicts with another parent-child relationship")
                conflict = check_conflicting_edges(edges_df,deeper_parents,right_recomb_parent,nd_idx)
                if conflict:
                    logging.error("Adding right recombinant parent conflicts with another parent-child relationship")
                
                # Add left/right recomb parent to all edges between parent and current node
                for parent_nd in deeper_parents: # should be nd_parent with recomb_parents removed
                    edges_df = split_edge_on_rec(edges_df,parent_nd,nd_idx,left_recomb_parent,right_recomb_parent,rec_breakpoint)
                    
            elif mrp_type == 0: # coalescent event
                
                # Find parents deeper in ARG than most recent parent (mrp)
                deeper_parents = copy.deepcopy(nd_parents)
                deeper_parents.remove(mrp)
                recent_parent = mrp #coal_parents[np.argmin(coal_parent_times)]
                
                # Check if adding more recent parent between node and deeper parent creates a conflicting parent-child edge
                conflict = check_conflicting_edges(edges_df,deeper_parents,recent_parent,nd_idx)
                if conflict:
                    
                    logging.debug("Found conflict in parent-child relations: need to add additional recombination node")
                    
                    # Find recomb breakpoint by extent of parent-child relationships  
                    child_recent_parent_edges = edges_df.loc[(edges_df['parent'] == recent_parent) & (edges_df['child'] == nd_idx)]
                    child_deeper_parent_edges = edges_df.loc[(edges_df['parent'] == deeper_parents[0]) & (edges_df['child'] == nd_idx)]
                    if child_recent_parent_edges.left.min() < child_deeper_parent_edges.left.min():
                        rec_breakpoint = child_recent_parent_edges.right.max() # if left of conflicting edges
                    elif child_recent_parent_edges.left.min() > child_deeper_parent_edges.left.min(): # replace with genome length param
                        rec_breakpoint = child_recent_parent_edges.left.min() # if right of conflicting edges
                    else:
                        # Not sure if this ever happens
                        logging.debug('Could not localize recomb event based on child-recent_parent start/end')
                    
                    # Add two addtional hidden rec nodes to nodes_df
                    time = nodes_df.at[nd_idx,'time'] + (nodes_df.at[recent_parent,'time'] - nodes_df.at[nd_idx,'time'])/2 # can just pick a random time between parent and child
                    new_rec_node = {'flags': 131072, 'population': -1, 'individual': -1, 'time':time, 'metadata':time} # add extra edge for old parent and split node
                    nodes_df = nodes_df.append(new_rec_node, ignore_index = True)
                    left_recomb_parent = len(nodes_df.index) - 1
                    nodes_df = nodes_df.append(new_rec_node, ignore_index = True)
                    right_recomb_parent = len(nodes_df.index) - 1
                    deeper_parents.append(recent_parent) # need to include recent parent
                    for parent_nd in deeper_parents:
                        edges_df = split_edge_on_rec(edges_df,parent_nd,nd_idx,left_recomb_parent,right_recomb_parent,rec_breakpoint)
                        
                    # Resort nodes_df in chronological order an re-index parents/children in edges_df
                    nodes_df, edges_df = resort_dfs(tip_count, nodes_df, edges_df)
                
                else:
                    
                    # Add more recent_parent to all edges between parent and child node
                    for parent_nd in deeper_parents:
                        # Split all edges with parent and child node into two edges with recent_parent between them
                        edges_df = split_edge(edges_df,parent_nd,recent_parent,nd_idx)
        
            # Check if we can assemble TalbeCollection and TreeSequence after each merger
            # edges_df.drop_duplicates(inplace=True) # need to drop dups or squash will throw conflicting child/parent error
            # updated_tables = df2TreeTables(edges_df,nodes_df,1000)
            # updated_tables.edges.squash() # squash sorts and combines "equivalent" adjacent edges
            # updated_ts = updated_tables.tree_sequence()            
            # for tree in updated_ts.trees():
            #     print("-" * 20)
            #     print("tree {}: interval = {}".format(tree.index, tree.interval))
            #     print(tree.draw(format="unicode"))
            
        nd_idx += 1
        
    return nodes_df, edges_df

def preorder_check(edges_df,nodes_df,tip_count):
    
    """
        Run pre-order ARG check making sure all nodes have no more than two children
        And that rec nodes have no more than one child.
        Adds additional "hidden" nodes required to have proper TreeSequence
        
        TODO: Never fully implemented. Works in "most" cases where we need to add a single node so that parent has single child across all local trees.
        But does not necessarily work in more complex cases where added node still has multiple children.
        Also in the future we should check if splitting edge on a new child creates a conflicting parent-child relationship in the ARG.
        If it does we likely need to add additional recombination nodes as we do for conflicting relationships in the post-order check.
        
    """
    
    logging.info("Pre-order ARG traversal checking parent-child relationships")   
    
    # Pre-order traversal: length of nodes_df can change while iterating due to adding nodes 
    nd_idx = len(nodes_df.index) - 1
    while nd_idx >= 0:
        
        #print("Current node: " + str(nd_idx))
        
        # Get children of current node
        nd_children = list(set(edges_df.loc[edges_df['parent'] == nd_idx, 'child'].tolist()))
        #print("Children: ", nd_children)
        
        # Get parent type (coalescent/recombination)
        parent_type = nodes_df.at[nd_idx,'flags']
        
        if parent_type == 131072: # recomb event
            
            # Not currently accounted for
            if len(nd_children) > 1:
                logging.error('Warning: recombination node has too many children')
        
        elif parent_type == 0:
        
            if len(nd_children) > 2:
                
                logging.debug("Coalescent node has multiple children: need to add additional recombination node")
                
                child_edges = edges_df.loc[edges_df['parent'] == nd_idx]
                
                # Get left and right limits of inteveral over which parent has multiple children
                breakpoints = list(set(child_edges.left.values.tolist() + child_edges.right.values.tolist()))
                breakpoints.sort()
                interval_lefts = []
                interval_rights = []
                for bp_idx in range(len(breakpoints)-1):
                    start_loc = breakpoints[bp_idx]
                    end_loc = breakpoints[bp_idx+1]
                    interval_children = list(set(child_edges.loc[(child_edges['left'] <= start_loc) & (child_edges['right'] >= end_loc), 'child'].tolist()))
                    if len(interval_children) > 1:
                        interval_lefts.append(start_loc)
                        interval_rights.append(end_loc)
                        if len(interval_children) > 2:
                            # This should not happen and not sure how to fix it if it does
                            logging.error("WARNING: Coalescent node has more than two children in the same genomic interval")
                left_limit = min(interval_lefts)
                right_limit = max(interval_rights)
                
                # Find children of node that change between intervals
                changing_children = []
                for child in nd_children:
                    child_recent_parent_edges = edges_df.loc[(edges_df['parent'] == nd_idx) & (edges_df['child'] == child)]
                    child_left_limit = child_recent_parent_edges.left.min()
                    child_right_limit = child_recent_parent_edges.right.max()
                    if child_left_limit > left_limit or child_right_limit < right_limit:
                        changing_children.append(child)
                #print("Changing children: " + str(len(changing_children)))
                
                child_times = [nodes_df.at[x,'time'] for x in changing_children] 
                min_coal_time = max(child_times) #max(nodes_df.at[changing_children[0],'time'],nodes_df.at[changing_children[1],'time'])
                time = nodes_df.at[nd_idx,'time'] - (nodes_df.at[nd_idx,'time'] - min_coal_time)/2 # can just pick a random time between parent and child
                #time = nodes_df.at[nd_idx,'time'] - (nodes_df.at[nd_idx,'time'] - min_coal_time) * 0.01
                new_coal_node = {'flags': 0, 'population': -1, 'individual': -1, 'time':time, 'metadata':time} # add extra edge for old parent and split node
                nodes_df = nodes_df.append(new_coal_node, ignore_index = True)
                new_node_index = len(nodes_df.index) - 1
                
                for child_nd in changing_children:
                    "Split edge with deeper_parent and child node into two edges"
                    edges_df = split_edge(edges_df,nd_idx,new_node_index,child_nd)
                
                # Add recombination event if necessary -- never implemented
                # Find recomb breakpoint by extent of parent-child relationships  
                # child1_edges = edges_df.loc[(edges_df['parent'] == nd_idx) & (edges_df['child'] == changing_children[0])]
                # child2_edges = edges_df.loc[(edges_df['parent'] == nd_idx) & (edges_df['child'] == changing_children[1])]
                # if child1_edges.left.min() < child2_edges.left.min():
                #     rec_breakpoint = child1_edges.right.max() # if left of conflicting edges
                # elif child1_edges.left.min() > child2_edges.left.min(): # replace with genome length param
                #     rec_breakpoint = child1_edges.left.min() # if right of conflicting edges
                # else:
                #     # Not sure if this can ever happen
                #     print('Could not localize recomb event based on child-recent_parent start/end')
                
                # add two rec nodes to nodes_df
                # min_rec_time = max(nodes_df.at[changing_children[0],'time'],nodes_df.at[changing_children[1],'time'])
                # time = nodes_df.at[nd_idx,'time'] - (nodes_df.at[nd_idx,'time'] - min_rec_time)/2 # can just pick a random time between parent and child
                # new_rec_node = {'flags': 131072, 'population': -1, 'individual': -1, 'time':time, 'metadata':time} # add extra edge for old parent and split node
                # nodes_df = nodes_df.append(new_rec_node, ignore_index = True)
                # left_recomb_parent = len(nodes_df.index) - 1
                # nodes_df = nodes_df.append(new_rec_node, ignore_index = True)
                # right_recomb_parent = len(nodes_df.index) - 1
                
                # for child_nd in changing_children:
                #     edges_df = split_edge_on_rec(edges_df,nd_idx,child_nd,left_recomb_parent,right_recomb_parent,rec_breakpoint)
                    
                # Resort nodes_df in chronological order an re-index parents/children in edges_df
                nodes_df, edges_df = resort_dfs(tip_count, nodes_df, edges_df)
                
        nd_idx -= 1
        
    return nodes_df, edges_df

def tree2nodesdf(tree, prev_node_df=pd.DataFrame()):

    """
        Converts tree into nodes dataframe compatible with ts.tables.nodes.
        Unique node ids are assigned to each node based on edge_bitmasks for coalescent nodes 
        or a cumulative negative counter for recombination events.
        Also checks if unique id already exists in previous nodes_df
        If it does, checks node height equality and assigns new unique id to nodes with different heights
    """
    
    # Create empty nodes_df compatible with ts.tables.nodes
    nodes_df = pd.DataFrame(columns=['flags', 'time', 'population', 'metadata','unique_id'])
    
    # Encode bipartitions -- used fo unique node ids
    tree.encode_bipartitions(suppress_unifurcations=False)
    
    # Calculate node ages in tree
    tree.calc_node_ages(ultrametricity_precision=False)
    
    # Set negative counter for unique recombination node ids
    if not prev_node_df.empty:
        min_id = prev_node_df.unique_id.min()
        if min_id < 0:
            recomb_id = min_id - 1    
        else:
            recomb_id = -1           
    else:
        recomb_id = -1
    
    # Process each node in reverse chronological order
    for node in tree.ageorder_node_iter():
        
        children = list(node.child_nodes())
        unique_id = node.edge.split_bitmask

        if len(children) == 0: # node is a sampled tip
            flags = tskit.NODE_IS_SAMPLE
            nodes_df = nodes_df.append({'flags':flags, 'time':node.age, 'population':-1, 'metadata':int(node.taxon.label), 'unique_id':unique_id},ignore_index=True)
            node.unique_id = unique_id
        elif len(children) == 1: # node is recombination event
            flags = 131072 #recomb_flag
            nodes_df = nodes_df.append({'flags':flags, 'time':node.age, 'population':-1, 'metadata':None, 'unique_id':recomb_id},ignore_index=True)
            node.unique_id = recomb_id
            recomb_id -= 1 # decrement recombinant node unique id counter
        else: # Coalescent node
            flags = 0
            if not prev_node_df.empty:
                
                # Check if duplicate node with same unique_id but different coal time exists
                duplicate = prev_node_df.loc[prev_node_df['unique_id'] == unique_id]
                if not duplicate.empty:
                    if duplicate.iloc[0].time != node.age:
                        # Assign new unique id for duplicate nodes with different coal times
                        max_unique_id = prev_node_df.unique_id.max()
                        unique_id = int(max_unique_id + unique_id)
                        
                # Check if coal event with same time exists
                # TODO: Check if this causes problem if a coal event with the same time exists but is not topologically equivalent
                coal_time_duplicates = prev_node_df[(prev_node_df['flags']==0) & (prev_node_df['time']==node.age)]
                if not coal_time_duplicates.empty:
                    unique_id = coal_time_duplicates.iloc[0].unique_id
                
            nodes_df = nodes_df.append({'flags':flags, 'time':node.age, 'population':-1, 'metadata':None, 'unique_id':unique_id},ignore_index=True)
            node.unique_id = unique_id
    
    return nodes_df, tree                 

def tree2edgedf(tree,left,right,id_dict):
    
    """
        Converts tree into edges dataframe compatible with ts.tables.edges.
    """   
    
    edges_df = pd.DataFrame(columns=['left', 'right', 'parent', 'child','parent_unique_id','child_unique_id'])
    for node in tree.ageorder_node_iter():
        children = list(node.child_nodes())
        parent_id = id_dict[node.unique_id]
        for child in children:
            child_id = id_dict[child.unique_id]
            edges_df = edges_df.append({'left':left, 'right':right, 'parent':parent_id, 'child':child_id, 'parent_unique_id':node.unique_id, 'child_unique_id':child.unique_id}, ignore_index=True)
                    
    return edges_df     


def reindex_edgedf(edge_df,id_dict):
    
    """
        Reindex parent/child ids in edge table based on new id_dict
        id_dict maps unique_ids -> updated node ids
    """

    for index, row in edge_df.iterrows():
        row.parent = id_dict[row.parent_unique_id]
        row.child = id_dict[row.child_unique_id]
        
    return edge_df


def df2TreeTables(edge_df,node_df,total_length):
    
    """
        Convert pandas dataframe to tskit tree tables
    """
    
    tables = tskit.TableCollection(total_length)
    for index, row in edge_df.iterrows():
        tables.edges.add_row(row['left'], row['right'], int(row['parent']), int(row['child']))  
    for index, row in node_df.iterrows():
        tables.nodes.add_row(flags=int(row['flags']), time=row['time'], population=int(row['population'])) 
        #tables.nodes.add_row(flags=int(row['flags']), time=row['time'], population=int(row['population']), metadata=row['metadata'])        
    tables.sort()
    
    return tables


def treeTables2df(tables):
    
    """
        Convert tskit tree tables to pandas dataframe
    """
    
    "As pandas dataframe"
    edges_dict = {'left':tables.edges.left,
                  'right':tables.edges.right,
                  'parent':tables.edges.parent,
                  'child':tables.edges.child}
    edges_df = pd.DataFrame.from_dict(edges_dict)
    
    nodes_dict = {'flags':tables.nodes.flags,
                  'population':tables.nodes.population,
                  'individual':tables.nodes.individual,
                  'time':tables.nodes.time,
                  'metadata':tables.nodes.time}
    nodes_df = pd.DataFrame.from_dict(nodes_dict)
    
    return edges_df, nodes_df

def check_time_constraints(edge_df,node_df):
    
    """
        Check to make sure parent/child nodes are in chronological order in node/edge dfs
        TODO: Not currently used -- can remove in future
    """
    
    for index, row in edge_df.iterrows():
        parent = int(row['parent'])
        child = int(row['child'])
        parent_time = node_df.iloc[parent].time
        child_time = node_df.iloc[child].time
        if child_time > parent_time:
            print("WARNING: child node time is greater than parent node time for child " + str(child) + " and parent " + str(parent))

def check_parent_child_intervals(edge_df,node_df):
    
    """
        Check that each child only has one parent over each genomic interval
    """

    breakpoints = list(set(edge_df.left.values.tolist() + edge_df.right.values.tolist()))
    breakpoints.sort()
    for bp_idx in range(len(breakpoints)-1):
        start_loc = breakpoints[bp_idx]
        end_loc = breakpoints[bp_idx+1]
        for node_idx, node in node_df.iterrows(): 
            interval_parents = list(set(edge_df.loc[(edge_df['left'] <= start_loc) & (edge_df['right'] >= end_loc) & (edge_df['child'] == node_idx), 'parent'].tolist()))
            if len(interval_parents) > 1:
                logging.debug("Node " + str(node_idx) + "has more than two parents over interval: " + str(start_loc) + ":" + str(end_loc))


def split_edge(edges_df,parent_nd,split_nd,child_nd):
    
    """
        Split all edges in edges_df between parent_nd and child_nd by adding split_nd between them
    """
    
    edges_to_split = edges_df.loc[(edges_df['parent'] == parent_nd) & (edges_df['child'] == child_nd)]
    for edge_idx, edge in edges_to_split.iterrows(): 
    
        edges_df.loc[edge_idx,'parent'] = split_nd # replace old parent with split node
        new_parent_edge = {'left': edge.left, 'right': edge.right, 'parent': parent_nd, 'child': split_nd} # add extra edge for old parent and split node
        edges_df = edges_df.append(new_parent_edge, ignore_index = True)
    
    edges_df = edges_df.astype({"parent": int, "child": int}) # wish we didn't have to do this!
    return edges_df

def split_edge_on_rec(edges_df,parent_nd,child_nd,left_recomb_parent,right_recomb_parent,rec_breakpoint):
    
    """
        Split all edges in edges_df between parent_nd and child_nd by adding additional recombinant node (split_nd) between them
        Whether left or right recombinant parent is used depends on whether edge is to left/right of breakpoint
    """
    
    edges_to_split = edges_df.loc[(edges_df['parent'] == parent_nd) & (edges_df['child'] == child_nd)]
    for edge_idx, edge in edges_to_split.iterrows(): 
    
        if edge.right <= rec_breakpoint:
            rec_node = left_recomb_parent
        else:
            rec_node = right_recomb_parent    
    
        edges_df.loc[edge_idx,'parent'] = rec_node # replace old parent with split node
        new_parent_edge = {'left': edge.left, 'right': edge.right, 'parent': parent_nd, 'child': rec_node} # add extra edge for old parent and split node
        edges_df = edges_df.append(new_parent_edge, ignore_index = True)
    
    edges_df = edges_df.astype({"parent": int, "child": int}) # wish we didn't have to do this!
    return edges_df

def check_conflicting_edges(edges_df,deeper_parents,recent_parent,child_nd):
    
    """
        Check if adding recent_parent creates a conflicting parent-child edge
        Such that the child has two parents in same genomic interval.
        If it does, then we need to add another recombination event to reconcile the difference
    """

    conflict = False
    for parent_nd in deeper_parents:
        edges_to_split = edges_df.loc[(edges_df['parent'] == parent_nd) & (edges_df['child'] == child_nd)]
        recent_parent_edges = edges_df.index[edges_df['child'] == recent_parent]
        for ed_idx, ed in edges_to_split.iterrows():
            for ped_idx in recent_parent_edges:
                overlap = test_overlap(edges_df,ed_idx,ped_idx)
                if overlap:
                    conflict = True
                    
    return conflict

def test_overlap(edges_df,idx1,idx2):
    
    """
        Test if the span of two edges overlap over any genomic interval
    """
    
    overlap = False
    if edges_df.at[idx1,'left'] <= edges_df.at[idx2,'left']:
        left_idx = idx1
        right_idx = idx2
    else:
        left_idx = idx2
        right_idx = idx1
    if edges_df.at[left_idx,'right'] > edges_df.at[right_idx,'left']: # if overlap
        overlap = True
        
    return overlap

def resort_dfs(tip_count, nodes_df, edges_df):
    
    """
        Resort nodes_df by time and then reindex parent/child nodes in edges_df 
        But do not sort tip nodes or else tip taxon labels will change
    """

    tip_nodes_df = nodes_df.iloc[0:tip_count]
    sorted_nodes_df = nodes_df.iloc[tip_count:len(nodes_df.index)]
    sorted_nodes_df = sorted_nodes_df.sort_values('time', ignore_index=False)
    old_indexes = list(range(tip_count)) + sorted_nodes_df.index.tolist()
    nodes_df = pd.concat([tip_nodes_df, sorted_nodes_df], ignore_index=True)
    id_dict = dict(zip(old_indexes, nodes_df.index))
    for index, row in edges_df.iterrows():
        edges_df.loc[index,'parent'] = id_dict[row.parent]
        edges_df.loc[index,'child'] = id_dict[row.child]
    
    return nodes_df, edges_df

def merged_overlapping_edges(edges_df):
    
    """
        Check for and merge overlapping edges
        This just sets overlapping edges to have same left/right
        Still need to drop duplicates to revome redundant info
        
        TODO: Not currently used -- can remove in future
    """
    
    for index, edge in edges_df.iterrows():
        parent = edge['parent']
        child = edge['child']
        equiv_edges = edges_df.index[(edges_df['parent'] == parent) & (edges_df['child'] == child)]
        for equiv_idx in equiv_edges:
            if equiv_idx != index:
                "Check for overlap"                     
                if edges_df.at[index,'left'] <= edges_df.at[equiv_idx,'left']:
                    left_idx = index
                    right_idx = equiv_idx
                else:
                    left_idx = equiv_idx
                    right_idx = index
                if edges_df.at[left_idx,'right'] > edges_df.at[right_idx,'left']: # if overlap
                    "Set left/right same for both edges so one will be dropped as a duplicate"
                    #edges_df.loc[left_idx,'left'] = edges_df.at[left_idx,'left']
                    edges_df.loc[right_idx,'left'] = edges_df.at[left_idx,'left']
                    max_right = max(edges_df.at[left_idx,'right'], edges_df.at[right_idx,'right'])
                    edges_df.loc[left_idx,'right'] = max_right
                    edges_df.loc[right_idx,'right'] = max_right
                    
    return edges_df
