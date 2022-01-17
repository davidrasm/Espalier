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
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import re
from Espalier.viz import balticmod as bt

def plot(tree_files,fig_name,**kwargs):
    
    """
        
        Plot trees in tree_files as a tanglegram
        
        Parameters:     
            tree_files (list): list of newick tree files containing local trees
            fig_name (str): figure name for output png file
              
        Optional keyword arguments:
            tree_labels (list) : Labels corresponding to trees files.
            numerical_taxa_names (boolean): set True if taxa names are numeric
            cmap (matplotlib colormap): colormap
            dispalce_scalar (float): scales displacement between neighboring trees based on max tree height 
            branch_width (float): branch thinkness/width in plotted trees
            xax_border (float): size of border on x-axis
            yax_border (float): size of border on y-axis
            tip_font_size (float): tip label font size
    
    """

    tree_labels = kwargs.get('tree_labels', None)
    numerical_taxa_names = kwargs.get('numerical_taxa_names', False) 

    # Can fiddle with these params to improve fig appearance    
    cmap = kwargs.get('cmap', mpl.cm.viridis) # mpl.cm.Spectral is nice too
    dispalce_scalar = kwargs.get('dispalce_scalar', 0.5)
    branch_width = kwargs.get('branch_width', 2) 
    xax_border = kwargs.get('xax_border', 0.1) # should be relative to cumulative-displace
    yax_border = kwargs.get('yax_border', 2) # should be relative to ySpan of trees
    tip_font_size = kwargs.get('tip_font_size', 12)
    
    # Load trees into tree dict
    trees={}
    segments = list(range(len(tree_files)))
    for idx,tr in enumerate(tree_files):
        output_tree = tr.replace('.tre','.nexus')
        convert2nexus(tr,output_tree,numerical_taxa_names)
        ll=bt.loadNexus(output_tree,absoluteTime=False)
        #ll.setAbsoluteTime(2020.0)
        trees[idx]=ll
    print('\nDone!')
    
    # Rescale tree heights so they are all equal
    tree_heights = []
    for t,tr in enumerate(trees.keys()): ## iterate over trees
        cur_tree=trees[tr] ## fetch tree object
        tree_heights.append(cur_tree.treeHeight)
    max_height_cap = max(tree_heights)
    for t,tr in enumerate(trees.keys()): ## iterate over trees
        cur_tree=trees[tr] ## fetch tree object
        for k in cur_tree.Objects: ## iterate over a flat list of branches
            k.length = k.length * (max_height_cap/cur_tree.treeHeight)
        cur_tree.traverse_tree() ## required to set heights
        cur_tree.treeStats() ## report stats about tree
    
    # Compute displaceAmount on the same scale as the tree heights
    displaceAmount= max_height_cap * dispalce_scalar
    
    # Extract tip positions
    tip_positions={x:{} for x in segments} ## remember the position of each tip in each tree
    for t,tr in enumerate(trees.keys()): ## iterate over trees
        cur_tree=trees[tr] ## fetch tree object
        for k in cur_tree.Objects:
            if k.branchType=='leaf':
                tip_positions[tr][k.name]=(k.height,k.y) ## remember (X, Y) position of tip
    
    for X in range(10): ## 10 untangling iterations
        print('iteration %d'%(X+1))
        for t,tr in enumerate(segments): ## iterate over each tree
            print(tr)
            ptr=segments[t-1] ## previous tree
            ntr=segments[t] ## next tree
            seg=trees[ptr] ## fetch appropriate tree
            nex_seg=trees[ntr]
            for k in sorted(nex_seg.Objects,key=lambda q:q.height): ## iterate over branches from most recent to oldest
                if k.branchType=='node': ## can only sort nodes
                    leaves=[[seg.tipMap[tip] for tip in w.leaves if tip in seg.tipMap] if w.branchType=='node' else [w.name] for w in k.children] ## descendent tips in current order
                    
                    for c in range(len(leaves)):
                        leaves[c]=sorted(leaves[c],key=lambda x:tip_positions[ntr][x][1] if x in tip_positions[ntr] else 0.0) ## sort leaves according to their positions in the next tree
                    
                    ys=[sorted([tip_positions[ntr][w][1] for w in cl if w in tip_positions[ntr]]) for cl in leaves] ## extract y positions of descendents
                    merge_ys=sum(ys,[]) ## flatten list of tip y coordinates
                    ypos=range(min(merge_ys),max(merge_ys)+1) ## get y positions of tips in current order
                    order={i:x for i,x in enumerate(leaves)} ## dict of tip order: tip name
                    
                    new_order=sorted(order.keys(),key=lambda x:-np.mean([(tip_positions[ptr][order[x][w]][1]-ypos[w]) for w in range(min([len(order[x]),len(ypos)])) if order[x][w] in tip_positions[ptr]])) ## get new order by sorting existing order based on y position differences
                    if new_order!=range(len(leaves)): ## if new order is not current order
                        k.children=[k.children[i] for i in new_order] ## assign new order of child branches
                        nex_seg.drawTree() ## update y positions
    
                        for w in nex_seg.Objects: ## iterate over objects in next tree
                            if w.branchType=='leaf':
                                tip_positions[ntr][w.name]=(w.height,w.y) ## remember new positions
                    
            if t==0: ## if first tree
                trees[segments[t]].drawTree() ## update positions
                lvs=sorted([w for w in trees[segments[t]].Objects if w.branchType=='leaf'],key=lambda x:x.y) ## get leaves in y position order
                
                norm=mpl.colors.Normalize(0,len(lvs))
                pos_colours={w.name:cmap(norm(w.y)) for w in lvs} ## assign colour
                
    
    # Plotting all trees
    fig,ax = plt.subplots(figsize=(12,8),facecolor='w')
    cumulative_displace=0 ## this tracks the "current" x position, so trees are plotted one after another
    tree_names = segments
    ref_tree = segments[0] # first tree
    
    for t,tr in enumerate(tree_names): ## iterate over trees
        cur_tree=trees[tr] ## fetch tree object
        
        x_attr=lambda k: k.height+cumulative_displace
        #x_attr=lambda k: (k.height*(max_height/cur_tree.treeHeight))+cumulative_displace
        
        b_func=lambda k: branch_width 
        s_func=lambda k: 30
        su_func=lambda k: 60
        ct_func=lambda k: cmap(tip_positions[ref_tree][k.name][1]/float(cur_tree.ySpan))
        cu_func=lambda k: 'k'
        z_func=lambda k: 100
        zu_func=lambda k: 99
        
        # For tip naming
        text_func = lambda k: k.name.replace('_',' ')
        target_func = lambda k: k.is_leaf()
        position_func = lambda k: (k.height+cumulative_displace+0.2, k.y)
        
        def colour_func(node):
            #if traitName in node.traits:
            #    return 'indianred' if node.traits[traitName]=='V' else 'steelblue'
            #else:
                return 'k'
            
        cn_func=colour_func
        
        cur_tree.plotTree(ax,x_attr=x_attr,branchWidth=b_func,colour_function=cn_func)
        cur_tree.plotPoints(ax,x_attr=x_attr,size_function=s_func,colour_function=ct_func,zorder_function=z_func)
        cur_tree.plotPoints(ax,x_attr=x_attr,size_function=su_func,colour_function=cu_func,zorder_function=zu_func)
        
        # Add tip label if at last tree
        if t == len(tree_names) - 1: # last_tree
            cur_tree.addText(ax, text=text_func, position=position_func,fontsize=tip_font_size)
        
        for k in cur_tree.Objects: ## iterate over branches
            if isinstance(k,bt.leaf): ## if leaf...
                y=k.y
                pos_in_first_tree=tip_positions[ref_tree][k.name][1] ## fetch y coordinate of same tip in the first tree
                frac_pos=pos_in_first_tree/float(cur_tree.ySpan) ## normalize coordinate to be within interval [0.0,1.0]
    
                if t!=len(tree_names)-1: ## as long as we're not at the last tree - connect tips with coloured lines
                    next_x,next_y=tip_positions[tree_names[t+1]][k.name] ## fetch coordinates of same tip in next tree
                    next_x+=cumulative_displace+cur_tree.treeHeight+displaceAmount ## adjust x coordinate by current displacement and future displacement
                    nextIncrement=cumulative_displace+cur_tree.treeHeight
                    ax.plot([x_attr(k),nextIncrement+0.05*displaceAmount,nextIncrement+0.95*displaceAmount,next_x],[y,y,next_y,next_y],lw=1,ls='-',color=cmap(frac_pos),zorder=0) ## connect current tip with same tip in the next tree
        
        if tree_labels:
            add_tree_label(ax,cur_tree,tree_labels[tr],cumulative_displace)
        
        cumulative_displace+=cur_tree.treeHeight+displaceAmount ## increment displacement by the height of the tree
    
    [ax.spines[loc].set_visible(False) for loc in ['top','right','left','bottom']]
    
    ax.tick_params(axis='x',size=0)
    ax.tick_params(axis='y',size=0)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    
    ax.set_ylim(-yax_border,cur_tree.ySpan+yax_border) ## set y limits
    ax.set_xlim(-xax_border,cumulative_displace+xax_border)
    
    plt.savefig(fig_name, dpi=300)
    
def convert2nexus(in_tree,out_tree,numerical_taxa_names=False):
    
    """
        Convert newick file to nexus for tanglegram plotting
    """
    
    myTree=bt.loadNewick(in_tree, absoluteTime=False)
    myTree.traverse_tree() ## required to set heights
    myTree.treeStats() ## report stats about tree
    names = []
    for idx,k in enumerate(myTree.Objects): ## iterate over a flat list of branches
        if k.branchType=='leaf':
            curr_name = k.numName
            names.append(curr_name)
    
    date_str = '' #'_2020.00'
    
    # Write taxa names
    nex=open(out_tree,"w")
    nex.write("#NEXUS\n")
    nex.write("Begin taxa;\n")
    nex.write("\tDimensions ntax=" + str(len(names)) + ";\n")	
    nex.write("\t\tTaxlabels\n")
    for n in names:
        nex.write("\t\t\t" + n + date_str + "\n")    
    nex.write("\t\t\t;\n")
    nex.write("End;\n")		
    
    # Write translation 	
    nex.write("Begin trees;\n")	
    nex.write("\tTranslate\n")	
    for idx,n in enumerate(names):
        if numerical_taxa_names:
            nex.write("\t\t" + n + ' ' + n + date_str + "\n") # if taxa names are numbers
        else:
            nex.write("\t\t" + str(idx+1) + ' ' + n + date_str + "\n") # if taxa names are non-numerical strings
       
    nex.write(";\n")
    
    # Write tree
    with open(in_tree, 'r') as file:
        tree_str = file.read().replace('\n', '')
    if not numerical_taxa_names:
        for idx,n in enumerate(names):
            tree_str = re.sub(n, str(idx+1), tree_str) # if taxa names are non-numerical strings    
    nex.write("tree TREE1 = " + tree_str + "\n")
    nex.write("End;\n")

def add_tree_label(ax,tree,label_str,cumulative_displace):
    
    """
        Add a label to tree
    """
    
    curr_min_x = np.Inf
    curr_max_x = -np.Inf
    curr_min_y = np.Inf
    curr_max_y = -np.Inf
    for k in tree.Objects:
        if k.x > curr_max_x:
            curr_max_x = k.x
        if k.x < curr_min_x:
            curr_min_x = k.x
        if k.y > curr_max_y:
            curr_max_y = k.y
        if k.y < curr_min_y:
            curr_min_y = k.y
    x_text_pos = cumulative_displace + (curr_max_x - curr_min_x) / 2
    y_text_pos = curr_max_y + (curr_max_y - curr_min_y) * 0.05
    ax.text(x_text_pos,y_text_pos,label_str,horizontalalignment='center',fontsize=8)

# def plot_pair(tree_file1,tree_file2,fig_name,tree_labels=None):
    
#     """
#         Never finished implementing
#         Seems like things have changed in new version of baltic
#         Maybe just use bit about inverting x_attributes
#     """

#     "Params to vary to improve figure"
#     dispalce_scalar = 0.5 #0.2 # sets displacement between neighboring trees based on fraction of max tree height
#     branch_width = 2 # branch thinkness/width in plotted trees (default is 4)
#     xax_border = 0.1 # was 0.2 -- should be relative to cumulative-displace
#     yax_border = 2 # was 0.2 -- should be relative to ySpan of trees
    
#     "Load trees into tree dict"
#     trees={} ## dict
#     segments = [0,1] #list(range(len(tree_files)))
    
#     output_tree = tree_file1.replace('.tre','.nexus')
#     convert2nexus(tree_file1,output_tree)
#     trees[0]=bt.loadNexus(output_tree,absoluteTime=False)
    
#     output_tree = tree_file2.replace('.tre','.nexus')
#     convert2nexus(tree_file2,output_tree)
#     trees[1]=bt.loadNexus(output_tree,absoluteTime=False)
    
#     tip_positions={x:{} for x in segments} ## remember the position of each tip in each tree
    
#     for t,tr in enumerate(trees.keys()): ## iterate over trees
#         cur_tree=trees[tr] ## fetch tree object
#         for k in cur_tree.Objects:
#             if k.branchType=='leaf':
#                 tip_positions[tr][k.name]=(k.height,k.y) ## remember (X, Y) position of tip
    
#     cmap=mpl.cm.Spectral
    
#     tip_positions={x:{} for x in trees} ## remember the position of each tip in each tree
    
#     for t,tr in enumerate(trees.keys()): ## iterate over trees
#         cur_tree=trees[tr] ## fetch tree object
#         for k in cur_tree.Objects:
#             if k.branchType=='leaf':
#                 tip_positions[tr][k.name]=(k.height,k.y) ## remember (X, Y) position of tip
    
#     cmap=mpl.cm.Spectral
    
#     for X in range(10): ## 10 untangling iterations
#         print('iteration %d'%(X+1))
#         for t,tr in enumerate(segments): ## iterate over each tree
#             print(tr)
#             ptr=segments[t-1] ## previous tree
#             ntr=segments[t] ## next tree
#             seg=trees[ptr] ## fetch appropriate tree
#             nex_seg=trees[ntr]
#             for k in sorted(nex_seg.Objects,key=lambda q:q.height): ## iterate over branches from most recent to oldest
#                 if k.branchType=='node': ## can only sort nodes
#                     leaves=[[seg.tipMap[tip] for tip in w.leaves if tip in seg.tipMap] if w.branchType=='node' else [w.name] for w in k.children] ## descendent tips in current order
                    
#     #                 leaves=[[seg.tipMap[tip] for tip in w.leaves] if w.branchType=='node' else [w.name] for w in k.children] ## descendent tips in current order
                    
#                     for c in range(len(leaves)):
#     #                     leaves[c]=sorted(leaves[c],key=lambda x:tip_positions[ntr][x][1]) ## sort leaves according to their positions in the next tree
#                         leaves[c]=sorted(leaves[c],key=lambda x:tip_positions[ntr][x][1] if x in tip_positions[ntr] else 0.0) ## sort leaves according to their positions in the next tree
                    
#                     ys=[sorted([tip_positions[ntr][w][1] for w in cl if w in tip_positions[ntr]]) for cl in leaves] ## extract y positions of descendents
#                     merge_ys=sum(ys,[]) ## flatten list of tip y coordinates
#                     ypos=range(min(merge_ys),max(merge_ys)+1) ## get y positions of tips in current order
#                     order={i:x for i,x in enumerate(leaves)} ## dict of tip order: tip name
                    
#                     new_order=sorted(order.keys(),key=lambda x:-np.mean([(tip_positions[ptr][order[x][w]][1]-ypos[w]) for w in range(min([len(order[x]),len(ypos)])) if order[x][w] in tip_positions[ptr]])) ## get new order by sorting existing order based on y position differences
                    
#     #                 new_order=sorted(order.keys(),key=lambda x:-np.mean([(tip_positions[ptr][order[x][w]][1]-ypos[w]) for w in range(len(order[x]))])) ## get new order by sorting existing order based on y position differences
                    
#                     if new_order!=range(len(leaves)): ## if new order is not current order
#                         k.children=[k.children[i] for i in new_order] ## assign new order of child branches
#                         nex_seg.drawTree() ## update y positions
    
#                         for w in nex_seg.Objects: ## iterate over objects in next tree
#                             if w.branchType=='leaf':
#                                 tip_positions[ntr][w.name]=(w.height,w.y) ## remember new positions
                    
#             if t==0: ## if first tree
#                 trees[segments[t]].drawTree() ## update positions
#                 lvs=sorted([w for w in trees[segments[t]].Objects if w.branchType=='leaf'],key=lambda x:x.y) ## get leaves in y position order
                
#                 norm=mpl.colors.Normalize(0,len(lvs))
#                 pos_colours={w.name:cmap(norm(w.y)) for w in lvs} ## assign colour
    
#     fig = plt.subplots(figsize=(20,10),facecolor='w')
    
#     gs = GridSpec(1, 1,hspace=0.0,wspace=0.0)
#     ax = plt.subplot(gs[0])
    
#     tree1=trees[0]
#     tree2=trees[1]
    
#     cmap=mpl.cm.cividis
#     x_attr=lambda k: k.height ## branch x position is determined by height
#     c_func=lambda k:'k' #'indianred' if k.traits[traitName]=='V' else 'steelblue'
#     ct_func=lambda k: cmap(k.y/float(tree1.ySpan)) ## call colour map with fraction that represents the y position of a tip (returns colour)
    
    
#     #cur_tree.plotTree(ax,x_attr=x_attr,branchWidth=b_func,colour_function=cn_func)
#     #cur_tree.plotPoints(ax,x_attr=x_attr,size_function=s_func,colour_function=ct_func,zorder_function=z_func)
#     #cur_tree.plotPoints(ax,x_attr=x_attr,size_function=su_func,colour_function=cu_func,zorder_function=zu_func)
    
#     b_func=lambda k: branch_width
#     tree1.plotTree(ax,x_attr=x_attr,branchWidth=b_func,colour_function=c_func) ## plot black tree
#     tree1.plotPoints(ax,x_attr=x_attr,size=50,colour_function=ct_func,zorder=100) ## plot circles at tips
    
    
#     skip=tree1.treeHeight*0.35 ## skip this many units between trees
#     x_attr=lambda k: tree1.treeHeight+skip+tree2.treeHeight-k.height ## for tree2 we'll offset x coordinates by the height of the tree and invert branches
#     tree2.plotTree(ax,x_attr=x_attr,branchWidth=b_func,colour_function=c_func) ## plot black tree
#     tree2.plotPoints(ax,x_attr=x_attr,size=50,colour_function=ct_func,zorder=100) ## plot circles at tips
    
#     for k in filter(lambda x: x.branchType=='leaf',tree1.Objects): ## grab leaf objects in tree1
#         x=k.height ## get height
#         y=k.y ## get y position
        
#         matching_tip=tree2.getBranches(lambda x: x.branchType=='leaf' and x.name==k.name) ## fetch corresponding branch in tree2
#         match_y=matching_tip.y
#         xs=[x,tree1.treeHeight+0.15*skip,tree1.treeHeight+skip-0.15*skip,x_attr(matching_tip)] ## x coordinates for tangleline
#         ys=[y,y,match_y,match_y] ## y coordinates for tangleline
#         ax.plot(xs,ys,color=cmap(y/float(tree1.ySpan))) ## plot tangleline
        
#     [ax.spines[loc].set_visible(False) for loc in ['top','right','left','bottom']]
    
#     ax.tick_params(axis='x',size=0)
#     ax.tick_params(axis='y',size=0)
#     ax.set_xticklabels([])
#     ax.set_yticklabels([])
    
#     ax.set_ylim(-1,tree1.ySpan+1) ## set y limits
#     ax.set_xlim(-5,tree1.treeHeight+skip+tree2.treeHeight+5)
    
#     plt.show()


if  __name__ == '__main__':
    
    "Test plotting pairs of tree face-to-face"
    #tree_file1 = './examples/reconciler_MLTree_r1.tre'
    #tree_file2 = './examples/reconciler_MLTree_r2.tre'
    #tanglegram_fig_name = 'tanglegram_face-to-face-trees.png'
    #plot_pair(tree_file1,tree_file2,tanglegram_fig_name,tree_labels=None)
    
    "Test with multiple trees"
    segments = 4
    path = "../sim_trees/"
    ML_tree_files = [path + "disentangler_test0_MLTree" + str(i) + ".tre" for i in range(segments)]
    tanglegram_fig_name = 'tanglegram-ML-trees.png' 
    
    tree_labels = ['Tree ' + str(s) for s in range(1,segments+1)]
    
    cmap=mpl.cm.Spectral
    plot(ML_tree_files, tanglegram_fig_name,tree_labels=tree_labels, numerical_taxa_names=True, cmap=cmap, tip_font_size=10)

