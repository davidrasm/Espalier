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

from Espalier import TreeOps
import numpy as np
from Bio.SeqRecord import SeqRecord
import logging
from scipy.stats import chi2
import pandas as pd

# Build an exact-id lookup and tolerate DendroPy's underscore/space conversions.
def _build_seq_lookup(records, seq_file):
    seq_dict = {}
    for record in records:
        key = record.id.strip()
        normalized_record = SeqRecord(record.seq, id=key, name=key, description=key)
        aliases = {
            key,
            key.replace("_", " "),
        }
        for alias in aliases:
            existing = seq_dict.get(alias)
            if existing is not None and existing.id != normalized_record.id:
                raise ValueError(
                    f"Duplicate sequence key '{alias}' while loading {seq_file}"
                )
            seq_dict[alias] = normalized_record
    return seq_dict

class SPRTester(object):
    
    """
        A callable class used to run likelihood ratio test to identify
        statistically well supported recombination/reassortment events
        from the SPR events identified in a MAF
    """
    
    def __init__(self,raxml,lower_bound_ratio=0.1,**kwargs):
        
        '''             
            Parameters: 
                raxml (RAxML): instance of RAxML runner used to compute seq likelihoods
                #lower_bound_ratio (float): lower bound ratio used by branch-and-bound algorithm
               
            Optional keyword arguments:
                temp_dir (str): directory to store temp output
        '''
        
        self.raxml = raxml
        self.lower_bound_ratio = lower_bound_ratio        
        #self.prior_gamma = prior_gamma
        self.temp_dir = kwargs.get('temp_dir', './')
    
    def __call__(self,ref_tree,ref_seq_file,alt_tree,alt_seq_file,maf):
    
       
    
        """
        
            Run likelihood ratio test to identify statistically well supported recombination/reassortment events
            from the SPR events identified in a MAF
        
            Parameters:     
               ref_tree (dendropy.Tree): input ref tree (T1)
               ref_seq_file (str): fasta sequence file containing alignment for ref tree (T1)
               alt_tree (dendropy.Tree): input alt tree (T2)
               alt_seq_file (str): fasta sequence file containing alignment for alt tree (T2)
               maf (dendropy.TreeList): maximum agreement forest between ref and alt
               
            Optional keyword arguments:
               
            Returns:     
               results_df (pandas.DataFrame): Dataframe with test results
        """
                        
        # Write T1/ref tree file
        T1_file = self.temp_dir + 'temp_T1.tre'
        ref_tree.write(path=T1_file,schema='newick',suppress_annotations=True,suppress_rooting=True)
            
        # Write T2/alt tree file
        T2_file = self.temp_dir + 'temp_T2.tre'
        alt_tree.write(path=T2_file,schema='newick',suppress_annotations=True,suppress_rooting=True)
        
        "Get and set best model found by RAxML"
        #self.raxml.model = self.raxml.get_best_model(T1_file,ref_seq_file)
            
        "Precompute likelihoods for L_S1_T1 and L_S2_T2"
        L_S1_T1 = self.raxml.get_tree_likelihood(T1_file,ref_seq_file,opt_model=True,opt_branches=True)
        L_S2_T2 = self.raxml.get_tree_likelihood(T2_file,alt_seq_file,opt_model=True,opt_branches=True)
        
        # Reverse order MAF subtrees so we start with final connected component
        maf = list(reversed(maf))
        
        subtree_count = len(maf)
        
        "Init lists for storing results"
        subtree_idx_list = []
        subtree_list = []
        logL_no_spr_list = []
        logL_spr_list = []
        logL_ratio_list = []
        lambda_list = []
        p_value_list = []
        
        "Outer loop over all subtrees V in MAF"
        for V_idx in range(1,subtree_count):
            
            logging.info("Processing next SPR in MAF : %s", str(V_idx))
            
            # Print next subtree V for debugging
            V = maf[V_idx]
            subtree_str = V.as_string(schema='newick',suppress_rooting=True,suppress_edge_lengths=True)
            print(subtree_str)
            
            #if len(V.nodes()) > 1:
            #    V.print_plot()
            #else:
            #    print('---' + V.nodes()[0].taxon.label)
            
            # Create trees for alt and ref onto which we'll regraft V
            R_ref = maf[0].clone(depth=2)
            R_alt = maf[0].clone(depth=2)
            
            "Regraft all subtrees U in MAF onto R_ref and R_alt"
            for U_idx in range(1,subtree_count):
                
                U = maf[U_idx]
                
                if U_idx == V_idx:
                    
                    # Find alternative attachment edge in R_ref and R_alt for U
                    logging.debug("Finding attachment edges")
                    attachment_edge_ref = find_attachment_edge(R_ref, alt_tree, U)
                    attachment_edge_alt = find_attachment_edge(R_alt, ref_tree, U)
                    
                    logging.debug("Attaching subtrees")
                    U_clone = U.clone(depth=2)
                    R_ref = attach_subtree(R_ref,U_clone,ref_tree,attachment_edge_ref,midpoint=False)
                    U_clone = U.clone(depth=2) #Should we also deep clone sub_tree?
                    R_alt = attach_subtree(R_alt,U_clone,alt_tree,attachment_edge_alt,midpoint=False)
                    
                else:
                
                    # Find attachment edge in R_ref and R_alt for U
                    logging.debug("Finding attachment edges")
                    attachment_edge_ref = find_attachment_edge(R_ref, ref_tree, U)
                    attachment_edge_alt = find_attachment_edge(R_alt, alt_tree, U)
                    
                    # Attach sub_tree U to R_ref and R_alt by adding a new "graft" node along attachment edge        
                    logging.debug("Attaching subtrees")
                    U_clone = U.clone(depth=2)
                    R_ref = attach_subtree(R_ref,U_clone,ref_tree,attachment_edge_ref,midpoint=False)
                    U_clone = U.clone(depth=2) #Should we also deep clone sub_tree?
                    R_alt = attach_subtree(R_alt,U_clone,alt_tree,attachment_edge_alt,midpoint=False)
        
            # Print regrafted trees for debugging
            #R_ref.print_plot()
            #R_alt.print_plot()
            #print()
                
            """
                Perform likelihood ratio test for SPR with subtree V
                Test is performed with four trees:
                    T1: Original ref/T1 tree
                    T1*: Tree T1 w/ V regrafted onto parent in T2 (R_ref)
                    T2: Original alt/T2 tree
                    T2*: Tree T2 w/ V regrafted onto parent in T1 (R_alt)
                And two alignment files:
                    ref_seq_file = seq file for S1 region
                    alt_seq_file = seq file for S2 region
            """
            
            # Write T1* tree file
            T1_star_file = self.temp_dir + 'temp_T1_star.tre'
            R_ref.write(path=T1_star_file,schema='newick',suppress_annotations=True,suppress_rooting=True)
            
            # Write T2* tree file
            T2_star_file = self.temp_dir + 'temp_T2_star.tre'
            R_alt.write(path=T2_star_file,schema='newick',suppress_annotations=True,suppress_rooting=True)
            
            # Compute likelihoods
            L_S1_T1_star = self.raxml.get_tree_likelihood(T1_star_file,ref_seq_file,opt_model=True,opt_branches=True)
            L_S2_T2_star = self.raxml.get_tree_likelihood(T2_star_file,alt_seq_file,opt_model=True,opt_branches=True)
            
            # Compute likelihood ratio
            logL_no_spr = (L_S1_T1 + L_S2_T2_star + L_S1_T1_star + L_S2_T2) / 2 # average of the two simpler models w/ no SPR move
            logL_spr = L_S1_T1 + L_S2_T2 # more complex model with SPR move
            logL_ratio = logL_no_spr - logL_spr
            
            logging.info("ln(L) no SPR = %s", str(logL_no_spr))
            logging.info("ln(L) SPR = %s", str(logL_spr))
            logging.info("LRT = %s", str(logL_ratio))
            
            # Compute AIC
            #num_params = 0 # num of model params
            #AIC_no_spr = 2*num_params - 2*logL_no_spr
            #AIC_spr =  2*num_params - 2*logL_spr
            #relative_like = np.exp((AIC_spr - AIC_no_spr)/2) # probability that no_spr model minimizes information loss
            #logging.info("Relative likelihood = %s", str(relative_like))
            
            # Test for significance -- reject no SPR model if p_value < alpha
            lambda_lr = -2*(logL_ratio)# likelihood ratio test statistic
            k = 1 # degrees of freedom i.e. the difference in the number of free params between models
            p_value = chi2.sf(lambda_lr,k) # the survival function: 1 - CDF
            logging.info("p-value = %s", str(p_value))
            
            # Append results to lists
            subtree_idx_list.append(V_idx)
            subtree_list.append(subtree_str)
            logL_no_spr_list.append(logL_no_spr)
            logL_spr_list.append(logL_spr)
            logL_ratio_list.append(logL_ratio)
            lambda_list.append(lambda_lr)
            p_value_list.append(p_value)
        
        # Pack results into dataframe
        data_dict = {'Subtree_index': subtree_idx_list,
                     'Subtree': subtree_list,
                     'logL_no_SPR': logL_no_spr_list,
                     'logL_SPR': logL_spr_list,
                     'logL_ratio': logL_ratio_list,
                     'lambda': lambda_list,
                     'p_value': p_value_list}
        results_df = pd.DataFrame(data_dict)
        
        return results_df
        
class ReconciledTree:
    """
        Simple class for storing reconciled trees along with their metadata used for post-processing
    """
    def __init__(self, parent=None, tree=None, index=-1, like=0.0, source=None):
        self.parent = parent
        self.tree = tree
        self.index = index
        self.like = like
        self.source = source

def find_attachment_edge(rec_tree,guide_tree,subtree):
    
    """
    
        Find edge in (reconciled) rec_tree where subtree attaches.
        Equivalent to the subtree's position in guide_tree.
            
        Proceeds in two steps:
            
        Step 1: Find attachment edge in guide_tree
        Attachment edge will be the first edge that includes all subtree_taxa in its leaf set
        
        Step 2: Find edge equivalent to the attachment edge in rec_tree
        Note: this edge will be the parent of the edge leading to the subtree after regrafting
        Hence we find this edge by the parent node's taxa set using edge.tail
        Equivalent edge will have all taxa in attachment_edge_taxa_set that are also taxa in the rec_tree
        But there may be no equivalent edge to attach to if none of the taxa in attachment_edge_taxa_set have been added to rec_tree yet
        So if taxa_set_intersect is null set move to parent of attachment edge until union is not null set
    
        Parameters:     
           rec_tree (dendropy.Tree): partially reconciled tree to which subtree is attached
           guide_tree (dendropy.Tree): full tree used to determine where to attach subtree
           subtree (dendropy.Tree): subtree to attach
           
        Returns:     
            attachment_edge (dendropy.Edge): edge to attach subtree to
    """
    
    # Get taxa set in subtree and rectree we are attaching"
    subtree_taxa = set([lf.taxon.label for lf in subtree.leaf_node_iter()])
    rectree_taxa = set([lf.taxon.label for lf in rec_tree.leaf_node_iter()])
    
    # Step 1: Find attachment edge in guide_tree
    for edge in guide_tree.postorder_edge_iter():
        edge_taxa = set([lf.taxon.label for lf in edge.head_node.leaf_iter()])
        if subtree_taxa.issubset(edge_taxa): # subtree_taxa are subset of all edge_taxa
            break
    if not edge.tail_node:
        print("Error finding attachment edge: edge has no parent/tail")
    
    # Step 2: Find equivalent edge in rec_tree to attach subtree 
    attachment_edge_taxa_set = set([lf.taxon.label for lf in edge.tail_node.leaf_iter()]) # tail_node because we need parent edge of node we will regraft to 
    taxa_set_intersection = attachment_edge_taxa_set & rectree_taxa
    while not taxa_set_intersection: # while intersect is null set
        edge = edge.tail_node.edge # get parent edge
        attachment_edge_taxa_set = set([lf.taxon.label for lf in edge.tail_node.leaf_iter()])
        taxa_set_intersection = attachment_edge_taxa_set & rectree_taxa
    for edge in rec_tree.postorder_edge_iter():
        edge_taxa = set([lf.taxon.label for lf in edge.head_node.leaf_iter()])
        if not taxa_set_intersection.difference(edge_taxa): 
            break # if there is no difference
    attachment_edge = edge # edge subtree will be attached to in rec_tree
    
    return attachment_edge

def attach_subtree(rec_tree,sub_tree,attachment_tree,attachment_edge,midpoint=False):
    
    """
    
        Attach sub_tree to rec_tree along attachment_edge by regrafting.
        Attachment tree (should be called guide tree) is only used to determine height of node subtree will be regrafted to.
        
        Parameters:     
           rec_tree (dendropy.Tree): partially reconciled tree to which subtree is attached
           sub_tree (dendropy.Tree): subtree to attach
           guide_tree (dendropy.Tree): full tree used to determine regraft node height
           attachment_edge (dendropy.Edge): edge to attach subtree to
           
        Optional keyword arguments:
            midpoint = If true, subtrees will be regrafted to node at midpoint of attachment_edge 
            
        Returns:     
            rec_tree (dendropy.Tree): Reconciled tree with subtree regrafted to attachment edge
    """
    
    min_edge_length = 0.0 # minimum branch length
    
    # Create a new 'graft' node to attach sub_tree and sister of sub_tree
    attachment_parent_node = attachment_edge.tail_node
    attachment_sister_node = attachment_edge.head_node
    deep_rooted = False
    if not attachment_parent_node:
        
        # Grafting above root of rec tree such that graft node will be new root
        # This can happen if cut one of the root's child edges
        
        logging.debug("Deep rooting tree")
        
        rec_tree = TreeOps.add_deep_root(rec_tree) # reroot == True by default now
        deep_rooted = True
        attachment_parent_node = attachment_edge.tail_node
        attachment_sister_node = attachment_edge.head_node
    
    # Get temporal constraints on graft node height:
    rec_tree.calc_node_ages(ultrametricity_precision=False)
    attachment_parent_height = attachment_parent_node.age
    max_height =  attachment_parent_height # Graft node must be below attrachment parent node
    sister_subtree_height = attachment_sister_node.age # But above the root node of the grafted sub_tree and sister_subtree
    sub_tree.calc_node_ages(ultrametricity_precision=False)
    subtree_root = sub_tree.seed_node # TreeOps.get_root_node(sub_tree)
    subtree_height = subtree_root.age
    min_height = max(sister_subtree_height,subtree_height)
    
    if attachment_parent_height < sister_subtree_height:
        logging.warning("Parent height is below sister subtree height")
    if attachment_parent_height < subtree_height:
        logging.warning("Parent height is below subtree height")
    
    #if attachment_parent_height < min_height:
        #label_fn = lambda nd: nd.taxon.label if nd.is_leaf() else str(f"{nd.age:0.2f}")
        #rec_tree.print_plot(show_internal_node_labels=True,node_label_compose_fn=label_fn)
        #sub_tree.print_plot(show_internal_node_labels=True,node_label_compose_fn=label_fn)
        #attachment_tree.calc_node_ages(ultrametricity_precision=False)
        #attachment_tree.print_plot(show_internal_node_labels=True,node_label_compose_fn=label_fn)
        #logging.warning("Parent height is below sister or subtree height")

    # Add new graft node to attach sub_tree to
    attachment_edge_length = attachment_sister_node.edge_length
    graft_node = attachment_parent_node.new_child()

    # Alternative way to exact - does not remove dummy 'new_child' node we've grafted above
    #sister_subtree = rec_tree.extract_tree_with_taxa(taxa=set([lf.taxon for lf in attachment_sister_node.leaf_iter()]))
    
    # Extract sister subtree descending from attachment edge and prune it from rec_tree
    sister_subtree = TreeOps.extract_subtree(rec_tree,attachment_sister_node)
    rec_tree.prune_subtree(attachment_sister_node, update_bipartitions=False, suppress_unifurcations=True)

    # Get and set length of graft edge
    if midpoint:
        graft_node.edge.length = attachment_edge_length / 2.0
    else:
        subtree_taxon = sub_tree.leaf_nodes()[0].taxon.label # get one leaf in subtree
        sister_subtree_taxon = sister_subtree.leaf_nodes()[0].taxon.label # get one leaf in sister subtree
        attachment_tree.calc_node_ages(ultrametricity_precision=False)
        mrca = attachment_tree.mrca(taxon_labels=[subtree_taxon,sister_subtree_taxon]) #Should this be tree specific (rec/alt?)
        graft_node_age = mrca.age
        
        #Check constraints on graft node height/age"
        if graft_node_age < min_height or graft_node_age > max_height:
            graft_node_age = min_height + (max_height - min_height)/2.0 # center between max and min
        graft_node.edge.length = attachment_parent_height - graft_node_age #mrca.edge_length

    # Reattach sister subtree to graft node
    if midpoint:
        sister_subtree_root = sister_subtree.seed_node
        sister_subtree_root.edge.length = attachment_edge_length / 2.0
        graft_node.add_child(sister_subtree_root)
    else:
        TreeOps.regraft_subtree(sister_subtree,graft_node,graft_node_age,min_edge_length=min_edge_length)
    
    
    # Attach sub_tree to graft node"
    if midpoint:
        subtree_root = sub_tree.seed_node
        subtree_root.edge.length = attachment_edge_length / 2.0
        graft_node.add_child(subtree_root)
    else:
        TreeOps.regraft_subtree(sub_tree,graft_node,graft_node_age,min_edge_length=min_edge_length)
    
    
    # Remove deep root if deep rooted
    if deep_rooted:
        rec_tree = TreeOps.remove_deep_root(rec_tree)
        
    # Double check that there are no negative branch lengths
    for edge in rec_tree.postorder_edge_iter():
        if edge.length:
            if edge.length < 0.0:
                edge.length = min_edge_length
                logging.debug("Negative branch length in reconciled tree")
    
    # Recheck rooting
    if not rec_tree.is_rooted:
        logging.warning("Reconciled tree has lost rooting")   

    return rec_tree
