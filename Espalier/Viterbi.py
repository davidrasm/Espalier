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
import logging

def viterbi(tree_trellis, trans_probs, like_array):
    
    """
        
        Viterbi algorithm to select a most likely tree path from trellis of candidate trees.
        Tree path is the maximum a posterior (MAP) path given the sequence data and trans probs provided.
        
        Parameters:     
           tree_trellis (list): list of TreeLists containting candidate tree states for each genomic region
           trans_probs (3D array-like): transition probabilities where entry n,i,j contains prob of transitioning from tree i in region n-1 to tree j in region n
           like_array (2D array-like): sequence likelihoods where entry n,i gives likelihood of sequence data in region n given tree i  
           
        Returns:     
           tree_path (list): selected tree path
           opt (list): integer-valued indexes of trees in tree path
    """
    
    logging.info("Running Viterbi algorithm")
    
    segments = len(tree_trellis) # i.e. genomic regions/intervals
    tree_states = [len(trees) for trees in tree_trellis]
    max_states = max(tree_states)
    state_probs = np.zeros((segments,max_states))
    path = np.zeros((segments,max_states)).astype(int) # pointers
    prior = [1.0]*tree_states[0]
    
    # Compute state probs for trees in region 0 based on priors 
    for st in range(tree_states[0]):
        state_probs[0][st] = np.log(prior[st]) + like_array[0][st] # always working on log scale
        path[0][st] = 0
        
    # Forward pass: computes prob of path ending in state st given trees observed up to region loc
    for loc in range(1, segments):
        for st in range(tree_states[loc]):
            probs = []
            for prev_st in range(tree_states[loc-1]):
                probs.append(state_probs[loc-1][prev_st] + np.log(trans_probs[loc][prev_st][st]))
            path[loc][st] = np.argmax(probs) 
            state_probs[loc][st] = max(probs) + like_array[loc][st]
    
    # Select most likely state/tree for final region 
    final_probs = state_probs[-1][:tree_states[-1]]
    best_st = np.argmax(final_probs) 
    opt = []
    opt.append(best_st)
    max_prob = final_probs[best_st]
 
    # Backwards pass to select tree path based on pointers in path
    previous = best_st
    tree_path = []
    tree_path.append(tree_trellis[-1][best_st])
    for loc in range(segments - 2, -1, -1):
        previous = path[loc+1][previous]
        opt.append(previous)
        tree_path.append(tree_trellis[loc][previous])
    
    # Reverse tree_path and opt
    tree_path = list(reversed(tree_path))
    opt = list(reversed(opt))
    
    logging.info('The path of trees is ' + ' '.join(str(opt)) + ' with a path likelihood of %s' % max_prob)
    
    return tree_path, opt


def _norm_log_probs(log_probs):
    
    """
        Normalizes array of log probabilities.
        Log probabilities are first recentered by substracting lowest log prob (to prevent underflow)
        then probs are converted to a linear scale before normalizing.
    """
    
    log_probs = log_probs - np.min(log_probs) # recenter so min log prob is at zero
    probs = np.exp(log_probs) # convert back to linear scale
    norm_probs = probs / np.sum(probs)
    return norm_probs

def forward_back(tree_trellis, trans_probs, like_array):
    
    """
        
        Run forward-backwards algorithm to sample a tree path from trellis of candidate trees.
        Unlike in Viterbi algorithm, trees are sampled probabilistically based on the path probs.
        Note: to prevent numerical underflow in the path probs we work with normalized probs.
        
        Parameters:     
           tree_trellis (list): list of TreeLists containting candidate tree states for each genomic region
           trans_probs (3D array-like): transition probabilities where entry n,i,j contains prob of transitioning from tree i in region n-1 to tree j in region n
           like_array (2D array-like): sequence likelihoods where entry n,i gives likelihood of sequence data in region n given tree i  
           
        Returns:     
           tree_path (list): sampled tree path
    """
    
    
    segments = len(tree_trellis)
    tree_states = [len(trees) for trees in tree_trellis]
    max_states = max(tree_states)
    state_probs = np.zeros((segments,max_states))
    prior = [1.0]*tree_states[0]
    
    # Compute state probs for tree states in region 0 based on priors
    log_state_probs = []
    for st in range(tree_states[0]):
        log_state_probs.append(np.log(prior[st]) + like_array[0][st]) # always working on log scale
    
    # Normalize (log) probs
    norm_probs = _norm_log_probs(log_state_probs)
    for st in range(tree_states[0]):
        state_probs[0][st] = norm_probs[st]
    
    # Run forward pass
    for loc in range(1, segments):
        
        # Compute forward probs (i.e. conditional likelihood of segments seen so far for T1, T2... TN"
        log_state_probs = []
        for st in range(tree_states[loc]):
            psum = 0 # sum over prev_states
            for prev_st in range(tree_states[loc-1]):
                psum += state_probs[loc-1][prev_st] * trans_probs[loc][prev_st][st] 
            log_state_probs.append(np.log(psum) + like_array[loc][st])
        
        # Normalize (log) probs 
        norm_probs = _norm_log_probs(log_state_probs)
        for st in range(tree_states[loc]):
            state_probs[loc][st] = norm_probs[st]
    
    # Randomly sample state/tree for final loc proportional to probs" 
    final_probs = state_probs[-1][:tree_states[-1]]
    chosen_st = np.random.choice(list(range(tree_states[-1])),1,p=final_probs)[0]
    opt = []
    opt.append(chosen_st)
 
    # Backward pass
    previous = chosen_st
    tree_path = []
    back_prob = like_array[-1][chosen_st] # these are on a log scale
    tree_path.append(tree_trellis[-1][chosen_st])
    for loc in range(segments - 2, -1, -1):
        weights = []
        for st in range(tree_states[loc]):
            weights.append(state_probs[loc][st] * trans_probs[loc+1][st][previous]) # did we ever check indexing here?
        W = weights / np.sum(weights)
        k = np.random.choice(list(range(tree_states[loc])),1,p=W)[0]
        tprob = trans_probs[loc+1][k][previous]
        previous = k
        back_prob += like_array[loc][previous] + np.log(tprob)
        opt.append(previous)
        tree_path.append(tree_trellis[loc][previous])
    
    logging.info('The path of sampled trees is ' + ' '.join(str(opt)) + ' with a max path likelihood of %s' % back_prob)
    
    return list(reversed(tree_path))

