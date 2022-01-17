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
from scipy.linalg import expm
from scipy.optimize import minimize_scalar

class SCAR(object):
    
    """
        Model class for the structured coalescent with ancestral recombination (SCAR)
        Ancestral states can be given (known) or marginalized over (unknown).
        This version relaxes rules about the ARG:
        Specifically coal nodes need not be strictly bifurcating and can have any number of children.
        But all children must have one and only one parent.
    """
    
    def __init__(self,rec_rate,M,Ne,genome_length,**kwargs):
        
        '''             
            Parameters: 
                rec_rate (float): recombination rate per lineage per site
                M (2D list/array): migration rates (forward-time) between subpopulations 
                Ne (1D list/array): effective population size of each subpopulation
                genome_length (int): total length of genome
               
            Optional keyword arguments:
                bounds (tuple): lower and upper bounds on estimated parameter given as (lower,upper)
                dt_step (str): integration time step used to compute lineage state probabilities
                known_ancestral_states (boolean): ancestral states need to be given in tree series (ts) if True
        '''
        
        # Model params
        self.rec_rate = rec_rate
        self.M = np.transpose(np.array(M)) # transpose to get reverse time matrix        
        self.Ne = np.array(Ne)
        self.genome_length = genome_length
        
        # Likelihood calculation params
        self.bounds = kwargs.get('bounds', (0,np.Inf))
        self.dt_step = kwargs.get('dt_step', 0.1)
        self.known_ancestral_states = kwargs.get('known_ancestral_states', False)
        
    def compute_neg_log_like(self,rec_rate,ts):
        
        """
            Compute the negative log likelihood of ARG in TreeSequence under the SCAR model
            Here we compute the negatve log likelihood because the scipy opt minimizes this func
            Note: ts, M, Ne, genome_length need to be passed through optimizer as a tuple
            
            Parameters:
                rec_rate (float): recombination rate per site
                ts (tskit.TreeSequence): ARG encoded as a TreeSequence object
                
        """
    
        pops = ts.num_populations
        if pops == 0: # we never assigned populations
            pops = 1
        #samples = ts.num_samples
        #states = [st.id for st in ts.populations()]
        
        # Get transition rate matrix
        Q = self.M - np.diag(np.sum(self.M,axis=1)) # set diagonals to negative row sums
        
        # Get ARG data from ts.tables
        populations = np.array(ts.tables.nodes.population,ndmin=1)
        children = np.array(ts.tables.edges.child,ndmin=1)
        parents = np.array(ts.tables.edges.parent,ndmin=1)
        lefts = np.array(ts.tables.edges.left,ndmin=1)
        rights = np.array(ts.tables.edges.right,ndmin=1)
        
        # Init lineage arrays
        active_lines = [] # active lines in ARG
        active_rec_links = [] # tracks ancestral material that can recomine on active_lines
        line_state_probs = [] # lineage state probabilities for active lines
        log_like = 0.0 # log likelihood of full tree
        
        # Iterate through each event/node in ARG TreeSequence working backwards through time
        for idx, event in enumerate(ts.tables.nodes):
            
            # Get time of event and time of next event"
            event_time = event.time
            if (idx+1 < len(ts.tables.nodes)): # if not at final event
                next_time = ts.tables.nodes[idx+1].time
            else:
                next_time = event.time
            t_elapsed = next_time - event_time # time elapsed between events
            

            # Determine event type from tskit event.flags
            event_type = None
            if event.flags == 1:
                event_type = 'sample'
            if event.flags == 0:
                event_type = 'coalescent'
            if event.flags == 131072:
                event_type = 'recombination'
            if event.flags == 262144:
                event_type = 'hidden_coalescent'
            if event.flags == 524288:
                event_type = 'migration'
            
            # Initialize prob of observing events or no events
            event_prob = 1.0
            prob_no_coal = 1.0
            prob_no_mig = 1.0
            prob_no_recomb = 1.0
            
            # Update active lineages based on event type: coalescent/sampling/migration events
            if 'sample' == event_type:
                
                # Add sampled lineage
                active_lines.append(idx)
                active_rec_links.append(self._get_line_links(idx,children,rights,lefts))
                state_probs = np.zeros(pops)
                if event.population == -1: # we never assigned populations
                    state_probs[0] = 1.0 # set prob to 1.0 for sampled state
                else:
                    state_probs[event.population] = 1.0 # set prob to 1.0 for sampled state
                line_state_probs.append(state_probs)            
            
            if 'coalescent' == event_type:
                
                # Get children of parent node at coalescent event
                coal_children = children[parents == idx] # parent has id == idx in parent column of edges table
                
                # Get uniique children b/c the same parent/child edge may occur more than once in the tree series if not in contiguous local trees
                coal_children = np.unique(coal_children)
                

                # Find coal_children in active_lines
                coal_children = [x for x in coal_children if x in active_lines]
                child_indexes = [active_lines.index(x) for x in coal_children]
                
                # Compute coalescent event prob for arbitrary number of children 
                coal_probs = np.ones(pops)
                for child_idx in child_indexes:
                    coal_probs *= line_state_probs[child_idx]
                coal_probs = coal_probs / self.Ne
                lambda_sum = sum(coal_probs)
                event_prob = lambda_sum
                
                # Compute new parent state probs
                if self.known_ancestral_states:
                    parent_probs = np.zeros(pops)
                    parent_probs[event.population] = 1.0
                else:
                    parent_probs = coal_probs / lambda_sum # renormalize probs
                    
                # Update lineage arrays - overwriting child1 with parent
                active_lines[child_indexes[0]] = idx # name of parent
                active_rec_links[child_indexes[0]] = self._get_line_links(idx,children,rights,lefts)
                line_state_probs[child_indexes[0]] = parent_probs
                child_indexes.pop(0) # remove first index given to parent
                for child_idx in sorted(child_indexes, reverse=True): # remove in reverse order so indexes don't change
                    del active_lines[child_idx]
                    del active_rec_links[child_idx]
                    del line_state_probs[child_idx]
            
            if 'hidden_coalescent' == event_type:
                
                # Hidden coalescent in ARG not observed in local trees - only need to update active_lines but nothing else
                
                coal_children = children[parents == idx]
                coal_children = np.unique(coal_children)
                child1 = coal_children[0]
                child2 = coal_children[1]
                child1_idx = active_lines.index(child1)
                child2_idx = active_lines.index(child2)
                
                # Compute likelihood of coalescent event
                p1 = line_state_probs[child1_idx]
                p2 = line_state_probs[child2_idx]
                coal_probs = (p1 * p2) / self.Ne
                lambda_sum = sum(coal_probs)
                event_prob = lambda_sum
                
                # Compute new parent state probs"
                if self.known_ancestral_states:
                    parent_probs = np.zeros(pops)
                    parent_probs[event.population] = 1.0
                else:
                    parent_probs = coal_probs / lambda_sum
                
                # Update lineage arrays - overwriting child1 with parent"
                active_lines[child1_idx] = idx # name of parent
                active_rec_links[child1_idx] = self._get_line_links(idx,children,rights,lefts)
                line_state_probs[child1_idx] = parent_probs
                del active_lines[child2_idx]
                del active_rec_links[child2_idx]
                del line_state_probs[child2_idx]
            
            if "recombination" == event_type:
                
                # At a recombination event a child node will have two different parents
                # We need to find the child shared among these two parents
                # Then replace child with left parent and add right parent
                
                # Find child of parent node
                child = children[parents == idx]
                child = np.unique(child)
                assert len(child) == 1
                
                # Remember that child may have already been removed from active_lines
                if child in active_lines:
                
                    # Get indexes of both (left and right) parent of child"
                    recomb_parents = parents[children == child]
                    
                    # Parents edges may occur more than once in the tree series if not in contiguous trees
                    recomb_parents = np.unique(recomb_parents)
                    
                    # Make sure recombination event results in a child splitting into two parents"
                    recomb_parents = [x for x in recomb_parents if ts.tables.nodes[x].flags == 131072] # have to be recombination event
                    assert len(recomb_parents) == 2
                    
                    # Get parents
                    left_parent = recomb_parents[0]
                    right_parent = recomb_parents[1]
        
                    child_idx = active_lines.index(child)
                    
                    # Compute recombination event prob
                    links = active_rec_links[child_idx] # Links gives num of sites at which lineage carries material ancestral to the sample
                    event_prob = rec_rate * links
                    
                    # Compute new parent state probs
                    if self.known_ancestral_states:
                        parent_probs = np.zeros(pops)
                        parent_probs[event.population] = 1.0
                    else:
                        parent_probs = line_state_probs[child_idx]
                    
                    # Update lineage arrays - overwriting child with left parent"
                    active_lines[child_idx] = left_parent # name of parent
                    active_rec_links[child_idx] = self._get_line_links(left_parent,children,rights,lefts)
                    line_state_probs[child_idx] = parent_probs
                    
                    # Add other recombining parent
                    active_lines.append(right_parent)
                    active_rec_links.append(self._get_line_links(right_parent,children,rights,lefts))
                    line_state_probs.append(parent_probs)
            
            if 'migration' == event_type:
                
                # Find migrating (child) lineage
                mig_child = children[parents == idx] # parent has id == idx in parent column of edges table
                mig_child = np.unique(mig_child)
                
                # Get migration info from nodes list
                curr_state = populations[mig_child[0]]
                new_state = populations[idx]
                
                migrant_idx = active_lines.index(mig_child) #change this for ts index
                
                # Update lineage arrays
                active_lines[migrant_idx] = idx # name of parent
                
                # Compute event prob
                if self.known_ancestral_states:
                    new_probs = np.zeros(pops)
                    new_probs[new_state] = 1.0 # event. population
                    line_state_probs[migrant_idx] = new_probs
                    event_prob = self.M[curr_state][new_state]
                else:
                    event_prob = 1.0 # pretend as if we don't see migration events
                            
            # Compute prob of no coalescent over time interval
            if not np.isclose(t_elapsed, 0):
                
                if self.known_ancestral_states:
                    
                    # Sum line probs to get total number of lines in each state A
                    A = np.zeros(pops)
                    for probs in line_state_probs: A += probs
                    
                    # Compute prob of no coalescent over time interval
                    pairs = (A * (A-1)) / 2 # number of pairs in each pop
                    lambdas =  pairs * (1/self.Ne) # coal rate in each pop   
                    prob_no_coal = np.exp(-np.sum(lambdas)*t_elapsed)
                
                    # Compute prob of no migration over the time interval
                    sam = 0
                    for i in range(pops):
                        for z in range(pops):
                            sam += (A[i])*(self.M[i][z])
                    prob_no_mig = np.exp(-sam*t_elapsed)
                    
                    # Compute prob of no recombination event over the time interval
                    # Links are computed per population b/c we are assuming recombination can only happen in same pop
                    line_prod = np.array(line_state_probs) * np.array(active_rec_links)[:, np.newaxis]
                    sum_links = np.sum(np.sum(line_prod))
                              
                    prob_no_recomb = np.exp(-sum_links * rec_rate * t_elapsed) # assumes rho / genome_length is constant across pops
                    
                else: # Unknown ancestral lineage states
                
                    # Integrate lineage prob equations backwards
                    dt_times = list(np.arange(event_time,next_time,self.dt_step)) # integration steps going backwards in time
                    for idx,tx in enumerate(dt_times):
                        
                        # Get time step
                        if (idx+1 < len(dt_times)):
                            dt = dt_times[idx+1] - tx # integration time step
                        else:
                            dt = next_time - tx
    
                        # Should not need to exponentiate transition matrix if dt is small enough
                        expQdt = expm(Q*dt) # exponentiate time-scaled transition rate matrix
    
                        # Update line state probs using Euler integration
                        for ldx,probs in enumerate(line_state_probs):
                            line_state_probs[ldx] = np.matmul(probs,expQdt)
                        
                        # Update total number of lines in each state A
                        A = np.zeros(pops)
                        for probs in line_state_probs: A += probs # sum line probs to get total number of lines in each state
                        
                        # Compute prob of no coalescent over time interval
                        pairs = (A * (A-1)) / 2 # number of pairs in each pop
                        pairs = pairs.clip(min=0) # make sure non are negative
                        lambdas = pairs * (1/self.Ne) # coal rate in each pop
                        prob_no_coal *= np.exp(-np.sum(lambdas)*dt)
                        
                        # Compute prob of no migration over the time interal"
                        prob_no_mig = 1.0
                        
                        # Compute prob of no recombination event over the time interval
                        # Links are computed per population b/c we are assuming recombination can only happen in same pop
                        line_prod = np.array(line_state_probs) * np.array(active_rec_links)[:, np.newaxis]
                        sum_links = np.sum(np.sum(line_prod))
                                            
                        prob_no_recomb *= np.exp(-sum_links * rec_rate * dt)
            
            log_like += np.log(event_prob) + np.log(prob_no_coal) + np.log(prob_no_mig) + np.log(prob_no_recomb)
            
        return -log_like


    def _get_line_links(self,line,children,rights,lefts):
        
        """
            Compute the number of sites (links) ancestral to the sample and thus eligible to undergo recombination
        """
        
        line_links = 0
        if len(children[children==line])==1:
            line_links = (rights[children==line] - lefts[children==line])[0] - 1
        elif len(children[children==line])>=2:
            line_links = max(rights[children==line]) - min(lefts[children==line]) - 1 
            
        return line_links
    
    def opt_MLE(self,ts):
        
        """
            Find MLE of single parameter (assumed to be rec_rate) using numerical optimization.
            TODO: Generalize to allow for other demographic parameters to be estimataed.
        """
        
        # Optimize likelihood by minimizing negative log likelihood
        res = minimize_scalar(self.compute_neg_log_like, args=(ts), bounds=self.bounds, method='bounded')
        mle = res.x
    
        return mle

if __name__ == '__main__':
       
    from Espalier.sim import ARGSimulator
    
    # Specify sim params
    samples = 10
    genome_length = 1e4
    rho = 2 # recombination rate per genome
    rec_rate = rho / genome_length # recombination rate per site
    Ne = 1.0  # effective pop sizes
    #M = [[0.0,0.25],[0.25,0.0]]  # migration rate matrix
    M = [[0]]
    
    # Run sim
    ts = ARGSimulator.sim_ARG(sample_size=samples,Ne=Ne,length=genome_length,recombination_rate=rec_rate,min_breakpoints=1)    
    breaks = len(ts.breakpoints(as_array=True)) - 2 # minus 2 b/c tskit counts ends as breakpoints
    print('recombination breakpoints = ' + str(breaks))
    
    #mig_events = ts.num_migrations
    #print('Num migrations: ', str(ts.num_migrations))
    
    # Initialize SCAR model class
    bounds = (0.0,0.1)
    scar_model = SCAR(rec_rate,M,Ne,genome_length,bounds=bounds)
    
    # Check numerical optimization for MLE of single param
    mle = scar_model.opt_MLE(ts)
    print(mle)
    
    # Check likelihood is valid
    #L = compute_like(ts,**params)

