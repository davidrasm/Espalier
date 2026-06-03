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

from Espalier.ARGNodeTypes import classify_node_flags, is_recombinant

class SCAR(object):
    
    """
        Model class for the structured coalescent with ancestral recombination (SCAR)
        Ancestral states can be given (known) or marginalized over (unknown).
        This version relaxes rules about the ARG:
        Specifically coal nodes need not be strictly bifurcating and can have any number of children.
        But all children must have one and only one parent.
    """
    
    rate_model_name = "uniform_site"
    rate_units = "per_site_per_generation"
    rate_display_units = "per site per generation"

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
        self.bounds = kwargs.get('bounds', (0, np.inf))
        self.dt_step = kwargs.get('dt_step', 0.1)
        self.known_ancestral_states = kwargs.get('known_ancestral_states', False)
        self.last_opt_result = None
        self.last_rate_estimate_status = None
        self.last_likelihood_diagnostics = None
        
    def compute_neg_log_like(self,rec_rate,ts):
        
        """
            Compute the negative log likelihood of ARG in TreeSequence under the SCAR model
            Here we compute the negatve log likelihood because the scipy opt minimizes this func
            Note: ts, M, Ne, genome_length need to be passed through optimizer as a tuple
            
            Parameters:
                rec_rate (float): recombination rate per site
                ts (tskit.TreeSequence): ARG encoded as a TreeSequence object
                
        """
    
        self.last_likelihood_diagnostics = {
            "status": "running",
            "rate": float(rec_rate),
            "observed_recombination_events": 0,
            "zero_recomb_opportunity_events": 0,
            "recomb_exposure": 0.0,
            "first_nonfinite_event": None,
        }

        if rec_rate < 0 or not np.isfinite(rec_rate):
            self._record_nonfinite_likelihood(
                reason="invalid_rate",
                order_idx=None,
                node_id=None,
                event=None,
                event_type=None,
                event_prob=None,
            )
            return np.inf

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
        
        node_times = np.array(ts.tables.nodes.time, ndmin=1)
        node_order = np.lexsort((np.arange(len(node_times)), node_times))

        # Iterate through each event/node in the ARG working backwards through time.
        # Node IDs are stable, but tskit does not require node-table rows to be time-sorted.
        for order_idx, idx in enumerate(node_order):
            idx = int(idx)
            event = ts.tables.nodes[idx]
            
            # Get time of event and time of next event"
            event_time = event.time
            if (order_idx + 1 < len(node_order)): # if not at final event
                next_time = ts.tables.nodes[int(node_order[order_idx + 1])].time
            else:
                next_time = event.time
            t_elapsed = next_time - event_time # time elapsed between events
            

            # Determine event type from tskit/msprime node flags.
            event_type = classify_node_flags(event.flags)
            
            # Initialize prob of observing events and log-probs of no events.
            event_prob = 1.0
            log_no_coal = 0.0
            log_no_mig = 0.0
            log_no_recomb = 0.0
            
            # Update active lineages based on event type: coalescent/sampling/migration events
            if 'sample' == event_type:
                
                # Add sampled lineage
                active_lines.append(idx)
                active_rec_links.append(self._get_line_recomb_opportunities(idx,children,rights,lefts))
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

                if len(child_indexes) == 0:
                    event_prob = 1.0
                    coal_children = []
                else:
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
                    active_rec_links[child_indexes[0]] = self._get_line_recomb_opportunities(idx,children,rights,lefts)
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
                coal_children = [x for x in coal_children if x in active_lines]
                if len(coal_children) < 2:
                    event_prob = 1.0
                    coal_children = []
                else:
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
                    active_rec_links[child1_idx] = self._get_line_recomb_opportunities(idx,children,rights,lefts)
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
                if len(child) != 1:
                    child = None
                else:
                    child = int(child[0])
                
                # Remember that child may have already been removed from active_lines
                if child in active_lines:
                
                    # Get indexes of both (left and right) parent of child"
                    recomb_parents = parents[children == child]
                    
                    # Parents edges may occur more than once in the tree series if not in contiguous trees
                    recomb_parents = np.unique(recomb_parents)
                    
                    # Make sure recombination event results in a child splitting into two parents"
                    recomb_parents = [x for x in recomb_parents if is_recombinant(ts.tables.nodes[x].flags)]
                    if len(recomb_parents) != 2:
                        event_prob = 1.0
                        recomb_parents = []
                    else:
                        # Get parents
                        left_parent = recomb_parents[0]
                        right_parent = recomb_parents[1]

                        child_idx = active_lines.index(child)

                        # Compute recombination event prob
                        links = active_rec_links[child_idx] # Links gives num of sites at which lineage carries material ancestral to the sample
                        self.last_likelihood_diagnostics["observed_recombination_events"] += 1
                        if links <= 0:
                            self.last_likelihood_diagnostics["zero_recomb_opportunity_events"] += 1
                        event_prob = rec_rate * links

                        # Compute new parent state probs
                        if self.known_ancestral_states:
                            parent_probs = np.zeros(pops)
                            parent_probs[event.population] = 1.0
                        else:
                            parent_probs = line_state_probs[child_idx]

                        # Update lineage arrays - overwriting child with left parent"
                        active_lines[child_idx] = left_parent # name of parent
                        active_rec_links[child_idx] = self._get_line_recomb_opportunities(left_parent,children,rights,lefts)
                        line_state_probs[child_idx] = parent_probs

                        # Add other recombining parent
                        active_lines.append(right_parent)
                        active_rec_links.append(self._get_line_recomb_opportunities(right_parent,children,rights,lefts))
                        line_state_probs.append(parent_probs)
            
            if 'migration' == event_type:
                
                # Find migrating (child) lineage
                mig_child = children[parents == idx] # parent has id == idx in parent column of edges table
                mig_child = np.unique(mig_child)
                if len(mig_child) == 0:
                    event_prob = 1.0
                else:
                    mig_child = int(mig_child[0])
                    if mig_child not in active_lines:
                        event_prob = 1.0
                    else:
                        # Get migration info from nodes list
                        curr_state = populations[mig_child]
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
                    log_no_coal = -np.sum(lambdas) * t_elapsed
                
                    # Compute prob of no migration over the time interval
                    sam = 0
                    for i in range(pops):
                        for z in range(pops):
                            sam += (A[i])*(self.M[i][z])
                    log_no_mig = -sam * t_elapsed
                    
                    # Compute prob of no recombination event over the time interval
                    # Links are computed per population b/c we are assuming recombination can only happen in same pop
                    sum_links = self._sum_active_recomb_opportunities(line_state_probs, active_rec_links)
                    self.last_likelihood_diagnostics["recomb_exposure"] += float(sum_links * t_elapsed)
                    log_no_recomb = -sum_links * rec_rate * t_elapsed # assumes rho / genome_length is constant across pops
                    
                else: # Unknown ancestral lineage states
                
                    # Integrate lineage prob equations backwards
                    dt_times = list(np.arange(event_time,next_time,self.dt_step)) # integration steps going backwards in time
                    for step_idx,tx in enumerate(dt_times):
                        
                        # Get time step
                        if (step_idx+1 < len(dt_times)):
                            dt = dt_times[step_idx+1] - tx # integration time step
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
                        log_no_coal += -np.sum(lambdas) * dt
                        
                        # Compute prob of no migration over the time interal"
                        log_no_mig = 0.0
                        
                        # Compute prob of no recombination event over the time interval
                        # Links are computed per population b/c we are assuming recombination can only happen in same pop
                        sum_links = self._sum_active_recomb_opportunities(line_state_probs, active_rec_links)
                        self.last_likelihood_diagnostics["recomb_exposure"] += float(sum_links * dt)
                        log_no_recomb += -sum_links * rec_rate * dt
            
            log_terms = [log_no_coal, log_no_mig, log_no_recomb]
            if event_prob <= 0 or not np.isfinite(event_prob):
                self._record_nonfinite_likelihood(
                    reason="nonpositive_event_probability",
                    order_idx=order_idx,
                    node_id=idx,
                    event=event,
                    event_type=event_type,
                    event_prob=event_prob,
                )
                return np.inf
            if not all(np.isfinite(term) for term in log_terms):
                self._record_nonfinite_likelihood(
                    reason="nonfinite_waiting_time_probability",
                    order_idx=order_idx,
                    node_id=idx,
                    event=event,
                    event_type=event_type,
                    event_prob=event_prob,
                    log_no_coal=log_no_coal,
                    log_no_mig=log_no_mig,
                    log_no_recomb=log_no_recomb,
                )
                return np.inf
            log_like += np.log(event_prob) + log_no_coal + log_no_mig + log_no_recomb
            if not np.isfinite(log_like):
                self._record_nonfinite_likelihood(
                    reason="nonfinite_log_likelihood",
                    order_idx=order_idx,
                    node_id=idx,
                    event=event,
                    event_type=event_type,
                    event_prob=event_prob,
                    log_no_coal=log_no_coal,
                    log_no_mig=log_no_mig,
                    log_no_recomb=log_no_recomb,
                )
                return np.inf
            
        self.last_likelihood_diagnostics["status"] = "finite"
        self.last_likelihood_diagnostics["neg_log_likelihood"] = float(-log_like)
        return -log_like

    def _sum_active_recomb_opportunities(self, line_state_probs, active_rec_links):
        if not active_rec_links:
            return 0.0
        line_prod = np.array(line_state_probs) * np.array(active_rec_links)[:, np.newaxis]
        return float(np.sum(np.sum(line_prod)))

    def _record_nonfinite_likelihood(
        self,
        reason,
        order_idx,
        node_id,
        event,
        event_type,
        event_prob,
        log_no_coal=None,
        log_no_mig=None,
        log_no_recomb=None,
    ):
        if self.last_likelihood_diagnostics is None:
            self.last_likelihood_diagnostics = {}
        self.last_likelihood_diagnostics["status"] = "nonfinite"
        self.last_likelihood_diagnostics["nonfinite_reason"] = reason
        if event is None:
            event_details = None
        else:
            event_details = {
                "order_idx": int(order_idx),
                "node_id": int(node_id),
                "time": float(event.time),
                "flags": int(event.flags),
                "event_type": event_type,
                "event_prob": float(event_prob) if event_prob is not None else None,
            }
            if log_no_coal is not None:
                event_details["log_no_coal"] = float(log_no_coal)
            if log_no_mig is not None:
                event_details["log_no_mig"] = float(log_no_mig)
            if log_no_recomb is not None:
                event_details["log_no_recomb"] = float(log_no_recomb)
        self.last_likelihood_diagnostics["first_nonfinite_event"] = event_details

    def _get_line_recomb_opportunities(self,line,children,rights,lefts):
        
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
        res.fun = self.compute_neg_log_like(mle, ts)
        self.last_opt_result = res
        if not np.isfinite(res.fun):
            self.last_rate_estimate_status = "nonfinite_likelihood"
        else:
            self.last_rate_estimate_status = self._classify_rate_estimate(mle)

        return mle

    def _classify_rate_estimate(self, rate):

        """
            Classify whether an optimized rate estimate is at a bound or in the interior.
        """

        lower, upper = self.bounds
        span = upper - lower if np.isfinite(lower) and np.isfinite(upper) else max(abs(rate), 1.0)
        tol = max(1e-12, span * 1e-3)

        if np.isfinite(lower) and rate <= lower + tol:
            return "hit_lower_bound"
        if np.isfinite(upper) and rate >= upper - tol:
            return "hit_upper_bound"
        return "interior"


class BoundarySCAR(SCAR):

    """
        SCAR variant for segmented reassortment where recombination is only possible
        at predefined segment boundaries rather than uniformly at every genomic link.
    """

    rate_model_name = "boundary_hotspot"
    rate_units = "per_boundary_per_generation"
    rate_display_units = "per boundary per generation"

    def __init__(self, rec_rate, M, Ne, genome_length, boundary_positions, **kwargs):
        self.boundary_positions = sorted(int(position) for position in boundary_positions)
        self.n_boundaries = len(self.boundary_positions)
        super().__init__(rec_rate, M, Ne, genome_length, **kwargs)

    def _get_line_recomb_opportunities(self, line, children, rights, lefts):

        """
            Count configured segment boundaries that remain eligible for reassortment
            on a lineage's ancestral material.
        """

        line_mask = children == line
        if not np.any(line_mask):
            return 0

        intervals = [
            (left, right)
            for left, right in zip(lefts[line_mask], rights[line_mask])
            if right > left
        ]
        eligible = 0
        for boundary in self.boundary_positions:
            carries_left_side = any(left < boundary <= right for left, right in intervals)
            carries_right_side = any(left <= boundary < right for left, right in intervals)
            if carries_left_side and carries_right_side:
                eligible += 1
        return eligible

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
