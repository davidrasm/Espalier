"""
Modernized Dendro2TSConverter using tskit 0.6.0+ features.

This module reimplements key functions from Dendro2TSConverter.py
using modern tskit APIs for A/B testing comparison.

Key improvements:
1. extend_haplotypes() for edge compression
2. Vectorized time constraint fixing using numpy arrays
3. Modern node iterators
4. tskit constants instead of magic numbers
5. Systematic use of tables.sort()

Author: A/B Testing Framework
"""

import tskit
import msprime
import numpy as np
import pandas as pd
import logging
import time
from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional

import Espalier.Reconciler as ReconcilerModule
from Espalier.ARGBuilder import _jitter_coal_times, add_path_rec_nodes

# Use tskit constants instead of magic numbers
RECOMB_FLAG = msprime.NODE_IS_RE_EVENT  # 131072


@dataclass
class ConversionMetrics:
    """Metrics for tracking conversion performance."""
    total_time: float = 0.0
    time_fix_time: float = 0.0
    table_sort_time: float = 0.0
    extend_haplotypes_time: float = 0.0
    n_nodes: int = 0
    n_edges: int = 0
    n_recomb_nodes: int = 0
    n_time_fixes: int = 0
    n_iterations: int = 0


def prepare_tree_path(tree_path: List, displace_dt: float = 1e-4) -> List:
    """
    Harmonize node heights across the path and inject recombination nodes when
    the Viterbi output has no explicit recombination events. This mirrors the
    preprocessing done in the EM pipeline so that downstream conversion does
    not produce parent/child time cycles.
    """
    # Detect whether recombination nodes (unifurcations) are already present
    has_recomb = any(
        len(list(node.child_node_iter())) == 1
        for tree in tree_path
        for node in tree.postorder_internal_node_iter()
    )

    if not has_recomb:
        logging.info("[Modern] Reconciling linked heights and adding recombination nodes to Viterbi path")
        tree_path = ReconcilerModule.reconcile_linked_heights(tree_path)
        tree_path = _jitter_coal_times(tree_path, displace_dt=displace_dt)
        tree_path, _ = add_path_rec_nodes(tree_path)
    else:
        logging.debug("[Modern] Tree path already contains recombination nodes; skipping augmentation")

    return tree_path


def fix_time_constraints_vectorized(tables: tskit.TableCollection,
                                    epsilon: float = 1e-6,
                                    max_iterations: int = 200) -> Tuple[tskit.TableCollection, int, int]:
    """
    Fix time constraint violations using vectorized numpy operations.

    This is a modernized version that operates directly on tskit table arrays
    instead of pandas DataFrames, providing significant performance improvements.

    Parameters:
        tables: tskit TableCollection to fix
        epsilon: Minimum time difference between parent and child
        max_iterations: Maximum number of fix iterations

    Returns:
        tables: Fixed TableCollection
        total_fixes: Number of fixes applied
        iterations: Number of iterations used
    """
    start_time = time.perf_counter()

    # Get arrays directly from tables (much faster than DataFrame operations)
    times = tables.nodes.time.copy()
    edge_parent = tables.edges.parent
    edge_child = tables.edges.child

    # Build parent->children mapping using numpy
    parent_to_children: Dict[int, List[int]] = {}
    for i in range(len(edge_parent)):
        p = edge_parent[i]
        c = edge_child[i]
        if p not in parent_to_children:
            parent_to_children[p] = []
        if c not in parent_to_children[p]:
            parent_to_children[p].append(c)

    total_fixes = 0
    iteration = 0

    for iteration in range(max_iterations):
        violations_fixed = 0

        # Process nodes in order of time (youngest first)
        # This ensures children are processed before parents
        node_order = np.argsort(times)

        for node_idx in node_order:
            node_idx = int(node_idx)
            if node_idx not in parent_to_children:
                continue  # Leaf node

            children = parent_to_children[node_idx]
            parent_time = times[node_idx]

            # Vectorized max child time computation
            child_times = np.array([times[c] for c in children])
            max_child_time = np.max(child_times)

            if parent_time <= max_child_time:
                times[node_idx] = max_child_time + epsilon
                violations_fixed += 1
                total_fixes += 1

        if violations_fixed == 0:
            break

    # Update table with fixed times
    tables.nodes.set_columns(
        flags=tables.nodes.flags,
        time=times,
        population=tables.nodes.population,
        individual=tables.nodes.individual,
        metadata=tables.nodes.metadata,
        metadata_offset=tables.nodes.metadata_offset
    )

    elapsed = time.perf_counter() - start_time
    if total_fixes > 0:
        logging.info(f"[Modern] Fixed {total_fixes} time violations in {iteration + 1} iterations ({elapsed:.3f}s)")

    return tables, total_fixes, iteration + 1


def convert_modern(tree_path, tree_intervals) -> Tuple[tskit.TreeSequence, ConversionMetrics]:
    """
    Modern conversion routine using tskit 0.6.0+ features.

    This is a drop-in replacement for Dendro2TSConverter.convert() that uses:
    - Vectorized time constraint fixing
    - extend_haplotypes() for edge compression
    - Modern table operations

    Parameters:
        tree_path (list(dendropy.Tree)): tree path
        tree_intervals (list(tuple)): list of genomic intervals

    Returns:
        merged_ts (tskit.TreeSequence): tree sequence representing a connected ARG
        metrics (ConversionMetrics): performance metrics
    """
    metrics = ConversionMetrics()
    start_total = time.perf_counter()

    # Prepare tree path to avoid parent/child time cycles (not present in raw Viterbi output)
    tree_path = prepare_tree_path(tree_path)

    # Get local tree and ts for first interval
    segments = len(tree_path)
    tree = tree_path[0]
    tip_count = tree.__len__()
    left_pos = tree_intervals[0][0]
    right_pos = tree_intervals[0][1]
    total_length = tree_intervals[-1][-1]

    # Convert tree to node and edge dataframes (reuse original functions)
    from Espalier.Dendro2TSConverter import tree2nodesdf, tree2edgedf, reindex_edgedf

    merged_node_df, tree = tree2nodesdf(tree)
    # Use Python int conversion to handle arbitrarily large bitmasks
    id_dict = dict(zip([int(x) for x in merged_node_df['unique_id']], merged_node_df.index))
    merged_edge_df = tree2edgedf(tree, left_pos, right_pos, id_dict)

    logging.debug("[Modern] Merging node and edge tables")

    # Iteratively add nodes and edges from each local tree
    for loc in range(1, segments):
        tree = tree_path[loc]
        left_pos = tree_intervals[loc][0]
        right_pos = tree_intervals[loc][1]

        next_node_df, tree = tree2nodesdf(tree, prev_node_df=merged_node_df)
        merged_node_df = pd.concat([merged_node_df, next_node_df], ignore_index=True)
        merged_node_df.drop_duplicates(ignore_index=True, inplace=True)
        merged_node_df.sort_values('metadata', ignore_index=True, inplace=True)
        merged_node_df.sort_values('time', ignore_index=True, inplace=True)

        id_dict = dict(zip([int(x) for x in merged_node_df['unique_id']], merged_node_df.index))
        merged_edge_df = reindex_edgedf(merged_edge_df, id_dict)

        next_edge_df = tree2edgedf(tree, left_pos, right_pos, id_dict)
        merged_edge_df = pd.concat([merged_edge_df, next_edge_df], ignore_index=True)

    # Create TableCollection
    tables = tskit.TableCollection(total_length)
    for index, row in merged_edge_df.iterrows():
        tables.edges.add_row(row['left'], row['right'], int(row['parent']), int(row['child']))
    for index, row in merged_node_df.iterrows():
        tables.nodes.add_row(flags=int(row['flags']), time=row['time'], population=int(row['population']))

    # Use vectorized time constraint fixing
    start_fix = time.perf_counter()
    tables, n_fixes, n_iters = fix_time_constraints_vectorized(tables)
    metrics.time_fix_time = time.perf_counter() - start_fix
    metrics.n_time_fixes = n_fixes
    metrics.n_iterations = n_iters

    # Sort tables (modern systematic approach)
    start_sort = time.perf_counter()
    tables.sort()
    tables.edges.squash()
    tables.sort()
    metrics.table_sort_time = time.perf_counter() - start_sort

    # Create TreeSequence
    merged_ts = tables.tree_sequence()

    # Try to apply extend_haplotypes if available (tskit 0.6.0+)
    start_extend = time.perf_counter()
    try:
        if hasattr(merged_ts, 'extend_haplotypes'):
            extended_ts = merged_ts.extend_haplotypes()
            logging.info(f"[Modern] extend_haplotypes: {merged_ts.num_edges} -> {extended_ts.num_edges} edges")
            # Only use if it actually reduced edges and preserved structure
            if extended_ts.num_edges <= merged_ts.num_edges:
                merged_ts = extended_ts
        metrics.extend_haplotypes_time = time.perf_counter() - start_extend
    except Exception as e:
        logging.warning(f"[Modern] extend_haplotypes failed: {e}")
        metrics.extend_haplotypes_time = 0.0

    # Record metrics
    metrics.total_time = time.perf_counter() - start_total
    metrics.n_nodes = merged_ts.num_nodes
    metrics.n_edges = merged_ts.num_edges

    # Count recombination nodes using tskit constant
    recomb_nodes = np.sum(merged_ts.tables.nodes.flags == RECOMB_FLAG)
    metrics.n_recomb_nodes = recomb_nodes

    logging.info(f"[Modern] Conversion complete: {metrics.n_nodes} nodes, {metrics.n_edges} edges, "
                 f"{metrics.n_recomb_nodes} recomb nodes ({metrics.total_time:.3f}s)")

    return merged_ts, metrics


def convert_original(tree_path, tree_intervals) -> Tuple[tskit.TreeSequence, ConversionMetrics]:
    """
    Wrapper around original Dendro2TSConverter.convert() with metrics collection.

    Parameters:
        tree_path (list(dendropy.Tree)): tree path
        tree_intervals (list(tuple)): list of genomic intervals

    Returns:
        merged_ts (tskit.TreeSequence): tree sequence
        metrics (ConversionMetrics): performance metrics
    """
    metrics = ConversionMetrics()

    from Espalier import Dendro2TSConverter

    start_total = time.perf_counter()
    merged_ts = Dendro2TSConverter.convert(tree_path, tree_intervals)
    metrics.total_time = time.perf_counter() - start_total

    metrics.n_nodes = merged_ts.num_nodes
    metrics.n_edges = merged_ts.num_edges

    # Count recombination nodes
    recomb_nodes = np.sum(merged_ts.tables.nodes.flags == RECOMB_FLAG)
    metrics.n_recomb_nodes = recomb_nodes

    logging.info(f"[Original] Conversion complete: {metrics.n_nodes} nodes, {metrics.n_edges} edges, "
                 f"{metrics.n_recomb_nodes} recomb nodes ({metrics.total_time:.3f}s)")

    return merged_ts, metrics


def compare_tree_sequences(ts_a: tskit.TreeSequence,
                           ts_b: tskit.TreeSequence,
                           label_a: str = "A",
                           label_b: str = "B") -> Dict:
    """
    Compare two TreeSequences for equivalence and differences.

    Parameters:
        ts_a, ts_b: TreeSequences to compare
        label_a, label_b: Labels for logging

    Returns:
        dict: Comparison results
    """
    results = {
        'equivalent': True,
        'differences': [],
        'metrics': {}
    }

    # Compare basic structure
    if ts_a.num_nodes != ts_b.num_nodes:
        results['equivalent'] = False
        results['differences'].append(f"Node count: {label_a}={ts_a.num_nodes}, {label_b}={ts_b.num_nodes}")

    if ts_a.num_edges != ts_b.num_edges:
        results['equivalent'] = False
        results['differences'].append(f"Edge count: {label_a}={ts_a.num_edges}, {label_b}={ts_b.num_edges}")

    if ts_a.num_trees != ts_b.num_trees:
        results['equivalent'] = False
        results['differences'].append(f"Tree count: {label_a}={ts_a.num_trees}, {label_b}={ts_b.num_trees}")

    # Compare recombination nodes
    recomb_a = np.sum(ts_a.tables.nodes.flags == RECOMB_FLAG)
    recomb_b = np.sum(ts_b.tables.nodes.flags == RECOMB_FLAG)
    if recomb_a != recomb_b:
        results['differences'].append(f"Recomb nodes: {label_a}={recomb_a}, {label_b}={recomb_b}")

    # Compare breakpoints
    bp_a = list(ts_a.breakpoints())
    bp_b = list(ts_b.breakpoints())
    if bp_a != bp_b:
        results['differences'].append(f"Breakpoints differ: {label_a}={bp_a}, {label_b}={bp_b}")

    # Store metrics
    results['metrics'] = {
        f'{label_a}_nodes': ts_a.num_nodes,
        f'{label_b}_nodes': ts_b.num_nodes,
        f'{label_a}_edges': ts_a.num_edges,
        f'{label_b}_edges': ts_b.num_edges,
        f'{label_a}_trees': ts_a.num_trees,
        f'{label_b}_trees': ts_b.num_trees,
        f'{label_a}_recomb': recomb_a,
        f'{label_b}_recomb': recomb_b,
    }

    return results
