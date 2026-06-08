##########################
##
## Espalier: A Python package for tree reconciliation and reconstructing ARGs using maximum agreement forests.
##
## Copyright 2021-2022 David A. Rasmussen (drasmus@ncsu.edu)
##
############################

"""
Modern TreeSequence conversion utilities.

This module keeps the public behavior of ``Dendro2TSConverter.convert`` but
uses table-array operations for time repairs and tskit's newer edge
compression where available.
"""

import logging
import time
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import tskit

import Espalier.Reconciler as ReconcilerModule
from Espalier.ARGBuilder import _jitter_coal_times, add_path_rec_nodes
from Espalier.ARGNodeTypes import RECOMBINANT_FLAG, count_recombinant_nodes


@dataclass
class ConversionMetrics:
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
    has_recomb = any(
        len(list(node.child_node_iter())) == 1
        for tree in tree_path
        for node in tree.postorder_internal_node_iter()
    )

    if not has_recomb:
        logging.info("Reconciling linked heights and adding recombination nodes before TreeSequence conversion")
        tree_path = ReconcilerModule.reconcile_linked_heights(tree_path)
        tree_path = _jitter_coal_times(tree_path, displace_dt=displace_dt)
        tree_path, _ = add_path_rec_nodes(tree_path)

    return tree_path


def fix_time_constraints_vectorized(
    tables: tskit.TableCollection,
    epsilon: float = 1e-6,
    max_iterations: int = 200,
) -> Tuple[tskit.TableCollection, int, int]:
    times = tables.nodes.time.copy()
    edge_parent = tables.edges.parent
    edge_child = tables.edges.child

    parent_to_children: Dict[int, List[int]] = {}
    for parent, child in zip(edge_parent, edge_child):
        parent = int(parent)
        child = int(child)
        parent_to_children.setdefault(parent, [])
        if child not in parent_to_children[parent]:
            parent_to_children[parent].append(child)

    total_fixes = 0
    iterations_used = 0
    for iteration in range(max_iterations):
        iterations_used = iteration + 1
        violations_fixed = 0
        for node_idx in np.argsort(times, kind="mergesort"):
            node_idx = int(node_idx)
            children = parent_to_children.get(node_idx)
            if not children:
                continue
            max_child_time = float(np.max(times[children]))
            if times[node_idx] <= max_child_time:
                times[node_idx] = max_child_time + epsilon
                violations_fixed += 1
                total_fixes += 1
        if violations_fixed == 0:
            break

    tables.nodes.set_columns(
        flags=tables.nodes.flags,
        time=times,
        population=tables.nodes.population,
        individual=tables.nodes.individual,
        metadata=tables.nodes.metadata,
        metadata_offset=tables.nodes.metadata_offset,
    )
    return tables, total_fixes, iterations_used


def convert_with_metrics(tree_path, tree_intervals) -> Tuple[tskit.TreeSequence, ConversionMetrics]:
    metrics = ConversionMetrics()
    start_total = time.perf_counter()

    tree_path = prepare_tree_path(tree_path)

    from Espalier.Dendro2TSConverter import (
        _sort_nodes_df,
        reindex_edgedf,
        tree2edgedf,
        tree2nodesdf,
    )

    segments = len(tree_path)
    tree = tree_path[0]
    left_pos, right_pos = tree_intervals[0]
    total_length = tree_intervals[-1][-1]

    merged_node_df, tree = tree2nodesdf(tree)
    id_dict = dict(zip([int(x) for x in merged_node_df['unique_id']], merged_node_df.index))
    merged_edge_df = tree2edgedf(tree, left_pos, right_pos, id_dict)

    for loc in range(1, segments):
        tree = tree_path[loc]
        left_pos, right_pos = tree_intervals[loc]

        next_node_df, tree = tree2nodesdf(tree, prev_node_df=merged_node_df)
        merged_node_df = pd.concat([merged_node_df, next_node_df], ignore_index=True)
        merged_node_df.drop_duplicates(ignore_index=True, inplace=True)
        merged_node_df = _sort_nodes_df(merged_node_df)

        id_dict = dict(zip([int(x) for x in merged_node_df['unique_id']], merged_node_df.index))
        merged_edge_df = reindex_edgedf(merged_edge_df, id_dict)

        next_edge_df = tree2edgedf(tree, left_pos, right_pos, id_dict)
        merged_edge_df = pd.concat([merged_edge_df, next_edge_df], ignore_index=True)

    tables = tskit.TableCollection(total_length)
    for _, row in merged_edge_df.iterrows():
        tables.edges.add_row(row['left'], row['right'], int(row['parent']), int(row['child']))
    for _, row in merged_node_df.iterrows():
        tables.nodes.add_row(flags=int(row['flags']), time=row['time'], population=int(row['population']))

    start_fix = time.perf_counter()
    tables, metrics.n_time_fixes, metrics.n_iterations = fix_time_constraints_vectorized(tables)
    metrics.time_fix_time = time.perf_counter() - start_fix

    start_sort = time.perf_counter()
    tables.sort()
    tables.edges.squash()
    tables.sort()
    metrics.table_sort_time = time.perf_counter() - start_sort

    ts = tables.tree_sequence()

    start_extend = time.perf_counter()
    try:
        if hasattr(ts, "extend_haplotypes"):
            extended_ts = ts.extend_haplotypes()
            if extended_ts.num_edges <= ts.num_edges:
                ts = extended_ts
    except Exception as e:
        logging.warning("extend_haplotypes failed: %s", e)
    metrics.extend_haplotypes_time = time.perf_counter() - start_extend

    metrics.total_time = time.perf_counter() - start_total
    metrics.n_nodes = ts.num_nodes
    metrics.n_edges = ts.num_edges
    metrics.n_recomb_nodes = count_recombinant_nodes(ts.tables.nodes.flags)
    logging.info(
        "Modern conversion complete: %s nodes, %s edges, %s recombinant nodes",
        metrics.n_nodes,
        metrics.n_edges,
        metrics.n_recomb_nodes,
    )
    return ts, metrics


def convert(tree_path, tree_intervals) -> tskit.TreeSequence:
    ts, _ = convert_with_metrics(tree_path, tree_intervals)
    return ts

