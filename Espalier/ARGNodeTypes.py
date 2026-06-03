##########################
##
## Espalier: A Python package for tree reconciliation and reconstructing ARGs using maximum agreement forests.
##
## Copyright 2021-2022 David A. Rasmussen (drasmus@ncsu.edu)
##
############################

"""
Helpers for working with ARG node flags across msprime versions.

msprime now exposes ARG node types through ``msprime.NodeType``. The older
integer constants are still present but deprecated, so this module centralizes
the compatibility handling and keeps the rest of Espalier on bitwise checks.
"""

import numpy as np
import tskit

try:  # pragma: no cover - exercised only in environments without msprime
    import msprime
except ImportError:  # pragma: no cover
    msprime = None


SAMPLE_FLAG = int(tskit.NODE_IS_SAMPLE)


def _node_type_flag(node_type_name, legacy_name, fallback):
    if msprime is not None:
        node_type = getattr(msprime, "NodeType", None)
        if node_type is not None and hasattr(node_type, node_type_name):
            return int(getattr(node_type, node_type_name).value)
        if hasattr(msprime, legacy_name):
            return int(getattr(msprime, legacy_name))
    return int(fallback)


RECOMBINANT_FLAG = _node_type_flag("RECOMBINANT", "NODE_IS_RE_EVENT", 131072)
COMMON_ANCESTOR_FLAG = _node_type_flag("COMMON_ANCESTOR", "NODE_IS_CA_EVENT", 262144)
MIGRANT_FLAG = _node_type_flag("MIGRANT", "NODE_IS_MIG_EVENT", 524288)


def has_flag(flags, flag):
    return (int(flags) & int(flag)) != 0


def is_sample(flags):
    return has_flag(flags, SAMPLE_FLAG)


def is_recombinant(flags):
    return has_flag(flags, RECOMBINANT_FLAG)


def is_common_ancestor(flags):
    return has_flag(flags, COMMON_ANCESTOR_FLAG)


def is_migrant(flags):
    return has_flag(flags, MIGRANT_FLAG)


def classify_node_flags(flags):
    if is_sample(flags):
        return "sample"
    if is_recombinant(flags):
        return "recombination"
    if is_common_ancestor(flags):
        return "hidden_coalescent"
    if is_migrant(flags):
        return "migration"
    return "coalescent"


def count_recombinant_nodes(flags):
    flags_array = np.asarray(flags, dtype=np.uint64)
    return int(np.count_nonzero(flags_array & np.uint64(RECOMBINANT_FLAG)))


def summarize_recombination_events(ts):
    """
    Summarize recombination events from TreeSequence edge structure.

    A valid recombination event is represented by a child lineage splitting into
    exactly two recombinant parents. Counting flagged nodes and dividing by two
    can overcount when unary recombination annotations survive conversion as
    unpaired nodes.
    """

    tables = getattr(ts, "tables", None)
    nodes = getattr(tables, "nodes", None)
    edges = getattr(tables, "edges", None)
    if nodes is None or not hasattr(nodes, "flags"):
        return {
            "events": [],
            "n_events": 0,
            "recomb_node_ids": [],
            "n_recomb_nodes": 0,
            "paired_recomb_node_ids": [],
            "unpaired_recomb_node_ids": [],
            "children_with_unpaired_recomb_parents": [],
        }

    flags = np.asarray(nodes.flags, dtype=np.uint64)
    recomb_node_ids = [
        int(idx)
        for idx in np.where((flags & np.uint64(RECOMBINANT_FLAG)) != 0)[0]
    ]

    if edges is None or len(edges.child) == 0:
        return {
            "events": [],
            "n_events": 0,
            "recomb_node_ids": recomb_node_ids,
            "n_recomb_nodes": len(recomb_node_ids),
            "paired_recomb_node_ids": [],
            "unpaired_recomb_node_ids": recomb_node_ids,
            "children_with_unpaired_recomb_parents": [],
        }

    child_to_edge_indexes = {}
    for edge_index, child in enumerate(edges.child):
        child_to_edge_indexes.setdefault(int(child), []).append(edge_index)

    events = []
    children_with_unpaired = []
    paired_recomb_node_ids = set()
    for child, edge_indexes in sorted(child_to_edge_indexes.items()):
        recomb_parent_ids = sorted({
            int(edges.parent[edge_index])
            for edge_index in edge_indexes
            if is_recombinant(flags[int(edges.parent[edge_index])])
        })
        if len(recomb_parent_ids) == 2:
            parent_intervals = []
            for parent_id in recomb_parent_ids:
                intervals = [
                    (float(edges.left[edge_index]), float(edges.right[edge_index]))
                    for edge_index in edge_indexes
                    if int(edges.parent[edge_index]) == parent_id
                ]
                parent_intervals.append(
                    {
                        "parent": parent_id,
                        "intervals": intervals,
                    }
                )
            events.append(
                {
                    "child": int(child),
                    "parents": recomb_parent_ids,
                    "parent_intervals": parent_intervals,
                }
            )
            paired_recomb_node_ids.update(recomb_parent_ids)
        elif recomb_parent_ids:
            children_with_unpaired.append(
                {
                    "child": int(child),
                    "parents": recomb_parent_ids,
                }
            )

    unpaired_recomb_node_ids = sorted(set(recomb_node_ids) - paired_recomb_node_ids)
    return {
        "events": events,
        "n_events": len(events),
        "recomb_node_ids": recomb_node_ids,
        "n_recomb_nodes": len(recomb_node_ids),
        "paired_recomb_node_ids": sorted(paired_recomb_node_ids),
        "unpaired_recomb_node_ids": unpaired_recomb_node_ids,
        "children_with_unpaired_recomb_parents": children_with_unpaired,
    }


def count_recombination_events(ts):
    return summarize_recombination_events(ts)["n_events"]
