#!/usr/bin/env python
"""
Convert a tskit TreeSequence (.trees) into an Extended Newick network
with segment annotations, suitable for IcyTree/Baltic.

Example:
  python scripts/tskit_to_extended_newick.py \
    --trees /path/to/arg_treesequence.trees \
    --output /path/to/arg_network.newick \
    --tips-from-tree /path/to/segment1_ARG_local.tre
"""
from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from typing import Dict, List, Set, Tuple, Optional, Iterable

import tskit

from Espalier.ARGNodeTypes import RECOMBINANT_FLAG


@dataclass(frozen=True)
class Edge:
    parent: int
    child: int


def _read_tip_labels_from_newick(path: str) -> List[str]:
    """Read tip labels from a Newick file and return sorted unique labels."""
    # A lightweight parser: match leaf labels by regex outside of bracketed metadata.
    # This is not a full parser but works for standard Newick tip names.
    data = open(path, "r").read()
    # Remove comments/metadata in square brackets
    data = re.sub(r"\[.*?\]", "", data)
    # Extract tip labels from patterns like '(A:0.1,B:0.2)'
    # Tip names can contain many chars; we take tokens before ':' or ',' or ')'.
    tokens = re.findall(r"[\(,]([^\(\),:;]+)", data)
    labels = [t.strip() for t in tokens if t.strip()]
    # De-quote if quoted
    cleaned = []
    for lab in labels:
        if (lab.startswith("'") and lab.endswith("'")) or (lab.startswith('"') and lab.endswith('"')):
            lab = lab[1:-1]
        cleaned.append(lab)
    # Preserve unique labels in order of appearance
    seen = set()
    unique = []
    for lab in cleaned:
        if lab in seen:
            continue
        seen.add(lab)
        unique.append(lab)
    return unique


def _format_label(label: Optional[str]) -> str:
    if label is None:
        return ""
    # Allow unquoted if only safe characters
    if re.match(r"^[A-Za-z0-9_\-\.|%/\+&]+$", label):
        return label
    # Otherwise quote with single quotes and escape existing single quotes
    safe = label.replace("'", "''")
    return f"'{safe}'"


def _format_segments(segments: Set[int]) -> str:
    return "{" + ",".join(str(i) for i in sorted(segments)) + "}"


def _format_meta_value(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return str(value)
    return str(value)


def _format_metadata(
    segments: Set[int],
    num_segments: int,
    posterior: Optional[float] = None,
    extras: Optional[Dict[str, object]] = None,
    seg_bools: bool = True,
) -> str:
    parts = [f"segments={_format_segments(segments)}", f"segsCarried={len(segments)}"]
    if posterior is not None:
        parts.append(f"posterior={posterior}")
    for i in range(num_segments):
        if seg_bools:
            parts.append(f"seg{i}={'true' if i in segments else 'false'}")
        else:
            parts.append(f"seg{i}={1 if i in segments else 0}")
    if extras:
        for key, value in extras.items():
            parts.append(f"{key}={_format_meta_value(value)}")
    return "[&" + ",".join(parts) + "]"


def _pair_recomb_nodes(recomb_nodes: List[int]) -> List[Tuple[int, int]]:
    """Pair recombination nodes deterministically (matches summary: nodes/2)."""
    recomb_nodes_sorted = sorted(recomb_nodes)
    pairs = []
    for i in range(0, len(recomb_nodes_sorted) - 1, 2):
        pairs.append((recomb_nodes_sorted[i], recomb_nodes_sorted[i + 1]))
    return pairs


def _choose_recomb_rep(a: int,
                       b: int,
                       presence: Dict[int, Set[int]],
                       follow_segment: int) -> int:
    """Choose a representative node for a recombination pair."""
    pa = presence.get(a, set())
    pb = presence.get(b, set())
    if follow_segment >= 0:
        if (follow_segment in pa) and (follow_segment not in pb):
            return a
        if (follow_segment in pb) and (follow_segment not in pa):
            return b
    if len(pa) != len(pb):
        return a if len(pa) > len(pb) else b
    # Fall back to deterministic choice
    return a


def _select_best_parents(parents: Iterable[int],
                         child: int,
                         edge_to_segs: Dict[Edge, Set[int]],
                         num_segments: int) -> List[int]:
    """Pick up to two parents that best cover segments for a hybrid node."""
    parent_list = list(parents)
    if len(parent_list) <= 2:
        return parent_list
    best_pair = None
    best_score = (-1, -1)  # (coverage, total_segments)
    for i in range(len(parent_list)):
        for j in range(i + 1, len(parent_list)):
            p1, p2 = parent_list[i], parent_list[j]
            s1 = edge_to_segs.get(Edge(p1, child), set())
            s2 = edge_to_segs.get(Edge(p2, child), set())
            coverage = len(s1.union(s2))
            total = len(s1) + len(s2)
            score = (coverage, total)
            if score > best_score:
                best_score = score
                best_pair = (p1, p2)
            if coverage == num_segments and total == coverage:
                best_pair = (p1, p2)
                best_score = score
                break
        if best_score[0] == num_segments and best_score[1] == best_score[0]:
            break
    return list(best_pair) if best_pair is not None else parent_list[:2]


def _choose_main_parent(parents: List[int],
                        child: int,
                        edge_to_segs: Dict[Edge, Set[int]],
                        follow_segment: int,
                        node_times: Dict[int, float]) -> int:
    """Choose a main parent for traversal (CoalRe-like logic)."""
    if len(parents) == 1:
        return parents[0]
    carrying = [p for p in parents if follow_segment in edge_to_segs.get(Edge(p, child), set())]
    candidates = carrying if carrying else parents

    def sort_key(p: int):
        segs = edge_to_segs.get(Edge(p, child), set())
        return (len(segs), node_times.get(p, 0.0), -p)

    return max(candidates, key=sort_key)


def convert(ts: tskit.TreeSequence,
            tip_labels: Optional[List[str]] = None,
            root_epsilon: float = 1e-6,
            posterior_value: Optional[float] = None,
            follow_segment: int = -1,
            collapse_root_edges: bool = False) -> str:
    num_segments = ts.num_trees

    # Map sample node id -> label
    sample_labels: Dict[int, str] = {}
    if tip_labels is not None:
        # Map sorted labels onto sample ids in order (ts samples are 0..n-1)
        labels_sorted = sorted(tip_labels)
        if len(labels_sorted) != ts.num_samples:
            raise ValueError(
                f"Tip label count ({len(labels_sorted)}) does not match number of samples ({ts.num_samples})"
            )
        for sid, label in zip(ts.samples(), labels_sorted):
            sample_labels[int(sid)] = label
    else:
        for sid in ts.samples():
            sample_labels[int(sid)] = str(int(sid))

    # Identify recombination nodes
    recomb_nodes = [i for i, n in enumerate(ts.nodes()) if n.flags & RECOMBINANT_FLAG]
    recomb_pairs = _pair_recomb_nodes(recomb_nodes)
    recomb_set = set(recomb_nodes)

    # Track segment presence for recombination nodes
    presence: Dict[int, Set[int]] = {n: set() for n in recomb_nodes}
    for seg_idx, tree in enumerate(ts.trees()):
        for u in tree.nodes():
            if u in recomb_set:
                presence[u].add(seg_idx)

    # Auto-select follow_segment if negative: maximize recomb-pair coverage
    if follow_segment < 0 and num_segments > 0:
        pair_counts = []
        for seg in range(num_segments):
            count = sum(1 for a, b in recomb_pairs if seg in presence.get(a, set()) or seg in presence.get(b, set()))
            pair_counts.append(count)
        follow_segment = max(range(num_segments), key=lambda s: pair_counts[s])
        print(
            f"Auto follow_segment={follow_segment} (pair coverage {pair_counts})",
            file=sys.stderr,
        )

    merge_map: Dict[int, int] = {}
    hybrid_nodes: List[int] = []
    dangling_pairs = 0
    for a, b in recomb_pairs:
        rep = _choose_recomb_rep(a, b, presence, follow_segment)
        other = b if rep == a else a
        merge_map[other] = rep
        hybrid_nodes.append(rep)
        if not presence.get(rep) and not presence.get(other):
            dangling_pairs += 1
    hybrid_nodes = sorted(set(hybrid_nodes))
    hybrid_id_by_node: Dict[int, int] = {n: i for i, n in enumerate(hybrid_nodes)}
    if len(recomb_nodes) % 2 == 1:
        print(f"Warning: odd number of recombination nodes ({len(recomb_nodes)}); last node ignored.", file=sys.stderr)
    if dangling_pairs > 0:
        print(
            f"Warning: {dangling_pairs} recombination pair(s) have no segment presence; "
            "they will be attached as dangling events.",
            file=sys.stderr,
        )

    # Build edges with segment sets
    edge_to_segs: Dict[Edge, Set[int]] = {}
    node_times: Dict[int, float] = {i: ts.node(i).time for i in range(ts.num_nodes)}
    for a, b in recomb_pairs:
        rep = merge_map.get(a, a)
        other = b if rep == a else a
        node_times[rep] = max(node_times.get(rep, 0.0), node_times.get(other, 0.0))

    for seg_idx, tree in enumerate(ts.trees()):
        for u in tree.nodes():
            parent = tree.parent(u)
            if parent != tskit.NULL:
                parent = merge_map.get(int(parent), int(parent))
                u = merge_map.get(int(u), int(u))
                if parent == u:
                    continue
                edge = Edge(int(parent), int(u))
                edge_to_segs.setdefault(edge, set()).add(seg_idx)

    # Add synthetic root
    synthetic_root = ts.num_nodes
    max_time = max(node_times.values()) if node_times else 0.0
    node_times[synthetic_root] = max_time + root_epsilon

    # Add edges from synthetic root to each segment root
    for seg_idx, tree in enumerate(ts.trees()):
        for root in tree.roots:
            edge = Edge(synthetic_root, int(root))
            edge_to_segs.setdefault(edge, set()).add(seg_idx)

    # Ensure each hybrid has two parents by adding synthetic parents if needed
    raw_parents_tmp: Dict[int, List[int]] = {}
    for edge in edge_to_segs.keys():
        raw_parents_tmp.setdefault(edge.child, []).append(edge.parent)

    all_segments = set(range(num_segments))
    next_node_id = synthetic_root + 1
    for h in hybrid_nodes:
        parents = raw_parents_tmp.get(h, [])
        if len(parents) >= 2:
            continue
        # Determine segments already covered by existing parents
        covered = set()
        for p in parents:
            covered.update(edge_to_segs.get(Edge(p, h), set()))
        missing = all_segments - covered
        if len(parents) == 0:
            # Create two synthetic parents to cover all segments
            if num_segments >= 2:
                seg_list = sorted(all_segments)
                seg_a = {seg_list[0]}
                seg_b = set(seg_list[1:]) if len(seg_list) > 1 else set()
                for segs in (seg_a, seg_b if seg_b else seg_a):
                    synth_parent = next_node_id
                    next_node_id += 1
                    node_times[synth_parent] = node_times[h] + root_epsilon
                    edge_to_segs[Edge(synth_parent, h)] = set(segs)
                    edge_to_segs.setdefault(Edge(synthetic_root, synth_parent), set(segs))
            else:
                synth_parent = next_node_id
                next_node_id += 1
                node_times[synth_parent] = node_times[h] + root_epsilon
                edge_to_segs[Edge(synth_parent, h)] = set(all_segments)
                edge_to_segs.setdefault(Edge(synthetic_root, synth_parent), set(all_segments))
        elif len(parents) == 1:
            # Add one synthetic parent to cover missing segments (or duplicate if none missing)
            segs = missing if missing else covered if covered else all_segments
            synth_parent = next_node_id
            next_node_id += 1
            node_times[synth_parent] = node_times[h] + root_epsilon
            edge_to_segs[Edge(synth_parent, h)] = set(segs)
            edge_to_segs.setdefault(Edge(synthetic_root, synth_parent), set(segs))

    # Build raw parent/child mappings
    raw_parents: Dict[int, List[int]] = {}
    raw_children: Dict[int, List[int]] = {}
    for edge in edge_to_segs.keys():
        raw_parents.setdefault(edge.child, []).append(edge.parent)
        raw_children.setdefault(edge.parent, []).append(edge.child)

    # CoalRe-style "realCoal" annotation: overlap of segments between child edges
    real_coal_by_node: Dict[int, bool] = {}
    for node, children in raw_children.items():
        if len(children) < 2:
            real_coal_by_node[node] = False
            continue
        real = False
        for i in range(len(children)):
            for j in range(i + 1, len(children)):
                s1 = edge_to_segs.get(Edge(node, children[i]), set())
                s2 = edge_to_segs.get(Edge(node, children[j]), set())
                if s1.intersection(s2):
                    real = True
                    break
            if real:
                break
        real_coal_by_node[node] = real

    # Choose a single main parent for non-hybrids to keep a tree backbone
    main_parent: Dict[int, int] = {}
    for node, parents in raw_parents.items():
        if node in hybrid_id_by_node:
            continue
        if len(parents) >= 1:
            main_parent[node] = _choose_main_parent(parents, node, edge_to_segs, follow_segment, node_times)

    # Build filtered adjacency: keep only main-parent edges for non-hybrids,
    # and up to two best parents for hybrids.
    node_children: Dict[int, List[int]] = {}
    node_parents: Dict[int, List[int]] = {}
    for edge in edge_to_segs.keys():
        parent, child = edge.parent, edge.child
        if child in hybrid_id_by_node:
            allowed_parents = _select_best_parents(raw_parents.get(child, []), child, edge_to_segs, num_segments)
            if parent not in allowed_parents:
                continue
        else:
            if main_parent.get(child, parent) != parent:
                continue
        node_children.setdefault(parent, []).append(child)
        node_parents.setdefault(child, []).append(parent)

    # Choose a main parent for hybrids based on follow_segment/cardinality
    hybrid_main_parent: Dict[int, int] = {}
    for node, parents in node_parents.items():
        if node not in hybrid_id_by_node:
            continue
        if len(parents) >= 1:
            hybrid_main_parent[node] = _choose_main_parent(parents, node, edge_to_segs, follow_segment, node_times)

    dummy_tips: Set[int] = set()
    dummy_for: Dict[int, int] = {}

    # Add dangling children for hybrid nodes with no descendants to ensure destinations render
    dangling_hybrids = [h for h in hybrid_nodes if h not in node_children]
    if dangling_hybrids:
        for h in dangling_hybrids:
            dummy = next_node_id
            next_node_id += 1
            node_times[dummy] = node_times.get(h, 0.0) - root_epsilon
            sample_labels[dummy] = f"_dangling_H{hybrid_id_by_node[h]}"
            dummy_tips.add(dummy)
            dummy_for[dummy] = hybrid_id_by_node[h]
            incoming = set()
            for p in raw_parents.get(h, []):
                incoming.update(edge_to_segs.get(Edge(p, h), set()))
            if not incoming:
                incoming = set(range(num_segments))
            edge_to_segs[Edge(h, dummy)] = incoming
            node_children.setdefault(h, []).append(dummy)
            node_parents.setdefault(dummy, []).append(h)
        print(
            f"Warning: added {len(dangling_hybrids)} dangling leaf node(s) to render hybrid destinations.",
            file=sys.stderr,
        )

    # Deterministic child order
    for parent in node_children:
        node_children[parent] = sorted(set(node_children[parent]))

    # Prune unlabeled dead-end nodes to avoid empty tips like "([&...]:len)"
    keep_seeds: Set[int] = set(sample_labels.keys()) | set(hybrid_id_by_node.keys()) | {synthetic_root}
    keep_cache: Dict[int, bool] = {}

    def keep_node(n: int) -> bool:
        if n in keep_cache:
            return keep_cache[n]
        kept = n in keep_seeds
        for c in node_children.get(n, []):
            if keep_node(c):
                kept = True
        keep_cache[n] = kept
        return kept

    keep_node(synthetic_root)
    for parent in list(node_children.keys()):
        kept_children = [c for c in node_children[parent] if keep_cache.get(c, False)]
        if kept_children:
            node_children[parent] = sorted(set(kept_children))
        else:
            del node_children[parent]

    # After pruning, ensure all hybrids still have at least one child
    post_prune_dangling = [h for h in hybrid_nodes if h not in node_children]
    if post_prune_dangling:
        for h in post_prune_dangling:
            dummy = next_node_id
            next_node_id += 1
            node_times[dummy] = node_times.get(h, 0.0) - root_epsilon
            sample_labels[dummy] = f"_dangling_H{hybrid_id_by_node[h]}"
            dummy_tips.add(dummy)
            dummy_for[dummy] = hybrid_id_by_node[h]
            incoming = set()
            for p in raw_parents.get(h, []):
                incoming.update(edge_to_segs.get(Edge(p, h), set()))
            if not incoming:
                incoming = set(range(num_segments))
            edge_to_segs[Edge(h, dummy)] = incoming
            node_children.setdefault(h, []).append(dummy)
            node_parents.setdefault(dummy, []).append(h)
        print(
            f"Warning: added {len(post_prune_dangling)} dangling leaf node(s) after pruning.",
            file=sys.stderr,
        )

    def render_edge(parent: int, child: int) -> str:
        edge = Edge(parent, child)
        segments = edge_to_segs[edge]
        is_hybrid = child in hybrid_id_by_node
        traverse = True
        if is_hybrid:
            main_hybrid_parent = hybrid_main_parent.get(child, parent)
            if parent != main_hybrid_parent:
                traverse = False
            else:
                if child in rendered_hybrids:
                    traverse = False
                else:
                    rendered_hybrids.add(child)
        else:
            if main_parent.get(child, parent) != parent:
                traverse = False

        parts: List[str] = []
        # Traverse children if allowed
        if traverse and child in node_children:
            parts.append("(")
            parts.append(",".join(render_edge(child, gc) for gc in node_children[child]))
            parts.append(")")

        # Label for leaves
        label = sample_labels.get(child)
        if is_hybrid:
            # Hybrid node marker
            parts.append(f"#H{hybrid_id_by_node[child]}")
        elif label is not None:
            parts.append(_format_label(label))

        # Metadata for this edge
        posterior = posterior_value if is_hybrid else None
        extras = {"realCoal": real_coal_by_node.get(child, False)}
        if child in dummy_tips:
            extras["dummy"] = True
            extras["dummyFor"] = f"H{dummy_for.get(child)}"
        parts.append(_format_metadata(segments, num_segments, posterior=posterior, extras=extras))

        # Branch length
        if collapse_root_edges and parent == synthetic_root:
            length = root_epsilon
        else:
            parent_time = node_times[parent]
            child_time = node_times[child]
            length = parent_time - child_time
        if length <= 0:
            length = root_epsilon
        parts.append(f":{length}")

        return "".join(parts)

    # Render from synthetic root; it should have children
    if synthetic_root not in node_children:
        raise ValueError("Synthetic root has no children; check tree sequence.")

    rendered_hybrids: Set[int] = set()
    newick_body = "(" + ",".join(render_edge(synthetic_root, c) for c in node_children[synthetic_root]) + ")" + ";"
    return newick_body


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert tskit .trees to Extended Newick network")
    parser.add_argument("--trees", required=True, help="Path to .trees file")
    parser.add_argument("--output", required=True, help="Output file (.nex or .newick)")
    parser.add_argument("--tips-from-tree", default=None, help="Optional Newick file to derive tip labels")
    parser.add_argument("--root-epsilon", type=float, default=1e-6, help="Small time offset for synthetic root")
    parser.add_argument("--posterior", type=float, default=None, help="Posterior value to attach to reticulation edges")
    parser.add_argument(
        "--follow-segment",
        type=int,
        default=-1,
        help="Segment index to prefer when traversing hybrids (-1 for auto)",
    )
    parser.add_argument(
        "--collapse-root-edges",
        action="store_true",
        help="Set synthetic-root edge lengths to epsilon to reduce long basal branches",
    )
    args = parser.parse_args()

    ts = tskit.load(args.trees)
    tip_labels = _read_tip_labels_from_newick(args.tips_from_tree) if args.tips_from_tree else None

    newick = convert(
        ts,
        tip_labels=tip_labels,
        root_epsilon=args.root_epsilon,
        posterior_value=args.posterior,
        follow_segment=args.follow_segment,
        collapse_root_edges=args.collapse_root_edges,
    )

    # Write extended Newick (no NEXUS wrapper)
    with open(args.output, "w") as fh:
        fh.write(newick)


if __name__ == "__main__":
    main()
