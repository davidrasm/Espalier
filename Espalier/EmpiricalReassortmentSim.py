import csv
import json
import math
import random
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import dendropy
from Espalier import MAF

try:
    import pyvolve
except ImportError:  # pragma: no cover - exercised in environments without pyvolve
    pyvolve = None


DEFAULT_MUTATION_RATE = 0.003
DEFAULT_KAPPA = 2.75
DEFAULT_STATE_FREQS = [0.25, 0.25, 0.25, 0.25]
DEFAULT_SEGMENT_NAMES = ("segment01", "segment02")
DEFAULT_GRAFT_FRACTION = 0.9
EMPIRICAL_SIMULATION_MODE = "empirical_reassortment"


@dataclass(frozen=True)
class EmpiricalReassortmentEvent:
    event_index: int
    moved_taxa: List[str]
    source_taxa: List[str]
    target_taxa: List[str]
    clade_size: int
    branch_distance: float
    normalized_distance: float

    def to_dict(self) -> Dict[str, object]:
        return {
            "event_index": self.event_index,
            "moved_taxa": self.moved_taxa,
            "source_taxa": self.source_taxa,
            "target_taxa": self.target_taxa,
            "clade_size": self.clade_size,
            "branch_distance": self.branch_distance,
            "normalized_distance": self.normalized_distance,
        }


@dataclass(frozen=True)
class _PlannedEmpiricalEvent:
    source_taxa: List[str]
    target_taxa: List[str]
    branch_distance: float
    normalized_distance: float


@dataclass(frozen=True)
class EmpiricalReassortmentDataset:
    mode: str
    source_tree_path: str
    reference_tree_newick: str
    reassorted_tree_newick: str
    segment_lengths: List[int]
    segment_names: List[str]
    reassorted_segment_indexes: List[int]
    selected_boundary_indexes: List[int]
    selected_boundary_positions: List[int]
    num_reassortments: int
    min_clade_size: int
    max_clade_size: int
    min_normalized_distance: float
    max_normalized_distance: float
    mutation_rate: float
    kappa: float
    random_seed: int
    sequence_state_frequencies: List[float]
    events: List[EmpiricalReassortmentEvent]

    def truth_manifest(self) -> Dict[str, object]:
        return {
            "mode": self.mode,
            "source_tree_path": self.source_tree_path,
            "segment_lengths": self.segment_lengths,
            "segment_names": self.segment_names,
            "reassorted_segment_indexes": self.reassorted_segment_indexes,
            "selected_boundary_indexes": self.selected_boundary_indexes,
            "selected_boundary_positions": self.selected_boundary_positions,
            "num_reassortments": self.num_reassortments,
            "min_clade_size": self.min_clade_size,
            "max_clade_size": self.max_clade_size,
            "min_normalized_distance": self.min_normalized_distance,
            "max_normalized_distance": self.max_normalized_distance,
            "mutation_rate": self.mutation_rate,
            "kappa": self.kappa,
            "random_seed": self.random_seed,
            "sequence_state_frequencies": self.sequence_state_frequencies,
            "events": [event.to_dict() for event in self.events],
        }


def _require_pyvolve() -> None:
    if pyvolve is None:
        raise ImportError(
            "pyvolve is required for empirical reassortment sequence simulation. "
            "Install pyvolve to use this workflow."
        )


def _normalize_segment_lengths(segment_lengths: Sequence[int]) -> List[int]:
    lengths = [int(length) for length in segment_lengths]
    if len(lengths) < 2:
        raise ValueError("Empirical reassortment simulation requires at least two segments")
    if any(length <= 0 for length in lengths):
        raise ValueError("Segment lengths must all be positive integers")
    return lengths


def _normalize_segment_names(
    segment_lengths: Sequence[int],
    segment_names: Optional[Sequence[str]],
) -> List[str]:
    if segment_names is None:
        if len(segment_lengths) == len(DEFAULT_SEGMENT_NAMES):
            return list(DEFAULT_SEGMENT_NAMES)
        return [f"segment{idx + 1:02d}" for idx in range(len(segment_lengths))]
    names = [str(name) for name in segment_names]
    if len(names) != len(segment_lengths):
        raise ValueError("segment_names must match the number of segment lengths")
    if len(set(names)) != len(names):
        raise ValueError("segment_names must be unique")
    return names


def _normalize_reassorted_segment_indexes(
    segment_lengths: Sequence[int],
    reassorted_segment_indexes: Optional[Sequence[int]],
) -> List[int]:
    segment_count = len(segment_lengths)
    if reassorted_segment_indexes is None:
        indexes = [segment_count - 1]
    else:
        indexes = sorted({int(index) for index in reassorted_segment_indexes})
    if not indexes:
        raise ValueError("At least one reassorted segment index is required")
    if any(index < 0 or index >= segment_count for index in indexes):
        raise ValueError("Reassorted segment index is out of range")
    if len(indexes) == segment_count:
        raise ValueError("At least one segment must remain on the reference tree")
    return indexes


def _boundary_positions(segment_lengths: Sequence[int]) -> List[int]:
    positions = []
    total = 0
    for length in segment_lengths[:-1]:
        total += int(length)
        positions.append(total)
    return positions


def _selected_boundaries(
    segment_lengths: Sequence[int],
    reassorted_segment_indexes: Sequence[int],
) -> Tuple[List[int], List[int]]:
    boundary_positions = _boundary_positions(segment_lengths)
    reassorted = set(reassorted_segment_indexes)
    boundary_indexes = []
    for idx in range(len(segment_lengths) - 1):
        left_state = idx in reassorted
        right_state = (idx + 1) in reassorted
        if left_state != right_state:
            boundary_indexes.append(idx)
    return boundary_indexes, [boundary_positions[idx] for idx in boundary_indexes]


def _read_tree(tree_path: str, taxon_namespace: Optional[dendropy.TaxonNamespace] = None) -> dendropy.Tree:
    return dendropy.Tree.get(
        path=tree_path,
        schema="newick",
        rooting="default-rooted",
        preserve_underscores=True,
        taxon_namespace=taxon_namespace,
    )


def _write_tree(tree: dendropy.Tree, tree_path: Path) -> None:
    for node in tree.preorder_node_iter():
        if node.edge_length is None:
            node.edge_length = 0.0
    newick = tree.as_string(
        schema="newick",
        suppress_annotations=True,
        suppress_rooting=True,
    ).strip()
    if newick.endswith(";"):
        newick = f"{newick[:-1]}:0.0;"
    tree_path.write_text(newick)


def _sorted_leaf_labels(node: dendropy.Node) -> List[str]:
    return sorted(leaf.taxon.label for leaf in node.leaf_iter())


def _iter_ancestors(node: dendropy.Node) -> Iterable[dendropy.Node]:
    current = node.parent_node
    while current is not None:
        yield current
        current = current.parent_node


def _subtree_and_ancestor_nodes(node: dendropy.Node) -> set[dendropy.Node]:
    blocked = set(node.preorder_iter())
    blocked.update(_iter_ancestors(node))
    return blocked


def _node_distance(node_a: dendropy.Node, node_b: dendropy.Node) -> float:
    distance_to_ancestor: Dict[dendropy.Node, float] = {node_a: 0.0}
    distance = 0.0
    current = node_a
    while current.parent_node is not None:
        distance += current.edge_length or 0.0
        current = current.parent_node
        distance_to_ancestor[current] = distance

    distance = 0.0
    current = node_b
    while current not in distance_to_ancestor:
        if current.parent_node is None:
            raise ValueError("Could not find a common ancestor while computing node distance")
        distance += current.edge_length or 0.0
        current = current.parent_node

    return distance + distance_to_ancestor[current]


def _collect_candidate_nodes(
    tree: dendropy.Tree,
    min_clade_size: int,
    max_clade_size: int,
    excluded_taxa: Sequence[str],
) -> List[dendropy.Node]:
    excluded = set(excluded_taxa)
    candidates = []
    for node in tree.preorder_node_iter():
        if node.parent_node is None:
            continue
        leaf_labels = _sorted_leaf_labels(node)
        if excluded.intersection(leaf_labels):
            continue
        clade_size = len(leaf_labels)
        if min_clade_size <= clade_size <= max_clade_size:
            candidates.append(node)
    return candidates


def _valid_target_choices(
    tree: dendropy.Tree,
    source_node: dendropy.Node,
    min_normalized_distance: float,
    max_normalized_distance: float,
    blocked_nodes: Optional[Iterable[dendropy.Node]] = None,
) -> List[Tuple[dendropy.Node, float, float]]:
    source_parent = source_node.parent_node
    blocked = _subtree_and_ancestor_nodes(source_node)
    if blocked_nodes is not None:
        blocked.update(blocked_nodes)

    target_candidates = []
    distances = []
    for node in tree.preorder_node_iter():
        if node.parent_node is None or node in blocked:
            continue
        distance = _node_distance(source_parent, node)
        target_candidates.append(node)
        distances.append(distance)

    if not target_candidates:
        return []

    max_distance = max(distances)
    if max_distance <= 0.0:
        return []

    valid_candidates = []
    for node, distance in zip(target_candidates, distances):
        normalized = distance / max_distance
        if min_normalized_distance <= normalized <= max_normalized_distance:
            valid_candidates.append((node, distance, normalized))
    return valid_candidates


def _choose_target_node(
    tree: dendropy.Tree,
    source_node: dendropy.Node,
    min_normalized_distance: float,
    max_normalized_distance: float,
    rng: random.Random,
    blocked_nodes: Optional[Iterable[dendropy.Node]] = None,
) -> Optional[Tuple[dendropy.Node, float, float]]:
    valid_candidates = _valid_target_choices(
        tree,
        source_node,
        min_normalized_distance=min_normalized_distance,
        max_normalized_distance=max_normalized_distance,
        blocked_nodes=blocked_nodes,
    )

    if not valid_candidates:
        return None
    return rng.choice(valid_candidates)


def _select_non_overlapping_sources(
    tree: dendropy.Tree,
    num_reassortments: int,
    min_clade_size: int,
    max_clade_size: int,
    rng: random.Random,
) -> List[dendropy.Node]:
    candidates = _collect_candidate_nodes(
        tree,
        min_clade_size=min_clade_size,
        max_clade_size=max_clade_size,
        excluded_taxa=[],
    )
    rng.shuffle(candidates)

    selected_sources: List[dendropy.Node] = []
    included_taxa = set()
    for candidate in candidates:
        candidate_taxa = set(_sorted_leaf_labels(candidate))
        if included_taxa.intersection(candidate_taxa):
            continue
        selected_sources.append(candidate)
        included_taxa.update(candidate_taxa)
        if len(selected_sources) == num_reassortments:
            break
    return selected_sources


def _plan_empirical_reassortment_events(
    planning_tree: dendropy.Tree,
    num_reassortments: int,
    min_clade_size: int,
    max_clade_size: int,
    min_normalized_distance: float,
    max_normalized_distance: float,
    rng: random.Random,
    max_attempts: int = 250,
) -> List[_PlannedEmpiricalEvent]:
    for _ in range(max_attempts):
        selected_sources = _select_non_overlapping_sources(
            planning_tree,
            num_reassortments=num_reassortments,
            min_clade_size=min_clade_size,
            max_clade_size=max_clade_size,
            rng=rng,
        )
        if len(selected_sources) != num_reassortments:
            continue

        blocked_source_nodes = set()
        for source_node in selected_sources:
            blocked_source_nodes.update(_subtree_and_ancestor_nodes(source_node))

        def _assign_targets(
            source_nodes: Sequence[dendropy.Node],
            blocked_target_nodes: set[dendropy.Node],
            planned_events: List[_PlannedEmpiricalEvent],
        ) -> Optional[List[_PlannedEmpiricalEvent]]:
            if not source_nodes:
                return planned_events

            source_node = source_nodes[0]
            target_choices = _valid_target_choices(
                planning_tree,
                source_node,
                min_normalized_distance=min_normalized_distance,
                max_normalized_distance=max_normalized_distance,
                blocked_nodes=blocked_source_nodes | blocked_target_nodes,
            )
            rng.shuffle(target_choices)
            for target_node, branch_distance, normalized_distance in target_choices:
                next_blocked_targets = set(blocked_target_nodes)
                next_blocked_targets.update(_subtree_and_ancestor_nodes(target_node))
                next_planned_events = planned_events + [
                    _PlannedEmpiricalEvent(
                        source_taxa=_sorted_leaf_labels(source_node),
                        target_taxa=_sorted_leaf_labels(target_node),
                        branch_distance=round(branch_distance, 6),
                        normalized_distance=round(normalized_distance, 6),
                    )
                ]
                assigned = _assign_targets(
                    source_nodes[1:],
                    next_blocked_targets,
                    next_planned_events,
                )
                if assigned is not None:
                    return assigned
            return None

        planned_events = _assign_targets(selected_sources, set(), [])
        if planned_events is not None:
            return planned_events

    raise RuntimeError(
        "Could not find a valid set of non-overlapping empirical reassortment events "
        f"after {max_attempts} planning attempts"
    )


def _find_node_with_taxa(tree: dendropy.Tree, taxa: Sequence[str]) -> dendropy.Node:
    target_taxa = tuple(sorted(taxa))
    for node in tree.preorder_node_iter():
        if tuple(_sorted_leaf_labels(node)) == target_taxa:
            return node
    raise ValueError(f"Could not find a node with the requested taxa set of size {len(target_taxa)}")


def _tree_spr_distance(reference_tree: dendropy.Tree, reassorted_tree: dendropy.Tree) -> int:
    return int(MAF.get_spr_dist(reference_tree.clone(depth=2), reassorted_tree.clone(depth=2)))


def _apply_subtree_move(
    tree: dendropy.Tree,
    source_node: dendropy.Node,
    target_node: dendropy.Node,
    graft_fraction: float = DEFAULT_GRAFT_FRACTION,
) -> dendropy.Tree:
    if not (0.0 < graft_fraction < 1.0):
        raise ValueError("graft_fraction must be between 0 and 1")

    source_parent = source_node.parent_node
    if source_parent is None:
        raise ValueError("Cannot move the root node")

    subtree_root = source_node.extract_subtree()
    source_parent.remove_child(source_node)
    tree.suppress_unifurcations()

    target_parent = target_node.parent_node
    if target_parent is None:
        raise ValueError("Cannot regraft onto the root edge")

    original_target_length = target_node.edge_length or 0.0
    parent_length = original_target_length * graft_fraction
    child_length = max(original_target_length * (1.0 - graft_fraction), 1e-8)

    target_parent.remove_child(target_node)
    graft_node = target_parent.new_child(edge_length=parent_length)
    target_node.edge_length = child_length
    subtree_root.edge_length = child_length
    graft_node.add_child(target_node)
    graft_node.add_child(subtree_root)
    tree.suppress_unifurcations()
    return tree


def simulate_empirical_reassortment(
    tree_path: str,
    segment_lengths: Sequence[int],
    *,
    segment_names: Optional[Sequence[str]] = None,
    reassorted_segment_indexes: Optional[Sequence[int]] = None,
    num_reassortments: int = 1,
    min_clade_size: int = 1,
    max_clade_size: int = 10,
    min_normalized_distance: float = 0.0,
    max_normalized_distance: float = 1.0,
    mutation_rate: float = DEFAULT_MUTATION_RATE,
    kappa: float = DEFAULT_KAPPA,
    state_frequencies: Sequence[float] = DEFAULT_STATE_FREQS,
    random_seed: int = 13,
) -> EmpiricalReassortmentDataset:
    lengths = _normalize_segment_lengths(segment_lengths)
    names = _normalize_segment_names(lengths, segment_names)
    reassorted_indexes = _normalize_reassorted_segment_indexes(lengths, reassorted_segment_indexes)
    selected_boundary_indexes, selected_boundary_positions = _selected_boundaries(lengths, reassorted_indexes)
    if not selected_boundary_indexes:
        raise ValueError(
            "Chosen reassorted segments do not induce any segment boundary changes; "
            "select a subset that differs from its neighbors"
        )

    rng = random.Random(random_seed)
    taxon_namespace = dendropy.TaxonNamespace()
    base_reference_tree = _read_tree(tree_path, taxon_namespace=taxon_namespace)

    if num_reassortments <= 0:
        raise ValueError("num_reassortments must be positive")

    max_realization_attempts = 100
    last_observed_spr = None
    for _ in range(max_realization_attempts):
        reference_tree = base_reference_tree.clone(depth=2)
        reassorted_tree = base_reference_tree.clone(depth=2)
        planning_tree = base_reference_tree.clone(depth=2)
        planned_events = _plan_empirical_reassortment_events(
            planning_tree,
            num_reassortments=num_reassortments,
            min_clade_size=min_clade_size,
            max_clade_size=max_clade_size,
            min_normalized_distance=min_normalized_distance,
            max_normalized_distance=max_normalized_distance,
            rng=rng,
        )

        events: List[EmpiricalReassortmentEvent] = []
        for event_index, planned_event in enumerate(planned_events, start=1):
            source_node = _find_node_with_taxa(reassorted_tree, planned_event.source_taxa)
            target_node = _find_node_with_taxa(reassorted_tree, planned_event.target_taxa)
            moved_taxa = planned_event.source_taxa.copy()
            _apply_subtree_move(reassorted_tree, source_node, target_node)
            events.append(
                EmpiricalReassortmentEvent(
                    event_index=event_index,
                    moved_taxa=moved_taxa,
                    source_taxa=planned_event.source_taxa,
                    target_taxa=planned_event.target_taxa,
                    clade_size=len(moved_taxa),
                    branch_distance=planned_event.branch_distance,
                    normalized_distance=planned_event.normalized_distance,
                )
            )

        last_observed_spr = _tree_spr_distance(reference_tree, reassorted_tree)
        if last_observed_spr == num_reassortments:
            break
    else:
        raise RuntimeError(
            "Could not realize an empirical reassortment dataset whose true rooted SPR "
            f"distance matched the requested {num_reassortments} moves after "
            f"{max_realization_attempts} attempts; last observed SPR distance was {last_observed_spr}"
        )

    return EmpiricalReassortmentDataset(
        mode=EMPIRICAL_SIMULATION_MODE,
        source_tree_path=str(Path(tree_path).resolve()),
        reference_tree_newick=reference_tree.as_string(schema="newick").strip(),
        reassorted_tree_newick=reassorted_tree.as_string(schema="newick").strip(),
        segment_lengths=lengths,
        segment_names=names,
        reassorted_segment_indexes=reassorted_indexes,
        selected_boundary_indexes=selected_boundary_indexes,
        selected_boundary_positions=selected_boundary_positions,
        num_reassortments=num_reassortments,
        min_clade_size=min_clade_size,
        max_clade_size=max_clade_size,
        min_normalized_distance=min_normalized_distance,
        max_normalized_distance=max_normalized_distance,
        mutation_rate=float(mutation_rate),
        kappa=float(kappa),
        random_seed=int(random_seed),
        sequence_state_frequencies=[float(value) for value in state_frequencies],
        events=events,
    )


def _write_tree_from_newick(newick: str, output_path: Path) -> None:
    tree = dendropy.Tree.get(
        data=newick,
        schema="newick",
        rooting="default-rooted",
        preserve_underscores=True,
    )
    _write_tree(tree, output_path)


def _simulate_alignment(
    tree_path: Path,
    output_path: Path,
    sequence_length: int,
    mutation_rate: float,
    kappa: float,
    state_frequencies: Sequence[float],
    random_seed: int,
) -> None:
    if pyvolve is not None:
        try:
            model = pyvolve.Model(
                "nucleotide",
                {
                    "kappa": float(kappa),
                    "state_freqs": [float(value) for value in state_frequencies],
                },
            )
            partition = pyvolve.Partition(models=model, size=int(sequence_length))
            phylogeny = pyvolve.read_tree(file=str(tree_path), scale_tree=float(mutation_rate))
            evolver = pyvolve.Evolver(partitions=partition, tree=phylogeny)
            evolver(seqfile=str(output_path), ratefile=None, infofile=None)
            return
        except Exception:
            pass

    _simulate_alignment_fallback(
        tree_path=tree_path,
        output_path=output_path,
        sequence_length=sequence_length,
        mutation_rate=mutation_rate,
        kappa=kappa,
        state_frequencies=state_frequencies,
        random_seed=random_seed,
    )


def _simulate_alignment_fallback(
    tree_path: Path,
    output_path: Path,
    sequence_length: int,
    mutation_rate: float,
    kappa: float,
    state_frequencies: Sequence[float],
    random_seed: int,
) -> None:
    tree = dendropy.Tree.get(
        path=str(tree_path),
        schema="newick",
        rooting="default-rooted",
        preserve_underscores=True,
    )
    rng = random.Random(random_seed)
    bases = ["A", "C", "G", "T"]
    probs = [float(value) for value in state_frequencies]
    total = sum(probs)
    probs = [value / total for value in probs]

    root_sequence = rng.choices(bases, weights=probs, k=int(sequence_length))
    sequences: Dict[dendropy.Node, List[str]] = {tree.seed_node: root_sequence}

    transitions = {"A": "G", "G": "A", "C": "T", "T": "C"}
    transversions = {
        "A": ["C", "T"],
        "G": ["C", "T"],
        "C": ["A", "G"],
        "T": ["A", "G"],
    }
    transition_prob = float(kappa) / (float(kappa) + 2.0)

    for node in tree.preorder_node_iter():
        parent_seq = sequences[node]
        for child in node.child_node_iter():
            branch_length = float(child.edge_length or 0.0)
            change_prob = 1.0 - math.exp(-float(mutation_rate) * branch_length)
            child_seq = parent_seq.copy()
            if change_prob > 0.0:
                for idx, base in enumerate(child_seq):
                    if rng.random() < change_prob:
                        if rng.random() < transition_prob:
                            child_seq[idx] = transitions[base]
                        else:
                            child_seq[idx] = rng.choice(transversions[base])
            sequences[child] = child_seq

    with open(output_path, "w") as handle:
        for leaf in tree.leaf_node_iter():
            label = leaf.taxon.label
            seq = "".join(sequences[leaf])
            handle.write(f">{label}\n{seq}\n")


def write_empirical_reassortment_dataset(
    dataset: EmpiricalReassortmentDataset,
    output_dir: str,
) -> Dict[str, object]:
    output_path = Path(output_dir)
    align_dir = output_path / "alignments"
    tree_dir = output_path / "trees"
    align_dir.mkdir(parents=True, exist_ok=True)
    tree_dir.mkdir(parents=True, exist_ok=True)

    tree_paths = []
    seq_paths = []
    reassorted_set = set(dataset.reassorted_segment_indexes)
    for idx, (segment_name, segment_length) in enumerate(zip(dataset.segment_names, dataset.segment_lengths)):
        tree_newick = (
            dataset.reassorted_tree_newick if idx in reassorted_set else dataset.reference_tree_newick
        )
        tree_path = tree_dir / f"{segment_name}.fasta.rooted.tre"
        _write_tree_from_newick(tree_newick, tree_path)
        tree_paths.append(str(tree_path))

        seq_path = align_dir / f"{segment_name}.fasta"
        _simulate_alignment(
            tree_path=tree_path,
            output_path=seq_path,
            sequence_length=segment_length,
            mutation_rate=dataset.mutation_rate,
            kappa=dataset.kappa,
            state_frequencies=dataset.sequence_state_frequencies,
            random_seed=dataset.random_seed + idx,
        )
        seq_paths.append(str(seq_path))

    truth_path = output_path / "truth.json"
    with open(truth_path, "w") as handle:
        json.dump(dataset.truth_manifest(), handle, indent=2, sort_keys=True)
        handle.write("\n")

    return {
        "output_dir": str(output_path),
        "alignments_dir": str(align_dir),
        "trees_dir": str(tree_dir),
        "truth_path": str(truth_path),
        "seq_files": seq_paths,
        "tree_files": tree_paths,
    }


def generate_empirical_reassortment_panel(
    tree_path: str,
    output_dir: str,
    *,
    repeats: int,
    start_index: int = 1,
    segment_lengths: Sequence[int],
    segment_names: Optional[Sequence[str]] = None,
    reassorted_segment_indexes: Optional[Sequence[int]] = None,
    num_reassortments: int = 1,
    min_clade_size: int = 1,
    max_clade_size: int = 10,
    min_normalized_distance: float = 0.0,
    max_normalized_distance: float = 1.0,
    mutation_rate: float = DEFAULT_MUTATION_RATE,
    kappa: float = DEFAULT_KAPPA,
    state_frequencies: Sequence[float] = DEFAULT_STATE_FREQS,
    random_seed: int = 13,
) -> List[Dict[str, object]]:
    written_datasets = []
    for replicate_index in range(repeats):
        sample_index = start_index + replicate_index
        dataset = simulate_empirical_reassortment(
            tree_path,
            segment_lengths,
            segment_names=segment_names,
            reassorted_segment_indexes=reassorted_segment_indexes,
            num_reassortments=num_reassortments,
            min_clade_size=min_clade_size,
            max_clade_size=max_clade_size,
            min_normalized_distance=min_normalized_distance,
            max_normalized_distance=max_normalized_distance,
            mutation_rate=mutation_rate,
            kappa=kappa,
            state_frequencies=state_frequencies,
            random_seed=random_seed + replicate_index,
        )
        sample_dir = Path(output_dir) / f"sample{sample_index}"
        written_datasets.append(write_empirical_reassortment_dataset(dataset, str(sample_dir)))
    return written_datasets


def parse_csv_rows(csv_path: str) -> List[Dict[str, str]]:
    with open(csv_path, newline="") as handle:
        return list(csv.DictReader(handle))


def parse_int_list(value: str) -> List[int]:
    return [int(item.strip()) for item in value.split(",") if item.strip()]


def parse_float_list(value: str) -> List[float]:
    return [float(item.strip()) for item in value.split(",") if item.strip()]


def maybe_parse_segment_names(value: Optional[str]) -> Optional[List[str]]:
    if not value:
        return None
    return [item.strip() for item in value.split(",") if item.strip()]


def resolve_python_executable(python_executable: Optional[str] = None) -> str:
    return python_executable or sys.executable


def run_analysis_pipeline(
    dataset_dir: str,
    *,
    output_base: Optional[str] = None,
    python_executable: Optional[str] = None,
    raxml_path: str = "raxml-ng",
    em_iters: int = 1,
    generation_time_days: float = 3.0,
    skip_em: bool = False,
) -> Dict[str, object]:
    dataset_path = Path(dataset_dir)
    output_root = Path(output_base) if output_base else dataset_path / "analysis"
    output_root.mkdir(parents=True, exist_ok=True)
    log_path = dataset_path / "analysis.log"
    command = [
        resolve_python_executable(python_executable),
        str(Path(__file__).resolve().parents[1] / "scripts" / "run_arg_analysis.py"),
        "--alignments",
        str(dataset_path / "alignments"),
        "--trees",
        str(dataset_path / "trees"),
        "--output",
        str(output_root),
        "--analysis-mode",
        "reassortment",
        "--raxml-path",
        raxml_path,
        "--em-iters",
        str(em_iters),
        "--generation-time-days",
        str(generation_time_days),
    ]
    if skip_em:
        command.append("--skip-em")

    with open(log_path, "w") as log_handle:
        subprocess.run(command, check=True, stdout=log_handle, stderr=subprocess.STDOUT)

    analysis_dirs = sorted(output_root.glob("arg_output_*"))
    if not analysis_dirs:
        raise RuntimeError(f"No analysis output directory was created under {output_root}")
    analysis_dir = analysis_dirs[-1]
    return {
        "analysis_dir": str(analysis_dir),
        "analysis_log": str(log_path),
        "summary_path": str(analysis_dir / "recombination_summary.txt"),
        "diagnostics_path": str(analysis_dir / "em_diagnostics.json"),
    }
