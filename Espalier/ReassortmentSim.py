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

import json
import math
import os
from collections import Counter
from dataclasses import dataclass
from datetime import date, timedelta
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

import dendropy
import tskit
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from Espalier import MAF
from Espalier.ARGNodeTypes import COMMON_ANCESTOR_FLAG, count_recombinant_nodes

try:
    import msprime
except ImportError:  # pragma: no cover - exercised by import-time environments only
    msprime = None


EXACT_FIXTURE_MODE = "exact_fixture"
STOCHASTIC_BENCHMARK_MODE = "stochastic_benchmark"
SIMULATION_PROFILE_SIMPLE = "simple"
SIMULATION_PROFILE_FLU_CALIBRATED = "flu_calibrated"
DEFAULT_BOUNDARY_RATE = math.log(2.0)
DEFAULT_BACKGROUND_RATE = 0.0
DEFAULT_MUTATION_RATE = 5e-3
DEFAULT_POPULATION_SIZE = 10_000
DEFAULT_SAMPLE_SIZE = 16
DEFAULT_START_DATE = "2020-01-01"
DEFAULT_GENERATION_TIME_DAYS = 3.0
DEFAULT_MUTATION_KAPPA = 3.0
DEFAULT_GAMMA_SHAPE = 0.5
DEFAULT_FLU_FIXTURE_DIR = str(Path(__file__).resolve().parents[1] / "TestFiles" / "FluTest_subsampled")


@dataclass(frozen=True)
class FluSimulationCalibration:
    fixture_dir: str
    generation_time_days: float
    segment_names: List[str]
    segment_lengths: List[int]
    empirical_sample_dates: List[str]
    sample_span_days: int
    mutation_rate_per_year: float
    mutation_rate_per_generation: float
    invariant_fraction_by_segment: List[float]
    pattern_fraction_by_segment: List[float]
    gap_fraction_by_segment: List[float]
    nucleotide_diversity_by_segment: List[float]
    population_size_proxy: float

    def summary(self) -> Dict[str, object]:
        return {
            "fixture_dir": self.fixture_dir,
            "generation_time_days": self.generation_time_days,
            "segment_names": self.segment_names,
            "segment_lengths": self.segment_lengths,
            "sample_span_days": self.sample_span_days,
            "mutation_rate_per_year": self.mutation_rate_per_year,
            "mutation_rate_per_generation": self.mutation_rate_per_generation,
            "invariant_fraction_by_segment": self.invariant_fraction_by_segment,
            "pattern_fraction_by_segment": self.pattern_fraction_by_segment,
            "gap_fraction_by_segment": self.gap_fraction_by_segment,
            "nucleotide_diversity_by_segment": self.nucleotide_diversity_by_segment,
            "population_size_proxy": self.population_size_proxy,
        }


@dataclass
class SimulatedDataset:
    mode: str
    simulation_profile: str
    ancestry_ts: tskit.TreeSequence
    mutated_ts: tskit.TreeSequence
    segment_lengths: List[int]
    segment_names: List[str]
    segment_intervals: List[Tuple[int, int]]
    selected_boundary_indexes: List[int]
    selected_boundary_positions: List[int]
    observed_internal_breakpoints: List[int]
    topology_changed_internal_breakpoints: List[int]
    selected_topology_changed_boundary_positions: List[int]
    rate_map_positions: List[float]
    rate_map_rates: List[float]
    sample_labels: List[str]
    sample_size: int
    population_size: float
    mutation_rate: float
    boundary_rate: float
    background_rate: float
    ancestry_seed: int
    mutation_seed: int
    sample_dates: List[str]
    sample_times_generations: List[float]
    generation_time_days: float
    mutation_model: str
    mutation_rate_per_year: Optional[float] = None
    invariant_fraction_by_segment: Optional[List[float]] = None
    flu_calibration: Optional[Dict[str, object]] = None
    ancestral_sequence: Optional[str] = None

    def truth_manifest(self) -> Dict[str, object]:
        manifest = {
            "mode": self.mode,
            "simulation_profile": self.simulation_profile,
            "segment_lengths": self.segment_lengths,
            "segment_names": self.segment_names,
            "segment_intervals": [list(interval) for interval in self.segment_intervals],
            "selected_boundary_indexes": self.selected_boundary_indexes,
            "selected_boundary_positions": self.selected_boundary_positions,
            "observed_internal_breakpoints": self.observed_internal_breakpoints,
            "topology_changed_internal_breakpoints": self.topology_changed_internal_breakpoints,
            "selected_topology_changed_boundary_positions": self.selected_topology_changed_boundary_positions,
            "sample_size": self.sample_size,
            "population_size": self.population_size,
            "mutation_rate": self.mutation_rate,
            "boundary_rate": self.boundary_rate,
            "background_rate": self.background_rate,
            "ancestry_seed": self.ancestry_seed,
            "mutation_seed": self.mutation_seed,
            "sample_dates": self.sample_dates,
            "sample_times_generations": self.sample_times_generations,
            "rate_map": {
                "position": self.rate_map_positions,
                "rate": self.rate_map_rates,
            },
            "sample_labels": self.sample_labels,
            "generation_time_days": self.generation_time_days,
            "mutation_model": self.mutation_model,
            "record_full_arg": True,
            "full_arg_node_counts": {
                "recombinant": count_recombinant_nodes(self.ancestry_ts.tables.nodes.flags),
                "common_ancestor": int(
                    np.count_nonzero(self.ancestry_ts.tables.nodes.flags & COMMON_ANCESTOR_FLAG)
                ),
            },
        }
        if self.mutation_rate_per_year is not None:
            manifest["mutation_rate_per_year"] = self.mutation_rate_per_year
        if self.invariant_fraction_by_segment is not None:
            manifest["invariant_fraction_by_segment"] = self.invariant_fraction_by_segment
        if self.flu_calibration is not None:
            manifest["flu_calibration"] = self.flu_calibration
        return manifest


@dataclass
class StochasticBenchmarkResult:
    mode: str
    replicates: List[SimulatedDataset]
    summary: Dict[str, object]


def _require_msprime() -> None:
    if msprime is None:
        raise ImportError(
            "msprime is required for reassortment simulation support. "
            "Install Espalier with the 'dev' extra or install msprime separately."
        )


def annual_rate_to_generation_rate(rate_per_year: float, generation_time_days: float) -> float:
    return float(rate_per_year) * (float(generation_time_days) / 365.25)


def generation_rate_to_annual_rate(rate_per_generation: float, generation_time_days: float) -> float:
    return float(rate_per_generation) * (365.25 / float(generation_time_days))


def _site_heterozygosity(seqs: Sequence[str]) -> float:
    total = 0.0
    for column in zip(*seqs):
        counts = Counter(column)
        n = sum(counts.values())
        total += 1.0 - sum((count / n) ** 2 for count in counts.values())
    return total / len(seqs[0])


def _estimate_root_to_tip_clock_rate(
    tree_path: Path,
    date_lookup: Dict[str, date],
) -> Optional[float]:
    taxa = dendropy.TaxonNamespace()
    tree = dendropy.Tree.get(
        path=str(tree_path),
        schema="newick",
        rooting="default-rooted",
        taxon_namespace=taxa,
        preserve_underscores=True,
    )
    tree.calc_node_root_distances(return_leaf_distances_only=False)

    xs = []
    ys = []
    for leaf in tree.leaf_node_iter():
        sample_date = date_lookup.get(leaf.taxon.label)
        if sample_date is None:
            continue
        xs.append(sample_date.toordinal())
        ys.append(leaf.root_distance)

    if len(xs) < 2:
        return None

    x_years = (np.array(xs, dtype=float) - np.mean(xs)) / 365.25
    y = np.array(ys, dtype=float)
    slope = float(np.polyfit(x_years, y, 1)[0])
    if slope <= 0:
        return None
    return slope


@lru_cache(maxsize=8)
def load_flu_subsampled_calibration(
    fixture_dir: str = DEFAULT_FLU_FIXTURE_DIR,
    generation_time_days: float = DEFAULT_GENERATION_TIME_DAYS,
) -> FluSimulationCalibration:
    fixture_path = Path(fixture_dir)
    alignment_dir = fixture_path / "alignments"
    tree_dir = fixture_path / "trees"
    alignment_files = sorted(path for path in alignment_dir.iterdir() if path.suffix in {".aln", ".fasta"})
    if len(alignment_files) < 2:
        raise ValueError(f"Expected at least two alignment files in {alignment_dir}")

    segment_names = []
    segment_lengths = []
    invariant_fraction_by_segment = []
    pattern_fraction_by_segment = []
    gap_fraction_by_segment = []
    nucleotide_diversity_by_segment = []
    clock_rates = []
    empirical_sample_dates = None

    for alignment_file in alignment_files:
        records = list(SeqIO.parse(str(alignment_file), "fasta"))
        seqs = [str(record.seq) for record in records]
        segment_names.append(alignment_file.name)
        segment_lengths.append(len(seqs[0]))
        columns = list(zip(*seqs))
        invariant_sites = sum(1 for column in columns if len(set(column)) == 1)
        pattern_fraction_by_segment.append(len(set(columns)) / len(columns))
        invariant_fraction_by_segment.append(invariant_sites / len(columns))
        gap_fraction_by_segment.append(sum(seq.count("-") for seq in seqs) / sum(len(seq) for seq in seqs))
        nucleotide_diversity_by_segment.append(_site_heterozygosity(seqs))

        date_lookup = {
            record.id: date.fromisoformat(record.id.split("|")[-1])
            for record in records
        }
        tree_path = tree_dir / f"{alignment_file.name}.rooted.tre"
        if tree_path.exists():
            clock_rate = _estimate_root_to_tip_clock_rate(tree_path, date_lookup)
            if clock_rate is not None:
                clock_rates.append(clock_rate)

        if empirical_sample_dates is None:
            empirical_sample_dates = sorted(value.isoformat() for value in date_lookup.values())

    if not empirical_sample_dates:
        raise ValueError(f"Could not extract sample dates from {alignment_dir}")
    if not clock_rates:
        raise ValueError(f"Could not estimate any positive root-to-tip clock rates from {tree_dir}")

    mutation_rate_per_year = float(np.mean(clock_rates))
    mutation_rate_per_generation = annual_rate_to_generation_rate(
        mutation_rate_per_year,
        generation_time_days,
    )
    mean_diversity = float(np.mean(nucleotide_diversity_by_segment))

    # A conservative single-population proxy derived from flu diversity under a simplified
    # constant-size model. This intentionally biases downward relative to the raw diversity
    # estimate so simulated trees do not become unrealistically deep.
    population_size_proxy = mean_diversity / (4.0 * mutation_rate_per_generation)

    sample_dates = [date.fromisoformat(value) for value in empirical_sample_dates]
    return FluSimulationCalibration(
        fixture_dir=str(fixture_path),
        generation_time_days=generation_time_days,
        segment_names=segment_names,
        segment_lengths=segment_lengths,
        empirical_sample_dates=empirical_sample_dates,
        sample_span_days=(sample_dates[-1] - sample_dates[0]).days,
        mutation_rate_per_year=mutation_rate_per_year,
        mutation_rate_per_generation=mutation_rate_per_generation,
        invariant_fraction_by_segment=invariant_fraction_by_segment,
        pattern_fraction_by_segment=pattern_fraction_by_segment,
        gap_fraction_by_segment=gap_fraction_by_segment,
        nucleotide_diversity_by_segment=nucleotide_diversity_by_segment,
        population_size_proxy=float(population_size_proxy),
    )


def _normalize_segment_lengths(segment_lengths: Sequence[int]) -> List[int]:
    lengths = [int(length) for length in segment_lengths]
    if len(lengths) < 2:
        raise ValueError("At least two segment lengths are required for reassortment simulation")
    if any(length <= 0 for length in lengths):
        raise ValueError("Segment lengths must all be positive integers")
    return lengths


def _segment_intervals(segment_lengths: Sequence[int]) -> List[Tuple[int, int]]:
    intervals = []
    start = 0
    for length in segment_lengths:
        end = start + int(length)
        intervals.append((start, end))
        start = end
    return intervals


def _internal_boundary_positions(segment_lengths: Sequence[int]) -> List[int]:
    return [interval[1] for interval in _segment_intervals(segment_lengths)[:-1]]


def _normalize_boundary_indexes(
    segment_lengths: Sequence[int],
    reassorting_boundary_indexes: Optional[Sequence[int]],
) -> List[int]:
    max_index = len(segment_lengths) - 2
    if reassorting_boundary_indexes is None:
        return list(range(max_index + 1))
    indexes = sorted({int(index) for index in reassorting_boundary_indexes})
    for index in indexes:
        if index < 0 or index > max_index:
            raise ValueError(
                f"Boundary index {index} is out of range for {len(segment_lengths)} segments"
            )
    return indexes


def _normalize_segment_names(
    segment_lengths: Sequence[int],
    segment_names: Optional[Sequence[str]] = None,
) -> List[str]:
    if segment_names is None:
        return [f"segment{idx + 1:02d}" for idx in range(len(segment_lengths))]
    names = [str(name) for name in segment_names]
    if len(names) != len(segment_lengths):
        raise ValueError("segment_names must match the number of segment lengths")
    if len(set(names)) != len(names):
        raise ValueError("segment_names must be unique")
    return names


def _compress_rate_map(position: List[float], rate: List[float]) -> Tuple[List[float], List[float]]:
    if not rate:
        return position, rate
    compressed_position = [position[0]]
    compressed_rate = []
    current_rate = rate[0]
    for idx, next_rate in enumerate(rate[1:], start=1):
        if next_rate != current_rate:
            compressed_position.append(position[idx])
            compressed_rate.append(current_rate)
            current_rate = next_rate
    compressed_position.append(position[-1])
    compressed_rate.append(current_rate)
    return compressed_position, compressed_rate


def build_boundary_rate_map(
    segment_lengths: Sequence[int],
    reassorting_boundary_indexes: Optional[Sequence[int]],
    boundary_rate: float = DEFAULT_BOUNDARY_RATE,
    background_rate: float = DEFAULT_BACKGROUND_RATE,
) -> "msprime.RateMap":
    """
        Build a boundary-only recombination map for a concatenated segmented genome.

        The map uses one-site hotspot intervals at selected segment boundaries and a
        background rate elsewhere. In the default configuration, the background rate is zero.
    """

    _require_msprime()
    lengths = _normalize_segment_lengths(segment_lengths)
    selected_indexes = set(_normalize_boundary_indexes(lengths, reassorting_boundary_indexes))
    boundary_positions = _internal_boundary_positions(lengths)
    total_length = sum(lengths)

    if not boundary_positions:
        return msprime.RateMap(position=[0, total_length], rate=[background_rate])

    position = [0.0]
    rate = []
    cursor = 0
    for index, boundary in enumerate(boundary_positions):
        position.append(float(boundary))
        rate.append(float(background_rate))

        hotspot_end = min(boundary + 1, total_length)
        position.append(float(hotspot_end))
        hotspot_rate = boundary_rate if index in selected_indexes else background_rate
        rate.append(float(hotspot_rate))
        cursor = hotspot_end

    if cursor < total_length:
        position.append(float(total_length))
        rate.append(float(background_rate))

    position, rate = _compress_rate_map(position, rate)
    return msprime.RateMap(position=position, rate=rate)


def _build_sample_labels(sample_size: int, start_date: str = DEFAULT_START_DATE) -> List[str]:
    base_date = date.fromisoformat(start_date)
    labels = []
    for idx in range(sample_size):
        sample_date = base_date + timedelta(days=idx)
        labels.append(f"sample{idx + 1:04d}|simulated|{sample_date.isoformat()}")
    return labels


def _build_sample_labels_from_dates(sample_dates: Sequence[str]) -> List[str]:
    return [
        f"sample{idx + 1:04d}|simulated|{sample_date}"
        for idx, sample_date in enumerate(sample_dates)
    ]


def _sample_empirical_dates(
    sample_size: int,
    empirical_sample_dates: Sequence[str],
    random_seed: int,
) -> List[str]:
    rng = np.random.default_rng(random_seed)
    ordinals = np.array(
        [date.fromisoformat(sample_date).toordinal() for sample_date in empirical_sample_dates],
        dtype=int,
    )
    sampled = np.sort(rng.choice(ordinals, size=sample_size, replace=True))
    return [date.fromordinal(int(ordinal)).isoformat() for ordinal in sampled]


def _sample_times_from_dates(
    sample_dates: Sequence[str],
    generation_time_days: float,
) -> List[float]:
    ordinals = [date.fromisoformat(sample_date).toordinal() for sample_date in sample_dates]
    latest = max(ordinals)
    return [
        (latest - ordinal) / generation_time_days
        for ordinal in ordinals
    ]


def _build_heterochronous_sample_sets(sample_times_generations: Sequence[float]) -> List["msprime.SampleSet"]:
    counts = Counter(float(time) for time in sample_times_generations)
    return [
        msprime.SampleSet(num_samples=count, time=time, ploidy=1)
        for time, count in sorted(counts.items())
    ]


def build_mutation_rate_map(
    segment_lengths: Sequence[int],
    mean_mutation_rate: float,
    *,
    invariant_fraction_by_segment: Optional[Sequence[float]] = None,
    gamma_shape: float = DEFAULT_GAMMA_SHAPE,
    random_seed: int = 1,
) -> "msprime.RateMap":
    _require_msprime()
    lengths = _normalize_segment_lengths(segment_lengths)
    if invariant_fraction_by_segment is None:
        invariant_fraction_by_segment = [0.0] * len(lengths)
    fractions = [float(value) for value in invariant_fraction_by_segment]
    if len(fractions) != len(lengths):
        raise ValueError("invariant_fraction_by_segment must match the number of segments")
    if gamma_shape <= 0:
        raise ValueError("gamma_shape must be positive")

    rng = np.random.default_rng(random_seed)
    rates = []
    for seg_length, invariant_fraction in zip(lengths, fractions):
        active_fraction = max(0.0, min(1.0, 1.0 - invariant_fraction))
        segment_rates = np.zeros(int(seg_length), dtype=float)
        if active_fraction > 0.0:
            n_active = max(1, int(round(seg_length * active_fraction)))
            active_idx = rng.choice(seg_length, size=n_active, replace=False)
            relative_rates = rng.gamma(shape=gamma_shape, scale=1.0 / gamma_shape, size=n_active)
            segment_rates[active_idx] = relative_rates * (mean_mutation_rate / active_fraction)
        rates.extend(segment_rates.tolist())
    positions = [float(idx) for idx in range(sum(lengths) + 1)]
    return msprime.RateMap(position=positions, rate=rates)


def _build_ancestral_sequence(
    sequence_length: int,
    random_seed: int,
    root_distribution: Optional[Sequence[float]] = None,
) -> str:
    rng = np.random.default_rng(random_seed)
    bases = np.array(list("ACGT"))
    if root_distribution is None:
        probs = np.full(4, 0.25)
    else:
        probs = np.array(root_distribution, dtype=float)
        probs = probs / probs.sum()
    return "".join(rng.choice(bases, size=sequence_length, p=probs))


def _render_full_haplotypes(dataset: SimulatedDataset) -> List[str]:
    sequence_length = sum(dataset.segment_lengths)
    if not dataset.ancestral_sequence or len(dataset.ancestral_sequence) != sequence_length:
        raise ValueError("SimulatedDataset is missing a full-length ancestral sequence")

    sequences = [list(dataset.ancestral_sequence) for _ in range(dataset.sample_size)]
    for variant in dataset.mutated_ts.variants():
        position = int(variant.site.position)
        ancestral_state = variant.site.ancestral_state
        for sample_index in range(dataset.sample_size):
            sequences[sample_index][position] = ancestral_state
        for sample_index, allele_index in enumerate(variant.genotypes):
            allele = variant.alleles[allele_index]
            if allele is not None:
                sequences[sample_index][position] = allele
    return ["".join(sequence) for sequence in sequences]


def _extract_internal_breakpoints(ts: tskit.TreeSequence) -> List[int]:
    breakpoints = list(ts.breakpoints())
    return [int(position) for position in breakpoints[1:-1]]


def _tskit_tree_to_dendropy_tree(
    tree: tskit.Tree,
    taxon_namespace: dendropy.TaxonNamespace,
) -> dendropy.Tree:
    return dendropy.Tree.get(
        data=tree.newick(root=tree.roots[0]),
        schema="newick",
        rooting="default-rooted",
        taxon_namespace=taxon_namespace,
        preserve_underscores=True,
    )


def _extract_topology_changed_breakpoints(ts: tskit.TreeSequence) -> List[int]:
    topology_changed = []
    for breakpoint in _extract_internal_breakpoints(ts):
        taxon_namespace = dendropy.TaxonNamespace()
        left_tree = _tskit_tree_to_dendropy_tree(ts.at(breakpoint - 1), taxon_namespace)
        right_tree = _tskit_tree_to_dendropy_tree(ts.at(breakpoint), taxon_namespace)
        if MAF.test_discordance(left_tree, right_tree):
            topology_changed.append(breakpoint)
    return topology_changed


def _simulate_dataset(
    mode: str,
    segment_lengths: Sequence[int],
    reassorting_boundary_indexes: Optional[Sequence[int]],
    segment_names: Optional[Sequence[str]],
    sample_size: int,
    population_size: float,
    mutation_rate: float,
    boundary_rate: float,
    background_rate: float,
    ancestry_seed: int,
    mutation_seed: int,
    start_date: str,
    simulation_profile: str = SIMULATION_PROFILE_SIMPLE,
    generation_time_days: float = DEFAULT_GENERATION_TIME_DAYS,
    flu_fixture_dir: str = DEFAULT_FLU_FIXTURE_DIR,
    mutation_rate_per_year: Optional[float] = None,
    mutation_kappa: float = DEFAULT_MUTATION_KAPPA,
    invariant_fraction_by_segment: Optional[Sequence[float]] = None,
    gamma_shape: float = DEFAULT_GAMMA_SHAPE,
) -> SimulatedDataset:
    lengths = _normalize_segment_lengths(segment_lengths)
    names = _normalize_segment_names(lengths, segment_names)
    boundary_indexes = _normalize_boundary_indexes(lengths, reassorting_boundary_indexes)
    rate_map = build_boundary_rate_map(
        lengths,
        boundary_indexes,
        boundary_rate=boundary_rate,
        background_rate=background_rate,
    )
    calibration_summary = None
    if simulation_profile == SIMULATION_PROFILE_FLU_CALIBRATED:
        calibration = load_flu_subsampled_calibration(
            fixture_dir=flu_fixture_dir,
            generation_time_days=generation_time_days,
        )
        sample_dates = _sample_empirical_dates(
            sample_size,
            calibration.empirical_sample_dates,
            random_seed=ancestry_seed,
        )
        sample_times_generations = _sample_times_from_dates(
            sample_dates,
            generation_time_days,
        )
        ancestry_samples = _build_heterochronous_sample_sets(sample_times_generations)
        mean_mutation_rate = (
            annual_rate_to_generation_rate(mutation_rate_per_year, generation_time_days)
            if mutation_rate_per_year is not None
            else calibration.mutation_rate_per_generation
        )
        effective_invariant_fractions = list(
            invariant_fraction_by_segment
            if invariant_fraction_by_segment is not None
            else (
                calibration.invariant_fraction_by_segment
                if len(calibration.invariant_fraction_by_segment) == len(lengths)
                else [float(np.mean(calibration.invariant_fraction_by_segment))] * len(lengths)
            )
        )
        mutation_rate_map = build_mutation_rate_map(
            lengths,
            mean_mutation_rate,
            invariant_fraction_by_segment=effective_invariant_fractions,
            gamma_shape=gamma_shape,
            random_seed=mutation_seed + 10_000,
        )
        mutation_model = msprime.HKY(kappa=mutation_kappa)
        mutation_model_name = f"HKY(kappa={mutation_kappa})"
        calibration_summary = calibration.summary()
        mutation_rate_year = (
            mutation_rate_per_year
            if mutation_rate_per_year is not None
            else calibration.mutation_rate_per_year
        )
    else:
        sample_dates = [
            date.fromisoformat(label.split("|")[-1]).isoformat()
            for label in _build_sample_labels(sample_size, start_date=start_date)
        ]
        sample_times_generations = [0.0] * sample_size
        ancestry_samples = sample_size
        effective_invariant_fractions = None
        mean_mutation_rate = mutation_rate
        mutation_rate_map = mutation_rate
        mutation_model = msprime.JC69()
        mutation_model_name = "JC69"
        mutation_rate_year = None

    ancestry_ts = msprime.sim_ancestry(
        samples=ancestry_samples,
        ploidy=1,
        recombination_rate=rate_map,
        sequence_length=sum(lengths),
        discrete_genome=True,
        population_size=population_size,
        random_seed=ancestry_seed,
        record_full_arg=True,
    )
    mutated_ts = msprime.sim_mutations(
        ancestry_ts,
        rate=mutation_rate_map,
        random_seed=mutation_seed,
        model=mutation_model,
        discrete_genome=True,
    )
    ancestral_sequence = _build_ancestral_sequence(
        sum(lengths),
        random_seed=mutation_seed + 20_000,
    )

    boundary_positions = _internal_boundary_positions(lengths)
    selected_positions = [boundary_positions[index] for index in boundary_indexes]
    topology_changed_breakpoints = _extract_topology_changed_breakpoints(ancestry_ts)
    selected_topology_changed_positions = [
        position for position in selected_positions
        if position in topology_changed_breakpoints
    ]
    sample_labels = (
        _build_sample_labels_from_dates(sample_dates)
        if simulation_profile == SIMULATION_PROFILE_FLU_CALIBRATED
        else _build_sample_labels(sample_size, start_date=start_date)
    )
    return SimulatedDataset(
        mode=mode,
        simulation_profile=simulation_profile,
        ancestry_ts=ancestry_ts,
        mutated_ts=mutated_ts,
        segment_lengths=list(lengths),
        segment_names=names,
        segment_intervals=_segment_intervals(lengths),
        selected_boundary_indexes=boundary_indexes,
        selected_boundary_positions=selected_positions,
        observed_internal_breakpoints=_extract_internal_breakpoints(ancestry_ts),
        topology_changed_internal_breakpoints=topology_changed_breakpoints,
        selected_topology_changed_boundary_positions=selected_topology_changed_positions,
        rate_map_positions=[float(position) for position in rate_map.position],
        rate_map_rates=[float(rate) for rate in rate_map.rate],
        sample_labels=sample_labels,
        sample_size=sample_size,
        population_size=population_size,
        mutation_rate=mean_mutation_rate,
        boundary_rate=boundary_rate,
        background_rate=background_rate,
        ancestry_seed=ancestry_seed,
        mutation_seed=mutation_seed,
        sample_dates=sample_dates,
        sample_times_generations=sample_times_generations,
        generation_time_days=generation_time_days,
        mutation_model=mutation_model_name,
        mutation_rate_per_year=mutation_rate_year,
        invariant_fraction_by_segment=effective_invariant_fractions,
        flu_calibration=calibration_summary,
        ancestral_sequence=ancestral_sequence,
    )


def simulate_exact_fixture(
    segment_lengths: Sequence[int],
    reassorting_boundary_indexes: Optional[Sequence[int]] = None,
    *,
    segment_names: Optional[Sequence[str]] = None,
    sample_size: int = DEFAULT_SAMPLE_SIZE,
    population_size: float = DEFAULT_POPULATION_SIZE,
    mutation_rate: float = DEFAULT_MUTATION_RATE,
    boundary_rate: float = DEFAULT_BOUNDARY_RATE,
    background_rate: float = DEFAULT_BACKGROUND_RATE,
    random_seed: int = 13,
    max_attempts: int = 100,
    start_date: str = DEFAULT_START_DATE,
    simulation_profile: str = SIMULATION_PROFILE_SIMPLE,
    generation_time_days: float = DEFAULT_GENERATION_TIME_DAYS,
    flu_fixture_dir: str = DEFAULT_FLU_FIXTURE_DIR,
    mutation_rate_per_year: Optional[float] = None,
    mutation_kappa: float = DEFAULT_MUTATION_KAPPA,
    invariant_fraction_by_segment: Optional[Sequence[float]] = None,
    gamma_shape: float = DEFAULT_GAMMA_SHAPE,
    require_topology_change: bool = False,
) -> SimulatedDataset:
    """
        Generate a deterministic fixture whose internal breakpoints exactly match the
        requested reassorting segment boundaries.
    """

    _require_msprime()
    lengths = _normalize_segment_lengths(segment_lengths)
    boundary_indexes = _normalize_boundary_indexes(lengths, reassorting_boundary_indexes)
    expected_positions = [
        _internal_boundary_positions(lengths)[index]
        for index in boundary_indexes
    ]

    for attempt in range(max_attempts):
        ancestry_seed = random_seed + (attempt * 2)
        mutation_seed = ancestry_seed + 1
        dataset = _simulate_dataset(
            EXACT_FIXTURE_MODE,
            lengths,
            boundary_indexes,
            segment_names,
            sample_size,
            population_size,
            mutation_rate,
            boundary_rate,
            background_rate,
            ancestry_seed,
            mutation_seed,
            start_date,
            simulation_profile=simulation_profile,
            generation_time_days=generation_time_days,
            flu_fixture_dir=flu_fixture_dir,
            mutation_rate_per_year=mutation_rate_per_year,
            mutation_kappa=mutation_kappa,
            invariant_fraction_by_segment=invariant_fraction_by_segment,
            gamma_shape=gamma_shape,
        )
        has_expected_breakpoints = dataset.observed_internal_breakpoints == expected_positions
        has_expected_topology_changes = (
            dataset.topology_changed_internal_breakpoints == expected_positions
        )
        if has_expected_breakpoints and (
            not require_topology_change or has_expected_topology_changes
        ):
            return dataset

    raise RuntimeError(
        "Could not simulate an exact reassortment fixture after "
        f"{max_attempts} attempts; expected breakpoints {expected_positions}"
        + (" with topology changes" if require_topology_change else "")
    )


def simulate_stochastic_benchmark(
    segment_lengths: Sequence[int],
    reassorting_boundary_indexes: Optional[Sequence[int]] = None,
    *,
    segment_names: Optional[Sequence[str]] = None,
    sample_size: int = DEFAULT_SAMPLE_SIZE,
    population_size: float = DEFAULT_POPULATION_SIZE,
    mutation_rate: float = DEFAULT_MUTATION_RATE,
    boundary_rate: float = DEFAULT_BOUNDARY_RATE,
    background_rate: float = DEFAULT_BACKGROUND_RATE,
    random_seed: int = 101,
    num_replicates: int = 5,
    start_date: str = DEFAULT_START_DATE,
    simulation_profile: str = SIMULATION_PROFILE_SIMPLE,
    generation_time_days: float = DEFAULT_GENERATION_TIME_DAYS,
    flu_fixture_dir: str = DEFAULT_FLU_FIXTURE_DIR,
    mutation_rate_per_year: Optional[float] = None,
    mutation_kappa: float = DEFAULT_MUTATION_KAPPA,
    invariant_fraction_by_segment: Optional[Sequence[float]] = None,
    gamma_shape: float = DEFAULT_GAMMA_SHAPE,
) -> StochasticBenchmarkResult:
    """
        Simulate a benchmark panel without exact-pattern conditioning.
    """

    _require_msprime()
    lengths = _normalize_segment_lengths(segment_lengths)
    boundary_indexes = _normalize_boundary_indexes(lengths, reassorting_boundary_indexes)
    boundary_positions = _internal_boundary_positions(lengths)
    selected_positions = [boundary_positions[index] for index in boundary_indexes]

    replicates = []
    for replicate_index in range(num_replicates):
        ancestry_seed = random_seed + (replicate_index * 2)
        mutation_seed = ancestry_seed + 1
        replicates.append(
            _simulate_dataset(
                STOCHASTIC_BENCHMARK_MODE,
                lengths,
                boundary_indexes,
                segment_names,
                sample_size,
                population_size,
                mutation_rate,
                boundary_rate,
                background_rate,
                ancestry_seed,
                mutation_seed,
                start_date,
                simulation_profile=simulation_profile,
                generation_time_days=generation_time_days,
                flu_fixture_dir=flu_fixture_dir,
                mutation_rate_per_year=mutation_rate_per_year,
                mutation_kappa=mutation_kappa,
                invariant_fraction_by_segment=invariant_fraction_by_segment,
                gamma_shape=gamma_shape,
            )
        )

    breakpoint_counts = Counter(
        len(dataset.observed_internal_breakpoints)
        for dataset in replicates
    )
    topology_changed_breakpoint_counts = Counter(
        len(dataset.topology_changed_internal_breakpoints)
        for dataset in replicates
    )
    realized_frequency = {}
    realized_topology_changed_frequency = {}
    for position in boundary_positions:
        realized_frequency[str(position)] = (
            sum(position in dataset.observed_internal_breakpoints for dataset in replicates)
            / len(replicates)
        )
        realized_topology_changed_frequency[str(position)] = (
            sum(position in dataset.topology_changed_internal_breakpoints for dataset in replicates)
            / len(replicates)
        )

    summary = {
        "mode": STOCHASTIC_BENCHMARK_MODE,
        "simulation_profile": simulation_profile,
        "segment_lengths": list(lengths),
        "segment_names": _normalize_segment_names(lengths, segment_names),
        "selected_boundary_indexes": boundary_indexes,
        "selected_boundary_positions": selected_positions,
        "sample_size": sample_size,
        "population_size": population_size,
        "mutation_rate": replicates[0].mutation_rate,
        "mutation_rate_per_year": replicates[0].mutation_rate_per_year,
        "boundary_rate": boundary_rate,
        "background_rate": background_rate,
        "generation_time_days": generation_time_days,
        "num_replicates": num_replicates,
        "realized_boundary_breakpoint_frequency": realized_frequency,
        "realized_topology_changed_boundary_frequency": realized_topology_changed_frequency,
        "breakpoint_count_distribution": {
            str(count): occurrences
            for count, occurrences in sorted(breakpoint_counts.items())
        },
        "topology_changed_breakpoint_count_distribution": {
            str(count): occurrences
            for count, occurrences in sorted(topology_changed_breakpoint_counts.items())
        },
    }
    if replicates and replicates[0].flu_calibration is not None:
        summary["flu_calibration"] = replicates[0].flu_calibration
    return StochasticBenchmarkResult(
        mode=STOCHASTIC_BENCHMARK_MODE,
        replicates=replicates,
        summary=summary,
    )


def _tree_to_labeled_dendropy_tree(
    tree: tskit.Tree,
    sample_labels: Sequence[str],
) -> dendropy.Tree:
    label_map = {
        str(int(sample_id) + 1): sample_labels[index]
        for index, sample_id in enumerate(tree.tree_sequence.samples())
    }
    rooted_newick = tree.newick(root=tree.roots[0])
    dendro_tree = dendropy.Tree.get(
        data=rooted_newick,
        schema="newick",
        rooting="default-rooted",
        preserve_underscores=True,
    )
    for leaf in dendro_tree.leaf_node_iter():
        leaf.taxon.label = label_map.get(leaf.taxon.label, leaf.taxon.label)
    return dendro_tree


def write_reassortment_dataset(dataset: SimulatedDataset, output_dir: str) -> Dict[str, object]:
    """
        Write a simulated dataset to a reassortment-ready directory structure.
    """

    output_path = Path(output_dir)
    align_dir = output_path / "alignments"
    tree_dir = output_path / "trees"
    align_dir.mkdir(parents=True, exist_ok=True)
    tree_dir.mkdir(parents=True, exist_ok=True)

    haplotypes = _render_full_haplotypes(dataset)
    seq_files = []
    tree_files = []
    for segment_name, (start, end) in zip(dataset.segment_names, dataset.segment_intervals):
        seq_file = align_dir / f"{segment_name}.fasta"
        records = []
        for label, haplotype in zip(dataset.sample_labels, haplotypes):
            records.append(
                SeqRecord(
                    Seq(haplotype[start:end]),
                    id=label,
                    name=label,
                    description=label,
                )
            )
        SeqIO.write(records, seq_file, "fasta")
        seq_files.append(str(seq_file))

        tree = dataset.ancestry_ts.at(start)
        dendro_tree = _tree_to_labeled_dendropy_tree(tree, dataset.sample_labels)
        tree_file = tree_dir / f"{segment_name}.fasta.rooted.tre"
        dendro_tree.write(
            path=str(tree_file),
            schema="newick",
            suppress_annotations=True,
            suppress_rooting=True,
        )
        tree_files.append(str(tree_file))

    truth_path = output_path / "truth.json"
    with open(truth_path, "w") as handle:
        json.dump(dataset.truth_manifest(), handle, indent=2, sort_keys=True)
        handle.write("\n")

    return {
        "output_dir": str(output_path),
        "alignments_dir": str(align_dir),
        "trees_dir": str(tree_dir),
        "truth_path": str(truth_path),
        "seq_files": seq_files,
        "tree_files": tree_files,
    }


def write_stochastic_benchmark(result: StochasticBenchmarkResult, output_dir: str) -> Dict[str, object]:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    replicate_dirs = []
    for replicate_index, dataset in enumerate(result.replicates, start=1):
        replicate_dir = output_path / f"replicate_{replicate_index:03d}"
        write_reassortment_dataset(dataset, str(replicate_dir))
        replicate_dirs.append(str(replicate_dir))

    summary_path = output_path / "benchmark_summary.json"
    with open(summary_path, "w") as handle:
        json.dump(result.summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    return {
        "output_dir": str(output_path),
        "replicate_dirs": replicate_dirs,
        "summary_path": str(summary_path),
    }
