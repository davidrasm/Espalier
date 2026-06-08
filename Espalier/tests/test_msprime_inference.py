import json
import math
from pathlib import Path

import pytest

msprime = pytest.importorskip("msprime")

from Espalier.ReassortmentSim import (
    simulate_exact_fixture,
    simulate_stochastic_benchmark,
    write_reassortment_dataset,
    write_stochastic_benchmark,
)
from Espalier.tests.conftest import (
    extract_pipeline_signature,
    has_raxml,
    load_inference_diagnostics,
    run_arg_analysis,
)


SIMULATION_GRID = [
    {
        "name": "flu_like_none_n16",
        "segment_lengths": [1719, 1464],
        "reassorting_boundary_indexes": [],
        "sample_size": 16,
        "boundary_rate": 0.0,
    },
    {
        "name": "flu_like_one_n16",
        "segment_lengths": [1719, 1464],
        "reassorting_boundary_indexes": [0],
        "sample_size": 16,
        "boundary_rate": 1e-5,
    },
    {
        "name": "flu_like_one_n32",
        "segment_lengths": [1719, 1464],
        "reassorting_boundary_indexes": [0],
        "sample_size": 32,
        "boundary_rate": 1e-4,
    },
    {
        "name": "toy_all_n16",
        "segment_lengths": [400, 300, 250],
        "reassorting_boundary_indexes": [0, 1],
        "sample_size": 16,
        "boundary_rate": 1e-4,
    },
    {
        "name": "toy_all_n32",
        "segment_lengths": [400, 300, 250],
        "reassorting_boundary_indexes": [0, 1],
        "sample_size": 32,
        "boundary_rate": 1e-3,
    },
]


def _boundary_metrics(true_boundaries, inferred_breakpoints):
    true_set = set(int(position) for position in true_boundaries)
    inferred_set = set(int(position) for position in inferred_breakpoints if position != 0)
    true_positives = len(true_set & inferred_set)
    precision = true_positives / len(inferred_set) if inferred_set else (1.0 if not true_set else 0.0)
    recall = true_positives / len(true_set) if true_set else (1.0 if not inferred_set else 0.0)
    return {
        "precision": precision,
        "recall": recall,
        "event_count_error": abs(len(inferred_set) - len(true_set)),
    }


def _rate_metrics(inferred_rate, generating_rate):
    if generating_rate <= 0:
        return {
            "rate_ratio": None,
            "log_error": None,
        }
    ratio = inferred_rate / generating_rate
    return {
        "rate_ratio": ratio,
        "log_error": abs(math.log(ratio)),
    }


@pytest.mark.simulation
def test_stochastic_benchmark_summary_tracks_boundary_frequency_increasing_with_rate():
    rates = [1e-5, 1e-4, 1e-3]
    realized = []
    for rate in rates:
        benchmark = simulate_stochastic_benchmark(
            [400, 300],
            [0],
            sample_size=16,
            population_size=10_000,
            mutation_rate=5e-3,
            boundary_rate=rate,
            random_seed=41,
            num_replicates=8,
        )
        realized.append(
            benchmark.summary["realized_boundary_breakpoint_frequency"]["400"]
        )
        assert "realized_topology_changed_boundary_frequency" in benchmark.summary
        assert "topology_changed_breakpoint_count_distribution" in benchmark.summary

    assert realized[0] <= realized[1] <= realized[2]


@pytest.mark.slow
@pytest.mark.requires_raxml
@pytest.mark.simulation
@pytest.mark.skipif(not has_raxml(), reason="raxml-ng is required for inference benchmarks")
def test_msprime_inference_recovers_detectable_boundary(tmp_path):
    dataset = simulate_exact_fixture(
        [80, 80],
        reassorting_boundary_indexes=[0],
        sample_size=8,
        population_size=5_000,
        mutation_rate=2e-2,
        boundary_rate=0.1,
        random_seed=1,
        max_attempts=1,
        require_topology_change=True,
    )
    written = write_reassortment_dataset(dataset, str(tmp_path / "detectable_fixture"))

    output_dir, completed = run_arg_analysis(
        tmp_path / "detectable_fixture_analysis",
        Path(written["alignments_dir"]),
        Path(written["trees_dir"]),
        extra_args=[
            "--em-iters",
            "1",
            "--rec-rate",
            "0.01",
        ],
    )

    signature = extract_pipeline_signature(output_dir)
    diagnostics = load_inference_diagnostics(output_dir)
    inferred_breakpoints = signature["summary"]["breakpoints"][1:-1]
    metrics = _boundary_metrics(
        dataset.selected_topology_changed_boundary_positions,
        inferred_breakpoints,
    )

    assert completed.returncode == 0
    assert dataset.selected_topology_changed_boundary_positions == [80]
    assert inferred_breakpoints == [80.0]
    assert metrics["precision"] == 1.0
    assert metrics["recall"] == 1.0
    assert diagnostics["rate_model_name"] == "boundary_hotspot"
    assert diagnostics["estimated_rate"] is not None


@pytest.mark.slow
@pytest.mark.requires_raxml
@pytest.mark.simulation
@pytest.mark.skipif(not has_raxml(), reason="raxml-ng is required for inference benchmarks")
@pytest.mark.parametrize("scenario", SIMULATION_GRID, ids=[scenario["name"] for scenario in SIMULATION_GRID])
def test_msprime_inference_scores_boundary_and_rate_metrics(tmp_path, scenario):
    benchmark = simulate_stochastic_benchmark(
        scenario["segment_lengths"],
        scenario["reassorting_boundary_indexes"],
        sample_size=scenario["sample_size"],
        population_size=10_000,
        mutation_rate=5e-3,
        boundary_rate=scenario["boundary_rate"],
        random_seed=97,
        num_replicates=3,
    )
    written = write_stochastic_benchmark(benchmark, str(tmp_path / scenario["name"]))
    first_replicate = Path(written["replicate_dirs"][0])
    truth = json.loads((first_replicate / "truth.json").read_text())
    benchmark_summary = json.loads(Path(written["summary_path"]).read_text())

    output_dir, completed = run_arg_analysis(
        tmp_path / f"{scenario['name']}_analysis",
        first_replicate / "alignments",
        first_replicate / "trees",
        extra_args=[
            "--em-iters",
            "1",
            "--rec-rate",
            str(max(scenario["boundary_rate"], 1e-6)),
        ],
    )

    signature = extract_pipeline_signature(output_dir)
    assert signature["has_label_map"]
    assert signature["summary"] is not None
    inferred_breakpoints = signature["summary"]["breakpoints"][1:-1]
    diagnostics = load_inference_diagnostics(output_dir)
    metrics = _boundary_metrics(truth["selected_boundary_positions"], inferred_breakpoints)
    metrics.update(_rate_metrics(diagnostics.get("estimated_rate"), truth["boundary_rate"]))

    assert 0.0 <= metrics["precision"] <= 1.0
    assert 0.0 <= metrics["recall"] <= 1.0
    assert metrics["event_count_error"] >= 0
    assert diagnostics["rate_model_name"] == "boundary_hotspot"
    assert diagnostics["rate_units"] == "per_boundary_per_generation"
    assert diagnostics["estimate_source"] in {"opt_mle", "fallback_tree_length", "em_failure"}

    realized_frequency = benchmark_summary["realized_boundary_breakpoint_frequency"]
    for value in realized_frequency.values():
        assert 0.0 <= value <= 1.0

    inferred_rate = diagnostics.get("estimated_rate")
    assert inferred_rate is not None
    if (
        truth["boundary_rate"] > 0
        and diagnostics["estimate_source"] == "opt_mle"
        and diagnostics["estimate_status"] == "interior"
    ):
        assert metrics["rate_ratio"] > 0
        assert metrics["log_error"] >= 0
