import json
from pathlib import Path

import pytest
from Bio import SeqIO

msprime = pytest.importorskip("msprime")

from Espalier.Reassortment import ANALYSIS_MODE_REASSORTMENT, prepare_analysis_inputs
from Espalier.ReassortmentSim import (
    DEFAULT_GENERATION_TIME_DAYS,
    EXACT_FIXTURE_MODE,
    SIMULATION_PROFILE_FLU_CALIBRATED,
    annual_rate_to_generation_rate,
    build_boundary_rate_map,
    load_flu_subsampled_calibration,
    simulate_exact_fixture,
    simulate_stochastic_benchmark,
    write_reassortment_dataset,
)
from Espalier.ARGNodeTypes import count_recombinant_nodes
from Espalier.tests.conftest import run_python_script


pytestmark = pytest.mark.simulation


def test_build_boundary_rate_map_creates_boundary_only_hotspots():
    rate_map = build_boundary_rate_map(
        [400, 300, 250],
        [0, 1],
        boundary_rate=0.1,
        background_rate=0.0,
    )

    assert list(rate_map.position) == [0.0, 400.0, 401.0, 700.0, 701.0, 950.0]
    assert list(rate_map.rate) == [0.0, 0.1, 0.0, 0.1, 0.0]


def test_exact_fixture_matches_requested_boundary_pattern(tmp_path):
    dataset = simulate_exact_fixture(
        [400, 300, 250],
        reassorting_boundary_indexes=[0, 1],
        sample_size=16,
        population_size=10_000,
        mutation_rate=5e-3,
        boundary_rate=0.1,
        random_seed=5,
    )

    assert dataset.mode == EXACT_FIXTURE_MODE
    assert dataset.selected_boundary_positions == [400, 700]
    assert dataset.observed_internal_breakpoints == [400, 700]
    assert dataset.topology_changed_internal_breakpoints == [400, 700]
    assert dataset.selected_topology_changed_boundary_positions == [400, 700]
    assert count_recombinant_nodes(dataset.ancestry_ts.tables.nodes.flags) > 0

    written = write_reassortment_dataset(dataset, str(tmp_path / "fixture"))
    truth = json.loads(Path(written["truth_path"]).read_text())
    assert truth["selected_boundary_positions"] == [400, 700]
    assert truth["observed_internal_breakpoints"] == [400, 700]
    assert truth["topology_changed_internal_breakpoints"] == [400, 700]
    assert truth["selected_topology_changed_boundary_positions"] == [400, 700]
    assert truth["record_full_arg"] is True
    assert truth["full_arg_node_counts"]["recombinant"] > 0

    prepared = prepare_analysis_inputs(
        written["seq_files"],
        written["tree_files"],
        str(tmp_path / "prepared"),
        analysis_mode=ANALYSIS_MODE_REASSORTMENT,
    )
    assert len(prepared.common_isolates) == 16


def test_exact_fixture_can_require_topology_changing_boundaries():
    dataset = simulate_exact_fixture(
        [400, 300, 250],
        reassorting_boundary_indexes=[0, 1],
        sample_size=16,
        population_size=10_000,
        mutation_rate=5e-3,
        boundary_rate=0.1,
        random_seed=5,
        require_topology_change=True,
    )

    assert dataset.observed_internal_breakpoints == [400, 700]
    assert dataset.topology_changed_internal_breakpoints == [400, 700]
    assert dataset.selected_topology_changed_boundary_positions == [400, 700]


def test_exact_fixture_writer_outputs_reassortment_ready_files(tmp_path):
    dataset = simulate_exact_fixture(
        [120, 80],
        reassorting_boundary_indexes=[0],
        sample_size=8,
        population_size=5_000,
        mutation_rate=1e-2,
        boundary_rate=0.05,
        random_seed=11,
    )
    written = write_reassortment_dataset(dataset, str(tmp_path / "dataset"))

    alignments_dir = Path(written["alignments_dir"])
    trees_dir = Path(written["trees_dir"])
    assert sorted(path.name for path in alignments_dir.iterdir()) == ["segment01.fasta", "segment02.fasta"]
    assert sorted(path.name for path in trees_dir.iterdir()) == [
        "segment01.fasta.rooted.tre",
        "segment02.fasta.rooted.tre",
    ]

    records = list(SeqIO.parse(str(alignments_dir / "segment01.fasta"), "fasta"))
    assert len(records) == 8
    assert records[0].id.startswith("sample0001|simulated|")


def test_simulation_cli_writes_truth_manifest_and_segments(tmp_path):
    output_dir = tmp_path / "cli_fixture"
    run_python_script(
        "scripts/simulate_reassortment_msprime.py",
        [
            "--mode",
            "exact_fixture",
            "--segment-lengths",
            "90,60",
            "--reassorting-boundaries",
            "0",
            "--sample-size",
            "8",
            "--population-size",
            "5000",
            "--boundary-rate",
            "0.05",
            "--mutation-rate",
            "0.01",
            "--output-dir",
            str(output_dir),
        ],
    )

    truth = json.loads((output_dir / "truth.json").read_text())
    assert truth["mode"] == "exact_fixture"
    assert truth["observed_internal_breakpoints"] == [90]
    assert "topology_changed_internal_breakpoints" in truth
    assert "selected_topology_changed_boundary_positions" in truth
    assert truth["record_full_arg"] is True
    assert truth["full_arg_node_counts"]["recombinant"] > 0
    assert (output_dir / "alignments" / "segment01.fasta").exists()
    assert (output_dir / "trees" / "segment01.fasta.rooted.tre").exists()


def test_flu_calibration_extracts_empirical_defaults():
    calibration = load_flu_subsampled_calibration()

    assert calibration.segment_lengths == [1719, 1464]
    assert calibration.sample_span_days == 4481
    assert calibration.mutation_rate_per_year > 0
    assert calibration.population_size_proxy > 0
    assert len(calibration.empirical_sample_dates) == 434
    assert calibration.invariant_fraction_by_segment[0] > 0.2
    assert calibration.invariant_fraction_by_segment[1] > 0.2


def test_flu_calibrated_profile_uses_heterochronous_dates_and_invariant_sites(tmp_path):
    calibration = load_flu_subsampled_calibration()
    benchmark = simulate_stochastic_benchmark(
        calibration.segment_lengths,
        [0],
        sample_size=16,
        population_size=calibration.population_size_proxy,
        boundary_rate=annual_rate_to_generation_rate(0.03, DEFAULT_GENERATION_TIME_DAYS),
        random_seed=19,
        num_replicates=1,
        simulation_profile=SIMULATION_PROFILE_FLU_CALIBRATED,
        generation_time_days=DEFAULT_GENERATION_TIME_DAYS,
        mutation_rate_per_year=calibration.mutation_rate_per_year,
    )

    dataset = benchmark.replicates[0]
    assert dataset.simulation_profile == SIMULATION_PROFILE_FLU_CALIBRATED
    assert max(dataset.sample_times_generations) > 0.0
    assert dataset.mutation_rate_per_year == pytest.approx(calibration.mutation_rate_per_year)

    written = write_reassortment_dataset(dataset, str(tmp_path / "flu_like"))
    seqs = [str(record.seq) for record in SeqIO.parse(Path(written["alignments_dir"]) / "segment01.fasta", "fasta")]
    columns = list(zip(*seqs))
    invariant_sites = sum(1 for column in columns if len(set(column)) == 1)
    assert invariant_sites > 0
