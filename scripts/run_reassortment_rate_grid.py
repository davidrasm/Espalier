#!/usr/bin/env python3
"""
Run a large simulated reassortment benchmark grid with temporary alignments only.

This driver:
1. Simulates flu-like segmented datasets with msprime.
2. Persists only the simulated input trees and truth metadata.
3. Runs the Espalier reassortment pipeline on each replicate.
4. Deletes all generated alignments and temp sequence artifacts.
5. Writes per-replicate and per-treatment summaries for later analysis.
"""

import argparse
import ast
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import tskit


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from Espalier.Reassortment import ANALYSIS_MODE_REASSORTMENT
from Espalier.ReassortmentSim import (
    DEFAULT_FLU_FIXTURE_DIR,
    DEFAULT_GENERATION_TIME_DAYS,
    SIMULATION_PROFILE_FLU_CALIBRATED,
    build_boundary_rate_map,
    annual_rate_to_generation_rate,
    generation_rate_to_annual_rate,
    load_flu_subsampled_calibration,
    simulate_stochastic_benchmark,
    write_reassortment_dataset,
)

try:
    from calibrate_reassortment_rate import calibrate_results_file
except ImportError:  # pragma: no cover - only used when the standalone script is unavailable
    calibrate_results_file = None


DEFAULT_SEGMENT_LENGTHS = [1719, 1464]
DEFAULT_TREE_SIZES = [10, 50, 100, 500, 1000, 1500]
DEFAULT_BOUNDARY_RATES_PER_YEAR = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3]
DEFAULT_REPLICATES = 10
DEFAULT_BACKGROUND_RATE = 0.0
DEFAULT_EM_ITERS = 1
DEFAULT_JOBS = 1
DEFAULT_RANDOM_SEED = 10_001


def _parse_int_list(text: str) -> List[int]:
    return [int(part.strip()) for part in text.split(",") if part.strip()]


def _parse_float_list(text: str) -> List[float]:
    return [float(part.strip()) for part in text.split(",") if part.strip()]


def _parse_boundaries(text: str, segment_lengths: Sequence[int]) -> List[int]:
    normalized = text.strip().lower()
    if normalized == "all":
        return list(range(len(segment_lengths) - 1))
    if normalized == "none":
        return []
    return _parse_int_list(text)


def _boundary_positions(segment_lengths: Sequence[int]) -> List[int]:
    cumulative = 0
    positions = []
    for length in segment_lengths[:-1]:
        cumulative += int(length)
        positions.append(cumulative)
    return positions


def _mean_genome_wide_rate_per_site(
    segment_lengths: Sequence[int],
    selected_boundary_indexes: Sequence[int],
    boundary_rate: float,
    background_rate: float,
) -> float:
    rate_map = build_boundary_rate_map(
        segment_lengths,
        selected_boundary_indexes,
        boundary_rate=boundary_rate,
        background_rate=background_rate,
    )
    position = list(rate_map.position)
    rates = list(rate_map.rate)
    total_length = position[-1] - position[0]
    total_mass = 0.0
    for idx, rate in enumerate(rates):
        total_mass += (position[idx + 1] - position[idx]) * rate
    return total_mass / total_length


def _boundary_equivalent_rate(
    inferred_rate_per_site: Optional[float],
    genome_length: int,
    selected_boundary_count: int,
) -> Optional[float]:
    if inferred_rate_per_site is None or selected_boundary_count <= 0:
        return None
    return inferred_rate_per_site * genome_length / selected_boundary_count


def _slug_float(value: float) -> str:
    formatted = f"{value:.6g}"
    return formatted.replace("-", "m").replace(".", "p")


def _parse_recombination_summary(summary_path: Path) -> Optional[Dict[str, object]]:
    if not summary_path.exists():
        return None
    result: Dict[str, object] = {}
    for line in summary_path.read_text().splitlines():
        if line.startswith("Total recombination events:"):
            result["n_recomb_events"] = int(line.split(": ", 1)[1])
        elif line.startswith("Breakpoints:"):
            result["breakpoints"] = ast.literal_eval(line.split(": ", 1)[1])
    return result


def _parse_logged_recombination_rate(log_text: str) -> Optional[float]:
    matches = re.findall(r"Estimated recombination rate: ([0-9.eE+-]+) per site", log_text)
    if not matches:
        return None
    return float(matches[-1])


def _extract_pipeline_signature(output_dir: Path) -> Dict[str, object]:
    signature = {
        "has_label_map": (output_dir / "derived" / "label_map.tsv").exists(),
        "has_treesequence": (output_dir / "arg_treesequence.trees").exists(),
        "summary": _parse_recombination_summary(output_dir / "recombination_summary.txt"),
        "ts_breakpoints": None,
        "ts_num_trees": None,
    }
    if signature["has_treesequence"]:
        ts = tskit.load(str(output_dir / "arg_treesequence.trees"))
        signature["ts_breakpoints"] = list(ts.breakpoints())
        signature["ts_num_trees"] = ts.num_trees
    return signature


def _load_inference_diagnostics(output_dir: Path) -> Dict[str, object]:
    diagnostics_path = output_dir / "em_diagnostics.json"
    if not diagnostics_path.exists():
        return {}
    with open(diagnostics_path) as handle:
        return json.load(handle)


def _boundary_metrics(true_boundaries: Iterable[int], inferred_breakpoints: Optional[Iterable[float]]) -> Dict[str, Optional[float]]:
    true_set = set(int(position) for position in true_boundaries)
    inferred_values = inferred_breakpoints or []
    inferred_set = set(int(position) for position in inferred_values if int(position) not in (0,))
    true_positives = len(true_set & inferred_set)

    if inferred_set:
        precision = true_positives / len(inferred_set)
    else:
        precision = 1.0 if not true_set else 0.0

    if true_set:
        recall = true_positives / len(true_set)
    else:
        recall = 1.0 if not inferred_set else 0.0

    return {
        "boundary_precision": precision,
        "boundary_recall": recall,
        "boundary_detected": float(recall > 0.0),
        "event_count_error": float(abs(len(inferred_set) - len(true_set))),
    }


def _rate_metrics(
    inferred_rate_native: Optional[float],
    expected_rate_native: float,
) -> Dict[str, Optional[float]]:
    if inferred_rate_native is None or expected_rate_native <= 0:
        return {
            "rate_ratio_vs_expected_native": None,
            "rate_log_error_vs_expected_native": None,
        }
    ratio = inferred_rate_native / expected_rate_native
    return {
        "rate_ratio_vs_expected_native": ratio,
        "rate_log_error_vs_expected_native": abs(math.log(ratio)),
    }


def _remove_sequence_artifacts(analysis_dir: Path) -> None:
    for removable_dir in [
        analysis_dir / "derived" / "normalized_alignments",
        analysis_dir / "temp",
    ]:
        if removable_dir.exists():
            shutil.rmtree(removable_dir)
    for removable_file in analysis_dir.rglob("*"):
        if removable_file.is_file() and removable_file.suffix in {".fasta", ".fa", ".aln", ".phy", ".phylip"}:
            removable_file.unlink()


def _copy_tree_inputs(source_tree_dir: Path, destination_tree_dir: Path) -> None:
    destination_tree_dir.mkdir(parents=True, exist_ok=True)
    for tree_file in sorted(source_tree_dir.glob("*.tre")):
        shutil.copy2(tree_file, destination_tree_dir / tree_file.name)


def _build_run_command(
    python_executable: str,
    alignments_dir: Path,
    trees_dir: Path,
    output_base: Path,
    initial_rec_rate: float,
    rec_rate_lower: float,
    rec_rate_upper: float,
    ne: float,
    em_iters: int,
    raxml_path: str,
    use_modern_converter: bool,
    generation_time_days: float,
) -> List[str]:
    command = [
        python_executable,
        str(REPO_ROOT / "scripts" / "run_arg_analysis.py"),
        "--alignments",
        str(alignments_dir),
        "--trees",
        str(trees_dir),
        "--output",
        str(output_base),
        "--analysis-mode",
        ANALYSIS_MODE_REASSORTMENT,
        "--rec-rate",
        str(initial_rec_rate),
        "--rec-rate-lower",
        str(rec_rate_lower),
        "--rec-rate-upper",
        str(rec_rate_upper),
        "--ne",
        str(ne),
        "--em-iters",
        str(em_iters),
        "--raxml-path",
        raxml_path,
        "--generation-time-days",
        str(generation_time_days),
    ]
    if use_modern_converter:
        command.append("--use-modern-converter")
    return command


def _write_json(path: Path, payload: Dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def _load_existing_result(result_path: Path) -> Optional[Dict[str, object]]:
    if not result_path.exists():
        return None
    with open(result_path) as handle:
        return json.load(handle)


def _select_analysis_output_dir(output_base: Path) -> Path:
    output_dirs = sorted(output_base.glob("arg_output_*"))
    if len(output_dirs) != 1:
        raise RuntimeError(
            f"Expected exactly one analysis output directory in {output_base}, found {len(output_dirs)}"
        )
    return output_dirs[0]


def _combination_seed(base_seed: int, sample_size: int, boundary_rate: float, replicate_index: int) -> int:
    rate_component = int(round(boundary_rate * 1_000_000))
    return base_seed + (sample_size * 100_000) + rate_component + replicate_index


def _replicate_root(output_dir: Path, sample_size: int, boundary_rate: float, replicate_index: int) -> Path:
    return (
        output_dir
        / f"n_{sample_size:04d}"
        / f"rate_{_slug_float(boundary_rate)}"
        / f"replicate_{replicate_index:03d}"
    )


def _run_single_replicate(
    *,
    output_dir: Path,
    sample_size: int,
    boundary_rate_annual: float,
    boundary_rate_generation: float,
    replicate_index: int,
    segment_lengths: Sequence[int],
    selected_boundary_indexes: Sequence[int],
    population_size: float,
    mutation_rate_per_year: Optional[float],
    background_rate: float,
    em_iters: int,
    ne: float,
    raxml_path: str,
    python_executable: str,
    random_seed: int,
    use_modern_converter: bool,
    timeout_seconds: Optional[int],
    rec_rate_upper_generation: Optional[float],
    force: bool,
    simulation_profile: str,
    generation_time_days: float,
    flu_fixture_dir: str,
) -> Dict[str, object]:
    replicate_root = _replicate_root(output_dir, sample_size, boundary_rate_annual, replicate_index)
    result_path = replicate_root / "result.json"
    if not force:
        existing = _load_existing_result(result_path)
        if existing is not None:
            existing["status"] = "skipped_existing"
            return existing

    replicate_root.mkdir(parents=True, exist_ok=True)
    truth_path = replicate_root / "truth.json"
    input_tree_dir = replicate_root / "input_trees"
    analysis_base = replicate_root / "analysis"
    analysis_log_path = replicate_root / "analysis.log"

    expected_rate_per_site = _mean_genome_wide_rate_per_site(
        segment_lengths,
        selected_boundary_indexes,
        boundary_rate_generation,
        background_rate,
    )
    rec_rate_lower = 0.0
    if rec_rate_upper_generation is not None:
        rec_rate_upper = rec_rate_upper_generation
    else:
        rec_rate_upper = max(
            boundary_rate_generation * 10.0,
            annual_rate_to_generation_rate(1.0, generation_time_days),
        )
    scratch_dir = Path(tempfile.mkdtemp(prefix="espalier-grid-"))
    start_time = time.perf_counter()

    try:
        benchmark = simulate_stochastic_benchmark(
            segment_lengths,
            selected_boundary_indexes,
            sample_size=sample_size,
            population_size=population_size,
            boundary_rate=boundary_rate_generation,
            background_rate=background_rate,
            random_seed=random_seed,
            num_replicates=1,
            simulation_profile=simulation_profile,
            generation_time_days=generation_time_days,
            flu_fixture_dir=flu_fixture_dir,
            mutation_rate_per_year=mutation_rate_per_year,
        )
        dataset = benchmark.replicates[0]
        written = write_reassortment_dataset(dataset, str(scratch_dir))

        shutil.copy2(written["truth_path"], truth_path)
        _copy_tree_inputs(Path(written["trees_dir"]), input_tree_dir)

        command = _build_run_command(
            python_executable,
            Path(written["alignments_dir"]),
            Path(written["trees_dir"]),
            analysis_base,
            boundary_rate_generation,
            rec_rate_lower,
            rec_rate_upper,
            ne,
            em_iters,
            raxml_path,
            use_modern_converter,
            generation_time_days,
        )

        completed = subprocess.run(
            command,
            cwd=str(REPO_ROOT),
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
            check=True,
        )
        runtime_seconds = time.perf_counter() - start_time
        combined_log = f"{completed.stdout}\n{completed.stderr}"
        analysis_log_path.write_text(combined_log)

        analysis_dir = _select_analysis_output_dir(analysis_base)
        _remove_sequence_artifacts(analysis_dir)
        signature = _extract_pipeline_signature(analysis_dir)
        diagnostics = _load_inference_diagnostics(analysis_dir)
        summary = signature["summary"] or {}
        inferred_breakpoints = summary.get("breakpoints")
        if inferred_breakpoints:
            inferred_breakpoints = inferred_breakpoints[1:-1]
        boundary_metrics = _boundary_metrics(dataset.selected_boundary_positions, inferred_breakpoints)
        topology_boundary_metrics = {
            f"topology_changed_{key}": value
            for key, value in _boundary_metrics(
                dataset.selected_topology_changed_boundary_positions,
                inferred_breakpoints,
            ).items()
        }
        inferred_rate_native = diagnostics.get("estimated_rate")
        rate_units = diagnostics.get("rate_units")
        estimate_source = diagnostics.get("estimate_source")
        estimate_status = diagnostics.get("estimate_status")
        calibration_usable = estimate_source == "opt_mle" and estimate_status == "interior"
        inferred_boundary_rate_generation = (
            inferred_rate_native if rate_units == "per_boundary_per_generation" else None
        )
        inferred_boundary_rate_year = (
            generation_rate_to_annual_rate(inferred_boundary_rate_generation, generation_time_days)
            if inferred_boundary_rate_generation is not None
            else None
        )
        result = {
            "status": "completed",
            "success": True,
            "sample_size": sample_size,
            "boundary_rate": boundary_rate_annual,
            "boundary_rate_per_year": boundary_rate_annual,
            "boundary_rate_per_generation": boundary_rate_generation,
            "replicate_index": replicate_index,
            "random_seed": random_seed,
            "segment_lengths": list(segment_lengths),
            "selected_boundary_indexes": list(selected_boundary_indexes),
            "selected_boundary_positions": list(dataset.selected_boundary_positions),
            "observed_internal_breakpoints": list(dataset.observed_internal_breakpoints),
            "topology_changed_internal_breakpoints": list(dataset.topology_changed_internal_breakpoints),
            "selected_topology_changed_boundary_positions": list(
                dataset.selected_topology_changed_boundary_positions
            ),
            "expected_mean_rate_per_site": expected_rate_per_site,
            "expected_boundary_rate_per_generation": boundary_rate_generation,
            "expected_boundary_rate_per_year": generation_rate_to_annual_rate(
                boundary_rate_generation,
                generation_time_days,
            ),
            "inferred_rate_native": inferred_rate_native,
            "inferred_rate_units": rate_units,
            "inferred_boundary_rate_per_generation": inferred_boundary_rate_generation,
            "inferred_boundary_rate_per_year": inferred_boundary_rate_year,
            "estimate_source": estimate_source,
            "estimate_status": estimate_status,
            "rate_estimate_usable_for_calibration": calibration_usable,
            "rate_model_name": diagnostics.get("rate_model_name"),
            "ts_conversion_failures": diagnostics.get("ts_conversion_failures"),
            "runtime_seconds": runtime_seconds,
            "analysis_dir": str(analysis_dir),
            "analysis_log": str(analysis_log_path),
            "truth_path": str(truth_path),
            "input_tree_dir": str(input_tree_dir),
            "simulation_profile": simulation_profile,
            "generation_time_days": generation_time_days,
            "population_size": population_size,
            "mutation_rate_per_year": dataset.mutation_rate_per_year,
            "mutation_rate_per_generation": dataset.mutation_rate,
            "rec_rate_lower_bound": rec_rate_lower,
            "rec_rate_upper_bound": rec_rate_upper,
            "has_label_map": signature["has_label_map"],
            "has_treesequence": signature["has_treesequence"],
            "ts_num_trees": signature["ts_num_trees"],
            "ts_breakpoints": signature["ts_breakpoints"],
            "summary_breakpoints": summary.get("breakpoints"),
            "n_recomb_events": summary.get("n_recomb_events"),
            **boundary_metrics,
            **topology_boundary_metrics,
            **_rate_metrics(
                inferred_boundary_rate_generation if calibration_usable else None,
                boundary_rate_generation,
            ),
        }
    except subprocess.CalledProcessError as exc:
        runtime_seconds = time.perf_counter() - start_time
        combined_log = f"{exc.stdout or ''}\n{exc.stderr or ''}"
        analysis_log_path.write_text(combined_log)
        if analysis_base.exists():
            try:
                analysis_dir = _select_analysis_output_dir(analysis_base)
                _remove_sequence_artifacts(analysis_dir)
            except Exception:
                pass
        result = {
            "status": "failed",
            "success": False,
            "sample_size": sample_size,
            "boundary_rate": boundary_rate_annual,
            "boundary_rate_per_year": boundary_rate_annual,
            "boundary_rate_per_generation": boundary_rate_generation,
            "replicate_index": replicate_index,
            "random_seed": random_seed,
            "runtime_seconds": runtime_seconds,
            "error": f"CalledProcessError: exit {exc.returncode}",
            "analysis_log": str(analysis_log_path),
            "truth_path": str(truth_path) if truth_path.exists() else None,
            "input_tree_dir": str(input_tree_dir) if input_tree_dir.exists() else None,
            "expected_mean_rate_per_site": expected_rate_per_site,
            "expected_boundary_rate_per_generation": boundary_rate_generation,
            "expected_boundary_rate_per_year": boundary_rate_annual,
        }
    except subprocess.TimeoutExpired as exc:
        runtime_seconds = time.perf_counter() - start_time
        combined_log = f"{exc.stdout or ''}\n{exc.stderr or ''}"
        analysis_log_path.write_text(combined_log)
        result = {
            "status": "timeout",
            "success": False,
            "sample_size": sample_size,
            "boundary_rate": boundary_rate_annual,
            "boundary_rate_per_year": boundary_rate_annual,
            "boundary_rate_per_generation": boundary_rate_generation,
            "replicate_index": replicate_index,
            "random_seed": random_seed,
            "runtime_seconds": runtime_seconds,
            "error": f"Timed out after {timeout_seconds} seconds",
            "analysis_log": str(analysis_log_path),
            "truth_path": str(truth_path) if truth_path.exists() else None,
            "input_tree_dir": str(input_tree_dir) if input_tree_dir.exists() else None,
            "expected_mean_rate_per_site": expected_rate_per_site,
            "expected_boundary_rate_per_generation": boundary_rate_generation,
            "expected_boundary_rate_per_year": boundary_rate_annual,
        }
    except Exception as exc:
        runtime_seconds = time.perf_counter() - start_time
        if not analysis_log_path.exists():
            analysis_log_path.write_text(f"{type(exc).__name__}: {exc}\n")
        result = {
            "status": "failed",
            "success": False,
            "sample_size": sample_size,
            "boundary_rate": boundary_rate_annual,
            "boundary_rate_per_year": boundary_rate_annual,
            "boundary_rate_per_generation": boundary_rate_generation,
            "replicate_index": replicate_index,
            "random_seed": random_seed,
            "runtime_seconds": runtime_seconds,
            "error": f"{type(exc).__name__}: {exc}",
            "analysis_log": str(analysis_log_path),
            "truth_path": str(truth_path) if truth_path.exists() else None,
            "input_tree_dir": str(input_tree_dir) if input_tree_dir.exists() else None,
            "expected_mean_rate_per_site": expected_rate_per_site,
            "expected_boundary_rate_per_generation": boundary_rate_generation,
            "expected_boundary_rate_per_year": boundary_rate_annual,
        }
    finally:
        shutil.rmtree(scratch_dir, ignore_errors=True)

    _write_json(result_path, result)
    return result


def _safe_mean(values: Sequence[Optional[float]]) -> Optional[float]:
    observed = [value for value in values if value is not None]
    if not observed:
        return None
    return sum(observed) / len(observed)


def _safe_median(values: Sequence[Optional[float]]) -> Optional[float]:
    observed = sorted(value for value in values if value is not None)
    if not observed:
        return None
    midpoint = len(observed) // 2
    if len(observed) % 2:
        return observed[midpoint]
    return 0.5 * (observed[midpoint - 1] + observed[midpoint])


def _boundary_frequency_rate_mle(
    inferred_frequency: Optional[float],
    truth_frequency: Optional[float],
    generating_rate: float,
) -> Optional[float]:
    if (
        inferred_frequency is None
        or truth_frequency is None
        or generating_rate <= 0
        or truth_frequency <= 0.0
        or truth_frequency >= 1.0
    ):
        return None
    if inferred_frequency <= 0.0:
        return 0.0
    if inferred_frequency >= 1.0:
        return None
    opportunity = -math.log(1.0 - truth_frequency) / generating_rate
    return -math.log(1.0 - inferred_frequency) / opportunity


def _has_inferred_selected_boundary(result: Dict[str, object]) -> float:
    selected = set(int(position) for position in result.get("selected_boundary_positions", []))
    breakpoints = result.get("summary_breakpoints") or []
    inferred = set(
        int(position)
        for position in breakpoints[1:-1]
        if int(position) not in (0,)
    )
    return float(bool(selected & inferred))


def _summarize_results(results: Sequence[Dict[str, object]]) -> Dict[str, object]:
    successful = [result for result in results if result.get("success")]
    grouped: Dict[Tuple[int, float], List[Dict[str, object]]] = {}
    for result in results:
        grouped.setdefault(
            (int(result["sample_size"]), float(result["boundary_rate"])),
            [],
        ).append(result)

    treatment_summaries = []
    for (sample_size, boundary_rate), group in sorted(grouped.items()):
        successes = [result for result in group if result.get("success")]
        calibration_successes = [
            result
            for result in successes
            if result.get("rate_estimate_usable_for_calibration")
        ]
        opt_successes = [
            result
            for result in successes
            if result.get("estimate_source") == "opt_mle"
        ]
        treatment_summary = {
            "sample_size": sample_size,
            "boundary_rate": boundary_rate,
            "n_runs": len(group),
            "n_success": len(successes),
            "success_rate": len(successes) / len(group),
            "n_opt_estimates": len(opt_successes),
            "n_rate_estimates_used_for_calibration": len(calibration_successes),
            "n_fallback_estimates": sum(
                1 for result in successes if result.get("estimate_source") == "fallback_tree_length"
            ),
            "n_upper_bound_estimates": sum(
                1 for result in successes if result.get("estimate_status") == "hit_upper_bound"
            ),
            "n_lower_bound_estimates": sum(
                1 for result in successes if result.get("estimate_status") == "hit_lower_bound"
            ),
            "mean_runtime_seconds": _safe_mean([result.get("runtime_seconds") for result in group]),
            "median_runtime_seconds": _safe_median([result.get("runtime_seconds") for result in group]),
            "truth_boundary_frequency": _safe_mean(
                [
                    float(
                        bool(set(result.get("selected_boundary_positions", []))
                             & set(result.get("observed_internal_breakpoints", [])))
                    )
                    for result in successes
                ]
            ),
            "truth_topology_changed_boundary_frequency": _safe_mean(
                [
                    float(
                        bool(set(result.get("selected_boundary_positions", []))
                             & set(result.get("topology_changed_internal_breakpoints", [])))
                    )
                    for result in successes
                ]
            ),
            "inference_boundary_detection_rate": _safe_mean(
                [result.get("boundary_detected") for result in successes]
            ),
            "inferred_selected_boundary_frequency": _safe_mean(
                [_has_inferred_selected_boundary(result) for result in successes]
            ),
            "mean_boundary_precision": _safe_mean(
                [result.get("boundary_precision") for result in successes]
            ),
            "mean_boundary_recall": _safe_mean(
                [result.get("boundary_recall") for result in successes]
            ),
            "mean_event_count_error": _safe_mean(
                [result.get("event_count_error") for result in successes]
            ),
            "topology_changed_inference_boundary_detection_rate": _safe_mean(
                [result.get("topology_changed_boundary_detected") for result in successes]
            ),
            "mean_topology_changed_boundary_precision": _safe_mean(
                [result.get("topology_changed_boundary_precision") for result in successes]
            ),
            "mean_topology_changed_boundary_recall": _safe_mean(
                [result.get("topology_changed_boundary_recall") for result in successes]
            ),
            "mean_topology_changed_event_count_error": _safe_mean(
                [result.get("topology_changed_event_count_error") for result in successes]
            ),
            "mean_expected_uniform_rate_per_site": _safe_mean(
                [result.get("expected_mean_rate_per_site") for result in successes]
            ),
            "mean_expected_boundary_rate_per_generation": _safe_mean(
                [result.get("expected_boundary_rate_per_generation") for result in successes]
            ),
            "mean_expected_boundary_rate_per_year": _safe_mean(
                [result.get("expected_boundary_rate_per_year") for result in successes]
            ),
            "mean_inferred_boundary_rate_per_generation": _safe_mean(
                [result.get("inferred_boundary_rate_per_generation") for result in calibration_successes]
            ),
            "mean_inferred_boundary_rate_per_year": _safe_mean(
                [result.get("inferred_boundary_rate_per_year") for result in calibration_successes]
            ),
            "mean_rate_ratio_vs_expected_native": _safe_mean(
                [result.get("rate_ratio_vs_expected_native") for result in calibration_successes]
            ),
            "mean_rate_log_error_vs_expected_native": _safe_mean(
                [result.get("rate_log_error_vs_expected_native") for result in calibration_successes]
            ),
        }
        frequency_mle_rate = _boundary_frequency_rate_mle(
            treatment_summary["inferred_selected_boundary_frequency"],
            treatment_summary["truth_topology_changed_boundary_frequency"],
            treatment_summary["mean_expected_boundary_rate_per_generation"] or 0.0,
        )
        treatment_summary["topology_changed_frequency_mle_rate_per_year"] = (
            generation_rate_to_annual_rate(
                frequency_mle_rate,
                results[0].get("generation_time_days", DEFAULT_GENERATION_TIME_DAYS),
            )
            if frequency_mle_rate is not None
            else None
        )
        treatment_summaries.append(treatment_summary)

    return {
        "n_runs": len(results),
        "n_success": len(successful),
        "success_rate": len(successful) / len(results) if results else 0.0,
        "treatments": treatment_summaries,
    }


def _write_results_table(path: Path, results: Sequence[Dict[str, object]]) -> None:
    rows = []
    for result in results:
        row = dict(result)
        for key, value in list(row.items()):
            if isinstance(value, (list, dict)):
                row[key] = json.dumps(value, sort_keys=True)
        rows.append(row)

    fieldnames = sorted({key for row in rows for key in row.keys()})
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _build_tasks(
    tree_sizes: Sequence[int],
    boundary_rates: Sequence[float],
    replicates: int,
) -> List[Tuple[int, float, int]]:
    tasks = []
    for sample_size in tree_sizes:
        for boundary_rate in boundary_rates:
            for replicate_index in range(1, replicates + 1):
                tasks.append((sample_size, boundary_rate, replicate_index))
    return tasks


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a large simulated reassortment benchmark grid with temporary alignments only."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "TestFiles" / "SimulationBenchmarks" / "flu_like_rate_grid",
        help="Directory for benchmark outputs",
    )
    parser.add_argument(
        "--segment-lengths",
        type=str,
        default=",".join(str(length) for length in DEFAULT_SEGMENT_LENGTHS),
        help="Comma-separated segment lengths for the concatenated genome",
    )
    parser.add_argument(
        "--reassorting-boundaries",
        type=str,
        default="0",
        help="Boundary indexes to activate ('all', 'none', or comma-separated indexes)",
    )
    parser.add_argument(
        "--tree-sizes",
        type=str,
        default=",".join(str(size) for size in DEFAULT_TREE_SIZES),
        help="Comma-separated sample sizes",
    )
    parser.add_argument(
        "--boundary-rates",
        type=str,
        default=",".join(str(rate) for rate in DEFAULT_BOUNDARY_RATES_PER_YEAR),
        help="Comma-separated boundary reassortment intensities per year",
    )
    parser.add_argument(
        "--replicates",
        type=int,
        default=DEFAULT_REPLICATES,
        help="Replicates per treatment combination",
    )
    parser.add_argument(
        "--population-size",
        type=float,
        default=None,
        help="Population size passed to msprime and SCAR Ne (defaults to the flu-derived proxy)",
    )
    parser.add_argument(
        "--mutation-rate-per-year",
        type=float,
        default=None,
        help="Per-site mutation rate per year (defaults to the flu-derived clock proxy)",
    )
    parser.add_argument(
        "--background-rate",
        type=float,
        default=DEFAULT_BACKGROUND_RATE,
        help="Background recombination rate outside selected boundaries",
    )
    parser.add_argument(
        "--em-iters",
        type=int,
        default=DEFAULT_EM_ITERS,
        help="EM iterations per inference run",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=DEFAULT_JOBS,
        help="Concurrent treatment jobs to run",
    )
    parser.add_argument(
        "--python-executable",
        type=str,
        default=sys.executable,
        help="Python interpreter used for the analysis subprocess",
    )
    parser.add_argument(
        "--raxml-path",
        type=str,
        default=shutil.which("raxml-ng") or "raxml-ng",
        help="Path to the raxml-ng executable",
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=DEFAULT_RANDOM_SEED,
        help="Base random seed for reproducible replicate seeds",
    )
    parser.add_argument(
        "--generation-time-days",
        type=float,
        default=DEFAULT_GENERATION_TIME_DAYS,
        help=f"Generation time in days for annual-to-generation rate conversion (default: {DEFAULT_GENERATION_TIME_DAYS})",
    )
    parser.add_argument(
        "--flu-fixture-dir",
        type=str,
        default=DEFAULT_FLU_FIXTURE_DIR,
        help="Flu fixture directory used to derive empirical calibration defaults",
    )
    parser.add_argument(
        "--timeout-seconds",
        type=int,
        default=None,
        help="Optional per-replicate timeout for the analysis subprocess",
    )
    parser.add_argument(
        "--rec-rate-upper-per-year",
        type=float,
        default=None,
        help="Optional fixed upper bound for reassortment rate estimates, in per-boundary per-year units",
    )
    parser.add_argument(
        "--no-modern-converter",
        action="store_true",
        help="Disable the modern TreeSequence converter for EM output",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rerun replicates even if result.json already exists",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    segment_lengths = _parse_int_list(args.segment_lengths)
    selected_boundary_indexes = _parse_boundaries(args.reassorting_boundaries, segment_lengths)
    tree_sizes = _parse_int_list(args.tree_sizes)
    boundary_rates_per_year = _parse_float_list(args.boundary_rates)
    calibration = load_flu_subsampled_calibration(
        fixture_dir=args.flu_fixture_dir,
        generation_time_days=args.generation_time_days,
    )
    population_size = (
        args.population_size
        if args.population_size is not None
        else calibration.population_size_proxy
    )
    mutation_rate_per_year = (
        args.mutation_rate_per_year
        if args.mutation_rate_per_year is not None
        else calibration.mutation_rate_per_year
    )
    boundary_rates_generation = [
        annual_rate_to_generation_rate(rate, args.generation_time_days)
        for rate in boundary_rates_per_year
    ]
    rec_rate_upper_generation = (
        annual_rate_to_generation_rate(args.rec_rate_upper_per_year, args.generation_time_days)
        if args.rec_rate_upper_per_year is not None
        else None
    )
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    benchmark_config = {
        "simulation_profile": SIMULATION_PROFILE_FLU_CALIBRATED,
        "segment_lengths": segment_lengths,
        "selected_boundary_indexes": selected_boundary_indexes,
        "selected_boundary_positions": [
            _boundary_positions(segment_lengths)[index]
            for index in selected_boundary_indexes
        ],
        "tree_sizes": tree_sizes,
        "boundary_rates_per_year": boundary_rates_per_year,
        "boundary_rates_per_generation": boundary_rates_generation,
        "replicates": args.replicates,
        "population_size": population_size,
        "mutation_rate_per_year": mutation_rate_per_year,
        "background_rate": args.background_rate,
        "em_iters": args.em_iters,
        "jobs": args.jobs,
        "python_executable": args.python_executable,
        "raxml_path": args.raxml_path,
        "random_seed": args.random_seed,
        "generation_time_days": args.generation_time_days,
        "flu_calibration": calibration.summary(),
        "timeout_seconds": args.timeout_seconds,
        "rec_rate_upper_per_year": args.rec_rate_upper_per_year,
        "rec_rate_upper_per_generation": rec_rate_upper_generation,
        "store_sequences": False,
    }
    _write_json(output_dir / "benchmark_config.json", benchmark_config)

    tasks = _build_tasks(tree_sizes, boundary_rates_per_year, args.replicates)
    results: List[Dict[str, object]] = []

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as executor:
        future_map = {}
        boundary_rate_generation_lookup = dict(zip(boundary_rates_per_year, boundary_rates_generation))
        for sample_size, boundary_rate_per_year, replicate_index in tasks:
            seed = _combination_seed(args.random_seed, sample_size, boundary_rate_per_year, replicate_index)
            future = executor.submit(
                _run_single_replicate,
                output_dir=output_dir,
                sample_size=sample_size,
                boundary_rate_annual=boundary_rate_per_year,
                boundary_rate_generation=boundary_rate_generation_lookup[boundary_rate_per_year],
                replicate_index=replicate_index,
                segment_lengths=segment_lengths,
                selected_boundary_indexes=selected_boundary_indexes,
                population_size=population_size,
                mutation_rate_per_year=mutation_rate_per_year,
                background_rate=args.background_rate,
                em_iters=args.em_iters,
                ne=population_size,
                raxml_path=args.raxml_path,
                python_executable=args.python_executable,
                random_seed=seed,
                use_modern_converter=not args.no_modern_converter,
                timeout_seconds=args.timeout_seconds,
                rec_rate_upper_generation=rec_rate_upper_generation,
                force=args.force,
                simulation_profile=SIMULATION_PROFILE_FLU_CALIBRATED,
                generation_time_days=args.generation_time_days,
                flu_fixture_dir=args.flu_fixture_dir,
            )
            future_map[future] = (sample_size, boundary_rate_per_year, replicate_index)

        completed_count = 0
        for future in as_completed(future_map):
            sample_size, boundary_rate, replicate_index = future_map[future]
            result = future.result()
            results.append(result)
            completed_count += 1
            _write_results_table(output_dir / "benchmark_results.csv", results)
            _write_json(output_dir / "benchmark_summary.json", _summarize_results(results))
            status = result.get("status")
            success = result.get("success")
            print(
                f"[{completed_count}/{len(tasks)}] "
                f"n={sample_size} rate={boundary_rate:g} rep={replicate_index:03d} "
                f"status={status} success={success}",
                flush=True,
            )

    if calibrate_results_file is not None:
        calibrate_results_file(
            output_dir / "benchmark_results.csv",
            output_dir / "rate_calibration.json",
            output_dir / "rate_calibration_summary.csv",
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
