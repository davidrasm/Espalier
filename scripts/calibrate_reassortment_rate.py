#!/usr/bin/env python3
"""
Calibrate reassortment-rate estimates from benchmark result panels.

The calibration targets inferred boundary frequency across replicate panels rather
than ARG recombination-node counts from individual reconstructions.
"""

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


DEFAULT_GENERATION_TIME_DAYS = 3.0


def generation_rate_to_annual_rate(rate_per_generation: float, generation_time_days: float) -> float:
    return float(rate_per_generation) * (365.25 / float(generation_time_days))


def _decode_value(value):
    if value is None or value == "" or value == "None":
        return None
    if value == "True":
        return True
    if value == "False":
        return False
    try:
        return json.loads(value)
    except Exception:
        pass
    try:
        return float(value)
    except Exception:
        return value


def _load_rows(results_csv: Path) -> List[Dict[str, object]]:
    with open(results_csv, newline="") as handle:
        rows = list(csv.DictReader(handle))
    return [
        {key: _decode_value(value) for key, value in row.items()}
        for row in rows
    ]


def _as_int_set(values) -> set:
    if not values:
        return set()
    return {int(value) for value in values}


def _inferred_internal_breakpoints(row: Dict[str, object]) -> set:
    breakpoints = row.get("summary_breakpoints") or []
    return {
        int(position)
        for position in breakpoints[1:-1]
        if int(position) != 0
    }


def _safe_mean(values: Sequence[Optional[float]]) -> Optional[float]:
    observed = [float(value) for value in values if value is not None]
    if not observed:
        return None
    return sum(observed) / len(observed)


def _safe_fraction(values: Sequence[bool]) -> Optional[float]:
    if not values:
        return None
    return sum(float(value) for value in values) / len(values)


def _boundary_frequency_rate_mle(
    inferred_frequency: Optional[float],
    truth_frequency: Optional[float],
    generating_rate: float,
    generation_time_days: float,
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
    return generation_rate_to_annual_rate(
        -math.log(1.0 - inferred_frequency) / opportunity,
        generation_time_days,
    )


def _fit_exponential_detection_curve(points: Sequence[Dict[str, object]]) -> Dict[str, object]:
    usable = [
        point
        for point in points
        if point["true_rate_per_year"] > 0
        and point["inferred_selected_boundary_frequency"] is not None
    ]
    if len(usable) < 2:
        return {
            "model": "one_minus_exp",
            "status": "insufficient_points",
            "beta_per_year": None,
            "rmse": None,
            "points": [],
        }
    observed_frequencies = [
        point["inferred_selected_boundary_frequency"]
        for point in usable
    ]
    if all(frequency <= 0.0 for frequency in observed_frequencies):
        return {
            "model": "one_minus_exp",
            "status": "no_detected_boundaries",
            "beta_per_year": None,
            "rmse": None,
            "points": usable,
        }
    if all(frequency >= 1.0 for frequency in observed_frequencies):
        return {
            "model": "one_minus_exp",
            "status": "saturated_detected_boundaries",
            "beta_per_year": None,
            "rmse": None,
            "points": usable,
        }

    candidate_betas = [
        10 ** (-6.0 + idx * (10.0 / 2000.0))
        for idx in range(2001)
    ]
    best_beta = None
    best_sse = math.inf
    for beta in candidate_betas:
        sse = 0.0
        for point in usable:
            expected = 1.0 - math.exp(-beta * point["true_rate_per_year"])
            observed = point["inferred_selected_boundary_frequency"]
            sse += (observed - expected) ** 2
        if sse < best_sse:
            best_sse = sse
            best_beta = beta

    fitted_points = []
    for point in usable:
        predicted_frequency = 1.0 - math.exp(-best_beta * point["true_rate_per_year"])
        fitted_points.append(
            {
                "true_rate_per_year": point["true_rate_per_year"],
                "observed_inferred_frequency": point["inferred_selected_boundary_frequency"],
                "predicted_inferred_frequency": predicted_frequency,
                "residual": point["inferred_selected_boundary_frequency"] - predicted_frequency,
            }
        )

    return {
        "model": "one_minus_exp",
        "status": "fit",
        "beta_per_year": best_beta,
        "rmse": math.sqrt(best_sse / len(usable)),
        "points": fitted_points,
    }


def estimate_rate_from_frequency(
    inferred_frequency: float,
    beta_per_year: float,
) -> Optional[float]:
    if beta_per_year <= 0:
        return None
    if inferred_frequency <= 0:
        return 0.0
    if inferred_frequency >= 1:
        return None
    return -math.log(1.0 - inferred_frequency) / beta_per_year


def _aggregate_rows(rows: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    grouped: Dict[Tuple[int, float], List[Dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(int(row["sample_size"]), float(row["boundary_rate_per_year"]))].append(row)

    treatments = []
    for (sample_size, true_rate), group in sorted(grouped.items()):
        successes = [row for row in group if row.get("success")]
        if not successes:
            treatments.append(
                {
                    "sample_size": sample_size,
                    "true_rate_per_year": true_rate,
                    "n_runs": len(group),
                    "n_success": 0,
                    "success_rate": 0.0,
                }
            )
            continue

        truth_topology_flags = []
        truth_raw_flags = []
        inferred_flags = []
        for row in successes:
            selected = _as_int_set(row.get("selected_boundary_positions"))
            truth_topology_flags.append(
                bool(selected & _as_int_set(row.get("topology_changed_internal_breakpoints")))
            )
            truth_raw_flags.append(
                bool(selected & _as_int_set(row.get("observed_internal_breakpoints")))
            )
            inferred_flags.append(bool(selected & _inferred_internal_breakpoints(row)))

        generation_time_days = _safe_mean(
            [row.get("generation_time_days") for row in successes]
        ) or DEFAULT_GENERATION_TIME_DAYS
        expected_generation_rate = _safe_mean(
            [row.get("expected_boundary_rate_per_generation") for row in successes]
        )
        truth_topology_frequency = _safe_fraction(truth_topology_flags)
        inferred_frequency = _safe_fraction(inferred_flags)
        frequency_mle_rate = (
            _boundary_frequency_rate_mle(
                inferred_frequency,
                truth_topology_frequency,
                expected_generation_rate,
                generation_time_days,
            )
            if expected_generation_rate is not None
            else None
        )
        scar_rates = [
            row.get("inferred_boundary_rate_per_year")
            for row in successes
            if row.get("inferred_boundary_rate_per_year") is not None
        ]
        treatments.append(
            {
                "sample_size": sample_size,
                "true_rate_per_year": true_rate,
                "n_runs": len(group),
                "n_success": len(successes),
                "success_rate": len(successes) / len(group),
                "truth_boundary_frequency": _safe_fraction(truth_raw_flags),
                "truth_topology_changed_boundary_frequency": truth_topology_frequency,
                "inferred_selected_boundary_frequency": inferred_frequency,
                "frequency_mle_rate_per_year": frequency_mle_rate,
                "mean_scar_rate_per_year": _safe_mean(scar_rates),
                "n_scar_upper_bound": sum(
                    1 for row in successes if row.get("estimate_status") == "hit_upper_bound"
                ),
                "n_scar_lower_bound": sum(
                    1 for row in successes if row.get("estimate_status") == "hit_lower_bound"
                ),
                "n_scar_interior": sum(
                    1 for row in successes if row.get("estimate_status") == "interior"
                ),
                "mean_n_recomb_events": _safe_mean(
                    [row.get("n_recomb_events") for row in successes]
                ),
                "mean_runtime_seconds": _safe_mean(
                    [row.get("runtime_seconds") for row in successes]
                ),
            }
        )
    return treatments


def calibrate_results_file(
    results_csv: Path,
    output_json: Path,
    output_csv: Optional[Path] = None,
) -> Dict[str, object]:
    rows = _load_rows(results_csv)
    treatments = _aggregate_rows(rows)

    by_sample_size: Dict[int, List[Dict[str, object]]] = defaultdict(list)
    for treatment in treatments:
        by_sample_size[int(treatment["sample_size"])].append(treatment)

    curves = {}
    for sample_size, points in sorted(by_sample_size.items()):
        fit_points = [
            {
                "true_rate_per_year": point["true_rate_per_year"],
                "inferred_selected_boundary_frequency": point.get(
                    "inferred_selected_boundary_frequency"
                ),
            }
            for point in points
        ]
        curve = _fit_exponential_detection_curve(fit_points)
        if curve.get("beta_per_year"):
            for point in points:
                point["curve_calibrated_rate_per_year"] = estimate_rate_from_frequency(
                    point.get("inferred_selected_boundary_frequency") or 0.0,
                    curve["beta_per_year"],
                )
        curves[str(sample_size)] = curve

    payload = {
        "results_csv": str(results_csv),
        "n_rows": len(rows),
        "treatments": treatments,
        "curves_by_sample_size": curves,
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    with open(output_json, "w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")

    if output_csv is not None:
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        fieldnames = sorted({key for row in treatments for key in row.keys()})
        with open(output_csv, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(treatments)
    return payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fit frequency-based reassortment-rate calibration from benchmark results."
    )
    parser.add_argument("results_csv", type=Path, help="Path to benchmark_results.csv")
    parser.add_argument(
        "--output-json",
        type=Path,
        default=None,
        help="Output JSON path (default: rate_calibration.json next to the CSV)",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=None,
        help="Output treatment summary CSV path (default: rate_calibration_summary.csv next to the CSV)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output_json = args.output_json or args.results_csv.with_name("rate_calibration.json")
    output_csv = args.output_csv or args.results_csv.with_name("rate_calibration_summary.csv")
    calibrate_results_file(args.results_csv, output_json, output_csv)
    print(f"Wrote calibration: {output_json}")
    print(f"Wrote summary: {output_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
