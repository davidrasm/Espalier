#!/usr/bin/env python3

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Espalier.ReassortmentSim import (
    EXACT_FIXTURE_MODE,
    STOCHASTIC_BENCHMARK_MODE,
    SIMULATION_PROFILE_FLU_CALIBRATED,
    SIMULATION_PROFILE_SIMPLE,
    DEFAULT_BACKGROUND_RATE,
    DEFAULT_BOUNDARY_RATE,
    DEFAULT_FLU_FIXTURE_DIR,
    DEFAULT_GENERATION_TIME_DAYS,
    DEFAULT_MUTATION_RATE,
    DEFAULT_SAMPLE_SIZE,
    annual_rate_to_generation_rate,
    load_flu_subsampled_calibration,
    simulate_exact_fixture,
    simulate_stochastic_benchmark,
    write_reassortment_dataset,
    write_stochastic_benchmark,
)


def _parse_int_list(raw_value):
    return [int(item.strip()) for item in raw_value.split(",") if item.strip()]


def _parse_boundary_indexes(raw_value, segment_lengths):
    value = raw_value.strip().lower()
    if value == "all":
        return list(range(len(segment_lengths) - 1))
    if value == "none":
        return []
    return _parse_int_list(raw_value)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate reassortment-ready fixtures and benchmarks with msprime"
    )
    parser.add_argument(
        "--mode",
        choices=[EXACT_FIXTURE_MODE, STOCHASTIC_BENCHMARK_MODE],
        default=EXACT_FIXTURE_MODE,
        help="Simulation workflow to run",
    )
    parser.add_argument(
        "--simulation-profile",
        choices=[SIMULATION_PROFILE_SIMPLE, SIMULATION_PROFILE_FLU_CALIBRATED],
        default=SIMULATION_PROFILE_SIMPLE,
        help="Choose the legacy simple simulator or the flu-calibrated natural-process simulator",
    )
    parser.add_argument(
        "--segment-lengths",
        required=True,
        help="Comma-delimited segment lengths, e.g. 1719,1464",
    )
    parser.add_argument(
        "--segment-names",
        default=None,
        help="Optional comma-delimited segment names",
    )
    parser.add_argument(
        "--reassorting-boundaries",
        default="all",
        help="Boundary indexes to activate: 'all', 'none', or a comma-delimited 0-based list",
    )
    parser.add_argument(
        "--sample-size",
        type=int,
        default=DEFAULT_SAMPLE_SIZE,
        help=f"Number of haploid samples (default: {DEFAULT_SAMPLE_SIZE})",
    )
    parser.add_argument(
        "--population-size",
        type=float,
        default=None,
        help="Effective population size in generations (defaults to the flu-derived proxy in flu_calibrated mode)",
    )
    parser.add_argument(
        "--mutation-rate",
        type=float,
        default=DEFAULT_MUTATION_RATE,
        help=f"Per-site mutation rate per generation for the simple profile (default: {DEFAULT_MUTATION_RATE})",
    )
    parser.add_argument(
        "--mutation-rate-per-year",
        type=float,
        default=None,
        help="Per-site mutation rate per year for the flu_calibrated profile (defaults to the flu-derived clock proxy)",
    )
    parser.add_argument(
        "--boundary-rate",
        type=float,
        default=DEFAULT_BOUNDARY_RATE,
        help=f"Boundary hotspot rate per generation (default: {DEFAULT_BOUNDARY_RATE})",
    )
    parser.add_argument(
        "--boundary-rate-per-year",
        type=float,
        default=None,
        help="Boundary hotspot rate per year for flu_calibrated simulations",
    )
    parser.add_argument(
        "--background-rate",
        type=float,
        default=DEFAULT_BACKGROUND_RATE,
        help=f"Background per-site recombination rate (default: {DEFAULT_BACKGROUND_RATE})",
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=13,
        help="Base random seed for ancestry simulation",
    )
    parser.add_argument(
        "--num-replicates",
        type=int,
        default=5,
        help="Replicate count for stochastic benchmarks",
    )
    parser.add_argument(
        "--max-attempts",
        type=int,
        default=100,
        help="Maximum accept/reject attempts for exact fixtures",
    )
    parser.add_argument(
        "--require-topology-change",
        action="store_true",
        help="Require requested exact-fixture breakpoints to change adjacent local-tree topology",
    )
    parser.add_argument(
        "--generation-time-days",
        type=float,
        default=DEFAULT_GENERATION_TIME_DAYS,
        help=f"Generation time in days for annual-to-generation rate conversion (default: {DEFAULT_GENERATION_TIME_DAYS})",
    )
    parser.add_argument(
        "--mutation-kappa",
        type=float,
        default=3.0,
        help="HKY kappa parameter used in flu_calibrated mode",
    )
    parser.add_argument(
        "--gamma-shape",
        type=float,
        default=0.5,
        help="Gamma shape for among-site rate heterogeneity in flu_calibrated mode",
    )
    parser.add_argument(
        "--flu-fixture-dir",
        default=DEFAULT_FLU_FIXTURE_DIR,
        help="Flu fixture directory used to derive empirical calibration defaults",
    )
    parser.add_argument(
        "--start-date",
        default="2020-01-01",
        help="ISO date used for the first sample label",
    )
    parser.add_argument(
        "--output-dir",
        "-o",
        required=True,
        help="Output directory for alignments, trees, and truth manifests",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    segment_lengths = _parse_int_list(args.segment_lengths)
    segment_names = None
    if args.segment_names:
        segment_names = [item.strip() for item in args.segment_names.split(",") if item.strip()]
    boundary_indexes = _parse_boundary_indexes(args.reassorting_boundaries, segment_lengths)
    population_size = args.population_size
    mutation_rate_per_year = args.mutation_rate_per_year
    boundary_rate = args.boundary_rate

    if args.simulation_profile == SIMULATION_PROFILE_FLU_CALIBRATED:
        calibration = load_flu_subsampled_calibration(
            fixture_dir=args.flu_fixture_dir,
            generation_time_days=args.generation_time_days,
        )
        if population_size is None:
            population_size = calibration.population_size_proxy
        if mutation_rate_per_year is None:
            mutation_rate_per_year = calibration.mutation_rate_per_year
        if args.boundary_rate_per_year is not None:
            boundary_rate = annual_rate_to_generation_rate(
                args.boundary_rate_per_year,
                args.generation_time_days,
            )
    elif population_size is None:
        population_size = 10_000

    if args.mode == EXACT_FIXTURE_MODE:
        dataset = simulate_exact_fixture(
            segment_lengths,
            boundary_indexes,
            segment_names=segment_names,
            sample_size=args.sample_size,
            population_size=population_size,
            mutation_rate=args.mutation_rate,
            boundary_rate=boundary_rate,
            background_rate=args.background_rate,
            random_seed=args.random_seed,
            max_attempts=args.max_attempts,
            start_date=args.start_date,
            simulation_profile=args.simulation_profile,
            generation_time_days=args.generation_time_days,
            flu_fixture_dir=args.flu_fixture_dir,
            mutation_rate_per_year=mutation_rate_per_year,
            mutation_kappa=args.mutation_kappa,
            gamma_shape=args.gamma_shape,
            require_topology_change=args.require_topology_change,
        )
        result = write_reassortment_dataset(dataset, args.output_dir)
        print(f"Wrote exact fixture to {result['output_dir']}")
        print(f"Truth manifest: {result['truth_path']}")
    else:
        benchmark = simulate_stochastic_benchmark(
            segment_lengths,
            boundary_indexes,
            segment_names=segment_names,
            sample_size=args.sample_size,
            population_size=population_size,
            mutation_rate=args.mutation_rate,
            boundary_rate=boundary_rate,
            background_rate=args.background_rate,
            random_seed=args.random_seed,
            num_replicates=args.num_replicates,
            start_date=args.start_date,
            simulation_profile=args.simulation_profile,
            generation_time_days=args.generation_time_days,
            flu_fixture_dir=args.flu_fixture_dir,
            mutation_rate_per_year=mutation_rate_per_year,
            mutation_kappa=args.mutation_kappa,
            gamma_shape=args.gamma_shape,
        )
        result = write_stochastic_benchmark(benchmark, args.output_dir)
        print(f"Wrote stochastic benchmark to {result['output_dir']}")
        print(f"Summary: {result['summary_path']}")


if __name__ == "__main__":
    main()
