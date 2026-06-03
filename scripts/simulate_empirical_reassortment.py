#!/usr/bin/env python3

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Espalier.EmpiricalReassortmentSim import (
    DEFAULT_KAPPA,
    DEFAULT_MUTATION_RATE,
    DEFAULT_STATE_FREQS,
    generate_empirical_reassortment_panel,
    maybe_parse_segment_names,
    parse_csv_rows,
    parse_float_list,
    parse_int_list,
)


def _parse_reassorted_segments(raw_value):
    normalized = raw_value.strip().lower()
    if normalized in {"last", ""}:
        return None
    return parse_int_list(raw_value)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate empirical reassortment simulation datasets from a fixed input tree"
    )
    parser.add_argument("--tree", help="Input rooted Newick tree")
    parser.add_argument("--csv", default="", help="Optional CSV file describing multiple runs")
    parser.add_argument("--output", required=True, help="Output directory for generated samples")
    parser.add_argument("--segment-lengths", default="1719,1464", help="Comma-delimited segment lengths")
    parser.add_argument("--segment-names", default="", help="Optional comma-delimited segment names")
    parser.add_argument(
        "--reassorted-segments",
        default="last",
        help="Comma-delimited reassorted segment indexes or 'last' (default)",
    )
    parser.add_argument("--num-reassortments", type=int, default=1, help="Number of subtree moves per dataset")
    parser.add_argument("--repeats", type=int, default=1, help="Number of replicate datasets to generate")
    parser.add_argument("--start", type=int, default=1, help="Starting replicate index")
    parser.add_argument("--min-size", type=int, default=1, help="Minimum moved clade size")
    parser.add_argument("--max-size", type=int, default=10, help="Maximum moved clade size")
    parser.add_argument("--min-distance", type=float, default=0.0, help="Minimum normalized reassortment distance")
    parser.add_argument("--max-distance", type=float, default=1.0, help="Maximum normalized reassortment distance")
    parser.add_argument("--mutation-rate", type=float, default=DEFAULT_MUTATION_RATE, help="Branch scaling factor for sequence simulation")
    parser.add_argument("--kappa", type=float, default=DEFAULT_KAPPA, help="HKY kappa parameter for pyvolve")
    parser.add_argument(
        "--state-frequencies",
        default="0.25,0.25,0.25,0.25",
        help="Comma-delimited equilibrium nucleotide frequencies",
    )
    parser.add_argument("--random-seed", type=int, default=13, help="Base random seed")
    return parser.parse_args()


def _run_configuration(row, args):
    tree_path = row["tree"] if row else args.tree
    if not tree_path:
        raise ValueError("Either --tree or --csv with a tree column is required")

    segment_lengths = parse_int_list(row["segment_lengths"]) if row and row.get("segment_lengths") else parse_int_list(args.segment_lengths)
    segment_names = maybe_parse_segment_names(row.get("segment_names")) if row else maybe_parse_segment_names(args.segment_names)
    if segment_names is None:
        segment_names = maybe_parse_segment_names(args.segment_names)
    reassorted_segments = _parse_reassorted_segments(row["reassorted_segments"]) if row and row.get("reassorted_segments") else _parse_reassorted_segments(args.reassorted_segments)
    state_frequencies = parse_float_list(row["state_frequencies"]) if row and row.get("state_frequencies") else parse_float_list(args.state_frequencies)

    repeats = int(row["repeats"]) if row and row.get("repeats") else args.repeats
    start = int(row["start"]) if row and row.get("start") else args.start
    num_reassortments = int(row["num_reassortments"]) if row and row.get("num_reassortments") else args.num_reassortments
    min_size = int(row["min_size"]) if row and row.get("min_size") else args.min_size
    max_size = int(row["max_size"]) if row and row.get("max_size") else args.max_size
    min_distance = float(row["min_distance"]) if row and row.get("min_distance") else args.min_distance
    max_distance = float(row["max_distance"]) if row and row.get("max_distance") else args.max_distance
    mutation_rate = float(row["mutation_rate"]) if row and row.get("mutation_rate") else args.mutation_rate
    kappa = float(row["kappa"]) if row and row.get("kappa") else args.kappa
    output_dir = row["output"] if row and row.get("output") else args.output
    random_seed = int(row["random_seed"]) if row and row.get("random_seed") else args.random_seed

    written = generate_empirical_reassortment_panel(
        tree_path=tree_path,
        output_dir=output_dir,
        repeats=repeats,
        start_index=start,
        segment_lengths=segment_lengths,
        segment_names=segment_names,
        reassorted_segment_indexes=reassorted_segments,
        num_reassortments=num_reassortments,
        min_clade_size=min_size,
        max_clade_size=max_size,
        min_normalized_distance=min_distance,
        max_normalized_distance=max_distance,
        mutation_rate=mutation_rate,
        kappa=kappa,
        state_frequencies=state_frequencies,
        random_seed=random_seed,
    )
    for record in written:
        print(f"Wrote empirical reassortment dataset to {record['output_dir']}")


def main():
    args = parse_args()
    if args.csv:
        for row in parse_csv_rows(args.csv):
            _run_configuration(row, args)
    else:
        _run_configuration(None, args)


if __name__ == "__main__":
    main()
