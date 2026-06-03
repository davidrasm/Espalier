#!/usr/bin/env python3

import argparse
import json
import os
import sys
from pathlib import Path

import dendropy

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Espalier.EmpiricalReassortmentSim import (
    generate_empirical_reassortment_panel,
    maybe_parse_segment_names,
    parse_float_list,
    parse_int_list,
    run_analysis_pipeline,
)
from Espalier import MAF


def _parse_reassorted_segments(raw_value):
    normalized = raw_value.strip().lower()
    if normalized in {"last", ""}:
        return None
    return parse_int_list(raw_value)


def _parse_breakpoints(summary_path: Path):
    if not summary_path.exists():
        return []
    for line in summary_path.read_text().splitlines():
        if line.startswith("Breakpoints:"):
            return json.loads(line.split(": ", 1)[1].replace("'", "\""))
    return []


def _load_tree(tree_path: Path, taxon_namespace: dendropy.TaxonNamespace) -> dendropy.Tree:
    return dendropy.Tree.get(
        path=str(tree_path),
        schema="newick",
        rooting="default-rooted",
        preserve_underscores=True,
        taxon_namespace=taxon_namespace,
    )


def _compute_spr_distance(tree_a_path: Path, tree_b_path: Path):
    if not tree_a_path.exists() or not tree_b_path.exists():
        return None
    taxon_namespace = dendropy.TaxonNamespace()
    tree_a = _load_tree(tree_a_path, taxon_namespace)
    tree_b = _load_tree(tree_b_path, taxon_namespace)
    return int(MAF.get_spr_dist(tree_a, tree_b))


def _get_reference_vs_reassorted_segment_names(truth):
    reassorted_indexes = set(truth["reassorted_segment_indexes"])
    reference_name = None
    reassorted_name = None
    for idx, name in enumerate(truth["segment_names"]):
        if idx in reassorted_indexes and reassorted_name is None:
            reassorted_name = name
        if idx not in reassorted_indexes and reference_name is None:
            reference_name = name
    if reference_name is None or reassorted_name is None:
        raise ValueError("Could not determine reference/reassorted segment pair from truth manifest")
    return reference_name, reassorted_name


def _compute_move_count_metrics(sample_dir: Path, analysis_dir: Path, truth):
    reference_name, reassorted_name = _get_reference_vs_reassorted_segment_names(truth)
    expected_spr = int(truth["num_reassortments"])
    true_tree_spr = _compute_spr_distance(
        sample_dir / "trees" / f"{reference_name}.fasta.rooted.tre",
        sample_dir / "trees" / f"{reassorted_name}.fasta.rooted.tre",
    )
    arg_local_spr = _compute_spr_distance(
        analysis_dir / f"{reference_name}_ARG_local.tre",
        analysis_dir / f"{reassorted_name}_ARG_local.tre",
    )
    arg_with_recomb_spr = _compute_spr_distance(
        analysis_dir / f"{reference_name}_ARG_with_recomb.tre",
        analysis_dir / f"{reassorted_name}_ARG_with_recomb.tre",
    )
    return {
        "expected_spr_count": expected_spr,
        "true_tree_spr_count": true_tree_spr,
        "arg_local_spr_count": arg_local_spr,
        "arg_with_recomb_spr_count": arg_with_recomb_spr,
        "true_tree_matches_expected": true_tree_spr == expected_spr,
        "arg_local_matches_expected": arg_local_spr == expected_spr,
        "arg_with_recomb_matches_expected": arg_with_recomb_spr == expected_spr,
    }


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate empirical reassortment replicates and optionally run the Espalier reassortment pipeline"
    )
    parser.add_argument("--tree", required=True, help="Input rooted Newick tree")
    parser.add_argument("--output", required=True, help="Benchmark output directory")
    parser.add_argument("--segment-lengths", default="1719,1464", help="Comma-delimited segment lengths")
    parser.add_argument("--segment-names", default="", help="Optional comma-delimited segment names")
    parser.add_argument(
        "--reassorted-segments",
        default="last",
        help="Comma-delimited reassorted segment indexes or 'last' (default)",
    )
    parser.add_argument("--num-reassortments", type=int, default=1, help="Number of subtree moves per dataset")
    parser.add_argument("--repeats", type=int, default=5, help="Number of replicate datasets")
    parser.add_argument("--start", type=int, default=1, help="Starting replicate index")
    parser.add_argument("--min-size", type=int, default=1, help="Minimum moved clade size")
    parser.add_argument("--max-size", type=int, default=10, help="Maximum moved clade size")
    parser.add_argument("--min-distance", type=float, default=0.0, help="Minimum normalized reassortment distance")
    parser.add_argument("--max-distance", type=float, default=1.0, help="Maximum normalized reassortment distance")
    parser.add_argument("--mutation-rate", type=float, default=0.003, help="Branch scaling factor for sequence simulation")
    parser.add_argument("--kappa", type=float, default=2.75, help="HKY kappa parameter")
    parser.add_argument(
        "--state-frequencies",
        default="0.25,0.25,0.25,0.25",
        help="Comma-delimited equilibrium nucleotide frequencies",
    )
    parser.add_argument("--random-seed", type=int, default=13, help="Base random seed")
    parser.add_argument("--skip-analysis", action="store_true", help="Only generate datasets; do not run Espalier")
    parser.add_argument("--skip-em", action="store_true", help="Pass --skip-em through to run_arg_analysis.py")
    parser.add_argument("--python-executable", default=sys.executable, help="Python executable used for pipeline runs")
    parser.add_argument("--raxml-path", default="raxml-ng", help="Path to raxml-ng")
    parser.add_argument("--em-iters", type=int, default=1, help="EM iterations for Espalier")
    parser.add_argument("--generation-time-days", type=float, default=3.0, help="Generation time for annualized reporting")
    return parser.parse_args()


def main():
    args = parse_args()
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    written = generate_empirical_reassortment_panel(
        tree_path=args.tree,
        output_dir=str(output_dir),
        repeats=args.repeats,
        start_index=args.start,
        segment_lengths=parse_int_list(args.segment_lengths),
        segment_names=maybe_parse_segment_names(args.segment_names),
        reassorted_segment_indexes=_parse_reassorted_segments(args.reassorted_segments),
        num_reassortments=args.num_reassortments,
        min_clade_size=args.min_size,
        max_clade_size=args.max_size,
        min_normalized_distance=args.min_distance,
        max_normalized_distance=args.max_distance,
        mutation_rate=args.mutation_rate,
        kappa=args.kappa,
        state_frequencies=parse_float_list(args.state_frequencies),
        random_seed=args.random_seed,
    )

    results = []
    for record in written:
        sample_dir = Path(record["output_dir"])
        truth = json.loads(Path(record["truth_path"]).read_text())
        result = {
            "dataset_dir": str(sample_dir),
            "truth_path": record["truth_path"],
            "status": "generated",
            "expected_spr_count": int(truth["num_reassortments"]),
        }
        if not args.skip_analysis:
            analysis_result = run_analysis_pipeline(
                str(sample_dir),
                python_executable=args.python_executable,
                raxml_path=args.raxml_path,
                em_iters=args.em_iters,
                generation_time_days=args.generation_time_days,
                skip_em=args.skip_em,
            )
            breakpoints = _parse_breakpoints(Path(analysis_result["summary_path"]))
            move_count_metrics = _compute_move_count_metrics(
                sample_dir=sample_dir,
                analysis_dir=Path(analysis_result["analysis_dir"]),
                truth=truth,
            )
            result.update(analysis_result)
            result.update(move_count_metrics)
            result["selected_boundary_positions"] = truth["selected_boundary_positions"]
            result["summary_breakpoints"] = breakpoints
            result["boundary_detected"] = any(
                boundary in breakpoints for boundary in truth["selected_boundary_positions"]
            )
            result["status"] = "analyzed"
        results.append(result)
        with open(sample_dir / "result.json", "w") as handle:
            json.dump(result, handle, indent=2, sort_keys=True)
            handle.write("\n")

    summary = {
        "n_runs": len(results),
        "n_analyzed": sum(result["status"] == "analyzed" for result in results),
        "n_boundary_detected": sum(bool(result.get("boundary_detected")) for result in results),
        "n_true_tree_exact_spr": sum(bool(result.get("true_tree_matches_expected")) for result in results),
        "n_arg_local_exact_spr": sum(bool(result.get("arg_local_matches_expected")) for result in results),
        "n_arg_with_recomb_exact_spr": sum(bool(result.get("arg_with_recomb_matches_expected")) for result in results),
        "results": results,
    }
    with open(output_dir / "benchmark_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(f"Wrote benchmark summary to {output_dir / 'benchmark_summary.json'}")


if __name__ == "__main__":
    main()
