#!/usr/bin/env python3

import argparse
import json
import re
import sys
from pathlib import Path

import dendropy
from dendropy.calculate import treecompare

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from Espalier import MAF


def _load_tree(tree_path: Path, taxon_namespace: dendropy.TaxonNamespace) -> dendropy.Tree:
    return dendropy.Tree.get(
        path=str(tree_path),
        schema="newick",
        rooting="default-rooted",
        preserve_underscores=True,
        taxon_namespace=taxon_namespace,
    )


def _spr_distance(tree_a_path: Path, tree_b_path: Path) -> int:
    taxon_namespace = dendropy.TaxonNamespace()
    tree_a = _load_tree(tree_a_path, taxon_namespace)
    tree_b = _load_tree(tree_b_path, taxon_namespace)
    if treecompare.symmetric_difference(tree_a, tree_b) == 0:
        return 0
    return int(MAF.get_spr_dist(tree_a, tree_b))


def _parse_tree_sources(log_path: Path):
    if not log_path.exists():
        return None
    match = re.search(r"Tree sources:\s+(.+)", log_path.read_text())
    return match.group(1).strip() if match else None


def _get_segment_names(truth):
    reassorted_indexes = set(truth["reassorted_segment_indexes"])
    reference_name = None
    reassorted_name = None
    for idx, name in enumerate(truth["segment_names"]):
        if idx in reassorted_indexes and reassorted_name is None:
            reassorted_name = name
        if idx not in reassorted_indexes and reference_name is None:
            reference_name = name
    if reference_name is None or reassorted_name is None:
        raise ValueError("Could not determine reference and reassorted segment names")
    return reference_name, reassorted_name


def analyze_sample(sample_dir: Path):
    truth = json.loads((sample_dir / "truth.json").read_text())
    result = json.loads((sample_dir / "result.json").read_text())
    analysis_dir = Path(result["analysis_dir"])
    reference_name, reassorted_name = _get_segment_names(truth)

    true_reference_tree = sample_dir / "trees" / f"{reference_name}.fasta.rooted.tre"
    true_reassorted_tree = sample_dir / "trees" / f"{reassorted_name}.fasta.rooted.tre"
    normalized_reference_tree = analysis_dir / "derived" / "normalized_trees" / f"{reference_name}.fasta.rooted.tre"
    normalized_reassorted_tree = analysis_dir / "derived" / "normalized_trees" / f"{reassorted_name}.fasta.rooted.tre"
    consensus_reference_tree = analysis_dir / "consensus_reference.tre"
    arg_local_reference_tree = analysis_dir / f"{reference_name}_ARG_local.tre"
    arg_local_reassorted_tree = analysis_dir / f"{reassorted_name}_ARG_local.tre"
    arg_recomb_reference_tree = analysis_dir / f"{reference_name}_ARG_with_recomb.tre"
    arg_recomb_reassorted_tree = analysis_dir / f"{reassorted_name}_ARG_with_recomb.tre"

    return {
        "sample": sample_dir.name,
        "tree_sources": _parse_tree_sources(sample_dir / "analysis.log"),
        "expected_spr_count": int(truth["num_reassortments"]),
        "true_tree_spr_count": _spr_distance(true_reference_tree, true_reassorted_tree),
        "normalized_tree_spr_count": _spr_distance(normalized_reference_tree, normalized_reassorted_tree),
        "true_to_normalized_reference_spr": _spr_distance(true_reference_tree, normalized_reference_tree),
        "true_to_normalized_reassorted_spr": _spr_distance(true_reassorted_tree, normalized_reassorted_tree),
        "consensus_to_reference_spr": _spr_distance(consensus_reference_tree, normalized_reference_tree),
        "consensus_to_reassorted_spr": _spr_distance(consensus_reference_tree, normalized_reassorted_tree),
        "arg_local_spr_count": _spr_distance(arg_local_reference_tree, arg_local_reassorted_tree),
        "arg_with_recomb_spr_count": _spr_distance(arg_recomb_reference_tree, arg_recomb_reassorted_tree),
        "arg_local_reference_to_reference_spr": _spr_distance(arg_local_reference_tree, normalized_reference_tree),
        "arg_local_reference_to_reassorted_spr": _spr_distance(arg_local_reference_tree, normalized_reassorted_tree),
        "arg_local_reassorted_to_reference_spr": _spr_distance(arg_local_reassorted_tree, normalized_reference_tree),
        "arg_local_reassorted_to_reassorted_spr": _spr_distance(arg_local_reassorted_tree, normalized_reassorted_tree),
    }


def parse_args():
    parser = argparse.ArgumentParser(description="Analyze stage-wise SPR collapse in empirical reassortment benchmarks")
    parser.add_argument("--benchmark-dir", required=True, help="Benchmark directory containing sample*/result.json")
    parser.add_argument("--output", default="", help="Optional output JSON path")
    return parser.parse_args()


def main():
    args = parse_args()
    benchmark_dir = Path(args.benchmark_dir)
    sample_dirs = sorted(
        sample_dir
        for sample_dir in benchmark_dir.glob("sample*")
        if (sample_dir / "result.json").exists()
    )
    results = [analyze_sample(sample_dir) for sample_dir in sample_dirs]
    summary = {
        "benchmark_dir": str(benchmark_dir),
        "n_samples": len(results),
        "n_true_tree_exact": sum(item["true_tree_spr_count"] == item["expected_spr_count"] for item in results),
        "n_normalized_tree_exact": sum(item["normalized_tree_spr_count"] == item["expected_spr_count"] for item in results),
        "n_arg_local_exact": sum(item["arg_local_spr_count"] == item["expected_spr_count"] for item in results),
        "n_arg_with_recomb_exact": sum(item["arg_with_recomb_spr_count"] == item["expected_spr_count"] for item in results),
        "results": results,
    }

    if args.output:
        output_path = Path(args.output)
    else:
        output_path = benchmark_dir / "stage_diagnostics.json"
    output_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"Wrote stage diagnostics to {output_path}")


if __name__ == "__main__":
    main()
