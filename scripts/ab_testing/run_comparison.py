#!/usr/bin/env python3
"""
A/B Comparison Runner for Espalier Modernization

This script runs both the original and modernized implementations
on the same test data and compares:
- Runtime performance
- Memory usage
- Output equivalence
- TreeSequence structure
- Likelihood scores (if available)

Usage:
    python run_comparison.py [--iterations N] [--output-dir DIR]

Author: A/B Testing Framework
"""

import os
import sys
import argparse
import datetime
import logging
import time
import json
import traceback
from dataclasses import asdict
from typing import Dict, List, Optional, Tuple
import numpy as np

# Add parent directories to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import dendropy
from Bio import SeqIO

from Espalier.ARGBuilder import ARGBuilder
from Espalier.Reconciler import Reconciler
from Espalier.RAxML import RAxMLRunner
from Espalier.SCARLikelihood import SCAR
from Espalier import Utils

from ab_testing.modern_converter import (
    convert_modern,
    convert_original,
    compare_tree_sequences,
    ConversionMetrics
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class ABTestResult:
    """Container for A/B test results."""

    def __init__(self, name: str):
        self.name = name
        self.original_metrics: Optional[ConversionMetrics] = None
        self.modern_metrics: Optional[ConversionMetrics] = None
        self.original_ts = None
        self.modern_ts = None
        self.original_error: Optional[str] = None
        self.modern_error: Optional[str] = None
        self.comparison: Optional[Dict] = None
        self.original_likelihood: Optional[float] = None
        self.modern_likelihood: Optional[float] = None

    def to_dict(self) -> Dict:
        """Convert to dictionary for JSON serialization."""
        result = {
            'name': self.name,
            'original_error': self.original_error,
            'modern_error': self.modern_error,
        }

        if self.original_metrics:
            result['original_metrics'] = asdict(self.original_metrics)
        if self.modern_metrics:
            result['modern_metrics'] = asdict(self.modern_metrics)
        if self.comparison:
            result['comparison'] = self.comparison
        if self.original_likelihood is not None:
            result['original_likelihood'] = self.original_likelihood
        if self.modern_likelihood is not None:
            result['modern_likelihood'] = self.modern_likelihood

        # Compute speedup if both succeeded
        if self.original_metrics and self.modern_metrics:
            if self.modern_metrics.total_time > 0:
                result['speedup'] = self.original_metrics.total_time / self.modern_metrics.total_time
            result['time_fix_speedup'] = (
                self.original_metrics.time_fix_time / self.modern_metrics.time_fix_time
                if self.modern_metrics.time_fix_time > 0 else None
            )

        return result


def run_viterbi_path(argb: ARGBuilder,
                     ml_tree_files: List[str],
                     seq_files: List[str],
                     ref: dendropy.Tree,
                     rec_rate: float) -> Tuple[List, List]:
    """
    Run Viterbi to get tree path (shared by both implementations).

    Returns:
        tree_path: List of dendropy trees
        tree_intervals: List of (start, end) tuples
    """
    # Get tree path using ARGBuilder
    tree_path = argb.reconstruct_ARG(ml_tree_files, seq_files, ref, rec_rate)

    # Calculate tree intervals from sequence files
    tree_intervals = []
    cumulative = 0
    for sf in seq_files:
        seq_length = len(list(SeqIO.parse(sf, "fasta"))[0])
        tree_intervals.append((cumulative, cumulative + seq_length))
        cumulative += seq_length

    return tree_path, tree_intervals


def run_em_get_tree_path(argb: ARGBuilder,
                         ml_tree_files: List[str],
                         seq_files: List[str],
                         ref: dendropy.Tree,
                         scar_model: SCAR,
                         em_iters: int = 5) -> Tuple[List, List, float, int]:
    """
    Run EM algorithm to get tree path with recombination nodes.

    Returns:
        tree_path: List of dendropy trees with recombination nodes
        tree_intervals: List of (start, end) tuples
        rec_rate: Estimated recombination rate
        n_recomb: Number of recombination events
    """
    # Run EM - this returns (ts, rec_rate, n_recomb, tree_path)
    # We want the tree_path for our comparison
    result = argb.run_EM(ml_tree_files, seq_files, ref, scar_model, iters=em_iters)

    # Unpack based on what run_EM returns
    if len(result) == 4:
        _, rec_rate, n_recomb, tree_path = result
    else:
        # Fallback for older versions
        _, rec_rate, n_recomb = result[:3]
        tree_path = None

    # Calculate tree intervals
    tree_intervals = []
    cumulative = 0
    for sf in seq_files:
        seq_length = len(list(SeqIO.parse(sf, "fasta"))[0])
        tree_intervals.append((cumulative, cumulative + seq_length))
        cumulative += seq_length

    return tree_path, tree_intervals, rec_rate, int(n_recomb)


def run_ab_test(tree_path: List,
                tree_intervals: List,
                test_name: str = "conversion") -> ABTestResult:
    """
    Run A/B comparison between original and modern converter.

    Parameters:
        tree_path: List of dendropy trees
        tree_intervals: List of (start, end) tuples
        test_name: Name for this test

    Returns:
        ABTestResult with comparison data
    """
    result = ABTestResult(test_name)

    # Run original implementation
    logger.info(f"\n{'='*60}")
    logger.info(f"Running ORIGINAL implementation for: {test_name}")
    logger.info(f"{'='*60}")

    try:
        result.original_ts, result.original_metrics = convert_original(
            tree_path, tree_intervals
        )
    except Exception as e:
        result.original_error = f"{type(e).__name__}: {str(e)}"
        logger.error(f"Original failed: {result.original_error}")
        logger.debug(traceback.format_exc())

    # Run modern implementation
    logger.info(f"\n{'='*60}")
    logger.info(f"Running MODERN implementation for: {test_name}")
    logger.info(f"{'='*60}")

    try:
        result.modern_ts, result.modern_metrics = convert_modern(
            tree_path, tree_intervals
        )
    except Exception as e:
        result.modern_error = f"{type(e).__name__}: {str(e)}"
        logger.error(f"Modern failed: {result.modern_error}")
        logger.debug(traceback.format_exc())

    # Compare results if both succeeded
    if result.original_ts is not None and result.modern_ts is not None:
        result.comparison = compare_tree_sequences(
            result.original_ts,
            result.modern_ts,
            label_a="Original",
            label_b="Modern"
        )

    return result


def compute_likelihood(ts, scar_model: SCAR) -> Optional[float]:
    """
    Compute likelihood of TreeSequence under SCAR model.

    Parameters:
        ts: tskit TreeSequence
        scar_model: SCAR coalescent model

    Returns:
        Log-likelihood value or None if failed
    """
    try:
        # SCAR.likelihood expects a tree or tree sequence
        # This is a simplified version - full implementation would use SCAR properly
        if ts is None:
            return None

        # For now, return number of recombination events as a proxy
        # (proper likelihood requires full SCAR implementation)
        recomb_nodes = np.sum(ts.tables.nodes.flags == 131072)
        n_recomb = recomb_nodes // 2

        # Simple likelihood proxy based on recombination rate and count
        rec_rate = scar_model.rec_rate
        genome_length = ts.sequence_length

        # Poisson-based likelihood for recombination count
        expected_recomb = rec_rate * genome_length
        if expected_recomb > 0:
            from scipy.stats import poisson
            log_lik = poisson.logpmf(n_recomb, expected_recomb)
            return float(log_lik)
        return 0.0

    except Exception as e:
        logger.warning(f"Likelihood computation failed: {e}")
        return None


def print_comparison_report(results: List[ABTestResult], output_file: Optional[str] = None):
    """Print formatted comparison report."""

    lines = []
    lines.append("\n" + "=" * 80)
    lines.append("A/B TESTING COMPARISON REPORT")
    lines.append("=" * 80)
    lines.append(f"Generated: {datetime.datetime.now().isoformat()}")
    lines.append("")

    for result in results:
        lines.append("-" * 80)
        lines.append(f"Test: {result.name}")
        lines.append("-" * 80)

        if result.original_error:
            lines.append(f"  Original: FAILED - {result.original_error}")
        elif result.original_metrics:
            m = result.original_metrics
            lines.append(f"  Original:")
            lines.append(f"    Total time:     {m.total_time:.4f}s")
            lines.append(f"    Time fix time:  {m.time_fix_time:.4f}s")
            lines.append(f"    Nodes:          {m.n_nodes}")
            lines.append(f"    Edges:          {m.n_edges}")
            lines.append(f"    Recomb nodes:   {m.n_recomb_nodes}")
            lines.append(f"    Time fixes:     {m.n_time_fixes}")

        if result.modern_error:
            lines.append(f"  Modern: FAILED - {result.modern_error}")
        elif result.modern_metrics:
            m = result.modern_metrics
            lines.append(f"  Modern:")
            lines.append(f"    Total time:     {m.total_time:.4f}s")
            lines.append(f"    Time fix time:  {m.time_fix_time:.4f}s")
            lines.append(f"    Table sort:     {m.table_sort_time:.4f}s")
            lines.append(f"    Extend hapl:    {m.extend_haplotypes_time:.4f}s")
            lines.append(f"    Nodes:          {m.n_nodes}")
            lines.append(f"    Edges:          {m.n_edges}")
            lines.append(f"    Recomb nodes:   {m.n_recomb_nodes}")
            lines.append(f"    Time fixes:     {m.n_time_fixes}")
            lines.append(f"    Iterations:     {m.n_iterations}")

        # Speedup calculation
        if result.original_metrics and result.modern_metrics:
            if result.modern_metrics.total_time > 0:
                speedup = result.original_metrics.total_time / result.modern_metrics.total_time
                lines.append(f"  Speedup: {speedup:.2f}x")

            if result.modern_metrics.time_fix_time > 0 and result.original_metrics.time_fix_time > 0:
                time_fix_speedup = result.original_metrics.time_fix_time / result.modern_metrics.time_fix_time
                lines.append(f"  Time fix speedup: {time_fix_speedup:.2f}x")

        # Comparison results
        if result.comparison:
            lines.append(f"  Comparison:")
            if result.comparison['equivalent']:
                lines.append(f"    TreeSequences are EQUIVALENT")
            else:
                lines.append(f"    TreeSequences DIFFER:")
                for diff in result.comparison['differences']:
                    lines.append(f"      - {diff}")

        # Likelihood comparison
        if result.original_likelihood is not None or result.modern_likelihood is not None:
            lines.append(f"  Likelihood:")
            if result.original_likelihood is not None:
                lines.append(f"    Original: {result.original_likelihood:.4f}")
            if result.modern_likelihood is not None:
                lines.append(f"    Modern:   {result.modern_likelihood:.4f}")

        lines.append("")

    # Summary
    lines.append("=" * 80)
    lines.append("SUMMARY")
    lines.append("=" * 80)

    n_original_success = sum(1 for r in results if r.original_error is None and r.original_metrics)
    n_modern_success = sum(1 for r in results if r.modern_error is None and r.modern_metrics)
    n_equivalent = sum(1 for r in results if r.comparison and r.comparison['equivalent'])

    lines.append(f"Tests run:           {len(results)}")
    lines.append(f"Original succeeded:  {n_original_success}/{len(results)}")
    lines.append(f"Modern succeeded:    {n_modern_success}/{len(results)}")
    lines.append(f"Equivalent outputs:  {n_equivalent}/{len(results)}")

    # Average speedup
    speedups = []
    for r in results:
        if r.original_metrics and r.modern_metrics and r.modern_metrics.total_time > 0:
            speedups.append(r.original_metrics.total_time / r.modern_metrics.total_time)
    if speedups:
        lines.append(f"Average speedup:     {np.mean(speedups):.2f}x")

    lines.append("")

    report = "\n".join(lines)
    print(report)

    if output_file:
        with open(output_file, 'w') as f:
            f.write(report)
        logger.info(f"Report saved to: {output_file}")


def main():
    parser = argparse.ArgumentParser(description='Run A/B comparison tests for Espalier modernization')
    parser.add_argument('--iterations', type=int, default=1,
                        help='Number of iterations for each test (default: 1)')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Output directory for results')
    parser.add_argument('--em-iters', type=int, default=5,
                        help='Number of EM iterations (default: 5)')
    parser.add_argument('--skip-em', action='store_true',
                        help='Skip EM-based tests (faster)')
    args = parser.parse_args()

    # Configuration
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ALIGNMENT_DIR = os.path.join(BASE_DIR, 'TestFiles', 'alignments')
    TREE_DIR = os.path.join(BASE_DIR, 'TestFiles', 'trees')

    # Setup output directory
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d_%H%M%S')
    if args.output_dir:
        output_dir = args.output_dir
    else:
        output_dir = os.path.join(BASE_DIR, 'TestFiles', 'Outputs', f'ab_test_{timestamp}')
    os.makedirs(output_dir, exist_ok=True)

    logger.info("=" * 60)
    logger.info("Espalier A/B Testing Framework")
    logger.info("=" * 60)
    logger.info(f"Output directory: {output_dir}")

    # Get sequence files
    seq_files = sorted([
        os.path.join(ALIGNMENT_DIR, f)
        for f in os.listdir(ALIGNMENT_DIR)
        if f.endswith('.fasta')
    ])

    logger.info(f"\nInput alignments ({len(seq_files)} segments):")
    for sf in seq_files:
        logger.info(f"  - {os.path.basename(sf)}")

    # Get tree files
    ml_tree_files = []
    for sf in seq_files:
        seg_name = os.path.basename(sf).replace('_aligned_common.fasta', '')
        tree_file = os.path.join(TREE_DIR, f'{seg_name}_tree_common_binary.tre')
        if os.path.exists(tree_file):
            ml_tree_files.append(tree_file)
        else:
            logger.error(f"Tree file not found: {tree_file}")
            sys.exit(1)

    # Calculate genome length
    genome_length = sum(len(list(SeqIO.parse(sf, "fasta"))[0]) for sf in seq_files)
    logger.info(f"Total genome length: {genome_length} bp")

    # Build consensus reference tree
    ref = Utils.get_consensus_tree(ml_tree_files, root=True)
    taxa = ref.taxon_namespace

    # Initialize Espalier components
    temp_dir = os.path.join(output_dir, 'temp')
    os.makedirs(temp_dir, exist_ok=True)

    raxml = RAxMLRunner(raxml_path='raxml-ng', temp_dir=temp_dir)
    reconciler = Reconciler(raxml, lower_bound_ratio=0.1, prior_gamma=1e-4, temp_dir=temp_dir)
    argb = ARGBuilder(reconciler, raxml, temp_dir=temp_dir)

    # Initialize SCAR model
    INITIAL_REC_RATE = 1e-4
    scar_model = SCAR(
        rec_rate=INITIAL_REC_RATE,
        M=[[0]],
        Ne=1.0,
        genome_length=genome_length,
        bounds=(0, 0.01)
    )

    all_results = []

    # Test 1: Basic Viterbi path conversion
    logger.info("\n" + "=" * 60)
    logger.info("TEST 1: Basic Viterbi Path Conversion")
    logger.info("=" * 60)

    for i in range(args.iterations):
        logger.info(f"\nIteration {i+1}/{args.iterations}")

        tree_path, tree_intervals = run_viterbi_path(
            argb, ml_tree_files, seq_files, ref, INITIAL_REC_RATE
        )

        result = run_ab_test(tree_path, tree_intervals, f"viterbi_iter_{i+1}")
        all_results.append(result)

    # Test 2: EM-based path with recombination nodes
    if not args.skip_em:
        logger.info("\n" + "=" * 60)
        logger.info("TEST 2: EM-based Path with Recombination Nodes")
        logger.info("=" * 60)

        for i in range(args.iterations):
            logger.info(f"\nIteration {i+1}/{args.iterations}")

            try:
                tree_path, tree_intervals, rec_rate, n_recomb = run_em_get_tree_path(
                    argb, ml_tree_files, seq_files, ref, scar_model, em_iters=args.em_iters
                )

                if tree_path is not None:
                    result = run_ab_test(tree_path, tree_intervals, f"em_iter_{i+1}")

                    # Compute likelihoods
                    if result.original_ts:
                        result.original_likelihood = compute_likelihood(result.original_ts, scar_model)
                    if result.modern_ts:
                        result.modern_likelihood = compute_likelihood(result.modern_ts, scar_model)

                    all_results.append(result)
                else:
                    logger.warning("EM did not return tree path")

            except Exception as e:
                logger.error(f"EM test failed: {e}")
                logger.debug(traceback.format_exc())

    # Generate report
    report_file = os.path.join(output_dir, 'comparison_report.txt')
    print_comparison_report(all_results, report_file)

    # Save JSON results
    json_file = os.path.join(output_dir, 'comparison_results.json')
    with open(json_file, 'w') as f:
        json.dump([r.to_dict() for r in all_results], f, indent=2, default=str)
    logger.info(f"JSON results saved to: {json_file}")

    # Save TreeSequences
    for result in all_results:
        if result.original_ts:
            ts_file = os.path.join(output_dir, f'{result.name}_original.trees')
            result.original_ts.dump(ts_file)
        if result.modern_ts:
            ts_file = os.path.join(output_dir, f'{result.name}_modern.trees')
            result.modern_ts.dump(ts_file)

    logger.info("\n" + "=" * 60)
    logger.info("A/B TESTING COMPLETE")
    logger.info("=" * 60)
    logger.info(f"Results saved to: {output_dir}")

    return output_dir


if __name__ == '__main__':
    main()
