#!/usr/bin/env python3
"""
ARG Reconstruction and Recombination Rate Estimation Script

This script demonstrates the full Espalier pipeline:
1. Infer ML trees for each genomic segment using RAxML-NG
2. Build a consensus reference tree
3. Reconstruct the ARG using the Viterbi algorithm
4. Estimate recombination rate using EM algorithm with SCAR model
5. Extract and visualize recombination nodes/breakpoints

For reassortment analysis (e.g., segmented viruses), each segment file
represents a breakpoint boundary.

Usage:
    # Default (Hantavirus test data):
    python run_arg_analysis.py

    # Custom dataset:
    python run_arg_analysis.py --alignments TestFiles/FluTest/alignments --trees TestFiles/FluTest/trees

    # With custom parameters:
    python run_arg_analysis.py --alignments /path/to/alignments --trees /path/to/trees --em-iters 10 --rec-rate 1e-5

Author: Generated for Espalier analysis
"""

import os
import sys
import argparse
import datetime
import logging
import numpy as np
import dendropy

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# Add scripts directory for ab_testing module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from Espalier.ARGBuilder import ARGBuilder
from Espalier.Reassortment import (
    ANALYSIS_MODE_REASSORTMENT,
    ANALYSIS_MODE_RECOMBINATION,
    find_matching_tree,
    prepare_analysis_inputs,
)
from Espalier.Reconciler import Reconciler
from Espalier.RAxML import RAxMLRunner
from Espalier.SCARLikelihood import SCAR
from Espalier import Utils

# Modern converter for improved TreeSequence conversion
try:
    from ab_testing.modern_converter import convert_modern, ConversionMetrics
    MODERN_CONVERTER_AVAILABLE = True
except ImportError:
    MODERN_CONVERTER_AVAILABLE = False

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def setup_directories(base_dir, output_base=None):
    """Create output and temp directories.

    Parameters:
        base_dir: Base directory of the project
        output_base: Base directory for outputs (default: TestFiles/Outputs)

    Returns:
        output_dir: Timestamped output directory
        temp_dir: Temp directory within output_dir
    """
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d_%H%M%S')

    if output_base is None:
        output_base = os.path.join(base_dir, 'TestFiles', 'Outputs')

    os.makedirs(output_base, exist_ok=True)

    output_dir = os.path.join(output_base, f'arg_output_{timestamp}')
    temp_dir = os.path.join(output_dir, 'temp')

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)

    return output_dir, temp_dir


def infer_ml_trees(seq_files, output_dir, temp_dir, raxml_path='raxml-ng'):
    """
    Infer ML trees for each segment using RAxML-NG.

    Parameters:
        seq_files: List of FASTA alignment files
        output_dir: Directory for output trees
        temp_dir: Directory for temp files
        raxml_path: Path to raxml-ng executable

    Returns:
        List of ML tree file paths
    """
    logger.info("=" * 60)
    logger.info("Step 1: Inferring ML trees for each segment")
    logger.info("=" * 60)

    raxml = RAxMLRunner(raxml_path=raxml_path, temp_dir=temp_dir)
    ml_tree_files = []

    for i, seq_file in enumerate(seq_files):
        segment_name = os.path.basename(seq_file).replace('.fasta', '').replace('_aligned_common', '')
        logger.info(f"  Inferring tree for segment {i+1}/{len(seq_files)}: {segment_name}")

        tree_file = os.path.join(output_dir, f'{segment_name}_ML.tre')

        try:
            # Run RAxML-NG to get ML tree
            raxml.get_raxml_tree(seq_file, tree_file)
            ml_tree_files.append(tree_file)
            logger.info(f"    -> Saved: {tree_file}")
        except Exception as e:
            logger.error(f"    -> Failed to infer tree: {e}")
            raise

    return ml_tree_files


def get_reference_tree(ml_tree_files, output_dir):
    """
    Compute a consensus reference tree from local ML trees.

    Parameters:
        ml_tree_files: List of ML tree file paths
        output_dir: Directory for output

    Returns:
        dendropy.Tree: Consensus reference tree
    """
    logger.info("=" * 60)
    logger.info("Step 2: Building consensus reference tree")
    logger.info("=" * 60)

    ref = Utils.get_consensus_tree(ml_tree_files, root=True)

    # Save reference tree
    ref_file = os.path.join(output_dir, 'consensus_reference.tre')
    ref.write(path=ref_file, schema='newick', suppress_annotations=True, suppress_rooting=True)
    logger.info(f"  Consensus tree saved: {ref_file}")

    return ref


def reconstruct_arg_basic(argb, ml_tree_files, seq_files, ref, rec_rate, output_dir):
    """
    Basic ARG reconstruction using fixed recombination rate.

    Parameters:
        argb: ARGBuilder instance
        ml_tree_files: List of ML tree files
        seq_files: List of sequence files
        ref: Reference tree
        rec_rate: Fixed recombination rate per site
        output_dir: Output directory

    Returns:
        tree_path: List of local trees in the ARG
    """
    logger.info("=" * 60)
    logger.info("Step 3a: Basic ARG reconstruction (fixed rec rate)")
    logger.info("=" * 60)
    logger.info(f"  Using recombination rate: {rec_rate:.2e} per site")

    tree_path = argb.reconstruct_ARG(ml_tree_files, seq_files, ref, rec_rate)

    # Save reconstructed local trees
    logger.info("  Saving ARG local trees...")
    arg_tree_files = []
    for idx, tree in enumerate(tree_path):
        seg_name = os.path.basename(seq_files[idx]).replace('.fasta', '').replace('_aligned_common', '')
        tree_file = os.path.join(output_dir, f'{seg_name}_ARG_local.tre')
        tree.write(path=tree_file, schema='newick', suppress_annotations=True, suppress_rooting=True)
        arg_tree_files.append(tree_file)
        logger.info(f"    -> {tree_file}")

    logger.info(f"  Tree sources: {' -> '.join(argb.tree_source_path)}")

    return tree_path, arg_tree_files


def run_em_estimation(argb, ml_tree_files, seq_files, ref, scar_model, output_dir, em_iters=10, use_modern_converter=False):
    """
    Run EM algorithm for joint ARG reconstruction and recombination rate estimation.

    Parameters:
        argb: ARGBuilder instance
        ml_tree_files: List of ML tree files
        seq_files: List of sequence files
        ref: Reference tree
        scar_model: SCAR coalescent model
        output_dir: Output directory
        em_iters: Number of EM iterations
        use_modern_converter: Use the modernized TreeSequence converter

    Returns:
        ts: TreeSequence with full ARG
        rec_rate_est: Estimated recombination rate
        n_recomb: Number of inferred recombination events
    """
    logger.info("=" * 60)
    logger.info("Step 3b: EM-based ARG reconstruction & rate estimation")
    logger.info("=" * 60)
    logger.info(f"  Initial recombination rate: {scar_model.rec_rate:.2e}")
    logger.info(f"  Running {em_iters} EM iterations...")
    if use_modern_converter:
        logger.info("  Using MODERN converter (3x faster, more robust)")

    ts, rec_rate_est, n_recomb, tree_path_with_recs = argb.run_EM(
        ml_tree_files,
        seq_files,
        ref,
        scar_model,
        iters=em_iters,
        return_tree_path=True,
    )

    logger.info(f"  Final estimated recombination rate: {rec_rate_est:.6f} per site")
    logger.info(f"  Inferred recombination events: {int(n_recomb)}")

    # Save trees with recombination nodes (even if TreeSequence conversion fails)
    if tree_path_with_recs is not None:
        logger.info("  Saving ARG trees with recombination nodes...")
        for idx, tree in enumerate(tree_path_with_recs):
            seg_name = os.path.basename(seq_files[idx]).replace('.fasta', '').replace('_aligned_common', '').replace('.aln', '')
            tree_file = os.path.join(output_dir, f'{seg_name}_ARG_with_recomb.tre')
            tree.write(path=tree_file, schema='newick', suppress_annotations=True, suppress_rooting=True)
            logger.info(f"    -> {tree_file}")

            # Count recombination nodes (unifurcations) in the tree
            n_rec_nodes = sum(1 for node in tree.preorder_node_iter() if len(node.child_nodes()) == 1 and node.parent_node is not None)
            logger.info(f"       Recombination nodes in {seg_name}: {n_rec_nodes}")

    # Try modern converter if requested and original failed or for better performance
    if use_modern_converter and MODERN_CONVERTER_AVAILABLE and tree_path_with_recs is not None:
        logger.info("  Converting tree path using modern converter...")
        from Bio import SeqIO
        tree_intervals = []
        cumulative = 0
        for sf in seq_files:
            seq_length = len(list(SeqIO.parse(sf, "fasta"))[0])
            tree_intervals.append((cumulative, cumulative + seq_length))
            cumulative += seq_length

        try:
            ts, metrics = convert_modern(tree_path_with_recs, tree_intervals)
            logger.info(f"  Modern converter: {metrics.n_nodes} nodes, {metrics.n_edges} edges, "
                       f"{metrics.n_recomb_nodes} recomb nodes in {metrics.total_time:.3f}s")
        except Exception as e:
            logger.warning(f"  Modern converter failed: {e}")
            logger.info("  Falling back to original TreeSequence if available")

    # Save TreeSequence if available
    if ts is not None:
        ts_file = os.path.join(output_dir, 'arg_treesequence.trees')
        ts.dump(ts_file)
        logger.info(f"  TreeSequence saved: {ts_file}")
    else:
        logger.warning("  TreeSequence conversion failed - no .trees file saved")
        # Don't raise - we still have the trees with recombination nodes
        if tree_path_with_recs is None:
            raise ValueError("TreeSequence conversion failed and no tree path available")

    return ts, rec_rate_est, n_recomb


def extract_recombination_info(ts, seq_files, output_dir):
    """
    Extract and summarize recombination nodes and breakpoints from TreeSequence.

    Parameters:
        ts: tskit TreeSequence
        seq_files: List of sequence files (for segment info)
        output_dir: Output directory

    Returns:
        dict: Recombination summary information
    """
    logger.info("=" * 60)
    logger.info("Step 4: Extracting recombination information")
    logger.info("=" * 60)

    # Get recombination nodes (flag = 131072)
    recomb_flag = 131072
    node_flags = ts.tables.nodes.flags
    node_times = ts.tables.nodes.time

    recomb_node_ids = np.where(node_flags == recomb_flag)[0]
    n_recomb_nodes = len(recomb_node_ids) // 2  # Each event creates 2 nodes

    logger.info(f"  Recombination nodes found: {len(recomb_node_ids)} ({n_recomb_nodes} events)")

    # Get breakpoints (where tree topology changes)
    breakpoints = list(ts.breakpoints())
    logger.info(f"  Genomic breakpoints: {len(breakpoints) - 2} internal breakpoints")
    logger.info(f"  Breakpoint positions: {breakpoints}")

    # Extract recombination node details
    recomb_info = []
    for node_id in recomb_node_ids:
        time = node_times[node_id]

        # Find edges involving this recombination node
        edges = ts.tables.edges
        parent_edges = edges.parent == node_id
        child_edges = edges.child == node_id

        if np.any(parent_edges):
            edge_idx = np.where(parent_edges)[0][0]
            left = edges.left[edge_idx]
            right = edges.right[edge_idx]
        else:
            left, right = None, None

        recomb_info.append({
            'node_id': node_id,
            'time': time,
            'genomic_left': left,
            'genomic_right': right
        })

    # Save recombination summary
    summary_file = os.path.join(output_dir, 'recombination_summary.txt')
    with open(summary_file, 'w') as f:
        f.write("ARG Recombination Summary\n")
        f.write("=" * 60 + "\n\n")

        f.write("Genomic Segments:\n")
        for i, sf in enumerate(seq_files):
            f.write(f"  Segment {i+1}: {os.path.basename(sf)}\n")
        f.write("\n")

        f.write(f"Total recombination events: {n_recomb_nodes}\n")
        f.write(f"Breakpoints: {breakpoints}\n\n")

        f.write("Recombination Node Details:\n")
        f.write("-" * 60 + "\n")
        for info in recomb_info:
            f.write(f"  Node {info['node_id']}: time={info['time']:.4f}, ")
            f.write(f"interval=[{info['genomic_left']}, {info['genomic_right']})\n")

    logger.info(f"  Summary saved: {summary_file}")

    # Print local tree topologies
    logger.info("\n  Local trees across genome:")
    for tree in ts.trees():
        interval = tree.interval
        logger.info(f"    [{interval.left:.0f}, {interval.right:.0f}): {tree.num_roots} root(s)")

    return {
        'n_recomb_events': n_recomb_nodes,
        'recomb_nodes': recomb_info,
        'breakpoints': breakpoints,
        'n_local_trees': ts.num_trees
    }


def map_recombinations_to_segments(recomb_info, seq_files):
    """
    Map recombination events to segment boundaries (for reassortment analysis).

    This helps identify which segment boundaries have recombination/reassortment.
    """
    logger.info("=" * 60)
    logger.info("Step 5: Mapping recombinations to segment boundaries")
    logger.info("=" * 60)

    # Calculate segment boundaries
    segment_boundaries = [0]
    cumulative = 0
    segment_names = []

    for sf in seq_files:
        from Bio import SeqIO
        seg_length = len(list(SeqIO.parse(sf, "fasta"))[0])
        cumulative += seg_length
        segment_boundaries.append(cumulative)
        segment_names.append(os.path.basename(sf).replace('.fasta', '').replace('_aligned_common', ''))

    logger.info("  Segment boundaries:")
    for i, (name, start, end) in enumerate(zip(segment_names, segment_boundaries[:-1], segment_boundaries[1:])):
        logger.info(f"    {name}: [{start}, {end})")

    # Check which boundaries have topology changes
    breakpoints = recomb_info['breakpoints']
    logger.info("\n  Topology changes at boundaries:")
    for i in range(len(segment_names) - 1):
        boundary = segment_boundaries[i + 1]
        has_break = boundary in breakpoints or any(
            abs(bp - boundary) < 10 for bp in breakpoints
        )
        status = "REASSORTMENT DETECTED" if has_break else "no change"
        logger.info(f"    {segment_names[i]} <-> {segment_names[i+1]}: {status}")

    return segment_boundaries, segment_names


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Espalier ARG Reconstruction and Recombination Rate Estimation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Default (Hantavirus test data):
  python run_arg_analysis.py

  # Flu dataset:
  python run_arg_analysis.py --alignments TestFiles/FluTest/alignments --trees TestFiles/FluTest/trees

  # Custom parameters:
  python run_arg_analysis.py --alignments /path/to/aln --trees /path/to/trees --em-iters 10
        """
    )

    parser.add_argument('--alignments', '-a', type=str, default=None,
                        help='Directory containing alignment files (.fasta or .aln)')
    parser.add_argument('--trees', '-t', type=str, default=None,
                        help='Directory containing pre-computed tree files (.tre)')
    parser.add_argument('--output', '-o', type=str, default=None,
                        help='Output directory (default: TestFiles/Outputs/arg_output_TIMESTAMP)')
    parser.add_argument('--rec-rate', type=float, default=1e-4,
                        help='Initial recombination rate per site (default: 1e-4)')
    parser.add_argument('--em-iters', type=int, default=5,
                        help='Number of EM iterations (default: 5)')
    parser.add_argument('--ne', type=float, default=1.0,
                        help='Effective population size in coalescent units (default: 1.0)')
    parser.add_argument('--raxml-path', type=str, default='raxml-ng',
                        help='Path to raxml-ng executable (default: raxml-ng)')
    parser.add_argument('--skip-em', action='store_true',
                        help='Skip EM estimation, only run basic ARG reconstruction')
    parser.add_argument('--use-modern-converter', action='store_true',
                        help='Use the modernized TreeSequence converter (3x faster, more robust)')
    parser.add_argument('--analysis-mode', choices=[ANALYSIS_MODE_RECOMBINATION, ANALYSIS_MODE_REASSORTMENT],
                        default=ANALYSIS_MODE_RECOMBINATION,
                        help='Choose exact-label recombination analysis or isolate-normalized reassortment mode')
    parser.add_argument('--isolate-field', type=int, default=0,
                        help='Pipe-delimited FASTA field used as the isolate key in reassortment mode (default: 0)')

    return parser.parse_args()


def main():
    """Main analysis pipeline."""

    args = parse_args()

    # Configuration
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    # Set alignment and tree directories
    if args.alignments:
        ALIGNMENT_DIR = args.alignments if os.path.isabs(args.alignments) else os.path.join(BASE_DIR, args.alignments)
    else:
        ALIGNMENT_DIR = os.path.join(BASE_DIR, 'TestFiles', 'alignments')

    if args.trees:
        TREE_DIR = args.trees if os.path.isabs(args.trees) else os.path.join(BASE_DIR, args.trees)
    else:
        TREE_DIR = os.path.join(BASE_DIR, 'TestFiles', 'trees')

    # Analysis parameters
    INITIAL_REC_RATE = args.rec_rate
    EM_ITERATIONS = args.em_iters
    NE = args.ne
    RAXML_PATH = args.raxml_path
    ANALYSIS_MODE = args.analysis_mode

    logger.info("=" * 60)
    logger.info("Espalier ARG Analysis Pipeline")
    logger.info("=" * 60)
    logger.info(f"Analysis mode: {ANALYSIS_MODE}")

    # Setup directories
    output_dir, temp_dir = setup_directories(BASE_DIR, args.output)
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Alignment directory: {ALIGNMENT_DIR}")
    logger.info(f"Tree directory: {TREE_DIR}")

    # Get sequence files - support both .fasta and .aln extensions
    raw_seq_files = sorted([
        os.path.join(ALIGNMENT_DIR, f)
        for f in os.listdir(ALIGNMENT_DIR)
        if f.endswith('.fasta') or f.endswith('.aln')
    ])

    logger.info(f"\nInput alignments ({len(raw_seq_files)} segments):")
    for sf in raw_seq_files:
        logger.info(f"  - {os.path.basename(sf)}")

    # Check if pre-computed trees exist
    use_precomputed = True
    precomputed_trees = []
    for sf in raw_seq_files:
        tree_file = find_matching_tree(sf, TREE_DIR)
        if tree_file:
            precomputed_trees.append(tree_file)
        else:
            use_precomputed = False
            logger.warning(f"No matching tree found for: {os.path.basename(sf)}")
            break

    if use_precomputed and len(precomputed_trees) == len(raw_seq_files):
        logger.info(f"\nUsing pre-computed trees from {TREE_DIR}")
        raw_tree_files = precomputed_trees
        for tf in raw_tree_files:
            logger.info(f"  - {os.path.basename(tf)}")
    else:
        raw_tree_files = None

    prepared_inputs = prepare_analysis_inputs(
        raw_seq_files,
        raw_tree_files,
        output_dir,
        analysis_mode=ANALYSIS_MODE,
        isolate_field=args.isolate_field,
    )
    seq_files = prepared_inputs.seq_files
    if prepared_inputs.label_map_path:
        logger.info(f"Derived label map: {prepared_inputs.label_map_path}")
        logger.info(f"Normalized alignments: {prepared_inputs.normalized_alignments_dir}")
        if prepared_inputs.tree_files:
            logger.info(f"Normalized trees: {prepared_inputs.normalized_trees_dir}")

    # Calculate total genome length after any reassortment normalization.
    from Bio import SeqIO
    genome_length = sum(
        len(list(SeqIO.parse(sf, "fasta"))[0])
        for sf in seq_files
    )
    logger.info(f"\nTotal genome length: {genome_length} bp")

    if prepared_inputs.tree_files:
        ml_tree_files = prepared_inputs.tree_files
    else:
        infer_tree_dir = (
            prepared_inputs.normalized_trees_dir
            if ANALYSIS_MODE == ANALYSIS_MODE_REASSORTMENT
            else output_dir
        )
        try:
            ml_tree_files = infer_ml_trees(seq_files, infer_tree_dir, temp_dir, RAXML_PATH)
        except Exception as e:
            logger.error(f"RAxML-NG failed: {e}")
            raise RuntimeError(f"Could not find or infer trees for all segments: {e}")

    # Build consensus reference tree
    ref = get_reference_tree(ml_tree_files, output_dir)

    # Initialize Espalier components
    logger.info("\nInitializing Espalier components...")
    raxml = RAxMLRunner(raxml_path=RAXML_PATH, temp_dir=temp_dir)
    reconciler = Reconciler(raxml, lower_bound_ratio=0.1, prior_gamma=INITIAL_REC_RATE, temp_dir=temp_dir)
    argb = ARGBuilder(reconciler, raxml, temp_dir=temp_dir)

    # Basic ARG reconstruction
    tree_path, arg_tree_files = reconstruct_arg_basic(
        argb, ml_tree_files, seq_files, ref, INITIAL_REC_RATE, output_dir
    )

    # EM-based rate estimation (optional)
    ts = None
    rec_rate_est = INITIAL_REC_RATE
    n_recomb = 0

    if args.skip_em:
        logger.info("\nSkipping EM estimation (--skip-em flag set)")
    else:
        logger.info("\nInitializing SCAR coalescent model...")
        scar_model = SCAR(
            rec_rate=INITIAL_REC_RATE,
            M=[[0]],  # Single population (no migration)
            Ne=NE,
            genome_length=genome_length,
            bounds=(0, 0.01)  # Reasonable bounds for rec rate
        )

        try:
            ts, rec_rate_est, n_recomb = run_em_estimation(
                argb, ml_tree_files, seq_files, ref, scar_model, output_dir, EM_ITERATIONS,
                use_modern_converter=args.use_modern_converter
            )

            # Extract recombination information
            recomb_info = extract_recombination_info(ts, seq_files, output_dir)

            if ANALYSIS_MODE == ANALYSIS_MODE_REASSORTMENT:
                map_recombinations_to_segments(recomb_info, seq_files)

        except Exception as e:
            logger.error(f"EM estimation failed: {e}")
            logger.info("Continuing with basic ARG reconstruction results...")
            import traceback
            logger.debug(traceback.format_exc())

    # Final summary
    logger.info("\n" + "=" * 60)
    logger.info("ANALYSIS COMPLETE")
    logger.info("=" * 60)
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Number of segments: {len(seq_files)}")
    logger.info(f"Estimated recombination rate: {rec_rate_est:.6f} per site")
    if ts:
        logger.info(f"Recombination events detected: {int(n_recomb)}")
        logger.info(f"TreeSequence file: {os.path.join(output_dir, 'arg_treesequence.trees')}")

    # Cleanup temp directory (optional)
    # shutil.rmtree(temp_dir)

    return output_dir


if __name__ == '__main__':
    output_dir = main()
