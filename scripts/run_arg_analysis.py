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
import json
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
from Espalier.SCARLikelihood import BoundarySCAR, SCAR
from Espalier import Utils
from Espalier.ARGNodeTypes import RECOMBINANT_FLAG, summarize_recombination_events

# Modern converter for improved TreeSequence conversion
try:
    from Espalier.ModernDendro2TSConverter import convert as modern_convert
    MODERN_CONVERTER_AVAILABLE = True
except ImportError:
    modern_convert = None
    MODERN_CONVERTER_AVAILABLE = False

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


DEFAULT_RECOMBINATION_REC_RATE = 1e-4
DEFAULT_REASSORTMENT_REC_RATE_PER_YEAR = 0.1
DEFAULT_REASSORTMENT_RATE_UPPER_PER_YEAR = 3.0
DEFAULT_REASSORTMENT_REC_RATE_PER_DIVERGENCE = 0.1
DEFAULT_REASSORTMENT_RATE_UPPER_PER_DIVERGENCE = 10.0

TREE_TIME_UNIT_DIVERGENCE = "divergence"
TREE_TIME_UNIT_TIME = "time"
TREE_TIME_UNIT_GENERATIONS = "generations"


def annual_rate_to_generation_rate(rate_per_year, generation_time_days):
    return float(rate_per_year) * (float(generation_time_days) / 365.25)


def generation_rate_to_annual_rate(rate_per_generation, generation_time_days):
    return float(rate_per_generation) * (365.25 / float(generation_time_days))


def get_reassortment_rate_units(tree_time_unit):
    if tree_time_unit == TREE_TIME_UNIT_DIVERGENCE:
        return (
            "per_boundary_per_substitution_site",
            "per boundary per substitution/site",
        )
    if tree_time_unit == TREE_TIME_UNIT_TIME:
        return (
            "per_boundary_per_year",
            "per boundary per year",
        )
    if tree_time_unit == TREE_TIME_UNIT_GENERATIONS:
        return (
            "per_boundary_per_generation",
            "per boundary per generation",
        )
    raise ValueError(f"Unsupported tree time unit '{tree_time_unit}'")


def annualize_reassortment_rate(rate, tree_time_unit, generation_time_days, clock_rate):
    if tree_time_unit == TREE_TIME_UNIT_TIME:
        return float(rate)
    if tree_time_unit == TREE_TIME_UNIT_GENERATIONS:
        return generation_rate_to_annual_rate(rate, generation_time_days)
    if tree_time_unit == TREE_TIME_UNIT_DIVERGENCE and clock_rate is not None:
        return float(rate) * float(clock_rate)
    return None


def get_segment_lengths(seq_files):
    from Bio import SeqIO

    return [
        len(list(SeqIO.parse(seq_file, "fasta"))[0])
        for seq_file in seq_files
    ]


def get_boundary_positions(segment_lengths):
    cumulative = 0
    boundary_positions = []
    for length in segment_lengths[:-1]:
        cumulative += int(length)
        boundary_positions.append(cumulative)
    return boundary_positions


def resolve_rate_configuration(args):
    if args.analysis_mode == ANALYSIS_MODE_REASSORTMENT:
        initial_rate = args.rec_rate
        if initial_rate is None:
            if args.tree_time_unit == TREE_TIME_UNIT_GENERATIONS:
                initial_rate = annual_rate_to_generation_rate(
                    DEFAULT_REASSORTMENT_REC_RATE_PER_YEAR,
                    args.generation_time_days,
                )
            elif args.tree_time_unit == TREE_TIME_UNIT_TIME:
                initial_rate = DEFAULT_REASSORTMENT_REC_RATE_PER_YEAR
            else:
                initial_rate = DEFAULT_REASSORTMENT_REC_RATE_PER_DIVERGENCE
        lower = args.rec_rate_lower if args.rec_rate_lower is not None else 0.0
        upper = args.rec_rate_upper
        if upper is None:
            if args.tree_time_unit == TREE_TIME_UNIT_GENERATIONS:
                upper = annual_rate_to_generation_rate(
                    DEFAULT_REASSORTMENT_RATE_UPPER_PER_YEAR,
                    args.generation_time_days,
                )
            elif args.tree_time_unit == TREE_TIME_UNIT_TIME:
                upper = DEFAULT_REASSORTMENT_RATE_UPPER_PER_YEAR
            else:
                upper = DEFAULT_REASSORTMENT_RATE_UPPER_PER_DIVERGENCE
        return float(initial_rate), (float(lower), float(upper))

    initial_rate = args.rec_rate if args.rec_rate is not None else DEFAULT_RECOMBINATION_REC_RATE
    lower = args.rec_rate_lower if args.rec_rate_lower is not None else 0.0
    upper = args.rec_rate_upper if args.rec_rate_upper is not None else 0.01
    return float(initial_rate), (float(lower), float(upper))


def write_inference_diagnostics(output_dir, diagnostics):
    path = os.path.join(output_dir, "em_diagnostics.json")
    with open(path, "w") as handle:
        json.dump(diagnostics, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return path


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


def run_em_estimation(
    argb,
    ml_tree_files,
    seq_files,
    ref,
    scar_model,
    output_dir,
    em_iters=10,
    use_modern_converter=False,
    generation_time_days=3.0,
):
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
        rec_rate_est: Estimated recombination/reassortment rate in model-native units
        n_recomb: Number of inferred recombination events
        diagnostics: Estimation diagnostics and provenance
    """
    logger.info("=" * 60)
    logger.info("Step 3b: EM-based ARG reconstruction & rate estimation")
    logger.info("=" * 60)
    logger.info(f"  Initial recombination rate: {scar_model.rec_rate:.2e}")
    logger.info(f"  Running {em_iters} EM iterations...")
    if use_modern_converter:
        if MODERN_CONVERTER_AVAILABLE:
            logger.info("  Using modern TreeSequence converter")
            ts_converter = modern_convert
            ts_converter_name = "modern"
        else:
            logger.warning("  Modern converter unavailable; using legacy converter")
            ts_converter = None
            ts_converter_name = None
    else:
        logger.info("  Using legacy TreeSequence converter")
        from Espalier import Dendro2TSConverter
        ts_converter = Dendro2TSConverter.convert
        ts_converter_name = "legacy"

    ts, rec_rate_est, n_recomb, tree_path_with_recs, diagnostics = argb.run_EM(
        ml_tree_files,
        seq_files,
        ref,
        scar_model,
        iters=em_iters,
        return_tree_path=True,
        return_diagnostics=True,
        ts_converter=ts_converter,
        ts_converter_name=ts_converter_name,
    )

    logger.info(
        f"  Final estimated rate: {rec_rate_est:.6f} {scar_model.rate_display_units}"
    )
    if scar_model.rate_units == "per_boundary_per_generation":
        annual_rate = generation_rate_to_annual_rate(rec_rate_est, generation_time_days)
        logger.info(
            f"  Annualized reassortment rate: {annual_rate:.6f} per boundary per year "
            f"(generation time {generation_time_days:.3f} days)"
        )
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

    return ts, rec_rate_est, n_recomb, diagnostics


def apply_rate_unit_metadata(scar_model, tree_time_unit):
    if getattr(scar_model, "rate_model_name", None) != "boundary_hotspot":
        return
    rate_units, rate_display_units = get_reassortment_rate_units(tree_time_unit)
    scar_model.rate_units = rate_units
    scar_model.rate_display_units = rate_display_units


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

    node_flags = ts.tables.nodes.flags
    node_times = ts.tables.nodes.time

    recomb_summary = summarize_recombination_events(ts)
    recomb_node_ids = recomb_summary["recomb_node_ids"]
    n_recomb_events = recomb_summary["n_events"]
    unpaired_recomb_node_ids = recomb_summary["unpaired_recomb_node_ids"]

    logger.info(
        "  Recombination nodes found: %s; structural events: %s; unpaired nodes: %s",
        recomb_summary["n_recomb_nodes"],
        n_recomb_events,
        len(unpaired_recomb_node_ids),
    )

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

        f.write(f"Total recombination events: {n_recomb_events}\n")
        f.write(f"Recombinant flagged nodes: {recomb_summary['n_recomb_nodes']}\n")
        f.write(f"Unpaired recombinant nodes: {unpaired_recomb_node_ids}\n")
        f.write(f"Breakpoints: {breakpoints}\n\n")

        f.write("Structural Recombination Events:\n")
        f.write("-" * 60 + "\n")
        for event in recomb_summary["events"]:
            parent_text = ", ".join(str(parent) for parent in event["parents"])
            f.write(f"  Child {event['child']}: parents=[{parent_text}]\n")
        f.write("\n")

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
        'n_recomb_events': n_recomb_events,
        'n_recomb_nodes': recomb_summary["n_recomb_nodes"],
        'unpaired_recomb_nodes': unpaired_recomb_node_ids,
        'unpaired_recomb_node_count': len(unpaired_recomb_node_ids),
        'structural_recombination_events': recomb_summary["events"],
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
    parser.add_argument('--rec-rate', type=float, default=None,
                        help='Initial rate in model-native units. Recombination mode uses per site per generation; reassortment mode uses per boundary per generation.')
    parser.add_argument('--rec-rate-lower', type=float, default=None,
                        help='Lower optimization bound in model-native units (defaults depend on analysis mode).')
    parser.add_argument('--rec-rate-upper', type=float, default=None,
                        help='Upper optimization bound in model-native units (defaults depend on analysis mode).')
    parser.add_argument('--em-iters', type=int, default=5,
                        help='Number of EM iterations (default: 5)')
    parser.add_argument('--ne', type=float, default=1.0,
                        help='Effective population size in coalescent units (default: 1.0)')
    parser.add_argument('--raxml-path', type=str, default='raxml-ng',
                        help='Path to raxml-ng executable (default: raxml-ng)')
    parser.add_argument('--skip-em', action='store_true',
                        help='Skip EM estimation, only run basic ARG reconstruction')
    parser.add_argument('--use-modern-converter', action='store_true', default=True,
                        help='Use the modernized TreeSequence converter (default)')
    parser.add_argument('--use-legacy-converter', action='store_false', dest='use_modern_converter',
                        help='Use the legacy TreeSequence converter')
    parser.add_argument('--analysis-mode', choices=[ANALYSIS_MODE_RECOMBINATION, ANALYSIS_MODE_REASSORTMENT],
                        default=ANALYSIS_MODE_RECOMBINATION,
                        help='Choose exact-label recombination analysis or isolate-normalized reassortment mode')
    parser.add_argument('--isolate-field', type=int, default=0,
                        help='Pipe-delimited FASTA field used as the isolate key in reassortment mode (default: 0)')
    parser.add_argument('--tree-time-unit',
                        choices=[
                            TREE_TIME_UNIT_DIVERGENCE,
                            TREE_TIME_UNIT_TIME,
                            TREE_TIME_UNIT_GENERATIONS,
                        ],
                        default=TREE_TIME_UNIT_DIVERGENCE,
                        help='Units of input tree branch lengths for reassortment rate reporting. '
                             'Use divergence for substitutions/site, time for calendar years, '
                             'or generations for coalescent-generation trees.')
    parser.add_argument('--clock-rate', type=float, default=None,
                        help='Clock rate in substitutions/site/year. Only used to annualize divergence-unit trees.')
    parser.add_argument('--generation-time-days', type=float, default=3.0,
                        help='Generation time used when annualizing generation-unit reassortment rates (default: 3.0)')

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
    INITIAL_REC_RATE, RATE_BOUNDS = resolve_rate_configuration(args)
    EM_ITERATIONS = args.em_iters
    NE = args.ne
    RAXML_PATH = args.raxml_path
    ANALYSIS_MODE = args.analysis_mode

    logger.info("=" * 60)
    logger.info("Espalier ARG Analysis Pipeline")
    logger.info("=" * 60)
    logger.info(f"Analysis mode: {ANALYSIS_MODE}")
    logger.info(f"Initial rate: {INITIAL_REC_RATE:.6g}")
    logger.info(f"Rate bounds: ({RATE_BOUNDS[0]:.6g}, {RATE_BOUNDS[1]:.6g})")
    if ANALYSIS_MODE == ANALYSIS_MODE_REASSORTMENT:
        logger.info(f"Tree time unit: {args.tree_time_unit}")
        if args.tree_time_unit == TREE_TIME_UNIT_DIVERGENCE and args.clock_rate is not None:
            logger.info(f"Clock rate for annualization: {args.clock_rate:.6g} substitutions/site/year")

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
    segment_lengths = get_segment_lengths(seq_files)
    genome_length = sum(segment_lengths)
    boundary_positions = get_boundary_positions(segment_lengths)
    logger.info(f"\nTotal genome length: {genome_length} bp")
    if ANALYSIS_MODE == ANALYSIS_MODE_REASSORTMENT:
        logger.info(f"Segment boundary positions: {boundary_positions}")

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
        em_diagnostics = {
            "estimate_source": "skip_em",
            "estimate_status": "skipped",
            "rate_model_name": "unestimated",
            "rate_units": "unestimated",
            "rate_display_units": "unestimated",
            "opt_bounds": list(RATE_BOUNDS),
        }
    else:
        logger.info("\nInitializing SCAR coalescent model...")
        if ANALYSIS_MODE == ANALYSIS_MODE_REASSORTMENT:
            scar_model = BoundarySCAR(
                rec_rate=INITIAL_REC_RATE,
                M=[[0]],  # Single population (no migration)
                Ne=NE,
                genome_length=genome_length,
                boundary_positions=boundary_positions,
                bounds=RATE_BOUNDS,
            )
            apply_rate_unit_metadata(scar_model, args.tree_time_unit)
        else:
            scar_model = SCAR(
                rec_rate=INITIAL_REC_RATE,
                M=[[0]],  # Single population (no migration)
                Ne=NE,
                genome_length=genome_length,
                bounds=RATE_BOUNDS,
            )

        try:
            ts, rec_rate_est, n_recomb, em_diagnostics = run_em_estimation(
                argb, ml_tree_files, seq_files, ref, scar_model, output_dir, EM_ITERATIONS,
                use_modern_converter=args.use_modern_converter,
                generation_time_days=args.generation_time_days,
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
            em_diagnostics = {
                "estimate_source": "em_failure",
                "estimate_status": "failed",
                "rate_model_name": getattr(scar_model, "rate_model_name", "unknown"),
                "rate_units": getattr(scar_model, "rate_units", "unknown"),
                "rate_display_units": getattr(scar_model, "rate_display_units", "unknown"),
                "opt_bounds": list(RATE_BOUNDS),
            }

    # Final summary
    logger.info("\n" + "=" * 60)
    logger.info("ANALYSIS COMPLETE")
    logger.info("=" * 60)
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Number of segments: {len(seq_files)}")
    if ANALYSIS_MODE == ANALYSIS_MODE_REASSORTMENT:
        rate_units, rate_display_units = get_reassortment_rate_units(args.tree_time_unit)
        annual_rate = annualize_reassortment_rate(
            rec_rate_est,
            args.tree_time_unit,
            args.generation_time_days,
            args.clock_rate,
        )
        logger.info(f"Estimated reassortment rate: {rec_rate_est:.6f} {rate_display_units}")
        if annual_rate is not None:
            logger.info(f"Estimated reassortment rate (annualized): {annual_rate:.6f} per boundary per year")
        else:
            logger.info("Estimated reassortment rate was not annualized; provide --clock-rate for divergence trees")
    else:
        logger.info(f"Estimated recombination rate: {rec_rate_est:.6f} per site")
    if ts:
        logger.info(f"Recombination events detected: {int(n_recomb)}")
        logger.info(f"TreeSequence file: {os.path.join(output_dir, 'arg_treesequence.trees')}")

    annual_rate = (
        annualize_reassortment_rate(
            rec_rate_est,
            args.tree_time_unit,
            args.generation_time_days,
            args.clock_rate,
        )
        if ANALYSIS_MODE == ANALYSIS_MODE_REASSORTMENT
        else None
    )
    rate_units, rate_display_units = (
        get_reassortment_rate_units(args.tree_time_unit)
        if ANALYSIS_MODE == ANALYSIS_MODE_REASSORTMENT
        else ("per_site_per_generation", "per site per generation")
    )
    em_diagnostics.update(
        {
            "analysis_mode": ANALYSIS_MODE,
            "generation_time_days": args.generation_time_days,
            "tree_time_unit": args.tree_time_unit if ANALYSIS_MODE == ANALYSIS_MODE_REASSORTMENT else None,
            "clock_rate_subs_per_site_per_year": args.clock_rate,
            "initial_rate": INITIAL_REC_RATE,
            "rate_bounds": list(RATE_BOUNDS),
            "estimated_rate": rec_rate_est,
            "estimated_rate_units": rate_units,
            "estimated_rate_display_units": rate_display_units,
            "estimated_rate_per_year": annual_rate,
            "segment_lengths": segment_lengths,
            "boundary_positions": boundary_positions,
            "n_recomb_events": int(n_recomb),
            "has_treesequence": bool(ts),
        }
    )
    diagnostics_path = write_inference_diagnostics(output_dir, em_diagnostics)
    logger.info(f"Inference diagnostics: {diagnostics_path}")

    # Cleanup temp directory (optional)
    # shutil.rmtree(temp_dir)

    return output_dir


if __name__ == '__main__':
    output_dir = main()
