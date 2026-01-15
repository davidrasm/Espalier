#!/usr/bin/env python3
"""
Subsample sequences by date, keeping one sequence per week.

Parses sequence headers with format:
    >A/swine/Michigan/A02635726/2021|1B.2.1|1998B|TTTPPT|2021-04-21

Extracts the date from the last field after "|" separator.
Groups sequences by ISO week and keeps one representative per week.

Usage:
    python subsample_by_date.py --input-dir TestFiles/FluTest --output-dir TestFiles/FluTest_subsampled
"""

import os
import sys
import argparse
import random
from datetime import datetime
from collections import defaultdict
from typing import Dict, List, Tuple

def parse_date_from_header(header: str) -> Tuple[str, datetime]:
    """
    Extract date from sequence header.

    Expected format: >name|field1|field2|...|YYYY-MM-DD
    Returns: (header, datetime) or (header, None) if parsing fails
    """
    try:
        # Remove '>' if present
        header = header.lstrip('>')

        # Split by '|' and get last field
        fields = header.split('|')
        date_str = fields[-1].strip()

        # Parse date (YYYY-MM-DD format)
        date = datetime.strptime(date_str, '%Y-%m-%d')
        return header, date
    except (ValueError, IndexError):
        return header, None


def get_iso_week(date: datetime) -> str:
    """Get ISO year-week string for grouping."""
    iso_cal = date.isocalendar()
    return f"{iso_cal[0]}-W{iso_cal[1]:02d}"


def read_fasta(filepath: str) -> Dict[str, str]:
    """Read FASTA file and return dict of header -> sequence."""
    sequences = {}
    current_header = None
    current_seq = []

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = ''.join(current_seq)
                current_header = line[1:]  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)

        # Don't forget last sequence
        if current_header:
            sequences[current_header] = ''.join(current_seq)

    return sequences


def write_fasta(sequences: Dict[str, str], filepath: str):
    """Write sequences to FASTA file."""
    with open(filepath, 'w') as f:
        for header, seq in sequences.items():
            f.write(f">{header}\n")
            # Write sequence in lines of 60 characters
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + '\n')


def subsample_by_week(sequences: Dict[str, str], seed: int = 42) -> Dict[str, str]:
    """
    Subsample sequences, keeping one per ISO week.

    Parameters:
        sequences: Dict of header -> sequence
        seed: Random seed for reproducibility

    Returns:
        Subsampled dict of header -> sequence
    """
    random.seed(seed)

    # Group sequences by week
    weeks: Dict[str, List[str]] = defaultdict(list)
    no_date_seqs = []

    for header in sequences:
        _, date = parse_date_from_header(header)
        if date:
            week = get_iso_week(date)
            weeks[week].append(header)
        else:
            no_date_seqs.append(header)

    # Select one sequence per week
    selected = {}
    for week in sorted(weeks.keys()):
        # Pick one random sequence from this week
        chosen = random.choice(weeks[week])
        selected[chosen] = sequences[chosen]

    print(f"  Original sequences: {len(sequences)}")
    print(f"  Sequences with valid dates: {len(sequences) - len(no_date_seqs)}")
    print(f"  Unique weeks: {len(weeks)}")
    print(f"  Selected sequences: {len(selected)}")
    if no_date_seqs:
        print(f"  Sequences without dates (excluded): {len(no_date_seqs)}")

    return selected


def prune_tree_to_taxa(tree_file: str, taxa_set: set, output_file: str):
    """
    Prune a Newick tree to only include specified taxa.

    Uses dendropy for tree manipulation.
    """
    import dendropy

    # Read tree
    tree = dendropy.Tree.get(path=tree_file, schema='newick', preserve_underscores=True)

    # Get taxa to remove
    taxa_to_remove = []
    for taxon in tree.taxon_namespace:
        if taxon.label not in taxa_set:
            taxa_to_remove.append(taxon)

    # Prune taxa
    tree.prune_taxa(taxa_to_remove)

    # Suppress unifurcations
    tree.suppress_unifurcations()

    # Write pruned tree
    tree.write(path=output_file, schema='newick', suppress_annotations=True)

    print(f"  Original taxa: {len(tree.taxon_namespace) + len(taxa_to_remove)}")
    print(f"  Pruned taxa: {len(tree.taxon_namespace)}")


def main():
    parser = argparse.ArgumentParser(
        description='Subsample sequences by date, keeping one per week'
    )
    parser.add_argument('--input-dir', '-i', type=str, required=True,
                        help='Input directory containing alignments/ and trees/ subdirs')
    parser.add_argument('--output-dir', '-o', type=str, required=True,
                        help='Output directory for subsampled data')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed for reproducibility (default: 42)')

    args = parser.parse_args()

    input_dir = args.input_dir
    output_dir = args.output_dir

    # Create output directories
    os.makedirs(os.path.join(output_dir, 'alignments'), exist_ok=True)
    os.makedirs(os.path.join(output_dir, 'trees'), exist_ok=True)

    print("=" * 60)
    print("Subsampling sequences by week")
    print("=" * 60)
    print(f"Input: {input_dir}")
    print(f"Output: {output_dir}")
    print()

    # Find alignment files
    aln_dir = os.path.join(input_dir, 'alignments')
    tree_dir = os.path.join(input_dir, 'trees')

    aln_files = sorted([f for f in os.listdir(aln_dir)
                        if f.endswith('.fasta') or f.endswith('.aln')])

    # Track common taxa across all segments
    all_selected_taxa = None
    segment_sequences = {}

    # First pass: subsample each segment and find common taxa
    for aln_file in aln_files:
        print(f"\nProcessing: {aln_file}")
        print("-" * 40)

        filepath = os.path.join(aln_dir, aln_file)
        sequences = read_fasta(filepath)

        selected = subsample_by_week(sequences, seed=args.seed)
        segment_sequences[aln_file] = (sequences, selected)

        selected_taxa = set(selected.keys())
        if all_selected_taxa is None:
            all_selected_taxa = selected_taxa
        else:
            all_selected_taxa = all_selected_taxa.intersection(selected_taxa)

    print(f"\n{'=' * 60}")
    print(f"Common taxa across all segments: {len(all_selected_taxa)}")
    print("=" * 60)

    # Second pass: write only common taxa
    for aln_file in aln_files:
        print(f"\nWriting: {aln_file}")

        sequences, selected = segment_sequences[aln_file]

        # Filter to common taxa only
        common_selected = {h: sequences[h] for h in all_selected_taxa if h in sequences}

        output_file = os.path.join(output_dir, 'alignments', aln_file)
        write_fasta(common_selected, output_file)
        print(f"  Wrote {len(common_selected)} sequences to {output_file}")

        # Find and prune corresponding tree
        # Try different naming patterns
        base = aln_file
        tree_patterns = [
            base + '.tre',
            base + '.rooted.tre',
            base.replace('.fasta', '.tre').replace('.aln', '.tre'),
            base.replace('.fasta', '.rooted.tre').replace('.aln', '.rooted.tre'),
        ]

        tree_file = None
        for pattern in tree_patterns:
            candidate = os.path.join(tree_dir, pattern)
            if os.path.exists(candidate):
                tree_file = candidate
                break

        if tree_file:
            output_tree = os.path.join(output_dir, 'trees', os.path.basename(tree_file))
            print(f"\nPruning tree: {os.path.basename(tree_file)}")
            try:
                prune_tree_to_taxa(tree_file, all_selected_taxa, output_tree)
                print(f"  Wrote pruned tree to {output_tree}")
            except Exception as e:
                print(f"  Error pruning tree: {e}")
        else:
            print(f"  No matching tree found for {aln_file}")

    print(f"\n{'=' * 60}")
    print("SUBSAMPLING COMPLETE")
    print("=" * 60)
    print(f"Output directory: {output_dir}")
    print(f"Sequences per segment: {len(all_selected_taxa)}")

    # Estimate time savings
    original_n = len(segment_sequences[aln_files[0]][0])
    subsampled_n = len(all_selected_taxa)
    speedup = (original_n / subsampled_n) ** 2
    print(f"\nEstimated speedup: ~{speedup:.0f}x faster than full dataset")


if __name__ == '__main__':
    main()
