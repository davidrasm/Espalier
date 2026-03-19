#!/usr/bin/env python3
"""
Subsample reassortment-ready segment alignments by date, keeping one isolate per week.

Usage:
    python subsample_by_date.py --input-dir TestFiles/FluTest --output-dir TestFiles/FluTest_subsampled
"""

import argparse
import os
import sys

import dendropy
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Espalier.Reassortment import (
    find_matching_tree,
    load_segment_data,
    resolve_tree_identity,
    select_isolates_by_week,
    validate_tree_matches_segment,
)


def write_fasta(records, filepath):
    SeqIO.write(records, filepath, "fasta")


def prune_tree_to_isolates(segment_data, keep_isolates, output_file):
    tree = dendropy.Tree.get(
        path=segment_data.tree_file,
        schema="newick",
        preserve_underscores=True,
    )
    keep_set = set(keep_isolates)
    taxa_to_remove = []
    for leaf in tree.leaf_node_iter():
        identity = resolve_tree_identity(segment_data, leaf.taxon.label)
        if identity.analysis_label not in keep_set:
            taxa_to_remove.append(leaf.taxon)
    original_taxa = len(tree.taxon_namespace)
    if taxa_to_remove:
        tree.prune_taxa(taxa_to_remove)
    tree.suppress_unifurcations()
    tree.write(path=output_file, schema="newick", suppress_annotations=True)
    print(f"  Original taxa: {original_taxa}")
    print(f"  Pruned taxa: {len(tree.taxon_namespace)}")


def main():
    parser = argparse.ArgumentParser(
        description="Subsample reassortment alignments by date, keeping one isolate per ISO week"
    )
    parser.add_argument(
        "--input-dir",
        "-i",
        type=str,
        required=True,
        help="Input directory containing alignments/ and trees/ subdirs",
    )
    parser.add_argument(
        "--output-dir",
        "-o",
        type=str,
        required=True,
        help="Output directory for subsampled data",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility (default: 42)",
    )
    parser.add_argument(
        "--isolate-field",
        type=int,
        default=0,
        help="Pipe-delimited FASTA field used as the isolate key (default: 0)",
    )
    args = parser.parse_args()

    input_dir = args.input_dir
    output_dir = args.output_dir
    aln_dir = os.path.join(input_dir, "alignments")
    tree_dir = os.path.join(input_dir, "trees")

    os.makedirs(os.path.join(output_dir, "alignments"), exist_ok=True)
    os.makedirs(os.path.join(output_dir, "trees"), exist_ok=True)

    print("=" * 60)
    print("Subsampling shared isolates by week")
    print("=" * 60)
    print(f"Input: {input_dir}")
    print(f"Output: {output_dir}")
    print()

    aln_files = sorted(
        f for f in os.listdir(aln_dir)
        if f.endswith(".fasta") or f.endswith(".aln")
    )

    segment_data = []
    for aln_file in aln_files:
        print(f"\nProcessing: {aln_file}")
        print("-" * 40)
        seq_file = os.path.join(aln_dir, aln_file)
        tree_file = find_matching_tree(seq_file, tree_dir)
        if tree_file is None:
            raise ValueError(f"No matching tree found for {aln_file}")
        data = load_segment_data(seq_file, isolate_field=args.isolate_field, tree_file=tree_file)
        validate_tree_matches_segment(tree_file, data)
        segment_data.append(data)
        print(f"  Sequences: {len(data.records_by_full_label)}")
        print(f"  Unique isolates: {len(data.records_by_analysis_label)}")

    selected_isolates = select_isolates_by_week(segment_data, seed=args.seed)
    print(f"\n{'=' * 60}")
    print(f"Selected shared isolates: {len(selected_isolates)}")
    print("=" * 60)

    for data in segment_data:
        print(f"\nWriting: {data.segment_name}")
        selected_records = []
        selected_set = set(selected_isolates)
        for full_label, record in data.records_by_full_label.items():
            identity = data.identities_by_full_label[full_label]
            if identity.analysis_label not in selected_set:
                continue
            selected_records.append(
                SeqRecord(
                    record.seq,
                    id=full_label,
                    name=full_label,
                    description=full_label,
                )
            )

        output_file = os.path.join(output_dir, "alignments", data.segment_name)
        write_fasta(selected_records, output_file)
        print(f"  Wrote {len(selected_records)} sequences to {output_file}")

        output_tree = os.path.join(output_dir, "trees", os.path.basename(data.tree_file))
        print(f"  Pruning tree: {os.path.basename(data.tree_file)}")
        prune_tree_to_isolates(data, selected_isolates, output_tree)
        print(f"  Wrote pruned tree to {output_tree}")

    print(f"\n{'=' * 60}")
    print("SUBSAMPLING COMPLETE")
    print("=" * 60)
    print(f"Output directory: {output_dir}")
    print(f"Sequences per segment: {len(selected_isolates)}")

    original_n = len(segment_data[0].records_by_full_label)
    subsampled_n = len(selected_isolates)
    speedup = (original_n / subsampled_n) ** 2
    print(f"\nEstimated speedup: ~{speedup:.0f}x faster than full dataset")


if __name__ == "__main__":
    main()
