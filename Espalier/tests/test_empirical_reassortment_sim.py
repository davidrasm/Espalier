import json
import subprocess
import sys
from pathlib import Path

import dendropy
import pytest
from Bio import SeqIO
from Espalier import MAF

pytest.importorskip("pyvolve")

from Espalier.EmpiricalReassortmentSim import (
    generate_empirical_reassortment_panel,
    simulate_empirical_reassortment,
    write_empirical_reassortment_dataset,
)


def _write_tree(path: Path, newick: str) -> None:
    path.write_text(newick)


def _balanced_newick(labels) -> str:
    if len(labels) == 1:
        return f"{labels[0]}:1"
    midpoint = len(labels) // 2
    left = _balanced_newick(labels[:midpoint])
    right = _balanced_newick(labels[midpoint:])
    return f"({left},{right}):1"


def _leaf_labels_from_tree(path: Path):
    tree = dendropy.Tree.get(
        path=str(path),
        schema="newick",
        rooting="default-rooted",
        preserve_underscores=True,
    )
    return sorted(leaf.taxon.label for leaf in tree.leaf_node_iter())


def _spr_distance_from_newicks(reference_newick: str, reassorted_newick: str) -> int:
    taxon_namespace = dendropy.TaxonNamespace()
    reference_tree = dendropy.Tree.get(
        data=reference_newick,
        schema="newick",
        rooting="default-rooted",
        preserve_underscores=True,
        taxon_namespace=taxon_namespace,
    )
    reassorted_tree = dendropy.Tree.get(
        data=reassorted_newick,
        schema="newick",
        rooting="default-rooted",
        preserve_underscores=True,
        taxon_namespace=taxon_namespace,
    )
    return int(MAF.get_spr_dist(reference_tree, reassorted_tree))


def test_simulate_empirical_reassortment_records_event_constraints(tmp_path):
    tree_path = tmp_path / "empirical.tre"
    _write_tree(
        tree_path,
        "((((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1):1,((I:1,J:1):1,(K:1,L:1):1):2);",
    )

    dataset = simulate_empirical_reassortment(
        str(tree_path),
        segment_lengths=[20, 20],
        num_reassortments=1,
        min_clade_size=2,
        max_clade_size=3,
        min_normalized_distance=0.2,
        max_normalized_distance=1.0,
        random_seed=7,
    )

    assert dataset.mode == "empirical_reassortment"
    assert len(dataset.events) == 1
    assert dataset.selected_boundary_positions == [20]
    assert dataset.reference_tree_newick != dataset.reassorted_tree_newick

    event = dataset.events[0]
    assert 2 <= event.clade_size <= 3
    assert 0.2 <= event.normalized_distance <= 1.0
    assert len(event.moved_taxa) == event.clade_size
    assert _spr_distance_from_newicks(dataset.reference_tree_newick, dataset.reassorted_tree_newick) == 1


def test_empirical_reassortment_events_are_non_overlapping(tmp_path):
    tree_path = tmp_path / "empirical_multi.tre"
    _write_tree(
        tree_path,
        "(((((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1):1,(((I:1,J:1):1,(K:1,L:1):1):1,((M:1,N:1):1,(O:1,P:1):1):1):1):1,(Q:1,R:1):4);",
    )

    dataset = simulate_empirical_reassortment(
        str(tree_path),
        segment_lengths=[18, 18],
        num_reassortments=2,
        min_clade_size=2,
        max_clade_size=2,
        min_normalized_distance=0.0,
        max_normalized_distance=1.0,
        random_seed=11,
    )

    moved = [set(event.moved_taxa) for event in dataset.events]
    assert len(moved) == 2
    assert moved[0].isdisjoint(moved[1])
    assert _spr_distance_from_newicks(dataset.reference_tree_newick, dataset.reassorted_tree_newick) == 2


def test_empirical_reassortment_supports_multiple_non_overlapping_events(tmp_path):
    tree_path = tmp_path / "empirical_k3.tre"
    labels = [f"T{idx:02d}" for idx in range(1, 33)]
    _write_tree(tree_path, _balanced_newick(labels) + ";")

    dataset = simulate_empirical_reassortment(
        str(tree_path),
        segment_lengths=[18, 18],
        num_reassortments=3,
        min_clade_size=2,
        max_clade_size=2,
        min_normalized_distance=0.1,
        max_normalized_distance=0.8,
        random_seed=19,
    )

    assert len(dataset.events) == 3
    moved_taxa = [set(event.moved_taxa) for event in dataset.events]
    for idx, taxa_i in enumerate(moved_taxa):
        for jdx, taxa_j in enumerate(moved_taxa):
            if idx == jdx:
                continue
            assert taxa_i.isdisjoint(taxa_j)
    assert len({tuple(event.target_taxa) for event in dataset.events}) == 3
    assert _spr_distance_from_newicks(dataset.reference_tree_newick, dataset.reassorted_tree_newick) == 3


def test_empirical_reassortment_supports_five_non_overlapping_events(tmp_path):
    tree_path = tmp_path / "empirical_k5.tre"
    labels = [f"T{idx:02d}" for idx in range(1, 65)]
    _write_tree(tree_path, _balanced_newick(labels) + ";")

    dataset = simulate_empirical_reassortment(
        str(tree_path),
        segment_lengths=[18, 18],
        num_reassortments=5,
        min_clade_size=2,
        max_clade_size=2,
        min_normalized_distance=0.1,
        max_normalized_distance=0.8,
        random_seed=23,
    )

    assert len(dataset.events) == 5
    moved_taxa = [set(event.moved_taxa) for event in dataset.events]
    for idx, taxa_i in enumerate(moved_taxa):
        for jdx, taxa_j in enumerate(moved_taxa):
            if idx == jdx:
                continue
            assert taxa_i.isdisjoint(taxa_j)
    assert len({tuple(event.target_taxa) for event in dataset.events}) == 5
    assert _spr_distance_from_newicks(dataset.reference_tree_newick, dataset.reassorted_tree_newick) == 5


def test_write_empirical_reassortment_dataset_creates_files(tmp_path):
    tree_path = tmp_path / "empirical_write.tre"
    _write_tree(
        tree_path,
        "((((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1):1,((I:1,J:1):1,(K:1,L:1):1):2);",
    )

    dataset = simulate_empirical_reassortment(
        str(tree_path),
        segment_lengths=[12, 14, 16],
        reassorted_segment_indexes=[1],
        num_reassortments=1,
        min_clade_size=2,
        max_clade_size=2,
        random_seed=5,
    )

    written = write_empirical_reassortment_dataset(dataset, str(tmp_path / "dataset"))
    truth = json.loads(Path(written["truth_path"]).read_text())

    assert truth["selected_boundary_positions"] == [12, 26]
    assert len(list((tmp_path / "dataset" / "alignments").glob("*.fasta"))) == 3
    assert len(list((tmp_path / "dataset" / "trees").glob("*.tre"))) == 3

    expected_lengths = {
        "segment01.fasta": 12,
        "segment02.fasta": 14,
        "segment03.fasta": 16,
    }
    for fasta_path in (tmp_path / "dataset" / "alignments").glob("*.fasta"):
        records = list(SeqIO.parse(str(fasta_path), "fasta"))
        assert records
        assert len(records[0].seq) == expected_lengths[fasta_path.name]

    tree_dir = tmp_path / "dataset" / "trees"
    seg1_tree = tree_dir / "segment01.fasta.rooted.tre"
    seg2_tree = tree_dir / "segment02.fasta.rooted.tre"
    seg3_tree = tree_dir / "segment03.fasta.rooted.tre"
    assert _leaf_labels_from_tree(seg1_tree) == _leaf_labels_from_tree(seg2_tree)
    assert seg1_tree.read_text() == seg3_tree.read_text()
    assert seg1_tree.read_text() != seg2_tree.read_text()


def test_generate_empirical_panel_creates_requested_replicates(tmp_path):
    tree_path = tmp_path / "empirical_panel.tre"
    _write_tree(
        tree_path,
        "((((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1):1,((I:1,J:1):1,(K:1,L:1):1):2);",
    )

    written = generate_empirical_reassortment_panel(
        str(tree_path),
        str(tmp_path / "panel"),
        repeats=2,
        start_index=4,
        segment_lengths=[10, 10],
        num_reassortments=1,
        min_clade_size=2,
        max_clade_size=2,
        random_seed=3,
    )

    assert [Path(record["output_dir"]).name for record in written] == ["sample4", "sample5"]


def test_empirical_benchmark_script_generates_summary(tmp_path):
    tree_path = tmp_path / "empirical_benchmark.tre"
    _write_tree(
        tree_path,
        "((((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1):1,((I:1,J:1):1,(K:1,L:1):1):2);",
    )

    output_dir = tmp_path / "benchmark"
    subprocess.run(
        [
            sys.executable,
            str(Path(__file__).resolve().parents[2] / "scripts" / "run_empirical_reassortment_benchmark.py"),
            "--tree",
            str(tree_path),
            "--output",
            str(output_dir),
            "--segment-lengths",
            "10,12",
            "--num-reassortments",
            "1",
            "--repeats",
            "1",
            "--min-size",
            "2",
            "--max-size",
            "2",
            "--skip-analysis",
        ],
        check=True,
    )

    summary = json.loads((output_dir / "benchmark_summary.json").read_text())
    assert summary["n_runs"] == 1
    assert summary["n_analyzed"] == 0
    sample_dir = output_dir / "sample1"
    assert (sample_dir / "truth.json").exists()
    assert (sample_dir / "result.json").exists()
