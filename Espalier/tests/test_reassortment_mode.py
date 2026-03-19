import csv
from pathlib import Path
from types import SimpleNamespace

import dendropy
import numpy as np
import pytest
from Bio import SeqIO

import Espalier.ARGBuilder as ARGBuilderModule
from Espalier import Utils
from Espalier.ARGBuilder import ARGBuilder
from Espalier.Dendro2TSConverter import convert
from Espalier.Reassortment import (
    ANALYSIS_MODE_REASSORTMENT,
    load_segment_data,
    prepare_analysis_inputs,
    select_isolates_by_week,
)


def _write_fasta(path: Path, records):
    with open(path, "w") as handle:
        for header, seq in records:
            handle.write(f">{header}\n{seq}\n")


def _write_tree(path: Path, labels):
    if len(labels) == 2:
        newick = f"({labels[0]}:1,{labels[1]}:1);"
    elif len(labels) == 3:
        newick = f"(({labels[0]}:1,{labels[1]}:1):1,{labels[2]}:2);"
    else:
        raise ValueError("Only 2-tip and 3-tip test trees are supported")
    path.write_text(newick)


def _read_fasta_records(path: Path):
    return {record.id: str(record.seq) for record in SeqIO.parse(path, "fasta")}


def test_prepare_analysis_inputs_normalizes_reassortment_inputs(tmp_path):
    align_dir = tmp_path / "alignments"
    tree_dir = tmp_path / "trees"
    output_dir = tmp_path / "output"
    align_dir.mkdir()
    tree_dir.mkdir()

    seg1 = align_dir / "segment1.fasta"
    seg2 = align_dir / "segment2.fasta"
    _write_fasta(
        seg1,
        [
            ("Iso A|HA|2020-01-01", "AAAA"),
            ("Iso B|HA|2020-01-08", "CCCC"),
            ("Iso C|HA|2020-01-15", "GGGG"),
        ],
    )
    _write_fasta(
        seg2,
        [
            ("Iso A|NA|2020-01-01", "TTTT"),
            ("Iso B|NA|2020-01-08", "GGGG"),
        ],
    )

    tree1 = tree_dir / "segment1.tre"
    tree2 = tree_dir / "segment2.tre"
    _write_tree(
        tree1,
        [
            "Iso_A|HA|2020-01-01",
            "Iso_B|HA|2020-01-08",
            "Iso_C|HA|2020-01-15",
        ],
    )
    _write_tree(
        tree2,
        [
            "Iso_A|NA|2020-01-01",
            "Iso_B|NA|2020-01-08",
        ],
    )

    prepared = prepare_analysis_inputs(
        [str(seg1), str(seg2)],
        [str(tree1), str(tree2)],
        str(output_dir),
        analysis_mode=ANALYSIS_MODE_REASSORTMENT,
    )

    assert prepared.common_isolates == ["Iso-A", "Iso-B"]
    assert Path(prepared.label_map_path).exists()

    for normalized_seq in prepared.seq_files:
        assert _read_fasta_records(Path(normalized_seq)).keys() == {"Iso-A", "Iso-B"}

    for normalized_tree in prepared.tree_files:
        tree = dendropy.Tree.get(
            path=normalized_tree,
            schema="newick",
            preserve_underscores=True,
        )
        assert {leaf.taxon.label for leaf in tree.leaf_node_iter()} == {"Iso-A", "Iso-B"}

    with open(prepared.label_map_path) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert any(
        row["original_full_label"] == "Iso C|HA|2020-01-15" and row["selected"] == "0"
        for row in rows
    )


def test_load_segment_data_rejects_duplicate_isolate_keys(tmp_path):
    seq_file = tmp_path / "segment.fasta"
    _write_fasta(
        seq_file,
        [
            ("Iso A|HA|2020-01-01", "AAAA"),
            ("Iso A|HA-dup|2020-01-01", "CCCC"),
        ],
    )

    with pytest.raises(ValueError, match="Duplicate isolate key"):
        load_segment_data(str(seq_file))


def test_prepare_analysis_inputs_rejects_tree_alignment_mismatch(tmp_path):
    seq_file = tmp_path / "segment.fasta"
    tree_file = tmp_path / "segment.tre"
    _write_fasta(
        seq_file,
        [
            ("Iso A|HA|2020-01-01", "AAAA"),
            ("Iso B|HA|2020-01-08", "CCCC"),
        ],
    )
    _write_tree(
        tree_file,
        [
            "Iso_A|HA|2020-01-01",
            "Iso_C|HA|2020-01-15",
        ],
    )

    with pytest.raises(ValueError, match="Tree/alignment mismatch"):
        prepare_analysis_inputs(
            [str(seq_file)],
            [str(tree_file)],
            str(tmp_path / "output"),
            analysis_mode=ANALYSIS_MODE_REASSORTMENT,
        )


def test_concate_aligns_matches_on_record_ids(tmp_path):
    seg1 = tmp_path / "seg1.fasta"
    seg2 = tmp_path / "seg2.fasta"
    concat = tmp_path / "concat.fasta"
    _write_fasta(
        seg1,
        [
            ("Alpha segmentA", "AAAA"),
            ("Beta segmentA", "CCCC"),
        ],
    )
    _write_fasta(
        seg2,
        [
            ("Alpha segmentB", "GG"),
            ("Beta segmentB", "TT"),
        ],
    )

    Utils.concate_aligns([str(seg1), str(seg2)], str(concat))

    assert _read_fasta_records(concat) == {
        "Alpha": "AAAAGG",
        "Beta": "CCCCTT",
    }


def test_concate_aligns_reports_missing_taxa(tmp_path):
    seg1 = tmp_path / "seg1.fasta"
    seg2 = tmp_path / "seg2.fasta"
    concat = tmp_path / "concat.fasta"
    _write_fasta(
        seg1,
        [
            ("Alpha segmentA", "AAAA"),
            ("Beta segmentA", "CCCC"),
        ],
    )
    _write_fasta(
        seg2,
        [
            ("Alpha segmentB", "GG"),
        ],
    )

    with pytest.raises(ValueError, match="missing taxa"):
        Utils.concate_aligns([str(seg1), str(seg2)], str(concat))


def test_convert_handles_mixed_numeric_and_string_tip_labels():
    taxa = dendropy.TaxonNamespace()
    tree_data = "((1:1,A:1):1,2:2);"
    tree1 = dendropy.Tree.get(
        data=tree_data,
        schema="newick",
        rooting="default-rooted",
        taxon_namespace=taxa,
        preserve_underscores=True,
    )
    tree2 = dendropy.Tree.get(
        data=tree_data,
        schema="newick",
        rooting="default-rooted",
        taxon_namespace=taxa,
        preserve_underscores=True,
    )

    ts = convert([tree1, tree2], [(0, 10), (10, 20)])

    assert ts.num_nodes > 0
    assert ts.sequence_length == 20


def test_select_isolates_by_week_is_shared_and_order_independent(tmp_path):
    seg1 = tmp_path / "segment1.fasta"
    seg2 = tmp_path / "segment2.fasta"
    _write_fasta(
        seg1,
        [
            ("Iso1|HA|2020-01-01", "AAAA"),
            ("Iso2|HA|2020-01-03", "CCCC"),
            ("Iso3|HA|2020-01-08", "GGGG"),
            ("Iso4|HA|2020-01-10", "TTTT"),
        ],
    )
    _write_fasta(
        seg2,
        [
            ("Iso2|NA|2020-01-03", "CCCC"),
            ("Iso1|NA|2020-01-01", "AAAA"),
            ("Iso4|NA|2020-01-10", "TTTT"),
            ("Iso3|NA|2020-01-08", "GGGG"),
        ],
    )

    data1 = load_segment_data(str(seg1))
    data2 = load_segment_data(str(seg2))

    selected = select_isolates_by_week([data1, data2], seed=7)
    selected_swapped = select_isolates_by_week([data2, data1], seed=7)

    assert selected == selected_swapped
    assert len(selected) == 2
    assert set(selected).issubset({"Iso1", "Iso2", "Iso3", "Iso4"})


class _FakeCoalModel:
    def __init__(self, rec_rate=0.1, genome_length=5):
        self.rec_rate = rec_rate
        self.genome_length = genome_length

    def opt_MLE(self, ts):
        return 0.25


def test_run_em_return_contract(monkeypatch):
    ref = dendropy.Tree.get(
        data="(A:1,B:1);",
        schema="newick",
        rooting="default-rooted",
        preserve_underscores=True,
    )
    builder = ARGBuilder(reconciler=object(), raxml=object(), temp_dir="/tmp/")
    fake_ts = SimpleNamespace(
        tables=SimpleNamespace(
            nodes=SimpleNamespace(flags=np.array([], dtype=np.int64))
        )
    )

    monkeypatch.setattr(
        ARGBuilder,
        "_get_genome_segments",
        lambda self, seq_files: setattr(
            self,
            "segments",
            [SimpleNamespace(start=0, end=5)],
        ),
    )
    monkeypatch.setattr(
        ARGBuilder,
        "_build_trellis",
        lambda self, local_tree_files, seq_files, ref: ([[]], [[]]),
    )
    monkeypatch.setattr(
        ARGBuilder,
        "_get_like_array",
        lambda self, seq_files, tree_trellis: np.zeros((1, 1)),
    )
    monkeypatch.setattr(
        ARGBuilder,
        "_get_recomb_trans_probs",
        lambda self, tree_trellis, rec_rate, rSPR_array=None: np.zeros((1, 1, 1)),
    )
    monkeypatch.setattr(ARGBuilderModule, "_get_rSPR_array", lambda tree_trellis: np.zeros((1, 1, 1)))
    monkeypatch.setattr(
        ARGBuilderModule,
        "viterbi",
        lambda tree_trellis, trans_probs, like_array: ([ref.clone(depth=2)], [0]),
    )
    monkeypatch.setattr(
        ARGBuilderModule.Reconciler,
        "reconcile_linked_heights",
        lambda tree_path: tree_path,
    )
    monkeypatch.setattr(ARGBuilderModule, "_jitter_coal_times", lambda tree_path, displace_dt=0.0001: tree_path)
    monkeypatch.setattr(ARGBuilderModule, "add_path_rec_nodes", lambda tree_path: (tree_path, 0))
    monkeypatch.setattr(ARGBuilderModule.Dendro2TSConverter, "convert", lambda tree_path, tree_intervals: fake_ts)

    default_result = builder.run_EM(["tree"], ["seq"], ref, _FakeCoalModel(), iters=1)
    tree_path_result = builder.run_EM(
        ["tree"],
        ["seq"],
        ref,
        _FakeCoalModel(),
        iters=1,
        return_tree_path=True,
    )

    assert len(default_result) == 3
    assert len(tree_path_result) == 4
    assert tree_path_result[3] is not None
    assert len(tree_path_result[3]) == 1
