import csv
from pathlib import Path
from types import SimpleNamespace

import dendropy
import numpy as np
import pandas as pd
import pytest
import tskit
from Bio import SeqIO

import Espalier.ARGBuilder as ARGBuilderModule
import scripts.run_arg_analysis as RunArgAnalysis
from Espalier import Dendro2TSConverter
from Espalier import Utils
from Espalier.ARGBuilder import ARGBuilder
from Espalier.ARGNodeTypes import (
    RECOMBINANT_FLAG,
    count_recombination_events,
    is_recombinant,
    summarize_recombination_events,
)
from Espalier.Dendro2TSConverter import convert
from Espalier.Reassortment import (
    ANALYSIS_MODE_REASSORTMENT,
    load_segment_data,
    prepare_analysis_inputs,
    select_isolates_by_week,
)
from Espalier.SCARLikelihood import BoundarySCAR


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


def test_convert_handles_mixed_numeric_and_string_tip_labels(tmp_path):
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
    ts_path = tmp_path / "mixed_labels.trees"
    ts.dump(str(ts_path))
    loaded = tskit.load(str(ts_path))
    assert loaded.sequence_length == ts.sequence_length


def test_reindex_edgedf_updates_parent_child_columns():
    edge_df = pd.DataFrame(
        [
            {
                "left": 0,
                "right": 10,
                "parent": 0,
                "child": 1,
                "parent_unique_id": 100,
                "child_unique_id": 200,
            }
        ]
    )

    reindexed = Dendro2TSConverter.reindex_edgedf(edge_df, {100: 3, 200: 4})

    assert reindexed.loc[0, "parent"] == 3
    assert reindexed.loc[0, "child"] == 4
    assert edge_df.loc[0, "parent"] == 0
    assert edge_df.loc[0, "child"] == 1


def test_tree_tables_to_df_does_not_alias_metadata_to_time():
    tables = tskit.TableCollection(10)
    sample_a = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=-1)
    sample_b = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=-1)
    parent = tables.nodes.add_row(flags=0, time=1, population=-1)
    tables.edges.add_row(0, 10, parent, sample_a)
    tables.edges.add_row(0, 10, parent, sample_b)
    tables.sort()

    _, nodes_df = Dendro2TSConverter.treeTables2df(tables)

    assert nodes_df["metadata"].tolist() == [None, None, None]


def test_recombinant_flag_checks_are_bitwise():
    assert is_recombinant(RECOMBINANT_FLAG)
    assert is_recombinant(RECOMBINANT_FLAG | tskit.NODE_IS_SAMPLE)


def test_boundary_scar_counts_adjacent_segment_boundary_opportunity():
    model = BoundarySCAR(
        rec_rate=0.1,
        M=[[0]],
        Ne=1.0,
        genome_length=20,
        boundary_positions=[10],
    )
    children = np.array([5, 5, 6])
    lefts = np.array([0, 10, 0])
    rights = np.array([10, 20, 5])

    assert model._get_line_recomb_opportunities(5, children, rights, lefts) == 1
    assert model._get_line_recomb_opportunities(6, children, rights, lefts) == 0


def test_boundary_scar_likelihood_is_finite_for_adjacent_segment_reassortment():
    tables = tskit.TableCollection(20)
    sample = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=-1)
    left_recomb = tables.nodes.add_row(flags=RECOMBINANT_FLAG, time=1, population=-1)
    right_recomb = tables.nodes.add_row(flags=RECOMBINANT_FLAG, time=1, population=-1)
    root = tables.nodes.add_row(flags=0, time=2, population=-1)
    tables.edges.add_row(0, 10, left_recomb, sample)
    tables.edges.add_row(10, 20, right_recomb, sample)
    tables.edges.add_row(0, 10, root, left_recomb)
    tables.edges.add_row(10, 20, root, right_recomb)
    tables.sort()
    ts = tables.tree_sequence()

    model = BoundarySCAR(
        rec_rate=0.1,
        M=[[0]],
        Ne=1.0,
        genome_length=20,
        boundary_positions=[10],
        bounds=(0, 1),
    )
    neg_log_like = model.compute_neg_log_like(0.1, ts)

    assert np.isfinite(neg_log_like)
    assert model.last_likelihood_diagnostics["status"] == "finite"
    assert model.last_likelihood_diagnostics["observed_recombination_events"] == 1
    assert model.last_likelihood_diagnostics["zero_recomb_opportunity_events"] == 0
    assert model.last_likelihood_diagnostics["recomb_exposure"] > 0


def test_structural_recombination_count_ignores_unpaired_recombinant_nodes(tmp_path):
    tables = tskit.TableCollection(20)
    sample = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=-1)
    singleton_sample = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=-1)
    left_recomb = tables.nodes.add_row(flags=RECOMBINANT_FLAG, time=1, population=-1)
    right_recomb = tables.nodes.add_row(flags=RECOMBINANT_FLAG, time=1, population=-1)
    unpaired_recomb = tables.nodes.add_row(flags=RECOMBINANT_FLAG, time=1, population=-1)
    root = tables.nodes.add_row(flags=0, time=2, population=-1)
    singleton_left_root = tables.nodes.add_row(flags=0, time=2, population=-1)
    singleton_right_root = tables.nodes.add_row(flags=0, time=2, population=-1)

    tables.edges.add_row(0, 10, left_recomb, sample)
    tables.edges.add_row(10, 20, right_recomb, sample)
    tables.edges.add_row(0, 10, root, left_recomb)
    tables.edges.add_row(10, 20, root, right_recomb)
    tables.edges.add_row(0, 10, unpaired_recomb, singleton_sample)
    tables.edges.add_row(0, 10, singleton_left_root, unpaired_recomb)
    tables.edges.add_row(10, 20, singleton_right_root, singleton_sample)
    tables.sort()
    ts = tables.tree_sequence()

    summary = summarize_recombination_events(ts)

    assert count_recombination_events(ts) == 1
    assert summary["n_recomb_nodes"] == 3
    assert summary["n_events"] == 1
    assert summary["unpaired_recomb_node_ids"] == [unpaired_recomb]

    seq_file = tmp_path / "segment.fasta"
    seq_file.write_text(">sample\nAAAA\n")
    recomb_info = RunArgAnalysis.extract_recombination_info(
        ts,
        [str(seq_file)],
        str(tmp_path),
    )

    assert recomb_info["n_recomb_events"] == 1
    assert recomb_info["n_recomb_nodes"] == 3
    assert recomb_info["unpaired_recomb_nodes"] == [unpaired_recomb]
    assert "Total recombination events: 1" in (tmp_path / "recombination_summary.txt").read_text()


def test_reassortment_rate_unit_helpers_distinguish_tree_time_units():
    assert RunArgAnalysis.get_reassortment_rate_units("divergence") == (
        "per_boundary_per_substitution_site",
        "per boundary per substitution/site",
    )
    assert RunArgAnalysis.annualize_reassortment_rate(
        2.0,
        "divergence",
        generation_time_days=3.0,
        clock_rate=0.005,
    ) == pytest.approx(0.01)
    assert RunArgAnalysis.annualize_reassortment_rate(
        2.0,
        "divergence",
        generation_time_days=3.0,
        clock_rate=None,
    ) is None
    assert RunArgAnalysis.annualize_reassortment_rate(
        2.0,
        "time",
        generation_time_days=3.0,
        clock_rate=None,
    ) == pytest.approx(2.0)
    assert RunArgAnalysis.annualize_reassortment_rate(
        2.0,
        "generations",
        generation_time_days=3.0,
        clock_rate=None,
    ) == pytest.approx(243.5)


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
        self.rate_model_name = "uniform_site"
        self.rate_units = "per_site_per_generation"
        self.rate_display_units = "per site per generation"
        self.bounds = (0.0, 1.0)
        self.last_opt_result = None

    def opt_MLE(self, ts):
        self.last_rate_estimate_status = "interior"
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
    fake_converter = lambda tree_path, tree_intervals: fake_ts

    default_result = builder.run_EM(
        ["tree"],
        ["seq"],
        ref,
        _FakeCoalModel(),
        iters=1,
        ts_converter=fake_converter,
    )
    tree_path_result = builder.run_EM(
        ["tree"],
        ["seq"],
        ref,
        _FakeCoalModel(),
        iters=1,
        return_tree_path=True,
        ts_converter=fake_converter,
    )

    assert len(default_result) == 3
    assert len(tree_path_result) == 4
    assert tree_path_result[3] is not None
    assert len(tree_path_result[3]) == 1

    diagnostics_result = builder.run_EM(
        ["tree"],
        ["seq"],
        ref,
        _FakeCoalModel(),
        iters=1,
        return_diagnostics=True,
        ts_converter=fake_converter,
        ts_converter_name="fake",
    )

    assert len(diagnostics_result) == 4
    diagnostics = diagnostics_result[3]
    assert diagnostics["estimate_source"] == "opt_mle"
    assert diagnostics["estimate_status"] == "interior"
    assert diagnostics["rate_units"] == "per_site_per_generation"
    assert diagnostics["ts_converter"] == "fake"
