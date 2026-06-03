from pathlib import Path

import pytest
from Bio import SeqIO

from Espalier.Reassortment import ANALYSIS_MODE_REASSORTMENT, prepare_analysis_inputs
from Espalier.tests.conftest import FLU_SUBSAMPLED_ALIGNMENTS, FLU_SUBSAMPLED_TREES


pytestmark = pytest.mark.realdata


def _collect_alignment_records(path: Path):
    return list(SeqIO.parse(str(path), "fasta"))


def _collect_isolate_dates(records):
    return {
        record.description.split("|")[0]: record.description.split("|")[-1]
        for record in records
    }


def test_flu_subsampled_fixture_has_shared_isolates_and_consistent_dates():
    ha_records = _collect_alignment_records(FLU_SUBSAMPLED_ALIGNMENTS / "HA-swine_H1_HANA.fasta.aln")
    na_records = _collect_alignment_records(FLU_SUBSAMPLED_ALIGNMENTS / "NA-swine_H1_HANA.fasta.aln")

    assert len(ha_records) == 434
    assert len(na_records) == 434
    assert len(ha_records[0]) == 1719
    assert len(na_records[0]) == 1464

    ha_dates = _collect_isolate_dates(ha_records)
    na_dates = _collect_isolate_dates(na_records)
    assert set(ha_dates) == set(na_dates)
    assert len(ha_dates) == 434
    assert ha_dates == na_dates


def test_prepare_analysis_inputs_normalizes_subsampled_flu_fixture(tmp_path):
    seq_files = [
        str(FLU_SUBSAMPLED_ALIGNMENTS / "HA-swine_H1_HANA.fasta.aln"),
        str(FLU_SUBSAMPLED_ALIGNMENTS / "NA-swine_H1_HANA.fasta.aln"),
    ]
    tree_files = [
        str(FLU_SUBSAMPLED_TREES / "HA-swine_H1_HANA.fasta.aln.rooted.tre"),
        str(FLU_SUBSAMPLED_TREES / "NA-swine_H1_HANA.fasta.aln.rooted.tre"),
    ]

    prepared = prepare_analysis_inputs(
        seq_files,
        tree_files,
        str(tmp_path / "output"),
        analysis_mode=ANALYSIS_MODE_REASSORTMENT,
        isolate_field=0,
    )

    assert prepared.label_map_path is not None
    assert Path(prepared.label_map_path).exists()
    assert prepared.common_isolates is not None
    assert len(prepared.common_isolates) == 434
    assert len(prepared.seq_files) == 2
    assert len(prepared.tree_files) == 2

    for normalized_seq in prepared.seq_files:
        records = list(SeqIO.parse(normalized_seq, "fasta"))
        assert len(records) == 434
        assert all("|" not in record.id for record in records)


def test_prepare_analysis_inputs_writes_complete_flu_tree_outputs(tmp_path):
    prepared = prepare_analysis_inputs(
        [
            str(FLU_SUBSAMPLED_ALIGNMENTS / "HA-swine_H1_HANA.fasta.aln"),
            str(FLU_SUBSAMPLED_ALIGNMENTS / "NA-swine_H1_HANA.fasta.aln"),
        ],
        [
            str(FLU_SUBSAMPLED_TREES / "HA-swine_H1_HANA.fasta.aln.rooted.tre"),
            str(FLU_SUBSAMPLED_TREES / "NA-swine_H1_HANA.fasta.aln.rooted.tre"),
        ],
        str(tmp_path / "output"),
        analysis_mode=ANALYSIS_MODE_REASSORTMENT,
        isolate_field=0,
    )

    for normalized_tree in prepared.tree_files:
        tree_text = Path(normalized_tree).read_text()
        assert "sample" not in tree_text
        assert "|" not in tree_text
