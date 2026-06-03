from pathlib import Path

import pytest
import tskit

from Espalier.tests.conftest import (
    FLU_SUBSAMPLED_ALIGNMENTS,
    FLU_SUBSAMPLED_TREES,
    extract_pipeline_signature,
    has_raxml,
    load_inference_diagnostics,
    run_arg_analysis,
)


pytestmark = [
    pytest.mark.slow,
    pytest.mark.requires_raxml,
    pytest.mark.realdata,
    pytest.mark.skipif(not has_raxml(), reason="raxml-ng is required for inference benchmarks"),
]


def test_flu_reassortment_pipeline_produces_expected_artifacts(tmp_path):
    output_dir, completed = run_arg_analysis(
        tmp_path / "flu_output",
        FLU_SUBSAMPLED_ALIGNMENTS,
        FLU_SUBSAMPLED_TREES,
        extra_args=["--em-iters", "1", "--rec-rate", "1e-4"],
    )

    assert (output_dir / "derived" / "label_map.tsv").exists()
    assert (output_dir / "derived" / "normalized_alignments").is_dir()
    assert (output_dir / "derived" / "normalized_trees").is_dir()
    assert (output_dir / "consensus_reference.tre").exists()
    assert any(path.name.endswith("_ARG_with_recomb.tre") for path in output_dir.iterdir())
    assert "Analysis mode: reassortment" in completed.stderr
    diagnostics = load_inference_diagnostics(output_dir)
    assert diagnostics["rate_model_name"] == "boundary_hotspot"
    assert diagnostics["rate_units"] == "per_boundary_per_generation"

    if (output_dir / "arg_treesequence.trees").exists():
        ts = tskit.load(str(output_dir / "arg_treesequence.trees"))
        assert ts.num_trees >= 1
        assert list(ts.breakpoints()) == sorted(ts.breakpoints())

    if (output_dir / "recombination_summary.txt").exists():
        summary_text = (output_dir / "recombination_summary.txt").read_text()
        assert "ARG Recombination Summary" in summary_text
        assert "Breakpoints:" in summary_text


def test_flu_reassortment_pipeline_is_repeatable_at_coarse_artifact_level(tmp_path):
    first_output, _ = run_arg_analysis(
        tmp_path / "flu_run_1",
        FLU_SUBSAMPLED_ALIGNMENTS,
        FLU_SUBSAMPLED_TREES,
        extra_args=["--em-iters", "1", "--rec-rate", "1e-4"],
    )
    second_output, _ = run_arg_analysis(
        tmp_path / "flu_run_2",
        FLU_SUBSAMPLED_ALIGNMENTS,
        FLU_SUBSAMPLED_TREES,
        extra_args=["--em-iters", "1", "--rec-rate", "1e-4"],
    )

    first_signature = extract_pipeline_signature(first_output)
    second_signature = extract_pipeline_signature(second_output)

    assert first_signature["has_label_map"]
    assert second_signature["has_label_map"]
    assert first_signature["summary"] == second_signature["summary"]
    if first_signature["has_trees"] and second_signature["has_trees"]:
        assert first_signature["ts_breakpoints"] == second_signature["ts_breakpoints"]
        assert first_signature["ts_num_trees"] == second_signature["ts_num_trees"]
