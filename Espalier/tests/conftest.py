import ast
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import tskit


REPO_ROOT = Path(__file__).resolve().parents[2]
FLU_SUBSAMPLED_DIR = REPO_ROOT / "TestFiles" / "FluTest_subsampled"
FLU_SUBSAMPLED_ALIGNMENTS = FLU_SUBSAMPLED_DIR / "alignments"
FLU_SUBSAMPLED_TREES = FLU_SUBSAMPLED_DIR / "trees"


def run_python_script(script_relative_path, args, cwd=None):
    script_path = REPO_ROOT / script_relative_path
    completed = subprocess.run(
        [sys.executable, str(script_path), *args],
        cwd=str(cwd or REPO_ROOT),
        capture_output=True,
        text=True,
        check=True,
    )
    return completed


def run_arg_analysis(output_base, alignments_dir, trees_dir, extra_args=None):
    extra_args = extra_args or []
    completed = run_python_script(
        "scripts/run_arg_analysis.py",
        [
            "--alignments",
            str(alignments_dir),
            "--trees",
            str(trees_dir),
            "--output",
            str(output_base),
            "--analysis-mode",
            "reassortment",
            *extra_args,
        ],
    )
    output_dirs = sorted(Path(output_base).glob("arg_output_*"))
    if len(output_dirs) != 1:
        raise AssertionError(
            f"Expected exactly one analysis directory in {output_base}, found {len(output_dirs)}"
        )
    return output_dirs[0], completed


def parse_recombination_summary(summary_path):
    lines = Path(summary_path).read_text().splitlines()
    result = {}
    for line in lines:
        if line.startswith("Total recombination events:"):
            result["n_recomb_events"] = int(line.split(": ", 1)[1])
        elif line.startswith("Breakpoints:"):
            result["breakpoints"] = ast.literal_eval(line.split(": ", 1)[1])
    return result


def extract_pipeline_signature(output_dir):
    summary_path = output_dir / "recombination_summary.txt"
    signature = {
        "has_label_map": (output_dir / "derived" / "label_map.tsv").exists(),
        "has_trees": (output_dir / "arg_treesequence.trees").exists(),
        "summary": None,
    }
    if summary_path.exists():
        signature["summary"] = parse_recombination_summary(summary_path)
    if signature["has_trees"]:
        ts = tskit.load(str(output_dir / "arg_treesequence.trees"))
        signature["ts_breakpoints"] = list(ts.breakpoints())
        signature["ts_num_trees"] = ts.num_trees
    return signature


def load_inference_diagnostics(output_dir):
    diagnostics_path = Path(output_dir) / "em_diagnostics.json"
    if not diagnostics_path.exists():
        return {}
    return json.loads(diagnostics_path.read_text())


def has_raxml():
    return shutil.which("raxml-ng") is not None


def parse_logged_recombination_rate(completed_process):
    combined_output = f"{completed_process.stdout}\n{completed_process.stderr}"
    matches = re.findall(
        r"Estimated (?:recombination rate|reassortment rate): ([0-9.eE+-]+)",
        combined_output,
    )
    if not matches:
        return None
    return float(matches[-1])
