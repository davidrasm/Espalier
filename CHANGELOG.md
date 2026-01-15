# Espalier Fork - Changelog

## Summary of Modifications

This fork addresses several compatibility issues with modern Python packages and extends Espalier's functionality for viral reassortment analysis. The core EM/Viterbi approach and SCAR likelihood model remain unchanged.

---

## Bug Fixes & Compatibility Updates

### 1. NumPy 2.0 Compatibility
**File:** `Espalier/SCARLikelihood.py`

- Replaced deprecated `np.Inf` with `np.inf`
- NumPy 2.0 removed the capitalized `Inf` constant

### 2. Pandas 2.0 Compatibility
**File:** `Espalier/Dendro2TSConverter.py`

- Replaced all `DataFrame.append()` calls with `pd.concat()`
- `append()` was deprecated in pandas 1.4 and removed in 2.0
- Affected functions: `tree2nodesdf()`, `split_edge()`, `split_edge_on_rec()`, `postorder_check()`, `preorder_check()`

### 3. RAxML Negative Branch Length Handling
**File:** `Espalier/RAxML.py`

- Added `_sanitize_tree_file()` method to fix negative branch lengths
- RAxML-NG occasionally produces small negative values (e.g., -1e-6) due to numerical precision
- These cause downstream parsing failures in dendropy
- Solution: Replace negative values with a small positive epsilon (1e-10)

### 4. Taxon Label Parsing
**File:** `Espalier/Dendro2TSConverter.py`

- Fixed `int()` conversion error for non-numeric taxon labels
- Labels like `MglDCr02|2012-01-01` now handled correctly
- Added try/except block to preserve string labels when numeric parsing fails

### 5. Large Bitmask Integer Overflow
**File:** `Espalier/Dendro2TSConverter.py`

- Fixed `Python int too large to convert to C long` error
- Occurs with large trees where bitmask unique_ids exceed 64-bit integer range
- Solution: Use Python list comprehension instead of `.astype('int')` for id_dict creation

---

## ARG Reconstruction Improvements

### 6. Time Constraint Violation Handling
**File:** `Espalier/ARGBuilder.py`

- Added relaxed placement strategy when `max_recomb_time < min_recomb_time`
- Previously, these violations caused recombination nodes to be skipped
- Now uses midpoint placement with a warning logged
- Ensures recombination events are captured even with edge-case timing conflicts

### 7. Robust TreeSequence Conversion
**File:** `Espalier/Dendro2TSConverter.py`

- Added `fix_time_constraints()` function with iterative bottom-up approach
- Ensures `time[parent] > time[child]` for all edges (tskit requirement)
- Added retry logic with increasing epsilon for stubborn constraint violations
- Added systematic `tables.sort()` calls before and after `edges.squash()`

### 8. Tree Path Return from EM
**File:** `Espalier/ARGBuilder.py`

- Modified `run_EM()` to return the tree path with recombination nodes
- Previously only returned TreeSequence, which loses tree structure if conversion fails
- New return signature: `(ts, rec_rate, n_recomb, tree_path)`

---

## New Scripts

### 9. ARG Analysis Pipeline
**File:** `scripts/run_arg_analysis.py`

Complete pipeline script demonstrating Espalier usage:
- ML tree inference per segment (RAxML-NG)
- Consensus reference tree construction
- Basic ARG reconstruction with fixed recombination rate
- EM-based joint ARG reconstruction and rate estimation
- Recombination node extraction and segment boundary mapping
- Outputs TreeSequence files and recombination summaries

### 10. A/B Testing Framework
**Directory:** `scripts/ab_testing/`

Framework for comparing original vs. modernized implementations:
- `modern_converter.py`: Reimplemented converter using tskit 0.6.0+ features
  - Vectorized time constraint fixing (3.3x faster)
  - `extend_haplotypes()` for edge compression
  - Proper tree path preparation for Viterbi outputs
- `run_comparison.py`: Automated A/B test runner with metrics collection
- `generate_results_pdf.py`: PDF report generation

**A/B Test Results:**
| Metric | Original | Modern |
|--------|----------|--------|
| Viterbi success rate | 0% | 100% |
| EM success rate | 100% | 100% |
| Average speedup | — | 3.31x |
| Likelihood | Identical | Identical |

---

## Configuration & Documentation

### 11. Environment Specification
**File:** `environment.yml`

- Added conda environment file with pinned dependencies
- Ensures reproducible installation

### 12. Updated Documentation Links
**Files:** `README.md`, `setup.py`, `docs/source/conf.py`, `docs/source/intro.rst`

- Updated repository URLs to fork location
- Updated paper reference

### 13. Gitignore Updates
**File:** `.gitignore`

- Added patterns for test outputs (`TestFiles/Outputs/`, `arg_output_*/`)
- Added `*.trees` for TreeSequence files
- Prevents large output files from being committed

---

## Files Modified

| File | Changes |
|------|---------|
| `Espalier/ARGBuilder.py` | Relaxed time placement, tree_path return |
| `Espalier/Dendro2TSConverter.py` | pandas 2.0 compat, time fixes, bitmask overflow |
| `Espalier/RAxML.py` | Negative branch length sanitization |
| `Espalier/SCARLikelihood.py` | NumPy 2.0 compat |
| `Espalier/Reconciler.py` | Minor fixes |
| `Espalier/viz/PlotTanglegrams.py` | Minor fixes |
| `scripts/run_arg_analysis.py` | New analysis pipeline |
| `scripts/ab_testing/*` | New A/B testing framework |

---

## Compatibility

Tested with:
- Python 3.10+
- NumPy 2.0+
- pandas 2.0+
- tskit 0.5.0+
- dendropy 4.5+
- RAxML-NG

---

## Notes for Future Development

The A/B testing revealed that the modern converter's `prepare_tree_path()` function (which reconciles linked heights and adds recombination nodes) is essential for Viterbi outputs. The original converter expects trees that have already been processed through the EM pipeline's preprocessing steps. This could be integrated into the main codebase to make `Dendro2TSConverter.convert()` more robust for direct Viterbi path conversion.
