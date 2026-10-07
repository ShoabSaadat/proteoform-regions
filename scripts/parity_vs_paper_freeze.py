#!/usr/bin/env python
"""Freeze-parity harness: package harmonize() vs the paper's frozen evidence table.

Runs the package's stage-1 harmonization over the source publication's 10-study
cohort (read-only against the paper repository's cached raw files) and compares
the result to the frozen ``peptide_evidence_table.csv`` (k=10, 70,020 rows,
freeze 2026-09-17) with semantic per-cell comparison:

- numeric-looking values compare as floats (relative tolerance 1e-6,
  absorbing cross-version pandas/numpy sum/median rounding drift);
- everything else compares as strings;
- ``extraction_date`` is excluded (run-time stamped by design).

Exit code 0 = parity; 1 = any divergence (printed per study/column with counts
and a first-diverging example).

Usage:
    python scripts/parity_vs_paper_freeze.py --paper-root /path/to/yinyan
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

REL_TOL = 1e-6  # absorbs cross-version float-repr/sum/median rounding drift
ABS_TOL = 1e-9  # near-zero values compare absolutely


def _as_float(v):
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return None
    try:
        f = float(str(v))
    except (ValueError, OverflowError):
        return None
    return f


def _cells_equal(a, b) -> bool:
    """Semantic cell equality: numeric within tolerance, else exact string.

    Numeric comparison (rel 1e-6) is the honest form of "6-significant-figure
    identity": formatting both sides to %.6g and string-comparing fails
    exactly AT rounding boundaries (0.00012175549 -> '0.000121755' vs
    0.00012175551 -> '0.000121756'), which is where pandas/numpy version
    drift (paper freeze: pandas 3.0.5/numpy 2.4.6) lands. NaN mismatches
    (one side numeric, other empty) stay SEMANTIC failures.
    """
    fa, fb = _as_float(a), _as_float(b)
    if fa is not None and fb is not None:
        return abs(fa - fb) <= max(REL_TOL * max(abs(fa), abs(fb)), ABS_TOL)
    return str(a) == str(b) or (a is None or pd.isna(a)) and (b is None or pd.isna(b))


def _count_diffs(got: pd.Series, want: pd.Series) -> int:
    return int(sum(not _cells_equal(a, b) for a, b in zip(got, want, strict=True)))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--paper-root", required=True, type=Path)
    ap.add_argument(
        "--data-dir",
        type=Path,
        default=None,
        help="raw peptide_tables dir (default: <paper>/analysis/data/raw/peptide_tables)",
    )
    ap.add_argument(
        "--cohort",
        type=Path,
        default=Path(__file__).parent.parent / "examples/studies/cohort-paper-tranche.yaml",
    )
    args = ap.parse_args()

    paper = args.paper_root
    data_dir = args.data_dir or (paper / "analysis/data/raw/peptide_tables")
    frozen_path = paper / "analysis/data/processed/peptide_evidence_table.csv"
    if not frozen_path.exists():
        print(f"FATAL: frozen table not found: {frozen_path}")
        return 2

    import proteoform_regions as pfr

    studies = pfr.load_cohort(args.cohort)
    print(f"cohort: {len(studies)} studies from {args.cohort.name}")
    print(f"raw dir: {data_dir}")

    evidence, stats = pfr.harmonize(studies, data_dir, download=False)
    frozen = pd.read_csv(frozen_path, low_memory=False)

    print(f"\npackage rows: {len(evidence)} | frozen rows: {len(frozen)}")
    print(stats.to_string(index=False))

    failures = 0
    if len(evidence) != len(frozen):
        failures += 1
        print("ROW-COUNT MISMATCH")

    skip_cols = {"extraction_date"}
    columns = [c for c in frozen.columns if c not in skip_cols]
    missing = [c for c in columns if c not in evidence.columns]
    if missing:
        print(f"MISSING COLUMNS in package output: {missing}")
        return 1

    # per-study comparison in cohort (frozen concatenation) order
    study_failed = []
    for study in studies:
        acc = study.dataset_accession
        got = evidence[evidence["dataset_accession"] == acc].reset_index(drop=True)
        want = frozen[frozen["dataset_accession"] == acc].reset_index(drop=True)
        if len(got) != len(want):
            print(f"  {acc}: ROW MISMATCH package={len(got)} frozen={len(want)}")
            study_failed.append(acc)
            failures += 1
            continue
        diffs = []
        for col in columns:
            n_diff = _count_diffs(got[col], want[col])
            if n_diff:
                first = next(i for i in range(len(got)) if not _cells_equal(got.at[i, col], want.at[i, col]))
                diffs.append(
                    f"{col} ({n_diff} cells, e.g. row {first}: "
                    f"{str(want.at[first, col])[:40]!r} != {str(got.at[first, col])[:40]!r})"
                )
        if diffs:
            print(f"  {acc}: {len(diffs)} diverging columns:")
            for d in diffs[:8]:
                print(f"      {d}")
            study_failed.append(acc)
            failures += 1
        else:
            print(f"  {acc}: PARITY ({len(got)} rows, {len(columns)} columns)")

    if failures == 0:
        print(
            "PARITY: PASSED (semantic equality on all studies, all columns;\n"
            f"numeric cells within rel {REL_TOL:g}, strings exact; extraction_date excluded)"
        )
        return 0
    print(f"PARITY: FAILED ({failures} divergences; studies: {study_failed})")
    return 1


if __name__ == "__main__":
    sys.exit(main())
