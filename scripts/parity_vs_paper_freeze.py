#!/usr/bin/env python
"""Freeze-parity harness: package harmonize() vs the paper's frozen evidence table.

Runs the package's stage-1 harmonization over the source publication's 10-study
cohort (read-only against the paper repository's cached raw files) and compares
the result to the frozen ``peptide_evidence_table.csv`` (k=10, 70,020 rows,
freeze 2026-09-17) with semantic per-cell comparison:

- numeric-looking values compare as floats;
- everything else compares as strings;
- ``extraction_date`` is excluded (run-time stamped by design).

Exit code 0 = parity; 1 = any divergence (printed per study/column with counts
and a first-diverging example).

Usage:
    python scripts/parity_vs_paper_freeze.py --paper-root /path/to/yinyan
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

import pandas as pd


def _norm(v) -> str:
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return ""
    s = str(v)
    try:
        f = float(s)
    except (ValueError, OverflowError):
        return s
    # semantic numeric identity: 6 significant figures absorbs cross-version
    # float-repr/summation artifacts (numpy pairwise-sum order, median pairing)
    if f == int(f) and abs(f) < 1e15:
        return str(int(f))
    return f"{f:.6g}"


def _canonical(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.apply(lambda col: col.map(_norm))


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
        got = _canonical(evidence[evidence["dataset_accession"] == acc].reset_index(drop=True))
        want = _canonical(frozen[frozen["dataset_accession"] == acc].reset_index(drop=True))
        if len(got) != len(want):
            print(f"  {acc}: ROW MISMATCH package={len(got)} frozen={len(want)}")
            study_failed.append(acc)
            failures += 1
            continue
        diffs = []
        for col in columns:
            neq = got[col] != want[col]
            if neq.any():
                i = int(neq.to_numpy().nonzero()[0][0])
                diffs.append(
                    f"{col} ({int(neq.sum())} cells, e.g. row {i}: {want.at[i, col][:40]!r} != {got.at[i, col][:40]!r})"  # noqa: E501
                )
        if diffs:
            print(f"  {acc}: {len(diffs)} diverging columns:")
            for d in diffs[:8]:
                print(f"      {d}")
            study_failed.append(acc)
            failures += 1
        else:
            print(f"  {acc}: PARITY ({len(got)} rows, {len(columns)} columns)")

    pkg_hash = hashlib.sha256(_canonical(evidence[columns]).to_csv(index=False).encode()).hexdigest()
    froz_hash = hashlib.sha256(_canonical(frozen[columns]).to_csv(index=False).encode()).hexdigest()
    print(f"\ncanonical-frame sha256: package={pkg_hash[:16]}... frozen={froz_hash[:16]}...")
    if failures == 0 and pkg_hash == froz_hash:
        print("PARITY: IDENTICAL (semantic equality on all studies, all columns)")
        return 0
    print(f"PARITY: FAILED ({failures} divergences; studies: {study_failed})")
    return 1


if __name__ == "__main__":
    sys.exit(main())
