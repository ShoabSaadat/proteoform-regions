"""Stage 1: harmonize heterogeneous peptide tables into the canonical evidence schema.

``harmonize`` is a pure function over study configs + cached raw files: it
downloads missing files (retry-once policy, failures documented never
fabricated - the paper's S1 ``ensure_downloads`` semantics), runs each study's
adapter, and concatenates canonical-schema rows. Write-out is the caller's
job (see ``pipeline`` / ``cli``).
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import pandas as pd

from .guard import sha256_of
from .parsers import parse_study
from .schema import CANONICAL_COLUMNS
from .study import StudyConfig


def ensure_downloads(
    studies: list[StudyConfig], data_dir: str | Path, timeout: float = 600.0
) -> pd.DataFrame:
    """Fetch missing raw files (skip-if-exists, one retry) and return a manifest."""
    import requests

    manifest_rows: list[dict] = []
    for study in studies:
        for spec in study.files:
            if spec.url is None:
                continue
            dest = Path(data_dir) / study.dataset_accession / spec.name
            dest.parent.mkdir(parents=True, exist_ok=True)
            if not dest.exists():
                _fetch(spec.url, dest, timeout)
            manifest_rows.append(
                {
                    "dataset_accession": study.dataset_accession,
                    "file_name": spec.name,
                    "source_url": spec.url,
                    "bytes_on_disk": dest.stat().st_size if dest.exists() else 0,
                    "sha256": sha256_of(dest) if dest.exists() else None,
                }
            )
    return pd.DataFrame(manifest_rows)


def _fetch(url: str, dest: Path, timeout: float, retries: int = 2) -> None:
    import time

    for attempt in range(1, retries + 1):
        try:
            with requests.get(url, timeout=timeout, stream=True) as response:
                response.raise_for_status()
                with open(dest, "wb") as handle:
                    for chunk in response.iter_content(chunk_size=8 << 20):
                        handle.write(chunk)
            return
        except Exception:
            if attempt == retries:
                # documented failure, never fabricated (paper policy)
                print(f"DOWNLOAD FAILED (documented, proceeding): {url}")
                return
            time.sleep(2.0)


def harmonize(
    studies: list[StudyConfig],
    data_dir: str | Path,
    download: bool = True,
    extraction_date: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Harmonize a cohort; returns (evidence DataFrame, run stats DataFrame).

    The evidence frame carries exactly the 44 canonical columns; the stats
    frame records per-study row counts and adapter-reported counters.
    """
    extraction_date = extraction_date or datetime.now().strftime("%Y-%m-%d")
    if download:
        ensure_downloads(studies, data_dir)
    all_rows: list[dict] = []
    stats_rows: list[dict] = []
    for study in studies:
        rows, stats = parse_study(study, data_dir)
        for row in rows:
            row["extraction_date"] = extraction_date
        all_rows.extend(rows)
        stats_rows.append(stats)
    evidence = pd.DataFrame(all_rows, columns=CANONICAL_COLUMNS)
    stats = pd.DataFrame(stats_rows)
    return evidence, stats
