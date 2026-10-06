"""Preview-input guards and manifest verification.

Preview rows are the FIRST N rows of each harmonization source file. Those files
are sorted (intensity, protein group, or run order), so preview rows are a
structurally NON-RANDOM sample: valid for schema probing, invalid for any
statistic. This module makes preview semantics machine-checkable: tables written
in preview mode carry ``_PREVIEW_ONLY = True``, and downstream code MUST call
:func:`assert_not_preview_inputs` on every input frame so preview-flagged data
hard-fails before reaching computation.

Extracted verbatim from the source publication's ``preview_guard.py``
(analysis-memory ANL-009 / Adjustment A3), extended with the sha256 manifest
API used for run provenance and freeze verification.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

PREVIEW_FLAG_COLUMNS = ("_PREVIEW_ONLY", "preview_only")


class PreviewInputError(RuntimeError):
    """Raised when a preview-flagged table is passed where full evidence is required."""


def table_is_preview_only(table) -> bool:
    """True if the frame carries a preview flag column with any True-ish value."""
    for column in PREVIEW_FLAG_COLUMNS:
        if column in getattr(table, "columns", []):
            series = table[column]
            if series.fillna(False).astype(bool).any():
                return True
    return False


def assert_not_preview_inputs(table, name: str = "input table") -> None:
    """Hard-fail if `table` is preview-flagged (ANL-009 / Adjustment A3).

    Raises PreviewInputError with an explanatory message when the frame carries
    `_PREVIEW_ONLY`/`preview_only` flags, i.e. it contains first-N preview rows
    that must never feed statistical or physicochemical computation.
    """
    if table_is_preview_only(table):
        raise PreviewInputError(
            f"'{name}' is a PREVIEW-STAGE table (preview flag column is set). "
            "Preview rows are the first N rows of sorted source files - a "
            "non-random sample used only for schema probing (ANL-009). Downstream "
            "computation must consume the full harmonized evidence tables instead."
        )


# --- full-evidence table validation ------------------------------------------------

FULL_EVIDENCE_REQUIRED_COLUMNS = (
    "peptide_sequence",
    "uniprot_accession",
    "protein_accession_raw",
    "confidence_tier",
    "source_file",
    "extraction_method",
    "extraction_date",
    "dataset_accession",
)


class FullEvidenceSchemaError(RuntimeError):
    """Raised when a table presented as full peptide evidence lacks the S1 schema."""


def assert_full_evidence_table(table, name: str = "peptide evidence table") -> None:
    """Validate that `table` is the FULL peptide-evidence table (not preview, schema present).

    Downstream consumers (uniprot mapping onward) must consume the full table and call
    this on load: it asserts (a) no preview flags, and (b) presence of the canonical
    evidence/provenance/tier columns.
    """
    assert_not_preview_inputs(table, name)
    missing = [c for c in FULL_EVIDENCE_REQUIRED_COLUMNS if c not in getattr(table, "columns", [])]
    if missing:
        raise FullEvidenceSchemaError(
            f"'{name}' is missing required full-evidence columns {missing}. "
            "Expected the canonical harmonized peptide_evidence_table produced by "
            "proteoform_regions.harmonize()."
        )


# --- sha256 manifest API -------------------------------------------------------------


def sha256_of(path: str | Path, chunk: int = 8 << 20) -> str:
    """Streaming sha256 of a file (8 MiB blocks; identical to the paper's S1/S8 helpers)."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(chunk)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def write_manifest(files: list[str | Path], path: str | Path) -> Path:
    """Write a ``<sha256>  <relative-path>`` manifest (freeze format) for `files`."""
    base = Path(path).parent
    lines = [
        f"{sha256_of(f)}  {Path(f).relative_to(base) if _is_relative_to(f, base) else Path(f)}" for f in files
    ]
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return Path(path)


class ManifestVerificationError(RuntimeError):
    """Raised when one or more manifest entries fail sha256 verification."""


def verify_manifest(path: str | Path, base_dir: str | Path | None = None) -> bool:
    """Verify every ``<sha256>  <rel-path>`` line of a freeze-format manifest.

    ``base_dir`` resolves relative paths (defaults to the manifest's parent directory).
    Returns True when all entries match; raises ManifestVerificationError naming each
    mismatch, mirroring the paper's freeze gate ("FREEZE INTEGRITY FAIL - inference
    aborted").
    """
    path = Path(path)
    base = Path(base_dir) if base_dir else path.parent
    failures = []
    n = 0
    for line in path.read_text(encoding="utf-8").strip().splitlines():
        if not line.strip():
            continue
        n += 1
        expected, rel = line.split("  ", 1)
        target = Path(rel)
        if not target.is_absolute():
            target = base / target
        if not target.exists():
            failures.append(f"{rel}: MISSING")
        elif sha256_of(target) != expected:
            failures.append(f"{rel}: SHA256 MISMATCH")
    if failures:
        raise ManifestVerificationError(
            f"MANIFEST VERIFICATION FAIL ({len(failures)}/{n} entries): " + "; ".join(failures)
        )
    return True


def _is_relative_to(path: str | Path, base: Path) -> bool:
    try:
        Path(path).relative_to(base)
        return True
    except ValueError:
        return False
