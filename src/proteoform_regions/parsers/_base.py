"""Shared parser plumbing: canonical-row seeding and coercion helpers."""

from __future__ import annotations

from typing import Any

import numpy as np

from ..schema import CANONICAL_COLUMNS
from ..study import StudyConfig


def base_row(config: StudyConfig, extraction_method: str, source_file: str) -> dict:
    """Seed a canonical-schema row from study metadata (S1/S6 ``base_row``)."""
    row: dict[str, Any] = {column: None for column in CANONICAL_COLUMNS}
    row.update(
        study_id=config.study_id,
        dataset_accession=config.dataset_accession,
        publication_doi=config.publication_doi,
        sample_type=config.sample_type,
        disease_context=config.disease_context,
        assay_type=config.assay_type,
        acquisition_mode=config.acquisition_mode,
        search_engine=config.search_engine,
        source_file=source_file,
        extraction_method=extraction_method,
        representative_file_arm=config.representative_file_arm,
        representative_file_species=config.representative_file_species,
        entity_level="peptide",
        tryptic_background_eligible=True,
        physicochemical_sequence_basis="unmodified_peptide",
        decoy_or_contaminant=False,
        run_count=1,
        source_parser_family=config.parser_family,
        confidence_tier="Tier C",
    )
    return row


def to_float(value) -> float | None:
    """Lenient float coercion (None for missing/unparseable; S1 ``_to_float``)."""
    try:
        if value is None or value == "" or (isinstance(value, float) and np.isnan(value)):
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def to_int(value) -> int | None:
    """Lenient int coercion (S1 ``_to_int``)."""
    try:
        if value is None or value == "":
            return None
        return int(float(value))
    except (TypeError, ValueError):
        return None
