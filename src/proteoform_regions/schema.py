"""Canonical evidence schema, confidence tiers, and pre-specified pI classification.

The 44-column canonical schema is the package's data contract: every parser
adapter emits rows in this schema, and every downstream stage consumes it.
Column list extracted verbatim from the source publication's S1 harmonization
(notebook ``peptide_table_harmonization`` §6.9).
"""

from __future__ import annotations

import math

# --- pre-specified pI classification (ANL-003; asserted BEFORE any outcome analysis) ---

ACIDIC_BOUND = 7.0  # pI < 7.0 -> acidic
NEUTRAL_BOUND = 7.5  # 7.0 <= pI < 7.5 -> near_neutral; pI >= 7.5 -> basic


def classify_pi(value) -> str | None:
    """Pre-specified pI classification: acidic < 7.0 <= near_neutral < 7.5 <= basic.

    Returns None for None/NaN input (the S3 variant of the paper's helper, which
    had the most defensive NaN handling of the four notebook copies).
    """
    if value is None:
        return None
    if isinstance(value, float) and math.isnan(value):
        return None
    if value < ACIDIC_BOUND:
        return "acidic"
    if value < NEUTRAL_BOUND:
        return "near_neutral"
    return "basic"


# --- confidence tiers (evidence semantics) ---

CONFIDENCE_TIERS = {
    "Tier A": "sequence + accession + q-value (FDR) + intensity",
    "Tier B": "sequence + accession + PEP/score; no FDR q exported by source format",
    "Tier C": "sequence + accession only; no q / intensity-FDR semantics",
    "Tier D": "no peptide sequences exposed (protein-group level); not peptide evidence",
}

ENTITY_LEVELS = (
    "peptide",  # unmodified sequence observation
    "peptidoform",  # sequence + modification state (e.g. glycopeptide backbone+glycan)
    "detected_region",  # stitched multi-peptide region of a precursor (this package's core output)
    "inferred_fragment",
    "proteoform",  # NEVER asserted from bottom-up evidence alone (ANL-002)
)

# --- canonical evidence columns (44; verbatim from S1) ---

CANONICAL_COLUMNS = [
    "study_id",
    "dataset_accession",
    "publication_doi",
    "sample_type",
    "disease_context",
    "assay_type",
    "acquisition_mode",
    "search_engine",
    "protein_accession_raw",
    "uniprot_accession",
    "protein_group",
    "gene_symbol",
    "protein_name_raw",
    "peptide_sequence",
    "modified_sequence",
    "proforma_sequence",
    "start_position",
    "end_position",
    "charge_state",
    "psm_count",
    "run_count",
    "peptide_intensity",
    "protein_intensity",
    "q_value",
    "posterior_error_probability",
    "score_value",
    "score_label",
    "score_source",
    "source_file",
    "extraction_method",
    "extraction_date",
    "confidence_tier",
    "entity_level",
    "glycan_annotation",
    "glycan_sialylation_flag",
    "naked_backbone_sequence",
    "proforma_style_sequence",
    "tryptic_background_eligible",
    "physicochemical_sequence_basis",
    "representative_file_arm",
    "representative_file_species",
    "usi_example",
    "decoy_or_contaminant",
    "source_parser_family",
]


def empty_evidence_row() -> dict:
    """A canonical-schema row with every column set to None (S1 ``base_row`` seed)."""
    return {column: None for column in CANONICAL_COLUMNS}
