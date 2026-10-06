"""Schema tests: 44 canonical columns, tier vocabulary, pI classification boundaries."""

import math

import pandas as pd

from proteoform_regions import schema


def test_canonical_column_count_and_order():
    assert len(schema.CANONICAL_COLUMNS) == 44
    # anchor positions to catch accidental reordering
    assert schema.CANONICAL_COLUMNS[0] == "study_id"
    assert schema.CANONICAL_COLUMNS[13] == "peptide_sequence"
    assert schema.CANONICAL_COLUMNS[-1] == "source_parser_family"  # appended by the S1 dispatch
    assert schema.CANONICAL_COLUMNS[-2] == "decoy_or_contaminant"
    assert "glycan_sialylation_flag" in schema.CANONICAL_COLUMNS
    assert "tryptic_background_eligible" in schema.CANONICAL_COLUMNS


def test_empty_row_has_all_columns_none():
    row = schema.empty_evidence_row()
    assert set(row) == set(schema.CANONICAL_COLUMNS)
    assert all(v is None for v in row.values())


def test_confidence_tiers_documented():
    assert set(schema.CONFIDENCE_TIERS) == {"Tier A", "Tier B", "Tier C", "Tier D"}


class TestClassifyPi:
    def test_acidic(self):
        assert schema.classify_pi(4.2) == "acidic"
        assert schema.classify_pi(6.999999) == "acidic"

    def test_boundary_7_0_is_near_neutral(self):
        assert schema.classify_pi(7.0) == "near_neutral"
        assert schema.classify_pi(7.499999) == "near_neutral"

    def test_boundary_7_5_is_basic(self):
        assert schema.classify_pi(7.5) == "basic"
        assert schema.classify_pi(12.3) == "basic"

    def test_none_and_nan(self):
        assert schema.classify_pi(None) is None
        assert schema.classify_pi(float("nan")) is None
        assert schema.classify_pi(math.nan) is None

    def test_pandas_na(self):
        frame = pd.DataFrame({"pi": [7.2, float("nan")]})
        assert schema.classify_pi(frame["pi"].iloc[1]) is None
