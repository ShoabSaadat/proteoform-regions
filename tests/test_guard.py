"""Guard tests: preview trip, full-evidence schema, sha256 manifest round-trip."""

import pandas as pd
import pytest

from proteoform_regions import guard


def _frame(**cols):
    return pd.DataFrame(cols)


class TestPreviewGuard:
    def test_preview_flag_trips(self):
        frame = _frame(peptide_sequence=["A", "B"], _PREVIEW_ONLY=[True, True])
        with pytest.raises(guard.PreviewInputError, match="PREVIEW-STAGE"):
            guard.assert_not_preview_inputs(frame, "preview evidence")

    def test_legacy_flag_trips(self):
        frame = _frame(peptide_sequence=["A"], preview_only=[True])
        with pytest.raises(guard.PreviewInputError):
            guard.assert_not_preview_inputs(frame)

    def test_clear_table_passes(self):
        frame = _frame(peptide_sequence=["A", "B"])
        guard.assert_not_preview_inputs(frame, "full evidence")

    def test_false_flag_passes(self):
        frame = _frame(peptide_sequence=["A"], _PREVIEW_ONLY=[False])
        guard.assert_not_preview_inputs(frame)


class TestFullEvidenceSchema:
    REQUIRED = list(guard.FULL_EVIDENCE_REQUIRED_COLUMNS)

    def test_valid_table_passes(self):
        frame = _frame(**{c: ["x"] for c in self.REQUIRED})
        guard.assert_full_evidence_table(frame)

    def test_missing_columns_raise(self):
        frame = _frame(peptide_sequence=["x"])
        with pytest.raises(guard.FullEvidenceSchemaError, match="missing required"):
            guard.assert_full_evidence_table(frame)

    def test_preview_fails_before_schema_check(self):
        frame = _frame(**{c: ["x"] for c in self.REQUIRED}, _PREVIEW_ONLY=[True])
        with pytest.raises(guard.PreviewInputError):
            guard.assert_full_evidence_table(frame)


class TestManifest:
    def test_roundtrip_and_tamper_detection(self, tmp_path):
        a = tmp_path / "a.csv"
        b = tmp_path / "sub" / "b.csv"
        b.parent.mkdir()
        a.write_text("hello")
        b.write_text("world")

        manifest = tmp_path / "manifest.txt"
        guard.write_manifest([a, b], manifest)
        assert guard.verify_manifest(manifest) is True

        # tamper -> mismatch named
        b.write_text("world!")
        with pytest.raises(guard.ManifestVerificationError, match="SHA256 MISMATCH"):
            guard.verify_manifest(manifest)

        # delete -> missing named
        b.unlink()
        with pytest.raises(guard.ManifestVerificationError, match="MISSING"):
            guard.verify_manifest(manifest)

    def test_sha256_matches_paper_format(self, tmp_path):
        f = tmp_path / "f.bin"
        f.write_bytes(b"proteoform-regions")
        import hashlib

        assert guard.sha256_of(f) == hashlib.sha256(b"proteoform-regions").hexdigest()
