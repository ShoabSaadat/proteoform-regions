"""Freeze-parity tests against the source publication's frozen artifacts.

Two layers:

1. **Committed fixtures** (always run, no paper repo needed): the frozen
   evidence schema and the frozen per-study row counts pin the k=10 tranche
   ground truth; the Sheet12 50-row sample pins the walkthrough dataset's
   structure (P6's fragpipe adapter consumes it).
2. **Full-parity tests** (``parity``-marked; require the paper repository at
   ``PROTEOFORM_REGIONS_PAPER_ROOT``, read-only): re-run the package's
   harmonization over all 10 studies and the S8 bias models over the frozen
   feature+background tables, asserting the published headline numbers.
"""

import os
from pathlib import Path

import pandas as pd
import pytest

from proteoform_regions import pipeline
from proteoform_regions.schema import CANONICAL_COLUMNS

DATA = Path(__file__).parent / "data"
PAPER_ROOT = Path(os.environ.get("PROTEOFORM_REGIONS_PAPER_ROOT", ""))
PARITY_AVAILABLE = PAPER_ROOT.exists() and (PAPER_ROOT / "analysis/data/processed").is_dir()

requires_paper = pytest.mark.skipif(
    not PARITY_AVAILABLE, reason="PROTEOFORM_REGIONS_PAPER_ROOT not set or incomplete"
)


# ---------------------------------------------------------------- fixtures layer


class TestFrozenSchemaFixture:
    def test_frozen_columns_are_the_canonical_schema(self):
        sample = pd.read_csv(DATA / "frozen_evidence_schema_sample.csv", low_memory=False)
        assert list(sample.columns) == list(CANONICAL_COLUMNS)

    def test_frozen_sample_passes_full_evidence_guard(self):
        # one row per study: shape-only columns the guard checks on a full table
        sample = pd.read_csv(DATA / "frozen_evidence_schema_sample.csv", low_memory=False)
        assert len(sample) == 10
        assert sample["dataset_accession"].nunique() == 10
        # glycopeptide rows are peptidoform-level (ANL-018); all others peptide-level
        glyco = sample["dataset_accession"] == "PXD057799"
        assert sample.loc[glyco, "entity_level"].eq("peptidoform").all()
        assert sample.loc[~glyco, "entity_level"].eq("peptide").all()
        assert set(sample["source_parser_family"]) == {
            "maxquant_or_text_peptide_table",
            "diann_report",
            "mzidentml",
            "specialized_glycopeptide_csv",
            "spectronaut_peptide_quant",
            "progenesis_peptide_csv",
        }

    def test_frozen_row_counts_total(self):
        counts = pd.read_csv(DATA / "frozen_study_row_counts.csv")
        assert int(counts["rows"].sum()) == 70_020
        assert len(counts) == 10
        # largest study is the S6 DIA-NN extension tranche
        top = counts.loc[counts["rows"].idxmax()]
        assert top["dataset_accession"] == "PXD069732" and top["rows"] == 20_416


class TestSheet12Fixture:
    """50-row sample of Sadeghi Sheet12 (FragPipe_Fractions_psm) for the P6 walkthrough."""

    def test_shape_and_key_columns(self):
        frame = pd.read_csv(DATA / "sheet12_sample50.csv")
        assert len(frame) == 50
        for col in (
            "Peptide",
            "Prev AA",
            "Next AA",
            "Protein Start",
            "Protein End",
            "Charge",
            "Number of Missed Cleavages",
            "Protein",
            "Protein ID",
        ):
            assert col in frame.columns, col

    def test_psm_content_sane(self):
        frame = pd.read_csv(DATA / "sheet12_sample50.csv")
        peptides = frame["Peptide"].dropna().astype(str)
        assert peptides.str.isalpha().all()  # unmodified stripped sequences
        assert frame["Protein Start"].notna().all()
        assert (frame["Number of Missed Cleavages"] >= 0).all()
        assert frame["Charge"].between(1, 7).all()


# ---------------------------------------------------------------- full-parity layer


@pytest.mark.parity
@requires_paper
class TestHarmonizeParity:
    def test_all_ten_studies_semantic_equality(self):
        """Package harmonize() == frozen peptide_evidence_table.csv (10 studies).

        Comparator (numeric rel 1e-6 / exact strings) imported from the
        gate-run harness so test and script cannot drift apart.
        """
        import sys

        from proteoform_regions import harmonize, load_cohort

        repo_root = Path(__file__).resolve().parents[1]
        sys.path.insert(0, str(repo_root / "scripts"))
        try:
            import parity_vs_paper_freeze as harness
        finally:
            sys.path.pop(0)

        cohort = repo_root / "examples/studies/cohort-paper-tranche.yaml"
        ev, _ = harmonize(
            load_cohort(cohort),
            PAPER_ROOT / "analysis/data/raw/peptide_tables",
            download=False,
        )
        frozen = pd.read_csv(
            PAPER_ROOT / "analysis/data/processed/peptide_evidence_table.csv", low_memory=False
        )
        assert len(ev) == len(frozen) == 70_020

        columns = [c for c in frozen.columns if c != "extraction_date"]
        worst = []
        for acc in sorted(frozen["dataset_accession"].unique()):
            got = ev[ev["dataset_accession"] == acc].reset_index(drop=True)
            want = frozen[frozen["dataset_accession"] == acc].reset_index(drop=True)
            assert len(got) == len(want), acc
            for col in columns:
                n_diff = harness._count_diffs(got[col], want[col])
                if n_diff:
                    worst.append((acc, col, n_diff))
        assert not worst, f"divergences: {worst[:10]}"


@pytest.fixture(scope="class")
def bias_result():
    processed = PAPER_ROOT / "analysis/data/processed"
    features = pd.read_csv(processed / "processed_feature_table.csv", low_memory=False)
    background = pd.read_csv(processed / "processed_tryptic_background.csv")
    return pipeline.quantify_bias(features, background, run_rq1=False)


@pytest.mark.parity
@requires_paper
class TestBiasParity:
    """S8 models on frozen feature+background inputs -> published headline numbers.

    Frozen (processed_inference_results.json, freeze 2026-09-17):
    FE OR 0.4519 [0.4431, 0.4608]; RE OR 0.461 p=7.14e-30; I2 97.6%;
    LOO [0.436, 0.494]; RQ1 T=-0.6712, null center +0.116, p=0.0005,
    CI [0.0814, 0.1455], 5388 proteins, 0/5388 BH q<0.05.
    """

    def test_rq2_fixed_effect_headline(self, bias_result):
        fe = bias_result["fe"]
        assert fe["or"] == pytest.approx(0.4519, abs=5e-5)
        assert fe["ci_low"] == pytest.approx(0.4431, abs=5e-5)
        assert fe["ci_high"] == pytest.approx(0.4608, abs=5e-5)
        assert fe["I2"] == pytest.approx(0.976, abs=5e-4)
        assert bias_result["k"] == 10

    def test_rq2_per_study_or_table(self, bias_result):
        frozen_or = {
            "PXD008583": 0.5559,
            "PXD052666": 0.4724,
            "PXD054594": 0.6086,
            "PXD055218": 0.4433,
            "PXD056620": 0.5284,
            "PXD057799": 0.4792,
            "PXD060933": 0.4186,
            "PXD068982": 0.3274,
            "PXD069732": 0.3605,
            "PXD071549": 0.4873,
        }
        table = bias_result["or_table"].set_index("dataset_accession")
        for study, or_ in frozen_or.items():
            assert table.loc[study, "or_"] == pytest.approx(or_, abs=5e-5), study

    def test_rq2_random_effects_and_loo(self, bias_result):
        random_effects = bias_result["random_effects"]
        assert random_effects["or_"] == pytest.approx(0.461, abs=5e-4)
        assert random_effects["p"] == pytest.approx(7.14e-30, rel=0.05)
        lo, hi = bias_result["loo_envelope"]
        assert (lo, hi) == pytest.approx((0.436, 0.494), abs=5e-4)

    def test_rq1_permutation_headline(self):
        processed = PAPER_ROOT / "analysis/data/processed"
        features = pd.read_csv(processed / "processed_feature_table.csv", low_memory=False)
        background = pd.read_csv(processed / "processed_tryptic_background.csv")
        rq1 = pipeline.quantify_bias(features, background, run_rq1=True)["rq1"]
        assert rq1["t_obs"] == pytest.approx(-0.6712, abs=5e-5)
        assert rq1["null_center"] == pytest.approx(0.116, abs=5e-4)
        assert rq1["p_global"] == pytest.approx(0.0005, abs=5e-5)
        assert rq1["ci_low"] == pytest.approx(0.0814, abs=5e-4)
        assert rq1["ci_high"] == pytest.approx(0.1455, abs=5e-4)
        assert rq1["n_proteins"] == 5_388
        assert rq1["n_significant"] == 0  # 0/5388 proteins at BH q<0.05 (published)
