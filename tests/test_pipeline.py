"""End-to-end pipeline tests on the synthetic two-protein cohort.

Exercises the full ``run()`` arc (harmonize -> map -> compute -> background ->
stitch -> bias) offline: schema-valid evidence, exact-substring mapping,
feature columns, digestion-matched background, region stitching, OR
meta-analysis, run report, and the sha256 run manifest.
"""

import json

import numpy as np
import pandas as pd
import pytest

from proteoform_regions import guard, pipeline
from tests.synthetic import PROTEINS, STUDY_PEPTIDES, write_cohort


@pytest.fixture(scope="module")
def cohort(tmp_path_factory):
    root = tmp_path_factory.mktemp("cohort")
    return write_cohort(root)


@pytest.fixture(scope="module")
def full_run(cohort, tmp_path_factory):
    cohort_yaml, data_dir, cache_dir = cohort
    out_dir = tmp_path_factory.mktemp("run")
    report = pipeline.run(
        cohort_yaml,
        out_dir,
        data_dir=data_dir,
        download=False,
        uniprot_cache=cache_dir,
        offline=True,
        b_perm=100,
        b_boot=100,
    )
    return out_dir, report


class TestRunAllStages:
    def test_all_artifacts_written(self, full_run):
        out_dir, report = full_run
        for name in (
            "peptide_evidence_table.csv",
            "processed_peptide_protein_map.csv",
            "processed_protein_properties.csv",
            "processed_feature_table.csv",
            "processed_tryptic_background.csv",
            "processed_region_inference.csv",
            "processed_inference_results.json",
            "run_report.json",
            "manifest.sha256",
        ):
            assert (out_dir / name).exists(), name

    def test_report_records_stages_and_counts(self, full_run):
        _out, report = full_run
        assert report.stages == list(pipeline.STAGES)
        assert report.counts["harmonize"] == sum(len(p) for p in STUDY_PEPTIDES.values())
        assert report.counts["bias"] == len(STUDY_PEPTIDES)  # k = studies
        assert report.parameters["seed"] == 42
        assert report.proteoform_regions_version

    def test_manifest_verifies(self, full_run):
        out_dir, _ = full_run
        assert guard.verify_manifest(out_dir / "manifest.sha256", base_dir=out_dir)

    def test_evidence_is_schema_valid(self, full_run):
        out_dir, _ = full_run
        evidence = pd.read_csv(out_dir / "peptide_evidence_table.csv", low_memory=False)
        guard.assert_full_evidence_table(evidence, "synthetic")
        assert set(evidence["dataset_accession"]) == set(STUDY_PEPTIDES)

    def test_mapping_exact_for_all_detected(self, full_run):
        out_dir, _ = full_run
        mapping = pd.read_csv(out_dir / "processed_peptide_protein_map.csv")
        assert mapping["mapping_tier"].isin(["exact_unique", "exact_multi_position"]).all()

    def test_features_carry_physchem_and_precursor(self, full_run):
        out_dir, _ = full_run
        features = pd.read_csv(out_dir / "processed_feature_table.csv", low_memory=False)
        for col in (
            "peptide_pi",
            "peptide_class_anl003",
            "net_charge_ph74",
            "charge_density",
            "precursor_pi",
            "pi_divergence_peptide_minus_precursor",
            "missed_cleavages",
        ):
            assert col in features.columns
        assert features["peptide_class_anl003"].notna().all()
        assert np.isfinite(features["net_charge_ph74"]).all()

    def test_background_is_digestion_matched(self, full_run):
        out_dir, _ = full_run
        bg = pd.read_csv(out_dir / "processed_tryptic_background.csv")
        assert set(bg["accession"]) == set(PROTEINS)
        assert (bg["peptide_length"].between(7, 60)).all()
        # every detected peptide is itself a tryptic peptide -> present in background
        detected = set()
        for peptides in STUDY_PEPTIDES.values():
            detected.update(seq for seq, _ in peptides)
        bg_seqs = set(bg["peptide_sequence"])
        assert detected <= bg_seqs

    def test_regions_stitch_and_tier(self, full_run):
        out_dir, _ = full_run
        regions = pd.read_csv(out_dir / "processed_region_inference.csv")
        supported = regions[regions["confidence_label"] == "region_supported"]
        assert len(supported) >= 1
        assert (supported["region_n_peptides"] >= 2).all()
        # study A sees 4 adjacent PTEST01 peptides -> one >=4-peptide region
        p1 = supported[
            (supported["dataset_accession"] == "PXDTESTA") & (supported["mapped_accession"] == "PTEST01")
        ]
        assert p1["region_n_peptides"].max() >= 4
        assert p1["coverage_fraction"].max() > 0.5

    def test_bias_json_structure(self, full_run):
        out_dir, _ = full_run
        payload = json.loads((out_dir / "processed_inference_results.json").read_text())
        assert payload["k"] == 2
        assert len(payload["or_table"]) == 2
        fe = payload["fe"]
        assert fe["or"] > 0 and fe["ci_low"] < fe["or"] < fe["ci_high"]
        assert len(payload["loo_envelope"]) == 2
        assert "rq1" in payload and payload["rq1"]["n_proteins"] >= 1


class TestStageSelection:
    def test_bias_only_resumes_from_disk(self, cohort, full_run, tmp_path):
        cohort_yaml, data_dir, cache_dir = cohort
        out_dir, _ = full_run
        report = pipeline.run(
            cohort_yaml,
            out_dir,
            stages="bias",
            data_dir=data_dir,
            download=False,
            uniprot_cache=cache_dir,
            offline=True,
            b_perm=50,
            b_boot=50,
        )
        assert report.stages == ["bias"]
        assert report.counts["bias"] == 2

    def test_unknown_stage_rejected(self, cohort, tmp_path):
        cohort_yaml, *_ = cohort
        with pytest.raises(ValueError, match="unknown stages"):
            pipeline.run(cohort_yaml, tmp_path / "nope", stages="harmonize,warp", download=False)

    def test_harmonize_only_writes_evidence(self, cohort, tmp_path):
        cohort_yaml, data_dir, _ = cohort
        report = pipeline.run(cohort_yaml, tmp_path, stages="harmonize", data_dir=data_dir, download=False)
        assert report.stages == ["harmonize"]
        assert (tmp_path / "peptide_evidence_table.csv").exists()
        assert not (tmp_path / "processed_feature_table.csv").exists()


class TestStageFunctions:
    def test_compute_features_calculator_switch(self, full_run, cohort):
        out_dir, _ = full_run
        _cache_dir = cohort[2]
        evidence = pd.read_csv(out_dir / "peptide_evidence_table.csv", low_memory=False)
        mapping = pd.read_csv(out_dir / "processed_peptide_protein_map.csv")
        proteins = pd.read_csv(out_dir / "processed_protein_properties.csv", low_memory=False)
        bj = pipeline.compute_features(evidence, mapping, proteins, calculator="bjellqvist")
        hh = pipeline.compute_features(evidence, mapping, proteins, calculator="hh_rodwell_emboss_lineage")
        joined = bj[["peptide_sequence", "peptide_pi"]].merge(
            hh[["peptide_sequence", "peptide_pi"]], on="peptide_sequence", suffixes=("_bj", "_hh")
        )
        assert (joined["peptide_pi_bj"] - joined["peptide_pi_hh"]).abs().max() > 0.01

    def test_unknown_calculator_rejected(self, full_run):
        out_dir, _ = full_run
        evidence = pd.read_csv(out_dir / "peptide_evidence_table.csv", low_memory=False)
        mapping = pd.read_csv(out_dir / "processed_peptide_protein_map.csv")
        proteins = pd.read_csv(out_dir / "processed_protein_properties.csv", low_memory=False)
        with pytest.raises(ValueError, match="unknown calculator"):
            pipeline.compute_features(evidence, mapping, proteins, calculator="nope")

    def test_stitch_empty_mapping_returns_empty(self, full_run):
        out_dir, _ = full_run
        features = pd.read_csv(out_dir / "processed_feature_table.csv", low_memory=False)
        empty_map = pd.read_csv(out_dir / "processed_peptide_protein_map.csv").iloc[:0]
        assert pipeline.stitch_regions(features, empty_map).empty

    def test_quantify_bias_without_rq1(self, full_run):
        out_dir, _ = full_run
        features = pd.read_csv(out_dir / "processed_feature_table.csv", low_memory=False)
        bg = pd.read_csv(out_dir / "processed_tryptic_background.csv")
        result = pipeline.quantify_bias(features, bg, run_rq1=False)
        assert "rq1" not in result
        assert result["k"] == 2
        assert np.isfinite(result["fe"]["or"])
