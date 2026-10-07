"""End-to-end walkthrough test: the committed Sheet12 sample, offline.

Executes exactly what ``examples/sadeghi_walkthrough.ipynb`` demonstrates
(the notebook itself additionally runs under nbmake in CI): harmonize the
50-PSM FragPipe sample through the fragpipe adapter, map against the
committed offline UniProt snapshot, compute features, build the
digestion-matched background, stitch regions - no network, no paper repo.
"""

import json
from pathlib import Path

import pandas as pd
import pytest

from proteoform_regions import run

REPO = Path(__file__).resolve().parents[1]
WALKTHROUGH = REPO / "examples/walkthrough"

STAGES = ("harmonize", "map", "compute", "background", "stitch")


@pytest.fixture(scope="module")
def report(tmp_path_factory):
    return run(
        WALKTHROUGH / "cohort-sheet12-sample.yaml",
        tmp_path_factory.mktemp("walkthrough"),
        stages=list(STAGES),
        data_dir=WALKTHROUGH,
        download=False,
        uniprot_cache=WALKTHROUGH,
        snapshot_date="20261007",
        offline=True,
    )


class TestWalkthroughEndToEnd:
    def test_all_stages_ran_offline(self, report):
        assert tuple(report.stages) == STAGES
        assert report.parameters["snapshot_date"] == "20261007"
        assert report.counts["harmonize"] == 23  # 50 PSMs -> 23 (peptide, protein)
        assert report.counts["map"] == 23
        assert report.counts["map_properties"] == 22  # one per parent protein
        assert report.counts["background"] > 1000  # digestion-matched pool
        assert report.counts["stitch"] >= 10

    def test_artifacts_use_paper_vocabulary(self, report):
        expected = {
            "harmonize": "peptide_evidence_table.csv",
            "map": "processed_peptide_protein_map.csv",
            "map_properties": "processed_protein_properties.csv",
            "compute": "processed_feature_table.csv",
            "background": "processed_tryptic_background.csv",
            "stitch": "processed_region_inference.csv",
        }
        assert report.artifacts.items() >= expected.items()
        assert report.artifacts["manifest"] == "manifest.sha256"
        for name in expected.values():
            assert (report.out_dir / name).exists(), name
        manifest = (report.out_dir / "manifest.sha256").read_text()
        for name in expected.values():
            assert name in manifest  # every artifact hash-pinned

    def test_region_table_carries_positions_and_classes(self, report):
        regions = pd.read_csv(report.out_dir / "processed_region_inference.csv")
        assert regions["dataset_accession"].eq("PXD077545").all()
        assert regions["region_start"].notna().all()
        assert (regions["region_end"] >= regions["region_start"]).all()
        assert set(regions["region_pi_class_anl003"]) <= {"acidic", "basic", "neutral"}
        # single-peptide regions are labeled peptide_only; >=2 peptides within
        # the gap are detected_region - the 50-PSM sample is mostly singletons
        assert set(regions["confidence_label"]) <= {"detected_region", "peptide_only"}

    def test_run_report_is_json_round_trippable(self, report):
        payload = json.loads((report.out_dir / "run_report.json").read_text())
        assert payload["stages"] == list(STAGES)
        assert payload["counts"]["harmonize"] == 23
