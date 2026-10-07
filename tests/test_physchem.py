"""Physicochemical engine tests.

Three layers of evidence:
1. GOLDEN  - the publication's frozen ExPASy benchmark (16 sequences; pI within
   0.005 of ExPASy, one documented extreme-composition outlier).
2. PARITY  - the package's parameterized refactor vs the verbatim notebook
   reference implementation (identical float results).
3. LITERAL - pK table values pinned to the Supplementary Methods tables.
"""

import csv
from importlib.resources import files

import numpy as np
import pytest

from proteoform_regions import digest, physchem
from tests import reference_implementations as ref

GOLDEN = files("proteoform_regions") / "data" / "physchem_golden.csv"

SEQS = [
    "ALPAPIEK",
    "AAAATGTIFTFR",
    "AAAGFNVSLTDYWGR",
    "EEEEEEEEEEEK",
    "KKKKKKKKKKK",
    "KPYEEELK",
    "GNLNEQVFLK",
    "THPGFQPSLK",
    "AVGDKLPDYHGK",
    "QEPERNECFLQHK",
]


class TestGoldenBenchmark:
    @pytest.mark.parametrize(
        "row",
        [pytest.param(r, id=r["name"]) for r in csv.DictReader(GOLDEN.open())],
    )
    def test_pi_and_mw_match_frozen_publication_values(self, row):
        pi = physchem.biopython_pi(row["sequence"])
        mw = physchem.biopython_mw(row["sequence"])
        assert pi == pytest.approx(float(row["golden_pi"]), abs=1e-9)
        assert mw == pytest.approx(float(row["golden_mw"]), rel=1e-12)

    def test_golden_pi_within_0p005_of_expasy_documented_outlier_excepted(self):
        rows = list(csv.DictReader(GOLDEN.open()))
        devs = [
            (abs(float(row["golden_pi"]) - float(row["expasy_pi"])), row["name"])
            for row in rows
            if row["expasy_pi"]
        ]
        outliers = [d for d in devs if d[0] > 0.006]
        # paper: agreement within 0.005 pI units (rounded claim; worst non-outlier is
        # 0.00502) on all but ONE extreme-composition poly-E control (~0.42)
        assert outliers == [(pytest.approx(0.4200, abs=1e-3), "synthetic/control peptide 1")]
        assert max(d[0] for d in devs if d[1] != "synthetic/control peptide 1") <= 0.0051

    def test_mw_agrees_with_uniprot_declared_within_documented_ppm(self):
        # NOTE: for fragment rows the frozen CSV's uniprot_mw column holds the PARENT
        # protein's molWeight (benchmark artifact), so ppm is only meaningful for
        # full-protein rows. Paper: median -11.5 ppm, p95 |deviation| 32 ppm (n=2,602).
        # 7 of the 8 benchmark proteins are in the UniProt snapshot (P00004 was
        # filtered out by the notebook's identical `if a in entries` gate).
        rows = [r for r in csv.DictReader(GOLDEN.open()) if r["uniprot_mw"] and r["kind"] == "protein"]
        assert len(rows) == 7
        for row in rows:
            mw = physchem.biopython_mw(row["sequence"])
            ppm = (mw - float(row["uniprot_mw"])) / float(row["uniprot_mw"]) * 1e6
            assert abs(ppm) < 100, row["name"]


class TestRefactorParity:
    @pytest.mark.parametrize("table", ["hh_rodwell_emboss_lineage", "hh_lehninger"])
    def test_hh_pi_identical_to_notebook_implementation(self, table):
        ref.hh_pi_configure(table)
        expected = ref.hh_pi_batch(SEQS)
        got = physchem.hh_pi_batch(SEQS, table=table)
        np.testing.assert_array_equal(got, expected)  # exact float identity, not approx

    def test_hh_charge_identical(self):
        ref.hh_pi_configure("hh_rodwell_emboss_lineage")
        counts = ref.count_matrix(SEQS)
        expected = ref.hh_charge(7.4, counts)
        got = physchem.hh_charge(7.4, counts, table="hh_rodwell_emboss_lineage")
        np.testing.assert_array_equal(got, expected)

    def test_charge_at_ph_maps_unique_sequences(self):
        mapping = physchem.charge_at_ph(["AAA", "AAA", "KRR"], ph=7.4)
        assert list(mapping) == ["AAA", "KRR"]
        assert mapping["KRR"] > 0 and mapping["AAA"] < 1.0


class TestPkTables:
    def test_rodwell_literals(self):
        side = physchem.PK_TABLES["hh_rodwell_emboss_lineage"]["side"]
        assert side["K"] == pytest.approx(11.5)
        assert side["R"] == pytest.approx(11.5)
        assert side["D"] == pytest.approx(3.86)
        assert side["C"] == pytest.approx(8.33)
        assert side["Y"] == pytest.approx(10.7)
        assert side["H"] == pytest.approx(6.0)
        assert side["E"] == pytest.approx(4.25)
        assert physchem.PK_TABLES["hh_rodwell_emboss_lineage"]["n_term"] == pytest.approx(8.0)
        assert physchem.PK_TABLES["hh_rodwell_emboss_lineage"]["c_term"] == pytest.approx(3.1)

    def test_lehninger_literals(self):
        t = physchem.PK_TABLES["hh_lehninger"]
        assert t["n_term"] == pytest.approx(9.69)
        assert t["c_term"] == pytest.approx(2.34)
        assert t["side"]["H"] == pytest.approx(6.0)
        assert t["side"]["R"] == pytest.approx(12.48)

    def test_pi_between_0_and_14(self):
        pis = physchem.hh_pi_batch(SEQS)
        assert ((pis > 0) & (pis < 14)).all()

    def test_extreme_controls_ordering(self):
        pis = physchem.hh_pi_batch(["EEEEEEEEEEEK", "KKKKKKKKKKK"])
        assert pis[0] < 7.0 < pis[1]


class TestCompositionFeatures:
    def test_peptide_features_values(self):
        feats = physchem.peptide_features("ALPAPIEK")
        assert feats["peptide_length"] == 8
        assert feats["acidic_frac"] == pytest.approx(1 / 8)
        assert feats["basic_frac"] == pytest.approx(1 / 8)
        assert feats["tryptic_terminus"] is True
        assert feats["missed_cleavages"] == 0

    def test_missed_cleavage_count_internal_kr(self):
        # internal R followed by A counts; terminal R is a normal cut site, not internal
        assert physchem.peptide_features("AAAKPAAARAAG")["missed_cleavages"] == 1
        # K immediately before P is a non-site for trypsin (never a missed cleavage)
        assert physchem.peptide_features("AAAKPAAAG")["missed_cleavages"] == 0
        assert physchem.peptide_features("KPEEEER")["missed_cleavages"] == 0

    def test_gravy_range(self):
        assert physchem.gravy("IIIIIIII") > 0
        assert physchem.gravy("DDDDDDDD") < 0
        assert physchem.gravy("") is None

    def test_x_replacement_in_pi(self):
        assert physchem.biopython_pi("XXXXXXXX") == physchem.biopython_pi("AAAAAAAA")


class TestDigest:
    def test_parity_with_notebook(self):
        seq = "MKTLLLTLAVL" + "KAAEEKAAAAKAAAAEKAAAAPKAAAAAARAAAARAAAAR"
        assert digest.tryptic_digest(seq) == ref.tryptic_digest(seq)

    def test_keil_rule_blocks_before_proline(self):
        peptides = [p for p, _ in digest.tryptic_digest("AAAKPAARAAPR")]
        # every retained peptide ends at a real cut site: R whose next residue is not P
        assert all(p[-1] == "R" for p in peptides)

    def test_length_bounds(self):
        seq = "K" + "A" * 100 + "R" + "A" * 3 + "R"
        for peptide, _ in digest.tryptic_digest(seq):
            assert 7 <= len(peptide) <= 60

    def test_missed_cleavage_bins_documented_enumeration(self):
        # DIVERGENCE #31 (pinned): the published enumeration with max_missed=2
        # admits mc in {0,1,2,3} - the frozen paper background contains 271,317
        # mc=3 rows. Supp Methods "0, 1, 2" describes the analysis strata
        # (mc0/mc1/mc2+), not this bound. max_missed=1 reproduces strict 0..2.
        seq = "K" + "A" * 8 + "R" + "A" * 8 + "R" + "A" * 8 + "R" + "A" * 8 + "R"
        mcs = {mc for _, mc in digest.tryptic_digest(seq)}
        assert mcs == {0, 1, 2, 3}
        strict = {mc for _, mc in digest.tryptic_digest(seq, max_missed=1)}
        assert strict == {0, 1, 2}

    def test_parameterized_bounds(self):
        seq = "K" + "A" * 4 + "R" + "A" * 8 + "R"
        assert any(len(p) == 5 for p, _ in digest.tryptic_digest(seq, min_len=5, max_len=60))

    def test_background_dedup_per_accession_peptide(self):
        rows = list(digest.digest_protein("P99999", "KAAAAAAAAR"))
        keys = [(r["accession"], r["peptide_sequence"]) for r in rows]
        assert len(keys) == len(set(keys))
