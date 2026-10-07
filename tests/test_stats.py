"""Detection-bias statistics tests (S8 extraction).

Golden values are HAND-COMPUTED (calculator algebra in the comments) so a
regression in Woolf/BH/FE algebra cannot hide behind notebook parity; the
frozen-paper parity lives in the ``parity``-marked tests and
``scripts/parity_vs_paper_freeze.py``.
"""

import numpy as np
import pandas as pd
import pytest
from scipy import stats as sps

from proteoform_regions import stats as st


class TestBhQ:
    def test_hand_computed(self):
        # BH step-up: sorted p [0.005, 0.01, 0.03, 0.04], m=4
        # q_raw = p*m/rank = [0.02, 0.02, 0.04, 0.04]; reverse cummin -> [0.02, 0.02, 0.04, 0.04]
        # (already monotone); remap to original order.
        p = [0.01, 0.04, 0.03, 0.005]
        assert st.bh_q(p) == pytest.approx([0.02, 0.04, 0.04, 0.02])

    def test_monotone_enforcement(self):
        # sorted [0.018, 0.019] -> raw q = [0.036, 0.019]; the step-up takes the
        # running min from the LARGEST rank, so rank 1 inherits rank 2's 0.019
        assert st.bh_q([0.019, 0.018]) == pytest.approx([0.019, 0.019])

    def test_capped_at_one_and_empty(self):
        assert (st.bh_q([0.9, 1.0]) <= 1.0).all()
        assert len(st.bh_q([])) == 0

    def test_all_identical(self):
        # raw q = [1.5, 0.75, 0.5]; step-up min -> all inherit the last rank's 0.5
        assert st.bh_q([0.5, 0.5, 0.5]) == pytest.approx([0.5, 0.5, 0.5])


class TestOr2x2:
    def test_hand_computed(self):
        # lor = ln((10/20)/(30/40)) = ln(2/3) = -0.405465
        # se  = sqrt(1/10 + 1/20 + 1/30 + 1/40) = sqrt(0.208333) = 0.456435
        lor, se = st.or_2x2(10, 20, 30, 40)
        assert lor == pytest.approx(np.log(2 / 3))
        assert se == pytest.approx(np.sqrt(0.2083333333))

    def test_unconditional_half_correction(self):
        # EVERY cell is max(x, 0.5) even when only one is zero (always-on policy)
        lor, se = st.or_2x2(0, 10, 20, 30)
        a, b, c, d = 0.5, 10, 20, 30
        assert lor == pytest.approx(np.log((a / b) / (c / d)))
        assert se == pytest.approx(np.sqrt(1 / a + 1 / b + 1 / c + 1 / d))

    def test_symmetry_or_one(self):
        lor, se = st.or_2x2(5, 5, 5, 5)
        assert lor == pytest.approx(0.0)
        assert se == pytest.approx(np.sqrt(4 / 5))


class TestStudyOr:
    def test_counts_and_inference(self):
        det = pd.Series(["basic"] * 12 + ["acidic"] * 8)
        bg = pd.Series(["basic"] * 30 + ["near_neutral"] * 70)
        r = st.study_or(det, bg, "basic")
        assert (r["a"], r["b"], r["c"], r["d"]) == (12, 8, 30, 70)
        lor, se = st.or_2x2(12, 8, 30, 70)
        assert r["log_or"] == pytest.approx(lor)
        assert r["se"] == pytest.approx(se)
        assert r["or_"] == pytest.approx(float(np.exp(lor)))
        assert r["ci_low"] == pytest.approx(float(np.exp(lor - 1.96 * se)))
        assert r["p"] == pytest.approx(float(2 * sps.norm.sf(abs(lor / se))))
        # basic ENRICHED in detections (60% vs 30%) -> OR > 1
        assert r["or_"] > 1.0

    def test_drops_nan_labels(self):
        det = pd.Series(["basic", None, "acidic"])
        bg = pd.Series(["basic", "basic"])
        r = st.study_or(det, bg)
        assert (r["a"], r["b"]) == (1, 1)


class TestFeSummary:
    def _table(self, rows):
        return pd.DataFrame(rows, columns=["dataset_accession", "log_or", "se"])

    def test_single_study_identity(self):
        t = self._table([("S1", -0.7, 0.1)])
        mu, se, p, q_stat, i2_stat = st.fe_summary(t)
        assert mu == pytest.approx(-0.7)
        assert se == pytest.approx(0.1)
        assert q_stat == pytest.approx(0.0)
        assert i2_stat == 0.0

    def test_two_study_hand_computed(self):
        # w1 = 1/0.1^2 = 100, w2 = 1/0.2^2 = 25
        # mu = (100*-0.7 + 25*-0.3)/125 = -77.5/125 = -0.62
        # se = sqrt(1/125) = 0.08944
        t = self._table([("S1", -0.7, 0.1), ("S2", -0.3, 0.2)])
        mu, se, p, q_stat, i2_stat = st.fe_summary(t)
        assert mu == pytest.approx(-0.62)
        assert se == pytest.approx(np.sqrt(1 / 125))
        z = -0.62 / np.sqrt(1 / 125)
        assert p == pytest.approx(float(2 * sps.norm.sf(abs(z))))
        # Q = 100*(-0.7+0.62)^2 + 25*(-0.3+0.62)^2 = 100*0.0064 + 25*0.1024 = 0.64 + 2.56 = 3.2
        assert q_stat == pytest.approx(3.2)
        # I2 = (Q - df)/Q with df=1 -> (3.2-1)/3.2 = 0.6875
        assert i2_stat == pytest.approx(0.6875)

    def test_i2_floored_at_zero(self):
        # homogeneous studies: Q ~ 0 < df -> I2 floors at 0 (never negative)
        t = self._table([("S1", -0.5, 0.1), ("S2", -0.5, 0.1), ("S3", -0.5, 0.1)])
        _, _, _, q_stat, i2_stat = st.fe_summary(t)
        assert q_stat == pytest.approx(0.0, abs=1e-12)
        assert i2_stat == 0.0


class TestDlRandomEffects:
    def test_tau2_zero_when_homogeneous(self):
        t = pd.DataFrame(
            [("S1", -0.5, 0.1), ("S2", -0.5, 0.1)], columns=["dataset_accession", "log_or", "se"]
        )
        dl = st.dl_random_effects(t)
        assert dl["tau2"] == 0.0
        mu_fe, *_ = st.fe_summary(t)
        assert dl["mu"] == pytest.approx(mu_fe)

    def test_tau2_positive_when_heterogeneous(self):
        t = pd.DataFrame(
            [("S1", -1.0, 0.05), ("S2", 0.0, 0.05), ("S3", -1.0, 0.05)],
            columns=["dataset_accession", "log_or", "se"],
        )
        dl = st.dl_random_effects(t)
        assert dl["tau2"] > 0
        # RE weights flatten -> mu_re pulled toward the unweighted mean (-2/3) vs FE
        _, se_fe, *_ = st.fe_summary(t)
        assert dl["se"] > se_fe

    def test_returns_or_scale(self):
        t = pd.DataFrame(
            [("S1", np.log(0.5), 0.1), ("S2", np.log(0.6), 0.1)],
            columns=["dataset_accession", "log_or", "se"],
        )
        dl = st.dl_random_effects(t)
        assert dl["or_"] == pytest.approx(float(np.exp(dl["mu"])))


class TestLooEnvelope:
    def test_two_study_envelope_is_pair(self):
        t = pd.DataFrame(
            [("S1", np.log(0.4), 0.1), ("S2", np.log(0.8), 0.1)],
            columns=["dataset_accession", "log_or", "se"],
        )
        lo, hi = st.loo_envelope(t)
        # dropping S1 -> FE over {S2} = 0.8 ; dropping S2 -> FE over {S1} = 0.4
        assert (lo, hi) == pytest.approx((0.4, 0.8))


class TestRq1Permutation:
    @staticmethod
    def _frames(acidic_shift: float):
        # one protein (60 aa, tryptic peptides ~10 residues) detected in one study
        features = pd.DataFrame(
            {
                "dataset_accession": ["PXDTEST"] * 3,
                "mapped_accession": ["P99999"] * 3,
                "peptide_sequence": ["AAAAAAAALK", "AAAAAAADK", "AAAAAAAEK"],
                "peptide_pi": [4.0 + acidic_shift, 4.2 + acidic_shift, 4.4 + acidic_shift],
                "precursor_pi": [7.5, 7.5, 7.5],
            }
        )
        background = pd.DataFrame(
            {
                "accession": ["P99999"] * 6,
                "peptide_pi": [6.0, 6.4, 6.8, 7.0, 7.3, 7.6],
            }
        )
        return features, background

    def test_seed_determinism(self):
        f, b = self._frames(-1.0)
        r1 = st.rq1_permutation(f, b, seed=42, b_perm=200, b_boot=200)
        r2 = st.rq1_permutation(f, b, b_perm=200, b_boot=200)  # default seed 42
        assert r1["p_global"] == r2["p_global"]
        assert r1["t_obs"] == r2["t_obs"]
        assert r1["ci_low"] == r2["ci_low"]
        assert np.array_equal(r1["per_protein"]["perm_p"].to_numpy(), r2["per_protein"]["perm_p"].to_numpy())

    def test_observed_statistic_and_direction(self):
        f, b = self._frames(-2.0)
        r = st.rq1_permutation(f, b, seed=42, b_perm=500, b_boot=200)
        # t_obs = median detected pI - precursor = median(2.2, 2.4, 2.6... shifted) - 7.5
        expected = float(np.median([2.0, 2.2, 2.4]) - 7.5)
        assert r["t_obs"] == pytest.approx(expected)
        # null: median of 3 with-replacement draws from the 6-value bg, centered near the
        # bg middle (with-replacement medians concentrate centrally; not exactly median(bg))
        assert -0.9 < r["null_center"] < -0.3
        # detected far more acidic than any null draw -> minimal p
        assert r["p_global"] == pytest.approx(1 / 201)

    @pytest.mark.filterwarnings("ignore:Mean of empty slice:RuntimeWarning")
    @pytest.mark.filterwarnings("ignore:invalid value encountered:RuntimeWarning")
    def test_skips_protein_with_short_background(self):
        # documented absence: bg has fewer peptides than detected n_unique
        features = pd.DataFrame(
            {
                "dataset_accession": ["PXDTEST"] * 2,
                "mapped_accession": ["P99999"] * 2,
                "peptide_sequence": ["AAAAAAAALK", "AAAAAAADK"],
                "peptide_pi": [4.0, 4.2],
                "precursor_pi": [7.5, 7.5],
            }
        )
        background = pd.DataFrame({"accession": ["P99999"], "peptide_pi": [6.0]})
        r = st.rq1_permutation(features, background, seed=42, b_perm=50, b_boot=50)
        assert r["n_proteins"] == 0
        assert np.isnan(r["t_obs"])

    def test_matches_verbatim_notebook_loop(self):
        """Package rq1_permutation == the notebook's inline loop, same seed."""
        f, b = self._frames(-1.5)
        result = st.rq1_permutation(f, b, seed=42, b_perm=300, b_boot=100)

        # verbatim S8 loop (B constants swapped for speed; rng stream otherwise identical)
        SEED, B_PERM, B_BOOT = 42, 300, 100
        rng = np.random.default_rng(SEED)
        per_protein = (
            f.dropna(subset=["precursor_pi"])
            .groupby(["dataset_accession", "mapped_accession"])
            .agg(
                n_unique=("peptide_sequence", "nunique"),
                median_peptide_pi=("peptide_pi", "median"),
                precursor_pi=("precursor_pi", "first"),
            )
            .reset_index()
        )
        per_protein = per_protein[per_protein["n_unique"] >= 2]
        bg_by_accession = {acc: g["peptide_pi"].dropna().to_numpy() for acc, g in b.groupby("accession")}
        obs_divs, null_divs_global, per_protein_p = [], [], []
        for record in per_protein.itertuples(index=False):
            bg = bg_by_accession.get(record.mapped_accession)
            if bg is None or len(bg) < record.n_unique:
                continue
            obs = record.median_peptide_pi - record.precursor_pi
            obs_divs.append(obs)
            draws = bg[rng.integers(0, len(bg), size=(B_PERM, record.n_unique))]
            null_medians = np.median(draws, axis=1)
            null_div = null_medians - record.precursor_pi
            null_divs_global.append(np.median(null_div))
            per_protein_p.append((1 + np.sum(np.abs(null_div) >= abs(obs))) / (B_PERM + 1))
        T_obs = float(np.median(np.array(obs_divs)))
        global_null = np.array(
            [
                np.median(rng.choice(null_divs_global, size=len(null_divs_global), replace=True))
                for _ in range(B_BOOT)
            ]
        )
        center = np.median(null_divs_global)
        exceed = np.abs(global_null - center) >= abs(T_obs - center)
        p_global = (1 + np.sum(exceed)) / (B_BOOT + 1)
        assert result["t_obs"] == pytest.approx(T_obs)
        assert result["p_global"] == pytest.approx(p_global)
        assert result["ci_low"] == pytest.approx(float(np.percentile(global_null, 2.5)))
        assert result["ci_high"] == pytest.approx(float(np.percentile(global_null, 97.5)))


class TestBootstrapAndEnrichment:
    def test_class_shares_exact_and_ci_brackets(self):
        labels = ["acidic"] * 50 + ["near_neutral"] * 30 + ["basic"] * 20
        shares, cis = st.bootstrap_class_shares(labels, n=500, seed=42)
        assert shares["acidic"] == pytest.approx(0.5)
        assert shares["near_neutral"] == pytest.approx(0.3)
        assert shares["basic"] == pytest.approx(0.2)
        for c in st.CLASSES:
            assert cis[c][0] <= shares[c] <= cis[c][1]

    def test_enrichment_ratio_exact(self):
        det = ["basic"] * 10 + ["acidic"] * 90
        bg = ["basic"] * 50 + ["acidic"] * 50
        out = st.enrichment_with_ci(det, bg, n=500, seed=42)
        # enrichment = share ratio (S5 definition): detected share / background share
        assert out["basic"]["enrichment"] == pytest.approx(0.10 / 0.50)
        assert out["acidic"]["enrichment"] == pytest.approx(0.90 / 0.50)
        assert out["basic"]["ci"][0] < 0.2 < out["basic"]["ci"][1]

    def test_empty_inputs(self):
        shares, cis = st.bootstrap_class_shares([], n=10)
        assert shares == {c: 0.0 for c in st.CLASSES}
        out = st.enrichment_with_ci([], ["acidic"], n=10)
        assert all(v["enrichment"] is None for v in out.values())


class TestTheoreticalCharge:
    def test_krh_plus_nterm(self):
        assert st.theoretical_charge_ph26("AAAA") == 1
        assert st.theoretical_charge_ph26("AAAK") == 2
        assert st.theoretical_charge_ph26("KRRHHAAA") == 6  # K+R+R+H+H... K,R,R,H,H = 5 + 1

    def test_verbatim_uppercase_only(self):
        # verbatim: counts literal 'KRH' characters; lowercase input counts as 0 (+ N-term)
        assert st.theoretical_charge_ph26("krh") == 1
