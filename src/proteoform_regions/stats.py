"""Detection-bias statistical models (Supplementary Methods Modules 4-6).

Pure, seed-aware functions extracted verbatim from the source publication's
S5 (descriptive), S8 (inferential), and charge-state (J4) layers:

- ``bh_q``                     Benjamini-Hochberg q-values
- ``or_2x2`` / ``study_or``    2x2 OR with Woolf SE (UNCONDITIONAL max(x, 0.5)
                               zero-cell correction on every cell)
- ``fe_summary``               inverse-variance fixed effect + Q + I2
- ``dl_random_effects``        DerSimonian-Laird tau2 / RE summary
- ``loo_envelope``             leave-one-study-out FE envelope
- ``rq1_permutation``          digestion-matched per-protein permutation null
                               (WITH replacement, B=2000, seed 42; two-stage
                               global p via resampling null medians)
- ``bootstrap_class_shares`` / ``enrichment_with_ci``   descriptive CIs
- ``theoretical_charge_ph26``  KRH+N-term protonation proxy (charge-state tests)
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

__all__ = [
    "bh_q",
    "or_2x2",
    "study_or",
    "fe_summary",
    "dl_random_effects",
    "loo_envelope",
    "rq1_permutation",
    "bootstrap_class_shares",
    "enrichment_with_ci",
    "theoretical_charge_ph26",
]

SEED = 42
B_PERM = 2000  # permutation iterations per protein (RQ1)
B_BOOT = 2000  # bootstrap iterations (CIs / envelopes)


def bh_q(pvals) -> np.ndarray:
    """Benjamini-Hochberg q-values (step-up, monotone via reverse cummin)."""
    pvals = np.asarray(pvals, dtype=float)
    order = np.argsort(pvals)
    ranked = pvals[order]
    m = len(ranked)
    q = ranked * m / (np.arange(m) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1.0)
    return out


def or_2x2(a: float, b: float, c: float, d: float) -> tuple[float, float]:
    """(log OR, Woolf SE) for detected-vs-background x basic-vs-rest.

    ``a``=detect&basic, ``b``=detect&rest, ``c``=bg&basic, ``d``=bg&rest.
    Zero-cell correction is UNCONDITIONAL: max(x, 0.5) on every cell always
    (always-on Haldane-Anscombe; simpler than the conditional variant).
    """
    a, b, c, d = max(a, 0.5), max(b, 0.5), max(c, 0.5), max(d, 0.5)
    lor = np.log((a / b) / (c / d))
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return lor, se


def study_or(det_classes, bg_classes, target: str = "basic") -> dict:
    """Per-study OR dictionary from class label arrays."""
    det_classes = pd.Series(det_classes).dropna()
    bg_classes = pd.Series(bg_classes).dropna()
    a = int((det_classes == target).sum())
    b = int((det_classes != target).sum())
    c = int((bg_classes == target).sum())
    d = int((bg_classes != target).sum())
    lor, se = or_2x2(a, b, c, d)
    return dict(
        a=a, b=b, c=c, d=d, log_or=lor, se=se, or_=float(np.exp(lor)),
        ci_low=float(np.exp(lor - 1.96 * se)), ci_high=float(np.exp(lor + 1.96 * se)),
        p=float(2 * stats.norm.sf(abs(lor / se))),
    )


def fe_summary(table: pd.DataFrame) -> tuple[float, float, float, float, float]:
    """Inverse-variance fixed effect over a per-study OR table.

    Requires columns ``log_or`` and ``se``. Returns (mu, se, p, Q, I2).
    """
    w = 1 / table["se"] ** 2
    mu = float(np.sum(w * table["log_or"]) / np.sum(w))
    se = float(np.sqrt(1 / np.sum(w)))
    z = mu / se
    p = float(2 * stats.norm.sf(abs(z)))
    Q = float(np.sum(w * (table["log_or"] - mu) ** 2))
    df = len(table) - 1
    I2 = max(0.0, (Q - df) / Q) if Q > 0 else 0.0
    return mu, se, p, Q, I2


def dl_random_effects(table: pd.DataFrame) -> dict:
    """DerSimonian-Laird random-effects summary (reported at k=10 per freeze spec)."""
    w_inv = 1 / table["se"] ** 2
    mu_fe, _, _, Q, _ = fe_summary(table)
    k = len(table)
    tau2 = (
        max(0.0, (Q - (k - 1)) / (np.sum(w_inv) - np.sum(w_inv**2) / np.sum(w_inv)))
        if Q > (k - 1)
        else 0.0
    )
    w_re = 1 / (table["se"] ** 2 + tau2)
    mu_re = float(np.sum(w_re * table["log_or"]) / np.sum(w_re))
    se_re = float(np.sqrt(1 / np.sum(w_re)))
    p_re = float(2 * stats.norm.sf(abs(mu_re / se_re)))
    return dict(tau2=float(tau2), mu=float(mu_re), se=se_re,
                or_=float(np.exp(mu_re)), p=p_re)


def loo_envelope(table: pd.DataFrame, by: str = "dataset_accession") -> tuple[float, float]:
    """Leave-one-study-out FE summary ORs -> (min, max) envelope (OR scale)."""
    mus = []
    for key in table[by]:
        sub = table[table[by] != key]
        mu_l, _, _, _, _ = fe_summary(sub)
        mus.append(float(np.exp(mu_l)))
    return float(min(mus)), float(max(mus))


def rq1_permutation(
    features: pd.DataFrame,
    background: pd.DataFrame,
    seed: int = SEED,
    b_perm: int = B_PERM,
    b_boot: int = B_BOOT,
    min_unique: int = 2,
) -> dict:
    """H1/RQ1: per-protein median peptide-pI divergence vs digestion-matched null.

    ``features`` needs (dataset_accession, mapped_accession, peptide_sequence,
    peptide_pi, precursor_pi); ``background`` needs (accession, peptide_pi).
    Sampling is WITH replacement from each protein's tryptic background
    (multinomial bootstrap); the global null resamples the SET of per-protein
    null medians. Returns observed statistic, two-sided global p, 95% null CI,
    and the per-protein table with BH q.

    Edge case (published): proteins with fewer background peptides than
    detected peptides are SKIPPED (documented absence, not corrected).
    """
    rng = np.random.default_rng(seed)
    per_protein = (
        features.dropna(subset=["precursor_pi"])
        .groupby(["dataset_accession", "mapped_accession"])
        .agg(n_unique=("peptide_sequence", "nunique"),
             median_peptide_pi=("peptide_pi", "median"),
             precursor_pi=("precursor_pi", "first"))
        .reset_index()
    )
    per_protein = per_protein[per_protein["n_unique"] >= min_unique]
    bg_by_accession = {acc: g["peptide_pi"].dropna().to_numpy() for acc, g in background.groupby("accession")}

    obs_divs, null_divs_global, per_protein_p = [], [], []
    for record in per_protein.itertuples(index=False):
        bg = bg_by_accession.get(record.mapped_accession)
        if bg is None or len(bg) < record.n_unique:
            continue
        obs = record.median_peptide_pi - record.precursor_pi
        obs_divs.append(obs)
        draws = bg[rng.integers(0, len(bg), size=(b_perm, record.n_unique))]
        null_medians = np.median(draws, axis=1)
        null_div = null_medians - record.precursor_pi
        null_divs_global.append(np.median(null_div))
        per_protein_p.append((1 + np.sum(np.abs(null_div) >= abs(obs))) / (b_perm + 1))

    obs_divs = np.array(obs_divs)
    T_obs = float(np.median(obs_divs))
    global_null = np.array([
        np.median(rng.choice(null_divs_global, size=len(null_divs_global), replace=True))
        for _ in range(b_boot)
    ])
    p_global = (
        (1 + np.sum(np.abs(global_null - np.median(null_divs_global)) >= abs(T_obs - np.median(null_divs_global))))
        / (b_boot + 1)
    )
    per_protein_valid = per_protein.iloc[: len(per_protein_p)].copy()
    per_protein_valid["perm_p"] = per_protein_p
    per_protein_valid["perm_q"] = bh_q(per_protein_p)
    return dict(
        t_obs=T_obs,
        p_global=float(p_global),
        ci_low=float(np.percentile(global_null, 2.5)),
        ci_high=float(np.percentile(global_null, 97.5)),
        null_center=float(np.median(null_divs_global)),
        n_proteins=int(len(obs_divs)),
        n_significant=int((per_protein_valid["perm_q"] < 0.05).sum()),
        per_protein=per_protein_valid,
    )


CLASSES = ("acidic", "near_neutral", "basic")


def bootstrap_class_shares(labels, n: int = 2000, seed: int = SEED) -> tuple[dict, dict]:
    """Bootstrap percentile CIs for class shares of a label array (S5)."""
    rng = np.random.default_rng(seed)
    labels = np.asarray(labels)
    n_total = len(labels)
    masks = {c: labels == c for c in CLASSES}
    shares = {c: float(masks[c].mean()) if n_total else 0.0 for c in CLASSES}
    cis = {c: [0.0, 0.0] for c in CLASSES}
    if n_total == 0:
        return shares, cis
    draws = np.empty((n, len(CLASSES)))
    for i in range(n):
        sample = labels[rng.integers(0, n_total, n_total)]
        draws[i] = [np.mean(sample == c) for c in CLASSES]
    for j, c in enumerate(CLASSES):
        cis[c] = [float(np.percentile(draws[:, j], 2.5)), float(np.percentile(draws[:, j], 97.5))]
    return shares, cis


def enrichment_with_ci(detected_labels, background_labels, n: int = 2000, seed: int = SEED) -> dict:
    """Enrichment (detected share / background share) with multinomial bootstrap CIs (S5)."""
    rng = np.random.default_rng(seed)
    d_arr = np.asarray(detected_labels)
    b_arr = np.asarray(background_labels)
    d_counts = np.array([(d_arr == c).sum() for c in CLASSES], dtype=float)
    b_counts = np.array([(b_arr == c).sum() for c in CLASSES], dtype=float)
    d_shares = d_counts / d_counts.sum() if d_counts.sum() else np.zeros(3)
    b_shares = b_counts / b_counts.sum() if b_counts.sum() else np.zeros(3)
    if d_counts.sum() == 0 or b_counts.sum() == 0:
        return {c: dict(enrichment=None, ci=[None, None]) for c in CLASSES}
    d_draws = rng.multinomial(int(d_counts.sum()), d_shares, size=n) / d_counts.sum()
    b_draws = rng.multinomial(int(b_counts.sum()), b_shares, size=n) / b_counts.sum()
    out = {}
    for j, c in enumerate(CLASSES):
        ratio = d_draws[:, j] / np.maximum(b_draws[:, j], 1e-12)
        out[c] = dict(
            enrichment=float(d_shares[j] / b_shares[j]) if b_shares[j] else None,
            ci=[float(np.percentile(ratio, 2.5)), float(np.percentile(ratio, 97.5))],
        )
    return out


def theoretical_charge_ph26(seq: str) -> int:
    """Dominant-charge proxy at acidic LC pH (~2.6): count(K)+count(R)+count(H)+1 (N-term).

    Upper-bound proxy; ESI charge-state distributions typically peak lower.
    """
    return sum(1 for aa in str(seq) if aa in "KRH") + 1
