"""Verbatim reference implementations from the source publication's notebooks.

These are COPIES OF THE ORIGINAL notebook code (global state and all), used by
the test suite to prove that the package's parameterized refactors are
numerically identical to the published implementations. Do not "improve" them:
their value is being byte-faithful to the ground truth.

Sources (paper repo, read-only):
- physicochemical_pipeline_validation.ipynb  -> hh_charge / hh_pi_configure /
  hh_pi_batch / tryptic_digest / peptide_features / gravy / classify_pi
- proteoform_region_inference.ipynb          -> stitching walk
"""

import numpy as np
from pyteomics import electrochem as ec

# ------------------------------------------------------------------ physchem (S3)

IONIZABLE = ["D", "E", "C", "Y", "H", "K", "R"]


def _table_pks(table):
    side = {res: table[res][0][0] for res in IONIZABLE if res in table}
    return {
        "n_term": table.get("H-", [("9.69", 1)])[0][0],
        "c_term": table.get("-OH", [("2.34", -1)])[0][0],
        "side": side,
    }


HH_TABLES = {
    "hh_rodwell_emboss_lineage": _table_pks(ec.pK_rodwell),
    "hh_lehninger": _table_pks(ec.pK_lehninger),
}


def count_matrix(sequences):
    seqs = [str(s).upper() for s in sequences]
    counts = {res: np.array([seq.count(res) for seq in seqs], dtype=float) for res in IONIZABLE}
    return counts


def hh_charge(ph, counts):
    charge = 1.0 / (1.0 + 10 ** (ph - _PK_N)) + -1.0 / (1.0 + 10 ** (_PK_C - ph))
    for res, p_k in _SIDE_PKS.items():
        n = counts[res]
        if p_k is None or not np.any(n):
            continue
        if res in ("D", "E", "C", "Y"):
            charge += -n / (1.0 + 10 ** (p_k - ph))
        else:
            charge += n / (1.0 + 10 ** (ph - p_k))
    return charge


_CALC_STATE = None


def hh_pi_configure(table_name):
    global _CALC_STATE, _PK_N, _PK_C, _SIDE_PKS
    table = HH_TABLES[table_name]
    _CALC_STATE = table_name
    _PK_N, _PK_C = table["n_term"], table["c_term"]
    _SIDE_PKS = table["side"]


def hh_pi_batch(sequences):
    """Vectorized bisection pI for all sequences with the configured table."""
    counts = count_matrix(sequences)
    lo = np.full(len(sequences), 0.0)
    hi = np.full(len(sequences), 14.0)
    for _ in range(60):
        mid = (lo + hi) / 2.0
        positive = hh_charge(mid, counts) > 0
        lo = np.where(positive, mid, lo)
        hi = np.where(positive, hi, mid)
    return (lo + hi) / 2.0


def tryptic_digest(sequence, max_missed=2):
    peptides, cuts = [], []
    for i, aa in enumerate(sequence[:-1]):
        if aa in "KR" and sequence[i + 1] != "P":
            cuts.append(i + 1)
    cuts.append(len(sequence))
    for start_index, start in enumerate(cuts):
        for end_index in range(start_index + 1, min(start_index + 2 + max_missed, len(cuts)) + 1):
            end = cuts[end_index] if end_index < len(cuts) else None
            if end is None:
                continue
            peptide = sequence[start:end]
            if 7 <= len(peptide) <= 60:
                peptides.append((peptide, end_index - start_index - 1))
    return peptides


# ------------------------------------------------------------------ stitching (S4)

REGION_GAP_AA = 25


def stitch_walk(spans, gap=REGION_GAP_AA):
    """The notebook's exact cluster walk over start-sorted spans."""
    spans = sorted(spans)
    clusters, current = [], [spans[0]]
    for span in spans[1:]:
        if span[0] - current[-1][1] <= gap:
            current.append(span)
        else:
            clusters.append(current)
            current = [span]
    clusters.append(current)
    return clusters
