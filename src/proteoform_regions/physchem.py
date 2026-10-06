"""pI / MW / net-charge engine with multiple pK models.

Primary calculator: Biopython ``IsoelectricPoint`` (Bjellqvist pK set, with
residue-specific terminal pKas). Sensitivity calculators: vectorized 60-iteration
bisection over Henderson-Hasselbalch charge with the pyteomics Rodwell/EMBOSS and
Lehninger pK tables.

Extraction note: the notebook version configured the Henderson-Hasselbalch model
through module globals (``hh_pi_configure``); here the pK table is an explicit
argument. Float operations are unchanged - equivalence against the original
global-state implementation is pinned by ``tests/test_physchem.py`` (reference
implementation comparison). Golden values from the publication's frozen ExPASy
benchmark (16 sequences) live in ``tests/data/physchem_golden.csv``.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence

import numpy as np
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from pyteomics import electrochem as ec

from .schema import classify_pi

__all__ = [
    "classify_pi",
    "biopython_pi",
    "biopython_mw",
    "hh_charge",
    "hh_pi_batch",
    "charge_at_ph",
    "count_matrix",
    "PK_TABLES",
    "gravy",
    "peptide_features",
    "fetch_expasy",
]

# ---------------------------------------------------------------- calculators

IONIZABLE = ["D", "E", "C", "Y", "H", "K", "R"]


def _table_pks(table) -> dict:
    """pyteomics pK table -> {'n_term', 'c_term', 'side': {res: pKa}} (verbatim helper)."""
    side = {res: table[res][0][0] for res in IONIZABLE if res in table}
    return {
        "n_term": table.get("H-", [("9.69", 1)])[0][0],
        "c_term": table.get("-OH", [("2.34", -1)])[0][0],
        "side": side,
    }


#: Henderson-Hasselbalch pK tables ('H-' = N-terminus, '-OH' = C-terminus).
#: Values resolved from pyteomics.electrochem at import (Rodwell: K=11.5, R=11.5,
#: D=3.86, C=8.33, Y=10.7, H=6.0, E=4.25, N=8.0, C-term=3.1; Lehninger: 10.53,
#: 12.48, 3.65, 8.18, 10.07, 6.0, 4.25, N=9.69, C-term=2.34). Pinned by tests.
PK_TABLES = {
    "hh_rodwell_emboss_lineage": _table_pks(ec.pK_rodwell),
    "hh_lehninger": _table_pks(ec.pK_lehninger),
}


def count_matrix(sequences: Sequence[str]) -> dict[str, np.ndarray]:
    """Per-ionizable-residue count vectors for a list of sequences (verbatim)."""
    seqs = [str(s).upper() for s in sequences]
    return {res: np.array([seq.count(res) for seq in seqs], dtype=float) for res in IONIZABLE}


def hh_charge(ph: float, counts: Mapping[str, np.ndarray], table: str | dict = "hh_rodwell_emboss_lineage"):
    """Net charge at `ph` from residue-count vectors under a pK table.

    `table` may be a PK_TABLES key or a {'n_term','c_term','side'} dict
    (the notebook's global-state configuration becomes an argument).
    """
    pks = PK_TABLES[table] if isinstance(table, str) else table
    p_k_n, p_k_c, side = pks["n_term"], pks["c_term"], pks["side"]
    charge = 1.0 / (1.0 + 10 ** (ph - p_k_n)) + -1.0 / (1.0 + 10 ** (p_k_c - ph))
    for res, p_k in side.items():
        n = counts[res]
        if p_k is None or not np.any(n):
            continue
        if res in ("D", "E", "C", "Y"):
            charge += -n / (1.0 + 10 ** (p_k - ph))
        else:
            charge += n / (1.0 + 10 ** (ph - p_k))
    return charge


def hh_pi_batch(
    sequences: Sequence[str],
    table: str = "hh_rodwell_emboss_lineage",
    iterations: int = 60,
) -> np.ndarray:
    """Vectorized bisection pI for all sequences with the given pK table.

    60 iterations of bisection on [0, 14]; convergence implicit (~14/2**60).
    Sign convention: pI is where charge crosses from positive (low pH) to
    non-positive (high pH).
    """
    pks = PK_TABLES[table] if isinstance(table, str) else table
    counts = count_matrix(sequences)
    lo = np.full(len(sequences), 0.0)
    hi = np.full(len(sequences), 14.0)
    for _ in range(iterations):
        mid = (lo + hi) / 2.0
        positive = hh_charge(mid, counts, pks) > 0
        lo = np.where(positive, mid, lo)
        hi = np.where(positive, hi, mid)
    return (lo + hi) / 2.0


def charge_at_ph(
    sequences: Sequence[str],
    ph: float = 7.4,
    table: str = "hh_rodwell_emboss_lineage",
) -> dict[str, float]:
    """Net charge per unique sequence at `ph` (paper default: pH 7.4, Rodwell lineage)."""
    uniq = list(dict.fromkeys(str(s) for s in sequences))
    pks = PK_TABLES[table] if isinstance(table, str) else table
    counts = count_matrix(uniq)
    charges = hh_charge(ph, counts, pks)
    return dict(zip(uniq, (float(c) for c in charges), strict=True))


def biopython_pi(sequence: str) -> float:
    """Primary pI calculator: Biopython IsoelectricPoint (Bjellqvist; X->A)."""
    return IsoelectricPoint(str(sequence).replace("X", "A")).pi()


def biopython_mw(sequence: str) -> float:
    """Molecular weight: Biopython average-isotopic mass (X->A replacement before calculation)."""
    return ProteinAnalysis(str(sequence).replace("X", "A")).molecular_weight()


# ---------------------------------------------------------------- composition features

_GRAVY_TABLE: Mapping[str, float] = ec.hydropathicity_KD


def gravy(seq: str) -> float | None:
    """Mean Kyte-Doolittle hydropathicity over recognized residues (None if empty)."""
    values = [_GRAVY_TABLE[aa] for aa in str(seq) if aa in _GRAVY_TABLE]
    return float(np.mean(values)) if values else None


def peptide_features(seq: str) -> dict:
    """Composition features (verbatim): length, acidic/basic fractions, Cys/Met,
    tryptic-terminus flag, missed-cleavage count from internal KR|not-P sites."""
    seq = str(seq)
    n = len(seq)
    acidic = sum(seq.count(r) for r in "DE")
    basic = sum(seq.count(r) for r in "KRH")
    return {
        "peptide_length": n,
        "acidic_frac": acidic / n,
        "basic_frac": basic / n,
        "cys_count": seq.count("C"),
        "met_count": seq.count("M"),
        "tryptic_terminus": seq[-1] in "KR",
        "missed_cleavages": sum(1 for i, aa in enumerate(seq[:-1]) if aa in "KR" and seq[i + 1] != "P"),
    }


# ---------------------------------------------------------------- live benchmark (network)


def fetch_expasy(sequence: str, timeout: float = 60.0, retries: int = 2):
    """Theoretical pI/MW from live ExPASy ProtParam (network; benchmark use only).

    Returns (pi, mw) floats or (None, None) on failure. Marked network-reliant:
    tests that use it are skipped in offline CI.
    """
    import time

    import requests

    url = "https://web.expasy.org/cgi-bin/protparam/protparam"
    for attempt in range(1, retries + 1):
        try:
            response = requests.get(
                url,
                params={"sequence": str(sequence)},
                timeout=timeout,
                headers={"User-Agent": "proteoform-regions benchmark (validation)"},
            )
            response.raise_for_status()
            import re

            text = re.sub(r"<[^>]+>", " ", response.text)
            p_i = re.search(r"Theoretical pI[/:]?\s*([\d.]+)", text)
            mw = re.search(r"Molecular weight[/:]?\s*([\d.]+)", text)
            return (
                float(p_i.group(1)) if p_i else None,
                float(mw.group(1)) if mw else None,
            )
        except Exception:
            if attempt == retries:
                return None, None
            time.sleep(2.0)


def pi_for_calculator(sequences: Iterable[str], calculator: str = "bjellqvist") -> dict[str, float]:
    """pI per sequence under a named calculator: 'bjellqvist' (Biopython primary)
    or any PK_TABLES key ('hh_rodwell_emboss_lineage', 'hh_lehninger')."""
    uniq = list(dict.fromkeys(str(s) for s in sequences))
    if calculator == "bjellqvist":
        return {s: biopython_pi(s) for s in uniq}
    values = hh_pi_batch(uniq, table=calculator)
    return dict(zip(uniq, (float(v) for v in values), strict=True))
