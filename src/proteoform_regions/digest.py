"""In-silico tryptic digestion and digestion-matched background construction.

Rules (Supplementary Methods Module 2 of the source publication):
- Keil cleavage rule: cut after K or R when the next residue is not P.
- Length filter 7 <= L <= 60 (parameterizable; paper values are the defaults).
- Background dedup is per (accession, peptide) pair - the same peptide from
  different proteins is a separate entry.
- Initiator Met: no special handling (digestion starts at canonical position 1).

**Exact enumeration semantics (verified against the frozen publication
background; divergence row #31).** With ``max_missed=2`` the published loop
enumerates peptides spanning up to THREE internal cut sites (missed_cleavages
in {0,1,2,3}), not two: the frozen ``processed_tryptic_background.csv``
contains 271,317 mc=3 rows (~28%). The Supplementary Methods statement
"missed cleavages 0, 1, 2" refers to the downstream ANALYSIS STRATA (mc0 /
mc1 / mc2+), into which mc=3 rows fold, not to the enumeration bound. This
package preserves the published enumeration verbatim; passing
``max_missed=1`` reproduces a strict 0..2 enumeration if that is intended.
"""

from __future__ import annotations

from collections.abc import Iterator

__all__ = ["tryptic_digest", "digest_protein", "background_peptides"]

MIN_LENGTH = 7
MAX_LENGTH = 60
MAX_MISSED = 2


def tryptic_digest(
    sequence: str,
    max_missed: int = MAX_MISSED,
    min_len: int = MIN_LENGTH,
    max_len: int = MAX_LENGTH,
) -> list[tuple[str, int]]:
    """Digest `sequence` in silico; returns [(peptide, missed_cleavages), ...].

    Verbatim from the publication's S3 engine, with the paper's inline length
    literals (7, 60) exposed as defaulted parameters. NOTE: the enumeration
    admits peptides with up to ``max_missed + 1`` internal cut sites (see
    module docstring) - preserved exactly to match the frozen background.
    """
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
            if min_len <= len(peptide) <= max_len:
                peptides.append((peptide, end_index - start_index - 1))
    return peptides


def digest_protein(
    accession: str,
    sequence: str,
    max_missed: int = MAX_MISSED,
    min_len: int = MIN_LENGTH,
    max_len: int = MAX_LENGTH,
) -> Iterator[dict]:
    """Canonical background rows for one protein: deduped per (accession, peptide)."""
    seen: set[tuple[str, str]] = set()
    for peptide, missed in tryptic_digest(sequence, max_missed, min_len, max_len):
        key = (accession, peptide)
        if key in seen:
            continue
        seen.add(key)
        yield {
            "accession": accession,
            "peptide_sequence": peptide,
            "missed_cleavages": missed,
            "peptide_length": len(peptide),
        }


def background_peptides(
    protein_sequences: dict[str, str],
    max_missed: int = MAX_MISSED,
    min_len: int = MIN_LENGTH,
    max_len: int = MAX_LENGTH,
) -> list[dict]:
    """Digestion-matched tryptic background over {accession: sequence}.

    Same-protocol digest as the detected peptides' parent proteins (RQ2
    comparator). Streams per protein; caller decides chunking/persistence.
    """
    rows: list[dict] = []
    for accession in sorted(protein_sequences):
        sequence = protein_sequences[accession]
        if not sequence:
            continue
        rows.extend(digest_protein(accession, sequence, max_missed, min_len, max_len))
    return rows
