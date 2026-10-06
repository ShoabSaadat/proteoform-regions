"""Single-linkage detected-region stitching from peptide evidence.

Algorithm (Supplementary Methods Module 3 of the source publication):

- Sort peptide spans by start position (tuple sort: start, end, then payload).
- Walk spans in order; join a span to the current cluster when
  ``next.start - last_member.end <= gap`` (paper default ``gap = 25``).
- Region bounds are ``min(member starts)`` .. ``max(member ends)``; the region
  sequence is the contiguous precursor slice between the bounds (unobserved
  intermediate residues included).
- A region with fewer than 2 distinct member peptide sequences stays
  ``peptide_only``; >= 2 is ``region_supported``; a numerically consistent MW
  constraint window upgrades to ``fragment_constrained`` (ANL-020).

**Exact-semantics note (pinned by tests).** The gap check uses the END OF THE
LAST-APPENDED MEMBER, not the running maximum end of the cluster. Because spans
are sorted by start, a later-starting nested span (smaller end) becomes the
reference end for the next comparison, which can SPLIT a span that textbook
running-max single linkage would join. This is the published algorithm's exact
behavior ("D1 last-member-end"); it is intentionally preserved, not normalized.
The join condition ``next.start - last.end <= 25`` corresponds to at most 24
unobserved residues between the last member and the next span.
"""

from __future__ import annotations

from collections.abc import Sequence

from .physchem import biopython_mw, biopython_pi
from .schema import classify_pi

__all__ = [
    "REGION_GAP_AA",
    "TERMINAL_FRACTION",
    "MIN_REGION_PEPTIDES",
    "stitch_spans",
    "region_bounds",
    "terminal_flags",
    "confidence_label",
    "region_sequence",
    "region_physchem",
    "region_record",
]

#: Peptides whose spans are separated by <= this many residues join one region.
REGION_GAP_AA = 25
#: Region touching the first/last 10% of the precursor counts as N/C-terminal-inclusive.
TERMINAL_FRACTION = 0.10
#: region_supported threshold (>= 2 distinct member sequences; ANL-002).
MIN_REGION_PEPTIDES = 2


def stitch_spans(spans: Sequence[tuple], gap: int = REGION_GAP_AA) -> list[list[tuple]]:
    """Cluster sorted-by-start spans into regions via the paper's exact loop.

    Each span is a tuple whose first two elements are (start, end); remaining
    elements are carried through untouched (the paper carries
    (sequence, peptide_pi, row_id)). Tuples are sorted exactly as the notebook
    does (lexicographic), so input order does not matter - determinism is total.
    """
    spans = sorted(spans)
    clusters: list[list[tuple]] = [[spans[0]]]
    current = clusters[0]
    for span in spans[1:]:
        if span[0] - current[-1][1] <= gap:
            current.append(span)
        else:
            current = [span]
            clusters.append(current)
    return clusters


def region_bounds(cluster: Sequence[tuple]) -> tuple[int, int]:
    """(min member start, max member end) -> 1-based inclusive region bounds."""
    return int(min(s[0] for s in cluster)), int(max(s[1] for s in cluster))


def terminal_flags(
    region_start: int,
    region_end: int,
    precursor_length: int,
    terminal_fraction: float = TERMINAL_FRACTION,
) -> tuple[bool, bool]:
    """(n_terminal_inclusive, c_terminal_inclusive) per the 10% rule."""
    return (
        region_start <= terminal_fraction * precursor_length,
        region_end >= (1 - terminal_fraction) * precursor_length,
    )


def confidence_label(
    n_unique_sequences: int,
    region_mw: float | None = None,
    constraint_window: tuple[float, float] | None = None,
    min_peptides: int = MIN_REGION_PEPTIDES,
) -> tuple[str, bool]:
    """Tiered region confidence with ANL-020 MW-window gating.

    Returns (label, mw_constraint_applied):

    - ``peptide_only``          : fewer than `min_peptides` distinct sequences, OR an
                                  available numeric window that the region MW violates
                                  (an out-of-window region cannot exceed the evidence
                                  of its peptides).
    - ``region_supported``      : >= `min_peptides` sequences and no numeric window.
    - ``fragment_constrained``  : >= `min_peptides` sequences and MW inside the window.

    The paper's pilots exposed no numeric windows, so every pilot region is
    peptide_only/region_supported with ``mw_constraint_applied=False`` - the
    unreachable path is exercised only by tests.
    """
    if n_unique_sequences < min_peptides:
        return "peptide_only", False
    if constraint_window is not None and region_mw is not None:
        low, high = constraint_window
        if not (low <= region_mw <= high):
            return "peptide_only", True
        return "fragment_constrained", True
    return "region_supported", False


def region_sequence(precursor_sequence: str, region_start: int, region_end: int) -> str | None:
    """Contiguous precursor slice covering the region (1-based inclusive bounds)."""
    if not precursor_sequence:
        return None
    return precursor_sequence[region_start - 1 : region_end]


def region_physchem(sequence: str | None) -> tuple[float | None, float | None]:
    """(pI, MW) of a region sequence (X->A; None inputs -> (None, None))."""
    if not sequence:
        return None, None
    cleaned = sequence.replace("X", "A")
    return biopython_pi(cleaned), biopython_mw(cleaned)


def region_record(
    study: str,
    accession: str,
    cluster: Sequence[tuple],
    precursor_length: int,
    precursor_pi: float | None = None,
    precursor_mw: float | None = None,
    gene_name: str | None = None,
    precursor_sequence: str = "",
    gap: int = REGION_GAP_AA,
    terminal_fraction: float = TERMINAL_FRACTION,
    min_peptides: int = MIN_REGION_PEPTIDES,
    constraint_window: tuple[float, float] | None = None,
    extra_flags: str | None = None,
) -> dict:
    """Build one canonical region row (S4 column set) from a stitched cluster.

    `cluster` spans are (start, end, sequence, peptide_pi, row_id) tuples, as
    produced by the paper's S4 walk. `extra_flags` carries RQ3 UniProt-feature
    flags when a feature index is supplied upstream.
    """
    region_start, region_end = region_bounds(cluster)
    n_peptides = len({s[2] for s in cluster})
    span_length = region_end - region_start + 1
    seq = region_sequence(precursor_sequence, region_start, region_end)
    region_pi, region_mw = region_physchem(seq)
    n_term, c_term = terminal_flags(region_start, region_end, precursor_length, terminal_fraction)
    label, constraint_applied = confidence_label(n_peptides, region_mw, constraint_window, min_peptides)
    return {
        "dataset_accession": study,
        "mapped_accession": accession,
        "gene_name": gene_name,
        "region_start": region_start,
        "region_end": region_end,
        "region_span_aa": span_length,
        "region_n_peptides": n_peptides,
        "region_unique_sequences": n_peptides,
        "member_peptide_rows": ";".join(str(int(s[4])) for s in cluster),
        "member_sequences": ";".join(sorted({s[2] for s in cluster})),
        "precursor_length": int(precursor_length),
        "coverage_fraction": span_length / precursor_length,
        "n_terminal_inclusive": n_term,
        "c_terminal_inclusive": c_term,
        "region_pi": region_pi,
        "region_mw_da": region_mw,
        "precursor_pi": precursor_pi,
        "precursor_mw_da": precursor_mw,
        "region_pi_divergence": (region_pi - precursor_pi) if region_pi is not None else None,
        "region_pi_class_anl003": classify_pi(region_pi),
        "precursor_pi_class_anl003": classify_pi(precursor_pi),
        "class_flip_region_vs_precursor": (
            classify_pi(region_pi) != classify_pi(precursor_pi)
            if region_pi is not None and precursor_pi is not None
            else None
        ),
        "confidence_label": label,
        "mw_constraint_applied": constraint_applied,
        "rq3_class_flags": extra_flags,
        "entity_level": "detected_region",
    }
