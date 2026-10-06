"""Region-stitching tests: gap boundary semantics, the D1 last-member-end
behavior, transitivity, terminal flags, and confidence gating.

These tests are the executable specification of Supplementary Methods Module 3.
"""

import pytest

from proteoform_regions import regions
from tests import reference_implementations as ref


def span(start, end, seq="PEPTIDE", row_id=0):
    """A paper-shaped span tuple: (start, end, sequence, peptide_pi, row_id)."""
    return (start, end, seq, 6.5, row_id)


class TestGapBoundary:
    """gap=25: join iff next.start - last.end <= 25 (i.e. <= 24 unobserved residues)."""

    def test_gap_of_25_joins(self):
        clusters = regions.stitch_spans([span(10, 20), span(45, 50)])
        # 45 - 20 = 25 <= 25 -> one region
        assert len(clusters) == 1

    def test_gap_of_26_splits(self):
        clusters = regions.stitch_spans([span(10, 20), span(46, 50)])
        # 46 - 20 = 26 > 25 -> two regions
        assert len(clusters) == 2

    def test_zero_gap_touching_joins(self):
        clusters = regions.stitch_spans([span(10, 20), span(20, 30)])
        assert len(clusters) == 1

    def test_overlap_joins(self):
        clusters = regions.stitch_spans([span(10, 30), span(25, 40)])
        assert len(clusters) == 1

    def test_parameterized_gap(self):
        assert len(regions.stitch_spans([span(10, 20), span(46, 50)], gap=26)) == 1
        assert len(regions.stitch_spans([span(10, 20), span(40, 50)], gap=10)) == 2


class TestLastMemberEndSemantics:
    """D1: the gap check references the LAST-APPENDED member's end, not the
    cluster's running maximum end. A nested span (later start, smaller end)
    becomes the reference and can SPLIT a span that max-end single linkage
    would join."""

    def test_nested_span_restricts_where_max_end_would_join(self):
        # sorted: A(10,40), B(30,44), C(41,42), D(70,80)
        # B joins A (30-40<=25); C joins B (41-44<=25); D: 70 - C.end(42) = 28 > 25 -> SPLIT
        clusters = regions.stitch_spans([span(10, 40), span(30, 44), span(41, 42), span(70, 80)])
        assert [len(c) for c in clusters] == [3, 1]
        assert regions.region_bounds(clusters[0]) == (10, 44)
        assert clusters[1][0][:2] == (70, 80)

    def test_running_max_equivalent_case_still_joins(self):
        # without the nested span, D would join the running max end 44: 70-44=26>25 -> split anyway;
        # use E(60,64): 60-44=16<=25 joins
        clusters = regions.stitch_spans([span(10, 40), span(30, 44), span(60, 64)])
        assert len(clusters) == 1

    def test_extending_member_becomes_reference(self):
        # B extends the cluster end; C is checked against B.end (the new max)
        clusters = regions.stitch_spans([span(10, 40), span(30, 60), span(85, 90)])
        # 85 - 60 = 25 <= 25 -> joins
        assert len(clusters) == 1


class TestTransitivityAndDeterminism:
    def test_chained_transitivity(self):
        # consecutive gaps of 10 each join; total span 10..80 far exceeds 25
        clusters = regions.stitch_spans([span(10, 20), span(30, 40), span(50, 60), span(70, 80)])
        assert len(clusters) == 1
        assert regions.region_bounds(clusters[0]) == (10, 80)

    def test_input_order_does_not_matter(self):
        ordered = [span(10, 20), span(45, 50), span(100, 110)]
        shuffled = [span(100, 110), span(45, 50), span(10, 20)]
        assert regions.stitch_spans(ordered) == regions.stitch_spans(shuffled)

    def test_parity_with_notebook_walk(self):
        import random

        rng = random.Random(42)
        spans = [(rng.randint(1, 300), 0, f"SEQ{i}", 6.0, i) for i in range(60)]
        spans = [(s, s + rng.randint(6, 30), seq, pi, rid) for s, _, seq, pi, rid in spans]
        assert regions.stitch_spans(spans) == ref.stitch_walk(spans)

    def test_single_span_is_its_own_cluster(self):
        clusters = regions.stitch_spans([span(5, 12)])
        assert clusters == [[(5, 12, "PEPTIDE", 6.5, 0)]]


class TestTerminalFlags:
    def test_n_and_c_terminal(self):
        flags = regions.terminal_flags(5, 95, 100)
        assert flags == (True, True)

    def test_interior_region(self):
        assert regions.terminal_flags(30, 40, 100) == (False, False)

    def test_exact_boundary_inclusive(self):
        assert regions.terminal_flags(10, 90, 100) == (True, True)
        assert regions.terminal_flags(11, 89, 100) == (False, False)


class TestConfidenceLabel:
    def test_single_sequence_peptide_only(self):
        assert regions.confidence_label(1) == ("peptide_only", False)

    def test_two_sequences_region_supported(self):
        assert regions.confidence_label(2) == ("region_supported", False)

    def test_numeric_window_consistent_upgrades(self):
        label, applied = regions.confidence_label(2, region_mw=15000.0, constraint_window=(10000.0, 20000.0))
        assert (label, applied) == ("fragment_constrained", True)

    def test_out_of_window_caps_to_peptide_only(self):
        label, applied = regions.confidence_label(3, region_mw=5000.0, constraint_window=(10000.0, 20000.0))
        assert (label, applied) == ("peptide_only", True)

    def test_window_without_mw_is_region_supported(self):
        # MW unknown + window present -> cannot gate; stays region_supported (conservative
        # upward direction only; documented)
        assert regions.confidence_label(2, constraint_window=(1.0, 2.0)) == ("region_supported", False)


class TestRegionRecord:
    PRECURSOR = "M" + "A" * 49 + "K" + "D" * 50  # length 102

    def test_record_shape_and_values(self):
        cluster = [span(52, 55, seq="KDDD", row_id=1), span(60, 63, seq="DDDD", row_id=2)]
        rec = regions.region_record(
            study="PXDTEST",
            accession="P00000",
            cluster=cluster,
            precursor_length=len(self.PRECURSOR),
            precursor_pi=5.0,
            precursor_sequence=self.PRECURSOR,
            gene_name="TEST",
        )
        assert rec["region_start"] == 52 and rec["region_end"] == 63
        assert rec["region_n_peptides"] == 2
        assert rec["confidence_label"] == "region_supported"
        assert rec["entity_level"] == "detected_region"
        assert rec["region_pi_divergence"] == pytest.approx(rec["region_pi"] - 5.0)
        assert rec["member_sequences"] == "DDDD;KDDD"

    def test_region_sequence_slice_is_contiguous_with_unobserved_residues(self):
        # 1-based inclusive: positions 12..16 of M A*10 K D*4 A*8 = K D D D D
        seq = regions.region_sequence("MAAAAAAAAAAKDDDDAAAAAAAA", 12, 16)
        assert seq == "KDDDD"
        assert regions.region_sequence("", 1, 5) is None

    def test_region_physchem_none_and_real(self):
        assert regions.region_physchem(None) == (None, None)
        pi, mw = regions.region_physchem("KDDD")
        assert 3.0 < pi < 7.0 and mw > 0
