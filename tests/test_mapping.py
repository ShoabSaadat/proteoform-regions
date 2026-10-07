"""UniProt mapping tests (S2 extraction) on synthetic entries - no network.

The chunked REST fetcher itself is exercised only by ``network``-marked
tests; everything here runs against snapshot files / synthetic entries.
"""

import json

import numpy as np
import pandas as pd
import pytest

from proteoform_regions import mapping


def make_entry(
    accession: str,
    sequence: str,
    gene: str = "TEST",
    name: str = "Test protein",
    entry_type: str = "UniProtKB reviewed (Swiss-Prot)",
) -> dict:
    return {
        "primaryAccession": accession,
        "sequence": {"value": sequence, "molWeight": 1000.0 * len(sequence)},
        "proteinDescription": {
            "recommendedName": {"fullName": {"value": name}},
        },
        "genes": [{"geneName": {"value": gene}}],
        "entryType": entry_type,
        "organism": {"scientificName": "Homo sapiens"},
        "entryAudit": {"sequenceVersion": 3},
    }


PROT_A = "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVK"
PROT_B = "MARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRK"


class TestNormalizeAccessions:
    def test_fasta_prefix_and_group(self):
        assert mapping.normalize_accessions("sp|P02768|ALBU_HUMAN;sp|P69999|CAH1_HUMAN") == [
            ("P02768", None),
            ("P69999", None),
        ]

    def test_comma_separator(self):
        assert mapping.normalize_accessions("P02768,P69999") == [("P02768", None), ("P69999", None)]

    def test_isoform_suffix_flagged(self):
        assert mapping.normalize_accessions("P02768-2") == [("P02768", "-2")]

    def test_isoform_only_for_pqoa_prefixes(self):
        # A0A0... isoform-like suffix kept; but e.g. 'ABC-2' is not an accession isoform
        assert mapping.normalize_accessions("Q9XYZ-1") == [("Q9XYZ", "-1")]

    def test_contaminant_and_decoy_dropped(self):
        assert mapping.normalize_accessions("CON__P12345;REV__P99999;P02768") == [("P02768", None)]

    def test_nan_and_none(self):
        assert mapping.normalize_accessions(None) == []
        assert mapping.normalize_accessions(float("nan")) == []

    def test_empty_tokens_skipped(self):
        assert mapping.normalize_accessions(" P02768 ;; ") == [("P02768", None)]


class TestEntryProperties:
    def test_recommended_name_and_genes(self):
        props = mapping.entry_properties(make_entry("P02768", PROT_A))
        assert props["sequence"] == PROT_A
        assert props["sequence_length"] == len(PROT_A)
        assert props["gene_name"] == "TEST"
        assert props["protein_name"] == "Test protein"
        # VERBATIM-BUG PIN (paper divergence #32): the publication's
        # endswith("Swiss-Prot") check is always False because real entryType is
        # "UniProtKB reviewed (Swiss-Prot)" (trailing paren); the frozen properties
        # table carries reviewed=False for all 5,777 proteins. Preserved for parity.
        assert props["reviewed"] is False
        assert props["sequence_version"] == 3

    def test_submission_name_fallback(self):
        e = make_entry("P02768", PROT_A, name=None)
        e["proteinDescription"] = {"submissionNames": [{"fullName": {"value": "Sub name"}}]}
        assert mapping.entry_properties(e)["protein_name"] == "Sub name"


class TestContaminantFlag:
    @pytest.mark.parametrize(
        ("accession", "name", "expected"),
        [
            ("P00760", "Trypsin-1", "TRYPSIN"),
            ("P35900", "Keratin, type II cytoskeletal 2 epidermal", "KERATIN"),
            ("P02769", "Bovine serum albumin", "BSA"),
            ("P01857", "Immunoglobulin heavy constant gamma 1", "IG_CHAIN"),
            ("P02768", "Albumin (human)", None),
        ],
    )
    def test_rules(self, accession, name, expected):
        assert mapping.contaminant_flag(accession, name) == expected


class TestMapPeptides:
    ENTRIES = {
        "P02768": make_entry("P02768", PROT_A, gene="ALB"),
        "P69999": make_entry("P69999", PROT_B.replace("ALPAPIEK", ""), gene="CAH1"),
        "P12345": make_entry("P12345", PROT_B, gene="H2A", name="Trypsin-like protein"),
    }

    def _evidence(self, rows):
        return pd.DataFrame(
            rows,
            columns=["peptide_row_id", "dataset_accession", "peptide_sequence", "protein_accession_raw"],
        )

    def test_exact_unique(self):
        ev = self._evidence([(0, "PXDTEST", "DAHKSEVAHR", "P02768")])
        m = mapping.map_peptides(ev, self.ENTRIES)
        assert m.loc[0, "mapping_tier"] == "exact_unique"
        assert m.loc[0, "mapped_accession"] == "P02768"
        assert m.loc[0, "positions"] == str(PROT_A.index("DAHKSEVAHR") + 1)

    def test_multi_position(self):
        seq = "ALPAPIEKGGGALPAPIEKR"
        entries = {"P1": make_entry("P1", seq)}
        ev = self._evidence([(0, "PXDTEST", "ALPAPIEK", "P1")])
        m = mapping.map_peptides(ev, entries)
        assert m.loc[0, "mapping_tier"] == "exact_multi_position"
        assert m.loc[0, "positions"] == "1;12"  # 8 aa + GGG -> 2nd copy at 1-based 12

    def test_failure_tiers(self):
        ev = self._evidence(
            [
                (0, "PXDTEST", "DAHKSEVAHR", None),  # no accession in row
                (1, "PXDTEST", "DAHKSEVAHR", "P00000"),  # accession not in UniProt
                (2, "PXDTEST", "WWWWWWWWWK", "P02768"),  # not a substring
            ]
        )
        m = mapping.map_peptides(ev, self.ENTRIES)
        assert m.loc[0, "mapping_tier"] == "no_accession_in_row"
        assert m.loc[1, "mapping_tier"] == "accession_unmapped_in_uniprot"
        assert m.loc[2, "mapping_tier"] == "not_found_in_group_sequences"

    def test_first_group_member_with_hit_wins(self):
        # group order preserved; first token whose sequence contains the peptide maps
        seq_a = "MKWVTFISLLFLFSSAYSR"
        entries = {"P1": make_entry("P1", "GGGGGGGGGG"), "P2": make_entry("P2", seq_a)}
        ev = self._evidence([(0, "PXDTEST", "MKWVTFISLL", "P1;P2")])
        m = mapping.map_peptides(ev, entries)
        assert m.loc[0, "mapped_accession"] == "P2"

    def test_degenerate_peptide_flag(self):
        entries = {
            "P1": make_entry("P1", "ALPAPIEKGG"),
            "P2": make_entry("P2", "GGALPAPIEK"),
        }
        ev = self._evidence(
            [
                (0, "PXDTEST", "ALPAPIEK", "P1"),
                (1, "PXDTEST", "ALPAPIEK", "P2"),
            ]
        )
        m = mapping.map_peptides(ev, entries)
        assert m["degenerate_peptide"].all()

    def test_contaminant_flag_propagates(self):
        ev = self._evidence([(0, "PXDTEST", "ARTKQTARK", "P12345")])  # trypsin-named protein
        m = mapping.map_peptides(ev, self.ENTRIES)
        assert m.loc[0, "contaminant_flag"] == "TRYPSIN"

    def test_isoform_row_flag(self):
        ev = self._evidence([(0, "PXDTEST", "DAHKSEVAHR", "P02768-2")])
        m = mapping.map_peptides(ev, self.ENTRIES)
        assert bool(m.loc[0, "row_had_isoform_suffix"]) is True
        assert m.loc[0, "mapping_tier"] == "exact_unique"


class TestSnapshotAndEnsure:
    def test_snapshot_roundtrip(self, tmp_path):
        entries = [make_entry("P02768", PROT_A)]
        (tmp_path / "uniprot_snapshot_20260101.json").write_text(json.dumps(entries))
        loaded = mapping.load_snapshot(tmp_path, "20260101")
        assert set(loaded) == {"P02768"}
        assert mapping.load_snapshot(tmp_path, "19990101") == {}

    def test_ensure_offline_uses_snapshot_only(self, tmp_path):
        entries = [make_entry("P02768", PROT_A), make_entry("P69999", PROT_B)]
        (tmp_path / "uniprot_snapshot_20260101.json").write_text(json.dumps(entries))
        got = mapping.ensure_uniprot_entries(
            {"P02768", "P69999", "P00000"}, tmp_path, "20260101", offline=True
        )
        assert set(got) == {"P02768", "P69999"}  # missing accession documented by absence


class TestProteinPropertiesTable:
    def test_columns_and_precursor_physchem(self):
        entries = {"P02768": make_entry("P02768", PROT_A, gene="ALB")}
        proteins = mapping.protein_properties_table(entries)
        row = proteins.iloc[0]
        assert row["accession"] == "P02768"
        assert row["sequence_length"] == len(PROT_A)
        assert row["precursor_pi"] == pytest.approx(
            float(
                __import__("Bio.SeqUtils.IsoelectricPoint", fromlist=["IsoelectricPoint"])
                .IsoelectricPoint(PROT_A)
                .pi()
            )
        )
        assert np.isfinite(row["precursor_mw_da"])
        assert "sequence" in proteins.columns  # kept in-memory, dropped only on write


class TestMapToUniprot:
    def test_end_to_end_offline(self, tmp_path):
        entries = [make_entry("P02768", PROT_A, gene="ALB"), make_entry("P69999", PROT_B, gene="CAH1")]
        (tmp_path / "uniprot_snapshot_20260917.json").write_text(json.dumps(entries))
        ev = pd.DataFrame(
            {
                "peptide_row_id": [0, 1],
                "dataset_accession": ["PXDTEST", "PXDTEST"],
                "peptide_sequence": ["DAHKSEVAHR", "ARTKQTARK"],
                "protein_accession_raw": ["P02768", "P69999"],
            }
        )
        mp, proteins = mapping.map_to_uniprot(ev, tmp_path, "20260917", offline=True)
        assert mp["mapping_tier"].tolist() == ["exact_unique", "exact_unique"]
        assert len(proteins) == 2

    def test_rejects_non_canonical_table(self):
        with pytest.raises(ValueError, match="protein_accession_raw"):
            mapping.map_to_uniprot(pd.DataFrame({"a": [1]}))
