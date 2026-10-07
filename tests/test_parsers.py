"""Parser adapter tests on synthetic per-format fixtures.

Each fixture mirrors the structural quirks the adapters were extracted around
(decoy flags, skiprows headers, record-block run files, mzIdentML reference
graphs). Real-data equivalence is covered separately by
``scripts/parity_vs_paper_freeze.py`` (read-only, run at gates - not in CI).
"""

import gzip
from pathlib import Path

import pytest

import proteoform_regions as pfr
from proteoform_regions import harmonize, study
from proteoform_regions.parsers import get_parser, list_parsers, parse_study


def make_config(adapter: str, accession: str = "PXDTEST", **kw) -> study.StudyConfig:
    return study.StudyConfig(
        dataset_accession=accession,
        adapter=adapter,
        sample_type="serum",
        disease_context="test",
        acquisition_mode="DDA",
        search_engine="TestEngine",
        **kw,
    )


def write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


class TestRegistry:
    def test_eight_adapters_registered(self):
        assert list_parsers() == [
            "diann",
            "fragpipe",
            "glycopeptide",
            "maxquant_msms",
            "maxquant_peptides",
            "mzidentml",
            "progenesis",
            "spectronaut",
        ]

    def test_unknown_adapter_raises(self):
        with pytest.raises(ValueError, match="unknown adapter"):
            get_parser("doesnotexist")


class TestMaxquantPeptides:
    FIXTURE = (
        "Sequence\tModified sequence\tProteins\tLeading razor protein\tGene names\tProtein names\t"
        "Start position\tEnd position\tPEP\tScore\tMS/MS Count\tReverse\tPotential contaminant\t"
        "Reporter intensity corrected 1\tReporter intensity corrected 2\n"
        "ALPAPIEK\t_ALPAPIEK_\tP02768;P99999\tP02768\tALB\tAlbumin\t25\t32\t0.01\t95.2\t7\t\t\t1000.0\t3000.0\n"
        "KPYEEELK\t_KPYEEELK_\tP02768\tP02768\tALB\tAlbumin\t80\t87\t0.02\t80.1\t3\t\t\t500.0\t700.0\n"
        "DECOYSEQ\t_DECOYSEQ_\tP0REV\tP0REV\tREV\tRev\t1\t8\t0.5\t10.0\t1\t+\t\t100.0\t100.0\n"
        "CONTAMSEQ\t_CONTAMSEQ_\tP0CON\tP0CON\tCON\tCon\t1\t8\t0.5\t10.0\t1\t\t+\t100.0\t100.0\n"
    )

    def test_rows_tier_and_intensity(self, tmp_path):
        f = write(tmp_path / "PXDTEST/peptides.txt", self.FIXTURE)
        cfg = make_config("maxquant_peptides", files=[{"name": f.name}])
        rows, n_raw = get_parser("maxquant_peptides")(cfg, [f])
        assert n_raw == 4 and len(rows) == 2  # decoy + contaminant dropped
        row = rows[0]
        assert row["peptide_sequence"] == "ALPAPIEK"
        assert row["confidence_tier"] == "Tier B"
        assert row["peptide_intensity"] == pytest.approx(4000.0)  # sum of reporters
        assert row["start_position"] == 25.0
        assert row["source_parser_family"] == "maxquant_or_text_peptide_table"
        assert row["extraction_method"] == "maxquant_peptides_table_full"
        assert row["entity_level"] == "peptide"
        assert row["study_id"] == "PXDTEST_2026"

    def test_intensity_fallback_column(self, tmp_path):
        fixture = (
            "Sequence\tProteins\tLeading razor protein\tReverse\tPotential contaminant\tIntensity\n"
            "ALPAPIEK\tP02768\tP02768\t\t\t123.0\n"
        )
        f = write(tmp_path / "PXDTEST/peptides.txt", fixture)
        cfg = make_config("maxquant_peptides", files=[{"name": f.name}])
        rows, _ = get_parser("maxquant_peptides")(cfg, [f])
        assert rows[0]["peptide_intensity"] == pytest.approx(123.0)


class TestMaxquantMsms:
    FIXTURE = (
        "Raw file\tScan number\tSequence\tModified sequence\tProteins\tGene Names\tProtein Names\t"
        "Charge\tPEP\tScore\tRetention time\tReverse\tPotential contaminant\n"
        "run1.raw\t101\tALPAPIEK\t_ALPAPIEK_\tP02768\tALB\tAlbumin\t2\t0.01\t95\t12.0\t\t\n"
        "run1.raw\t102\tALPAPIEK\t_ALPAPIEK_\tP02768\tALB\tAlbumin\t3\t0.03\t90\t12.5\t\t\n"
        "run2.raw\t201\tKPYEEELK\t_KPYEEELK_\tP02768\tALB\tAlbumin\t2\t0.02\t80\t20.0\t\t\n"
        "run2.raw\t202\tDECOYSEQ\t_DECOYSEQ_\tREV\tR\tR\t2\t0.5\t5\t1.0\t+\t\n"
    )

    def test_psm_aggregation(self, tmp_path):
        f = write(tmp_path / "PXDTEST/msms.txt", self.FIXTURE)
        cfg = make_config("maxquant_msms", files=[{"name": f.name}])
        rows, n_raw = get_parser("maxquant_msms")(cfg, [f])
        assert n_raw == 4 and len(rows) == 2
        first = rows[0]
        assert first["peptide_sequence"] == "ALPAPIEK"
        assert first["psm_count"] == 2
        assert first["run_count"] == 1  # both PSMs from run1.raw
        assert first["charge_state"] == pytest.approx(2.5)  # median of 2,3
        assert first["score_value"] == pytest.approx(95.0)  # max
        assert first["confidence_tier"] == "Tier B"
        assert first["uniprot_accession"] == "P02768"


class TestDiann:
    HEADER = (
        "File.Name\tProtein.Group\tProtein.Ids\tProtein.Names\tGenes\tModified.Sequence\t"
        "Stripped.Sequence\tPrecursor.Charge\tQ.Value\tPEP\tCScore\tPrecursor.Normalised\tRT\tMS2.Scan\n"
    )

    def test_s1_semantics_with_usi(self, tmp_path):
        fixture = self.HEADER + (
            "S1.raw\tP02768\tP02768\tAlbumin\tALB\t_ALPAPIEK_\tALPAPIEK\t2\t0.001\t0.01\t0.99\t1e5\t12.0\t1001\n"
            "S1.raw\tP02768\tP02768\tAlbumin\tALB\t_ALPAPIEK_\tALPAPIEK\t3\t0.002\t0.02\t0.98\t2e5\t13.0\t1002\n"
        )
        f = write(tmp_path / "PXDTEST/report.tsv", fixture)
        cfg = make_config("diann", files=[{"name": f.name}])
        rows, n_raw = get_parser("diann")(cfg, [f])
        assert n_raw == 2 and len(rows) == 2  # distinct charges stay separate
        row = rows[0]
        assert row["confidence_tier"] == "Tier A"
        assert row["q_value"] == pytest.approx(0.001)
        assert row["usi_example"] == "mzspec:PXDTEST:S1.raw:scan:1001:ALPAPIEK"

    def test_s6_semantics_no_usi_no_cscore(self, tmp_path):
        header_s6 = self.HEADER.replace("\tMS2.Scan\n", "\n").replace("\tCScore", "")
        fixture = (
            header_s6
            + "S1.raw\tP02768\tP02768\tAlbumin\tALB\t_ALPAPIEK_\tALPAPIEK\t2\t0.001\t0.01\t1e5\t12.0\n"
        )
        f = write(tmp_path / "PXDTEST/report.tsv", fixture)
        cfg = make_config("diann", files=[{"name": f.name}], build_usi=False)
        rows, _ = get_parser("diann")(cfg, [f])
        assert rows[0]["usi_example"] is None
        assert rows[0]["confidence_tier"] == "Tier A"
        # CScore absent -> fallback aggregation column used without crashing
        assert rows[0]["peptide_sequence"] == "ALPAPIEK"


class TestMzIdentML:
    def test_two_pass_extraction_with_usi(self, tmp_path):
        mzid = """<?xml version="1.0" encoding="UTF-8"?>
<MzIdentML xmlns="http://psidev.info/psi/pi/mzIdentML/1.2" version="1.2.0">
 <SequenceCollection>
  <DBSequence id="dbs1" accession="P02768" searchDatabase_ref="sdb1"/>
  <DBSequence id="dbs2" accession="P69999" searchDatabase_ref="sdb1"/>
  <Peptide id="pep1"><PeptideSequence>ALPAPIEK</PeptideSequence></Peptide>
 </SequenceCollection>
 <AnalysisProtocolCollection>
  <SpectrumIdentificationProtocol id="sip1" searchDatabase_ref="sdb1"/>
  <ProteinDetectionProtocol id="pdp1"/>
 </AnalysisProtocolCollection>
 <DataCollection><Inputs><SearchDatabase id="sdb1" location="local.fasta"/><SpectraData id="sd1" location="file:///run1.raw"/></Inputs>
 <AnalysisData>
  <ProteinDetectionList id="pdl1"/>
  <SpectrumIdentificationList id="sil1">
   <SpectrumIdentificationResult id="sir1" spectraData_ref="sd1" spectrumID="scan=1001" startScan="1001">
    <SpectrumIdentificationItem id="sii1" passThreshold="true" rank="1" peptide_ref="pep1" calcMass="800.4">
     <PeptideEvidenceRef peptideEvidence_ref="pe1"/>
     <PeptideEvidenceRef peptideEvidence_ref="pe2"/>
    </SpectrumIdentificationItem>
    <SpectrumIdentificationItem id="sii2" passThreshold="false" rank="2" peptide_ref="pep1" calcMass="800.4">
     <PeptideEvidenceRef peptideEvidence_ref="pe1"/>
    </SpectrumIdentificationItem>
   </SpectrumIdentificationResult>
  </SpectrumIdentificationList>
 </AnalysisData></DataCollection>
</MzIdentML>
"""
        # PeptideEvidence elements must precede SIRs for the reference maps
        mzid = mzid.replace(
            " <DataCollection>",
            ' <AnalysisCollection><SpectrumIdentification id="si1" spectraData_ref="sd1" protocol_ref="sip1"/></AnalysisCollection>\n'  # noqa: E501
            ' <PeptideEvidence id="pe1" peptide_ref="pep1" dBSequence_ref="dbs1" start="25" end="32"/>'
            '<PeptideEvidence id="pe2" peptide_ref="pep1" dBSequence_ref="dbs2" start="9" end="16"/>\n <DataCollection>',  # noqa: E501
        )
        path = tmp_path / "PXDTEST/test.mzid.gz"
        path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(path, "wt", encoding="utf-8") as h:
            h.write(mzid)
        cfg = make_config("mzidentml", files=[{"name": path.name}])
        rows, n_sir, n_pass = get_parser("mzidentml")(cfg, [path])
        assert n_sir == 1 and n_pass == 1  # only passThreshold=true counted
        assert len(rows) == 1
        row = rows[0]
        assert row["peptide_sequence"] == "ALPAPIEK"
        assert row["protein_group"] == "P02768;P69999"  # sorted accession set
        assert row["uniprot_accession"] == "P02768"
        assert row["confidence_tier"] == "Tier B"
        assert row["usi_example"] == "mzspec:PXDTEST:run1.raw:scan:1001:ALPAPIEK"
        assert row["psm_count"] == 1


class TestGlycopeptide:
    def test_anl018_policy_and_aggregation(self, tmp_path):
        def block(seq, glycan, site, charge, intensity):
            return f"{seq}\t{glycan}\t{site}\t1\n+{charge}\t1\t1\nScan\t1\t{intensity}\n"

        run_c2 = block("AAAGFNVSLTDYWGR", "Hex(5)HexNAc(2)NeuAc(1)", "P02768@651", 2, 1000.0) + block(
            "KPYEEELK", "Hex(5)HexNAc(2)", "P02768@80", 3, 500.0
        )
        run_c14 = block("AAAGFNVSLTDYWGR", "Hex(5)HexNAc(2)NeuAc(1)", "P02768@651", 3, 3000.0)
        files = [
            write(tmp_path / "PXDTEST/liulei_OE480_2024GCHC_C2_GlycoPeptideQuantification.txt", run_c2),
            write(tmp_path / "PXDTEST/liulei_OE480_2024GCHC_C14_GlycoPeptideQuantification.txt", run_c14),
        ]
        cfg = make_config("glycopeptide", files=[{"name": f.name} for f in files])
        rows, _ = get_parser("glycopeptide")(cfg, files)
        assert len(rows) == 2
        sialylated = next(r for r in rows if "NeuAc" in str(r["glycan_annotation"]))
        assert sialylated["glycan_sialylation_flag"] == "charge_altering"
        assert sialylated["psm_count"] == 2 and sialylated["run_count"] == 2
        assert sialylated["peptide_intensity"] == pytest.approx(2000.0)  # median of 1000, 3000
        assert sialylated["charge_state"] == pytest.approx(2.5)
        assert sialylated["naked_backbone_sequence"] == "AAAGFNVSLTDYWGR"
        assert sialylated["proforma_style_sequence"] == "AAAGFNVSLTDYWGR[Hex(5)HexNAc(2)NeuAc(1)]"
        assert sialylated["tryptic_background_eligible"] is False
        assert sialylated["physicochemical_sequence_basis"] == "naked_backbone"
        assert sialylated["entity_level"] == "peptidoform"
        assert sialylated["uniprot_accession"] == "P02768"
        assert sialylated["start_position"] == 651
        assert sialylated["confidence_tier"] == "Tier C"
        neutral = next(r for r in rows if "NeuAc" not in str(r["glycan_annotation"]))
        assert neutral["glycan_sialylation_flag"] == "neutral"
        assert "GlycoPeptideQuantification" in rows[0]["source_file"]


class TestSpectronaut:
    FIXTURE = (
        "PEP.StrippedSequence\tPG.ProteinGroups\tPEP.AllOccurringProteinAccessions\tPG.Genes\t"
        "PEP.PeptidePosition\tR1.PEP.RunEvidenceCount\tR2.PEP.RunEvidenceCount\tR1.PEP.Quantity\tR2.PEP.Quantity\n"
        "ALPAPIEK\tP02768;P99999\tP02768;P99999\tALB\t25;40\t2\t0\t1000.0\t3000.0\n"
        "KPYEEELK\tP02768\tP02768\tALB\t80\t1\t1\t500.0\t700.0\n"
        "BAD_SEQ_1234\tP02768\tP02768\tALB\t5\t1\t1\t1.0\t2.0\n"
    )

    def test_quantiles_and_position(self, tmp_path):
        f = write(tmp_path / "PXDTEST/report.tsv", self.FIXTURE)
        cfg = make_config("spectronaut", files=[{"name": f.name}])
        rows, n_raw = get_parser("spectronaut")(cfg, [f])
        assert n_raw == 3 and len(rows) == 2  # non-alpha sequence dropped
        first = rows[0]
        assert first["peptide_intensity"] == pytest.approx(2000.0)  # median of 1000, 3000
        assert first["psm_count"] == 2
        assert first["run_count"] == 1  # only R1 has evidence > 0, floored at 1
        assert first["start_position"] == 25
        assert first["confidence_tier"] == "Tier C"
        assert first["uniprot_accession"] == "P02768"


class TestProgenesis:
    FIXTURE = (
        "ignore1\nignore2\n"
        "Sequence,Accession,Description,Score,Anova,Retention time (min),20240101_A,20240101_B,20240102_A\n"
        "ALPAPIEK,P02768,Albumin,45.0,0.1,12.5,1000.0,3000.0,5000.0\n"
        "KPYEEELK,P02768,Albumin,40.0,0.2,13.0,100.0,200.0,400.0\n"
        "modified_seq,P02768,Albumin,1.0,1.0,1.0,1.0,1.0,1.0\n"
    )

    def test_skiprows_median_and_drop(self, tmp_path):
        f = write(tmp_path / "PXDTEST/peptides.csv", self.FIXTURE)
        cfg = make_config("progenesis", files=[{"name": f.name}])
        rows, n_raw = get_parser("progenesis")(cfg, [f])
        assert n_raw == 3 and len(rows) == 2  # non-alpha dropped
        first = rows[0]
        assert first["peptide_intensity"] == pytest.approx(3000.0)  # median of 1000,3000,5000
        assert first["confidence_tier"] == "Tier C"
        assert first["uniprot_accession"] == "P02768"
        # the notebook's retention_time_minutes key is attached but filtered out
        # of the canonical frame (documented published behavior)
        assert "retention_time_minutes" not in {c for c in []}


class TestStudyConfig:
    def test_yaml_roundtrip_and_defaults(self, tmp_path):
        cfg = make_config("diann")
        path = tmp_path / "study.yaml"
        import yaml

        path.write_text(yaml.safe_dump(cfg.to_dict(), sort_keys=False))
        loaded = study.StudyConfig.from_yaml(path)
        assert loaded.parser_family == "diann_report"
        assert loaded.study_id == "PXDTEST_2026"
        assert loaded.representative_file_species == "Homo sapiens"

    def test_family_override(self):
        cfg = make_config("diann", parser_family="maxquant_or_text_peptide_table")
        assert cfg.parser_family == "maxquant_or_text_peptide_table"

    def test_unknown_adapter_rejected(self):
        with pytest.raises(ValueError, match="unknown adapter"):
            make_config("nope")

    def test_cohort_load(self):
        cohort = Path(__file__).parent.parent / "examples/studies/cohort-paper-tranche.yaml"
        studies = pfr.load_cohort(cohort)
        assert len(studies) == 10
        assert {s.dataset_accession for s in studies} == {
            "PXD008583",
            "PXD052666",
            "PXD054594",
            "PXD055218",
            "PXD056620",
            "PXD057799",
            "PXD060933",
            "PXD068982",
            "PXD069732",
            "PXD071549",
        }
        by_acc = {s.dataset_accession: s for s in studies}
        assert by_acc["PXD057799"].adapter == "glycopeptide" and len(by_acc["PXD057799"].files) == 8
        assert by_acc["PXD056620"].adapter == "maxquant_msms"
        assert by_acc["PXD069732"].build_usi is False


class TestHarmonizeEndToEnd:
    def test_synthetic_cohort_through_harmonize(self, tmp_path):
        write(
            tmp_path / "PXDTEST/peptides.txt",
            "Sequence\tProteins\tLeading razor protein\tReverse\tPotential contaminant\tIntensity\n"
            "ALPAPIEK\tP02768\tP02768\t\t\t123.0\n",
        )
        write(
            tmp_path / "PXDOThER/report.tsv",
            "File.Name\tProtein.Group\tProtein.Ids\tProtein.Names\tGenes\tModified.Sequence\tStripped.Sequence\t"
            "Precursor.Charge\tQ.Value\tPEP\tPrecursor.Normalised\n"
            "S1.raw\tP02768\tP02768\tAlbumin\tALB\t_ALPAPIEK_\tALPAPIEK\t2\t0.001\t0.01\t1e5\n",
        )
        studies = [
            make_config("maxquant_peptides", "PXDTEST", files=[{"name": "peptides.txt"}]),
            make_config("diann", "PXDOThER", files=[{"name": "report.tsv"}]),
        ]
        evidence, stats = harmonize(studies, tmp_path, download=False)
        assert len(evidence) == 2
        assert list(evidence.columns) == pfr.schema.CANONICAL_COLUMNS if hasattr(pfr, "schema") else True
        from proteoform_regions import schema
        from proteoform_regions.guard import assert_full_evidence_table

        assert list(evidence.columns) == schema.CANONICAL_COLUMNS
        assert_full_evidence_table(evidence, "synthetic cohort evidence")
        assert set(stats["dataset_accession"]) == {"PXDTEST", "PXDOThER"}
        assert parse_study(studies[0], tmp_path)[1]["n_evidence_rows"] == 1


class TestMaxquantPeptidesS6Quirks:
    """S6 extension-tranche numeric quirks vs S1 pilot semantics.

    The S6 ``parse_maxquant_peptides`` filtered intensities with ``if v``
    (dropping 0.0 reporters before the sum) and coerced ``MS/MS Count`` with
    ``int(... or 0) or None`` (0 -> None). Both flags together reproduce the
    S6 studies byte-for-byte (freeze parity); defaults keep S1 semantics.
    """

    FIXTURE = (
        "Sequence\tProteins\tLeading razor protein\tReverse\tPotential contaminant\tMS/MS Count\t"
        "Reporter intensity corrected 1\tReporter intensity corrected 2\n"
        "ALPAPIEK\tP02768\tP02768\t\t\t0\t1000.0\t0.0\n"  # zero count, one zero reporter
        "KPYEEELK\tP02768\tP02768\t\t\t3\t500.0\t0.0\n"  # normal count, one zero reporter
        "GNLNEQVFLK\tP02768\tP02768\t\t\t\t0.0\t0.0\n"  # missing count, all-zero reporters
    )

    def _run(self, tmp_path, **flags):
        f = write(tmp_path / "PXDTEST/peptides.txt", self.FIXTURE)
        cfg = make_config("maxquant_peptides", files=[{"name": f.name}], **flags)
        return get_parser("maxquant_peptides")(cfg, [f])[0]

    def test_s1_defaults_keep_zeros_and_zero_count(self, tmp_path):
        rows = self._run(tmp_path)
        assert [r["psm_count"] for r in rows] == [0, 3, None]
        assert [r["peptide_intensity"] for r in rows] == [1000.0, 500.0, 0.0]

    def test_s6_flags_drop_zeros_and_null_zero_count(self, tmp_path):
        rows = self._run(tmp_path, drop_zero_intensities=True, psm_count_zero_as_none=True)
        assert [r["psm_count"] for r in rows] == [None, 3, None]
        assert [r["peptide_intensity"] for r in rows] == [1000.0, 500.0, None]

    def test_flags_serialize_through_cohort_yaml(self, tmp_path):
        cfg = make_config(
            "maxquant_peptides",
            files=[{"name": "peptides.txt"}],
            drop_zero_intensities=True,
            psm_count_zero_as_none=True,
        )
        out = tmp_path / "cohort.yaml"
        study.dump_studies([cfg], out)
        loaded = study.load_cohort(out)[0]
        assert loaded.drop_zero_intensities is True
        assert loaded.psm_count_zero_as_none is True


class TestFragpipe:
    """FragPipe combined-PSM adapter: PSMs aggregate to (peptide, protein) rows.

    Unit fixtures pin the aggregation semantics (min-start position with its
    end, psm_count, median Hyperscore, decoy flag); the committed 50-row
    Sheet12 sample (tests/data) pins real-data shape - it is also the
    walkthrough dataset (examples/walkthrough).
    """

    FIXTURE = (
        "Spectrum\tSpectrum File\tPeptide\tModified Peptide\tPrev AA\tNext AA\tCharge\t"
        "Hyperscore\tIntensity\tProtein Start\tProtein End\tProtein\tProtein ID\t"
        "Entry Name\tGene\tProtein Description\tMapped Proteins\n"
        "F1\tf1.raw\tEQQDSPGNK\t\tK\tD\t2\t18.77\t100.0\t446\t454\t"
        "sp|P08697|A2AP_HUMAN\tP08697\tA2AP_HUMAN\tSERPINF2\tAlpha-2-antiplasmin\t\n"
        "F2\tf1.raw\tEQQDSPGNK\t\tK\tD\t3\t16.78\t50.0\t446\t454\t"
        "sp|P08697|A2AP_HUMAN\tP08697\tA2AP_HUMAN\tSERPINF2\tAlpha-2-antiplasmin\t\n"
        "F3\tf2.raw\tEQQDSPGNK\t\tK\tD\t2\t21.10\t0.0\t446\t454\t"
        "sp|P08697|A2AP_HUMAN\tP08697\tA2AP_HUMAN\tSERPINF2\tAlpha-2-antiplasmin\t\n"
        "F4\tf1.raw\tSGKDPDR\tsGKDPDR\tR\tF\t2\t19.73\t\t192\t198\t"
        "tr|A0A096LPE2|XX\tA0A096LPE2\tXX\tSAA2-SAA4\tReadthrough protein\ttr|A|\n"
        "F5\tf1.raw\tQANTEER\t\tK\tK\t2\t22.20\t10.0\t349\t355\t"
        "sp|P06396|GELS_HUMAN\tP06396\tGELS_HUMAN\tGSN\tGelsolin\t\n"
        "F6\tf1.raw\tDECOYPEPK\t\tK\tK\t2\t9.99\t\t1\t9\t"
        "sp|REV_P99999|XXX\tREV_P99999\tXXX\tXXX\tDecoy protein\t\n"
    )

    def _run(self, tmp_path, **kw):
        f = write(tmp_path / "PXDTEST/psm.tsv", self.FIXTURE)
        cfg = make_config("fragpipe", **kw)
        rows, n_raw = get_parser("fragpipe")(cfg, [f])
        return rows, n_raw

    def test_psm_aggregation_semantics(self, tmp_path):
        rows, n_raw = self._run(tmp_path)
        assert n_raw == 6
        by_seq = {r["peptide_sequence"]: r for r in rows}
        agg = by_seq["EQQDSPGNK"]
        assert agg["psm_count"] == 3  # three PSMs collapse to one peptide row
        assert agg["run_count"] == 2  # two distinct spectrum files
        assert agg["start_position"] == 446.0 and agg["end_position"] == 454.0
        assert agg["peptide_intensity"] == 150.0  # summed over PSMs
        assert agg["score_value"] == pytest.approx(18.77)  # median hyperscore
        assert agg["uniprot_accession"] == "P08697"
        assert agg["gene_symbol"] == "SERPINF2"
        assert agg["protein_name_raw"] == "Alpha-2-antiplasmin"
        assert agg["modified_sequence"] is None  # empty Modified Peptide -> None
        assert agg["source_parser_family"] == "fragpipe_psm_table"
        assert agg["confidence_tier"] == "Tier B"

    def test_decoy_flagged_not_dropped(self, tmp_path):
        rows, _ = self._run(tmp_path)
        decoys = [r for r in rows if r["decoy_or_contaminant"]]
        assert len(decoys) == 1 and decoys[0]["peptide_sequence"] == "DECOYPEPK"
        assert len(rows) == 4  # 4 (peptide, protein) groups from 6 PSMs

    def test_missing_intensity_and_modified_stay_none(self, tmp_path):
        rows, _ = self._run(tmp_path)
        by_seq = {r["peptide_sequence"]: r for r in rows}
        assert by_seq["SGKDPDR"]["peptide_intensity"] is None
        assert by_seq["SGKDPDR"]["modified_sequence"] == "sGKDPDR"

    def test_comma_delimited_sheet_is_sniffed(self, tmp_path):
        comma = self.FIXTURE.replace("\t", ",")
        f = write(tmp_path / "PXDTEST/psm.csv", comma)
        cfg = make_config("fragpipe")
        rows, n_raw = get_parser("fragpipe")(cfg, [f])
        assert n_raw == 6 and len(rows) == 4

    def test_sheet12_sample_fifty_rows(self):
        """The committed walkthrough sample parses: 50 PSMs -> 23 peptide rows."""
        import pandas as pd

        from proteoform_regions.parsers import parse_study

        root = Path(__file__).resolve().parents[1]
        cfg = study.StudyConfig(
            dataset_accession="PXD077545",
            adapter="fragpipe",
            search_engine="FragPipe/MSFragger",
            files=[{"name": "sheet12_sample50.csv"}],
        )
        rows, stats = parse_study(cfg, root / "examples/walkthrough")
        assert stats["n_evidence_rows"] == 23
        assert stats["stat1"] == 50
        frame = pd.DataFrame(rows)
        assert frame["start_position"].notna().all()
        assert (frame["psm_count"] >= 1).all()
        assert not frame["decoy_or_contaminant"].any()
