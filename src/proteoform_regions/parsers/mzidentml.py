"""mzIdentML adapter: passThreshold PSM aggregates via a fast two-pass lxml scan (Tier B)."""

from __future__ import annotations

import gzip
import re
from pathlib import Path

from lxml import etree

from ..study import StudyConfig
from ._base import base_row


def parse_mzidentml(config: StudyConfig, paths: list[Path]) -> tuple[list[dict], int, int]:
    """mzIdentML (.mzid / .mzid.gz) -> one row per (sequence, accession set).

    Verbatim strategy from the paper's S1 ``parse_mzid``: two lxml passes
    (pass 1 builds DBSequence/Peptide/PeptideEvidence/SpectraData reference
    maps; pass 2 aggregates SpectrumIdentificationItems at passThreshold=true).
    pyteomics ``retrieve_refs=True`` is orders of magnitude slower on these
    files, hence the hand-rolled scan. Returns (rows, n_spectrum_results,
    n_pass_threshold_items).
    """
    scan_re = re.compile(r"scan=(\d+)", re.IGNORECASE)
    rows: list[dict] = []
    n_sir_total = 0
    n_pass_total = 0
    for path in paths:
        opener = gzip.open if path.suffix == ".gz" else open
        db_map: dict[str, str | None] = {}
        peptide_map: dict[str, str | None] = {}
        evidence_map: dict[str, str | None] = {}
        spectra_map: dict[str, str] = {}
        with opener(path, "rb") as handle:  # type: ignore[operator]
            for _, elem in etree.iterparse(handle, events=("end",)):
                tag = etree.QName(elem).localname
                if tag == "DBSequence":
                    db_map[elem.get("id")] = elem.get("accession")
                    elem.clear()
                elif tag == "Peptide":
                    sequence_elem = elem.find("{*}PeptideSequence")
                    peptide_map[elem.get("id")] = sequence_elem.text if sequence_elem is not None else None
                    elem.clear()
                elif tag == "PeptideEvidence":
                    evidence_map[elem.get("id")] = elem.get("dBSequence_ref")
                    elem.clear()
                elif tag == "SpectraData":
                    location = elem.get("location") or ""
                    spectra_map[elem.get("id")] = location.replace("\\", "/").split("/")[-1]
                    elem.clear()
        counter: dict[tuple[str, str], dict] = {}
        n_sir = 0
        n_pass = 0
        with opener(path, "rb") as handle:  # type: ignore[operator]
            for _, elem in etree.iterparse(handle, events=("end",)):
                if etree.QName(elem).localname != "SpectrumIdentificationResult":
                    continue
                n_sir += 1
                spectrum_id = elem.get("spectrumID") or ""
                spectra_ref = elem.get("spectraData_ref")
                scan_number = scan_re.search(spectrum_id)
                for item in elem.findall("{*}SpectrumIdentificationItem"):
                    if str(item.get("passThreshold", "")).lower() != "true":
                        continue
                    n_pass += 1
                    sequence = peptide_map.get(item.get("peptide_ref"))
                    if not sequence:
                        continue
                    accessions = sorted(
                        {
                            db_map.get(evidence_map.get(ref.get("peptideEvidence_ref") or ref.get("ref")))
                            for ref in item.findall("{*}PeptideEvidenceRef")
                            if evidence_map.get(ref.get("peptideEvidence_ref") or ref.get("ref")) in db_map
                        }
                    )
                    key = (sequence, ";".join(accessions) if accessions else "")
                    entry = counter.setdefault(key, {"psm_count": 0, "example_usi": None})
                    entry["psm_count"] += 1
                    if entry["example_usi"] is None and scan_number and spectra_ref in spectra_map:
                        entry["example_usi"] = (
                            f"mzspec:{config.dataset_accession}:{spectra_map[spectra_ref]}:"
                            f"scan:{scan_number.group(1)}:{sequence}"
                        )
                elem.clear()
        n_sir_total += n_sir
        n_pass_total += n_pass
        for (sequence, accession_set), entry in counter.items():
            row = base_row(config, "mzidentml_psm_aggregate_lxml", path.name)
            first_accession = accession_set.split(";")[0] if accession_set else None
            row.update(
                peptide_sequence=sequence,
                protein_group=accession_set or None,
                protein_accession_raw=first_accession,
                uniprot_accession=first_accession,
                psm_count=int(entry["psm_count"]),
                usi_example=entry["example_usi"],
                score_value=float(entry["psm_count"]),
                score_label="passThreshold PSM count",
                score_source="mzidentml",
                confidence_tier="Tier B",
            )
            rows.append(row)
    return rows, n_sir_total, n_pass_total
