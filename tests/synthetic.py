"""Synthetic two-protein, two-study cohort for e2e pipeline/CLI tests.

The proteins are designed so their tryptic digests yield in-range (7-60 aa)
peptides and the detected subsets stitch into multi-peptide regions. No
network: the UniProt "snapshot" is written from these sequences.
"""

from __future__ import annotations

import json
from pathlib import Path

PROTEINS = {
    "PTEST01": ("MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDK"),
    "PTEST02": ("MARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRK"),
}

#: detected peptides per study (rows in each study's peptides.txt)
STUDY_PEPTIDES = {
    "PXDTESTA": [
        ("DAHKSEVAHR", "PTEST01"),
        ("FKDLGEENFK", "PTEST01"),
        ("ALVLIAFAQYLQQCPFEDHVK", "PTEST01"),
        ("LVNEVTEFAK", "PTEST01"),
        ("STGGKAPRK", "PTEST02"),
        ("QLATKAARK", "PTEST02"),
        ("SAPATGGVK", "PTEST02"),
    ],
    "PXDTESTB": [
        ("DAHKSEVAHR", "PTEST01"),
        ("TCVADESAENCDK", "PTEST01"),
        ("STELLIRK", "PTEST02"),
        # mc1 tryptic peptide (the full 22-mer spans 4 internal cut sites = mc4,
        # beyond the enumeration bound; see digest module divergence note)
        ("KPHRYRPGTVALR", "PTEST02"),
    ],
}

_PEPTIDES_HEADER = "Sequence\tProteins\tLeading razor protein\tReverse\tPotential contaminant\tIntensity\n"


def write_cohort(root: Path) -> tuple[Path, Path, Path]:
    """Materialize raw files + snapshot + cohort YAML under ``root``.

    Returns (cohort_yaml, data_dir, uniprot_cache_dir).
    """
    data_dir = root / "raw"
    cache_dir = root / "uniprot"
    for accession, peptides in STUDY_PEPTIDES.items():
        study_dir = data_dir / accession
        study_dir.mkdir(parents=True, exist_ok=True)
        lines = [_PEPTIDES_HEADER]
        for seq, protein in peptides:
            lines.append(f"{seq}\t{protein}\t{protein}\t\t\t1000.0\n")
        (study_dir / "peptides.txt").write_text("".join(lines), encoding="utf-8")

    cache_dir.mkdir(parents=True, exist_ok=True)
    entries = []
    for i, (accession, sequence) in enumerate(sorted(PROTEINS.items()), start=1):
        entries.append(
            {
                "primaryAccession": accession,
                "sequence": {"value": sequence, "molWeight": 110.0 * len(sequence)},
                "proteinDescription": {"recommendedName": {"fullName": {"value": f"Test protein {i}"}}},
                "genes": [{"geneName": {"value": f"TST{i}"}}],
                "entryType": "UniProtKB reviewed (Swiss-Prot)",
                "organism": {"scientificName": "Homo sapiens"},
                "entryAudit": {"sequenceVersion": 1},
            }
        )
    (cache_dir / "uniprot_snapshot_20260917.json").write_text(json.dumps(entries), encoding="utf-8")

    cohort_yaml = root / "cohort.yaml"
    studies = [
        {
            "dataset_accession": accession,
            "adapter": "maxquant_peptides",
            "sample_type": "serum",
            "disease_context": "synthetic",
            "acquisition_mode": "DDA",
            "search_engine": "MaxQuant",
            "files": [{"name": "peptides.txt"}],
        }
        for accession in STUDY_PEPTIDES
    ]
    import yaml

    cohort_yaml.write_text(yaml.safe_dump({"studies": studies}, sort_keys=False), encoding="utf-8")
    return cohort_yaml, data_dir, cache_dir
