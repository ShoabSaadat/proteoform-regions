"""DIA-NN report.tsv adapter: precursor aggregates (Tier A)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from ..study import StudyConfig
from ._base import base_row, to_float


def parse_diann(config: StudyConfig, paths: list[Path]) -> tuple[list[dict], int]:
    """DIA-NN report -> one row per (seq, modseq, charge, protein group) (Tier A).

    Union of the paper's S1 ``parse_diann`` (USI construction from MS2.Scan)
    and S6 ``parse_diann_report`` (usecols filtered to header; PEP fallbacks
    when CScore/Precursor.Normalised are absent). USIs are built only when
    ``config.build_usi`` is not False AND MS2.Scan is available (S6's
    PXD069732 extraction never requested MS2.Scan -> USI None).
    """
    usecols_wanted = [
        "File.Name", "Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "Modified.Sequence",
        "Stripped.Sequence", "Precursor.Charge", "Q.Value", "PEP", "CScore", "Precursor.Normalised",
        "RT", "MS2.Scan",
    ]
    rows: list[dict] = []
    n_raw = 0
    for path in paths:
        header = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
        usecols = [c for c in usecols_wanted if c in header]
        frame = pd.read_csv(path, sep="\t", usecols=usecols, low_memory=False)
        n_raw += len(frame)
        agg = (
            frame.groupby(
                ["Stripped.Sequence", "Modified.Sequence", "Precursor.Charge", "Protein.Group", "Protein.Ids", "Genes"],
                dropna=False,
            )
            .agg(
                psm_count=("Stripped.Sequence", "size"),
                n_runs=("File.Name", "nunique"),
                q_min=("Q.Value", "min"),
                pep_median=("PEP", "median"),
                cscore_median=("CScore", "median") if "CScore" in usecols else ("PEP", "median"),
                intensity_median=("Precursor.Normalised", "median")
                if "Precursor.Normalised" in usecols
                else ("PEP", "median"),
                protein_name=("Protein.Names", "first") if "Protein.Names" in usecols else ("Protein.Group", "first"),
                example_file=("File.Name", "first"),
                example_scan=("MS2.Scan", "first") if "MS2.Scan" in usecols else ("PEP", "median"),
            )
            .reset_index()
        )
        build_usi = config.build_usi is not False and "MS2.Scan" in usecols
        for record in agg.to_dict(orient="records"):
            row = base_row(config, "diann_report_precursor_aggregate", path.name)
            usi = None
            if build_usi and record.get("example_file") and pd.notna(record.get("example_scan")):
                run_name = str(record["example_file"]).replace("\\", "/").split("/")[-1]
                usi = (
                    f"mzspec:{config.dataset_accession}:{run_name}:"
                    f"scan:{int(record['example_scan'])}:{record['Stripped.Sequence']}"
                )
            row.update(
                peptide_sequence=record["Stripped.Sequence"],
                modified_sequence=record["Modified.Sequence"],
                protein_group=record["Protein.Group"],
                protein_accession_raw=record["Protein.Ids"],
                uniprot_accession=str(record["Protein.Ids"]).split(";")[0] if record["Protein.Ids"] else None,
                gene_symbol=record["Genes"],
                protein_name_raw=record.get("protein_name"),
                charge_state=to_float(record["Precursor.Charge"]),
                psm_count=int(record["psm_count"]),
                run_count=int(record["n_runs"]),
                q_value=to_float(record["q_min"]),
                posterior_error_probability=to_float(record["pep_median"]),
                score_value=to_float(record["cscore_median"]),
                score_label="DIA-NN CScore (median per precursor)",
                score_source="diann_report",
                peptide_intensity=to_float(record["intensity_median"]),
                usi_example=usi,
                confidence_tier="Tier A",
            )
            rows.append(row)
    return rows, n_raw
