"""MaxQuant adapters: peptides.txt tables and msms.txt PSM aggregates."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..study import StudyConfig
from ._base import base_row, to_float, to_int


def parse_maxquant_peptides(config: StudyConfig, paths: list[Path]) -> tuple[list[dict], int]:
    """MaxQuant peptides.txt -> one evidence row per peptide (Tier B).

    Union of the paper's S1 ``parse_peptides_tissue`` and S6
    ``parse_maxquant_peptides``: reporter-corrected intensities (fallback:
    ``Intensity``) are summed over finite non-missing values (zeros KEPT - S1
    semantics; the frozen table keeps 0.0 sums); Reverse/contaminant rows
    dropped; ``psm_count`` 0 is preserved (S1 ``_to_int`` semantics).

    ``config.extract_protein_names=False`` replicates the S6 extension-tranche
    parsers, which did not map the ``Protein names`` column (the frozen table
    carries empty ``protein_name_raw`` for those studies).
    """
    rows: list[dict] = []
    n_raw = 0
    for path in paths:
        frame = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
        n_raw += len(frame)
        intensity_cols = [c for c in frame.columns if c.startswith("Reporter intensity corrected")] or [
            c for c in frame.columns if c == "Intensity"
        ]
        for record in frame.to_dict(orient="records"):
            if str(record.get("Reverse", "")) == "+" or str(record.get("Potential contaminant", "")) == "+":
                continue
            row = base_row(config, "maxquant_peptides_table_full", path.name)
            row.update(
                peptide_sequence=record.get("Sequence"),
                modified_sequence=record.get("Modified sequence"),
                protein_group=record.get("Proteins"),
                uniprot_accession=record.get("Leading razor protein"),
                protein_accession_raw=record.get("Leading razor protein"),
                gene_symbol=record.get("Gene names"),
                start_position=to_float(record.get("Start position")),
                end_position=to_float(record.get("End position")),
                posterior_error_probability=to_float(record.get("PEP")),
                score_value=to_float(record.get("Score")),
                score_label="MaxQuant Score",
                score_source="maxquant_peptides",
                psm_count=to_int(record.get("MS/MS Count")),
                confidence_tier="Tier B",
            )
            if getattr(config, "extract_protein_names", True):
                row["protein_name_raw"] = record.get("Protein names")
            values = [to_float(record.get(c)) for c in intensity_cols]
            values = [v for v in values if v is not None and np.isfinite(v)]
            if values:
                row["peptide_intensity"] = float(np.sum(values))
            rows.append(row)
    return rows, n_raw


def parse_maxquant_msms(config: StudyConfig, paths: list[Path]) -> tuple[list[dict], int]:
    """MaxQuant msms.txt -> PSM aggregates per (sequence, modseq, proteins) (Tier B).

    Verbatim from the paper's S1 ``parse_msms``: missing wanted columns are
    skipped with a note; Reverse/contaminant PSMs dropped before aggregation;
    medians/max over PSMs per peptide group.
    """
    wanted = [
        "Raw file", "Scan number", "Sequence", "Modified sequence", "Proteins", "Gene Names",
        "Protein Names", "Charge", "PEP", "Score", "Retention time", "Reverse",
        "Potential contaminant",
    ]
    rows: list[dict] = []
    n_raw = 0
    for path in paths:
        header = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
        usecols = [c for c in wanted if c in header]
        frame = pd.read_csv(path, sep="\t", usecols=usecols, low_memory=False)
        raw_len = len(frame)
        if "Reverse" in frame.columns:
            frame = frame[~frame["Reverse"].astype(str).str.strip().isin(["+"])]
        if "Potential contaminant" in frame.columns:
            frame = frame[~frame["Potential contaminant"].astype(str).str.strip().isin(["+"])]
        n_raw += raw_len
        agg = (
            frame.groupby(["Sequence", "Modified sequence", "Proteins"], dropna=False)
            .agg(
                psm_count=("Sequence", "size"),
                n_runs=("Raw file", "nunique"),
                pep_median=("PEP", "median"),
                score_max=("Score", "max"),
                charge_median=("Charge", "median"),
                gene=("Gene Names", "first"),
                protein_name=("Protein Names", "first"),
            )
            .reset_index()
        )
        for record in agg.to_dict(orient="records"):
            row = base_row(config, "maxquant_msms_psm_aggregate", path.name)
            row.update(
                peptide_sequence=record["Sequence"],
                modified_sequence=record["Modified sequence"],
                protein_group=record["Proteins"],
                protein_accession_raw=record["Proteins"].split(";")[0] if record["Proteins"] else None,
                uniprot_accession=record["Proteins"].split(";")[0] if record["Proteins"] else None,
                gene_symbol=record["gene"],
                protein_name_raw=record["protein_name"],
                psm_count=int(record["psm_count"]),
                run_count=int(record["n_runs"]),
                posterior_error_probability=float(record["pep_median"]) if pd.notna(record["pep_median"]) else None,
                score_value=float(record["score_max"]) if pd.notna(record["score_max"]) else None,
                score_label="MaxQuant Andromeda score (max per peptide)",
                score_source="maxquant_msms",
                charge_state=float(record["charge_median"]) if pd.notna(record["charge_median"]) else None,
                confidence_tier="Tier B",
            )
            rows.append(row)
    return rows, n_raw
