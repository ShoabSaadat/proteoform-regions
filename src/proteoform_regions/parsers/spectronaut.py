"""Spectronaut facility peptide-quant adapter (Tier C)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..study import StudyConfig
from ._base import base_row, to_float


def parse_spectronaut(config: StudyConfig, paths: list[Path]) -> tuple[list[dict], int]:
    """Spectronaut facility peptide report -> one row per alpha sequence.

    Verbatim from the paper's S6 ``parse_spectronaut``: run evidence counts
    sum into psm_count; median over peptide quantities; first parseable
    PeptidePosition token becomes start_position.
    """
    rows: list[dict] = []
    n_raw = 0
    for path in paths:
        frame = pd.read_csv(path, sep="\t", low_memory=False)
        n_raw += len(frame)
        evidence_count_cols = [c for c in frame.columns if c.endswith(".PEP.RunEvidenceCount")]
        quantity_cols = [c for c in frame.columns if ".PEP.Quantity" in c or ("[]" not in c and c.endswith(".AvgQuantity"))]
        if not quantity_cols:
            quantity_cols = [c for c in frame.columns if "Quantity" in c and "RunEvidence" not in c]
        for record in frame.to_dict(orient="records"):
            seq = record.get("PEP.StrippedSequence")
            if not isinstance(seq, str) or not seq.isalpha():
                continue
            ev_counts = [to_float(record.get(c)) for c in evidence_count_cols]
            run_count = sum(1 for v in ev_counts if v is not None and v > 0)
            quants = [to_float(record.get(c)) for c in quantity_cols]
            quants = [v for v in quants if v is not None and np.isfinite(v)]
            positions = str(record.get("PEP.PeptidePosition", "") or "")
            first_pos = None
            for token in positions.split(";"):
                try:
                    first_pos = int(float(token.strip()))
                    break
                except ValueError:
                    continue
            row = base_row(config, "spectronaut_facility_peptide_quant", path.name)
            row.update(
                peptide_sequence=seq,
                protein_accession_raw=record.get("PG.ProteinGroups"),
                uniprot_accession=str(record.get("PG.ProteinGroups")).split(";")[0]
                if record.get("PG.ProteinGroups")
                else None,
                protein_group=record.get("PEP.AllOccurringProteinAccessions") or record.get("PG.ProteinGroups"),
                gene_symbol=record.get("PG.Genes"),
                start_position=first_pos,
                run_count=max(run_count, 1),
                peptide_intensity=float(np.median(quants)) if quants else None,
                score_label="Spectronaut facility peptide quant (no per-peptide q exported)",
                score_source="spectronaut_facility",
                confidence_tier="Tier C",
                psm_count=int(sum(v for v in ev_counts if v is not None)),
            )
            rows.append(row)
    return rows, n_raw
