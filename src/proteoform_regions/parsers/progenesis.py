"""Progenesis QI peptide CSV adapter (Tier C)."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from ..study import StudyConfig
from ._base import base_row, to_float


def parse_progenesis(config: StudyConfig, paths: list[Path]) -> tuple[list[dict], int]:
    """Progenesis QI CSV export -> one row per alpha sequence (Tier C).

    Verbatim from the paper's S6 ``parse_progenesis_csv``: rows 1-2 of the
    export are abundance-block headers (header on line 3); abundance columns
    match ``^\\d{8}_``; intensities summarize as the median. NOTE: the notebook
    attached a ``retention_time_minutes`` key to each row dict, which the
    canonical column filter silently DROPPED (never landed in the frozen
    table) - the same key is attached here and dropped by the same mechanism,
    preserving published behavior rather than silently "fixing" it.
    """
    rows: list[dict] = []
    n_raw = 0
    for path in paths:
        frame = pd.read_csv(path, skiprows=2, low_memory=False)
        n_raw += len(frame)
        abundance_cols = [c for c in frame.columns if re.match(r"^\d{8}_", str(c))]
        for record in frame.to_dict(orient="records"):
            seq = record.get("Sequence")
            if not isinstance(seq, str) or not seq.isalpha():
                continue
            accession = record.get("Accession")
            row = base_row(config, "progenesis_peptide_csv", path.name)
            row.update(
                peptide_sequence=seq,
                protein_accession_raw=accession,
                uniprot_accession=accession,
                protein_name_raw=record.get("Description"),
                score_value=to_float(record.get("Score")),
                score_label="Progenesis/Mascot score",
                score_source="progenesis_csv",
                posterior_error_probability=to_float(record.get("Anova")) and None or None,
                confidence_tier="Tier C",
                retention_time_minutes=to_float(record.get("Retention time (min)")),
            )
            values = [to_float(record.get(c)) for c in abundance_cols]
            values = [v for v in values if v]
            if values:
                row["peptide_intensity"] = float(np.median(values))
            rows.append(row)
    return rows, n_raw
