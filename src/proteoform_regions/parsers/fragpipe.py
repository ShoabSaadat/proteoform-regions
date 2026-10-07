"""FragPipe/MSFragger PSM-table adapter (Tier B).

Consumes FragPipe's combined PSM export (the IonQuant ``combined_psm.tsv``
layout: ``Peptide`` / ``Modified Peptide`` / ``Protein Start`` / ``Protein
End`` / ``Charge`` / ``Hyperscore`` / ``Protein ID`` / ...) and aggregates
PSMs to one canonical row per (stripped peptide, protein ID) - the same
entity level the mzIdentML adapter produces. Positions come from the
FragPipe protein-coordinate columns (first occurrence: minimum start, with
its end), so region inference can run directly on the evidence table.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..study import StudyConfig
from ._base import base_row, to_float


def _text(value) -> str | None:
    """Non-empty string or None (pandas gives np.nan for missing cells)."""
    return value if isinstance(value, str) and value else None


def parse_fragpipe(config: StudyConfig, paths: list[Path]) -> tuple[list[dict], int]:
    """FragPipe PSM table -> one row per (peptide, protein), PSMs aggregated.

    Tab- or comma-delimited (separator sniffed). Rows keep degenerate
    multi-protein evidence via the raw ``Protein ID`` group; decoy/contaminant
    ``REV_``/``CON_`` prefixed accessions are flagged, never dropped.
    """
    rows: list[dict] = []
    n_raw = 0
    for path in paths:
        frame = pd.read_csv(path, sep=None, engine="python")
        n_raw += len(frame)
        for (peptide, protein_id), group in frame.groupby(["Peptide", "Protein ID"], sort=False):
            if not isinstance(peptide, str) or not peptide.isalpha():
                continue
            by_start = group.sort_values(["Protein Start", "Protein End"])
            first = by_start.iloc[0]
            hyperscores = [to_float(v) for v in group.get("Hyperscore")]
            hyperscores = [v for v in hyperscores if v is not None and np.isfinite(v)]
            intensities = [to_float(v) for v in group.get("Intensity")]
            intensities = [v for v in intensities if v is not None and np.isfinite(v)]
            spectra_files = group.get("Spectrum File")
            run_count = int(spectra_files.nunique()) if spectra_files is not None else 1
            raw_accessions = str(protein_id)
            row = base_row(config, "fragpipe_psm_table", path.name)
            row.update(
                peptide_sequence=peptide,
                modified_sequence=_text(first.get("Modified Peptide")),
                protein_accession_raw=raw_accessions,
                uniprot_accession=raw_accessions.replace(",", ";").split(";")[0].strip(),
                protein_group=_text(first.get("Mapped Proteins")) or raw_accessions,
                gene_symbol=_text(first.get("Gene")) or _text(first.get("Mapped Genes")),
                protein_name_raw=_text(first.get("Protein Description")),
                start_position=to_float(first.get("Protein Start")),
                end_position=to_float(first.get("Protein End")),
                run_count=max(run_count, 1),
                psm_count=int(len(group)),
                peptide_intensity=float(np.sum(intensities)) if intensities else None,
                score_value=float(np.median(hyperscores)) if hyperscores else None,
                score_label="FragPipe/MSFragger Hyperscore (median per peptide)",
                score_source="fragpipe_psm",
                confidence_tier="Tier B",
                decoy_or_contaminant=bool(raw_accessions.upper().startswith(("REV_", "CON_"))),
            )
            rows.append(row)
    return rows, n_raw
