"""GlyPep-Quant glycopeptide run-record adapter (Tier C; ANL-018 policy)."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from ..study import StudyConfig
from ._base import base_row, to_float, to_int

#: GlyPep-Quant per-run files interleave record blocks: a header line starting
#: with the (glycosylated-looking) sequence, optional '+' charge line, then a
#: 'Scan...' line carrying the summed intensity.
_RECORD_START = re.compile(r"^[A-Z][A-Z0-9]+\tHex")


def parse_glycopeptide(config: StudyConfig, paths: list[Path]) -> tuple[list[dict], int]:
    """pGlyco3/GlyPep-Quant per-run records -> peptidoform rows (Tier C).

    Verbatim semantics from the paper's S1 ``parse_glyco_human`` (decision D1:
    human clinical arm run files), including the ANL-018 physicochemical
    policy: pI/MW basis is the NAKED backbone; the glycan is kept in a
    ProForma-style column; NeuAc/NeuGc compositions are flagged
    ``charge_altering``; glycopeptide rows are ``tryptic_background_eligible =
    False``.
    """
    records: dict[tuple[str, str | None, str | None], dict] = {}
    total_lines = 0
    for path in paths:
        if not path.exists():
            continue
        run_id = path.stem.replace("_GlycoPeptideQuantification", "").split("_")[-1]
        with open(path, encoding="utf-8", errors="replace") as handle:
            current: dict | None = None
            for line in handle:
                total_lines += 1
                line = line.rstrip("\r\n")
                if _RECORD_START.match(line):
                    parts = line.split("\t")
                    current = {
                        "sequence": parts[0],
                        "glycan": parts[1] if len(parts) > 1 else None,
                        "sites": parts[2] if len(parts) > 2 else None,
                        "charge": None,
                        "sum_intensity": None,
                        "runs": set(),
                    }
                    continue
                if current is None:
                    continue
                if line.startswith("+") and current["charge"] is None:
                    current["charge"] = to_float(line.split("\t")[0].lstrip("+"))
                    continue
                if line.startswith("Scan"):
                    fields = [f for f in line.split("\t") if f != ""]
                    try:
                        current["sum_intensity"] = to_float(fields[-1])
                    except IndexError:
                        pass
                    key = (current["sequence"], current["glycan"], current["sites"])
                    entry = records.setdefault(
                        key, {"charges": [], "sum_intensities": [], "psm_count": 0, "runs": set()}
                    )
                    entry["psm_count"] += 1
                    entry["runs"].add(run_id)
                    if current["charge"] is not None:
                        entry["charges"].append(current["charge"])
                    if current["sum_intensity"] is not None:
                        entry["sum_intensities"].append(current["sum_intensity"])
                    current = None
    rows: list[dict] = []
    source_label = ""
    if paths:
        name0 = paths[0].name
        core = name0.rsplit("_", 2)[0]  # strip '<run>_GlycoPeptideQuantification.txt'
        source_label = f"{core}_*_{name0.rsplit('_', 1)[1]} ({len(paths)} runs)"
    for (sequence, glycan, sites), entry in records.items():
        row = base_row(config, "glypepquant_human_run_records", source_label)
        sialylation = "charge_altering" if ("NeuAc" in str(glycan) or "NeuGc" in str(glycan)) else "neutral"
        first_accession = None
        site_position = None
        if sites:
            first_site = sites.rstrip(";").split(";")[0]
            if "@" in first_site:
                first_accession, position = first_site.split("@", 1)
                site_position = to_int(position)
        row.update(
            peptide_sequence=sequence,
            protein_accession_raw=first_accession,
            uniprot_accession=first_accession,
            protein_group=sites or None,
            glycan_annotation=glycan,
            glycan_sialylation_flag=sialylation,
            entity_level="peptidoform",
            psm_count=int(entry["psm_count"]),
            run_count=len(entry["runs"]),
            charge_state=float(np.median(entry["charges"])) if entry["charges"] else None,
            peptide_intensity=float(np.median(entry["sum_intensities"])) if entry["sum_intensities"] else None,
            score_label="GlyPep-Quant record (no FDR q exported)",
            score_source="glypepquant_runs",
            start_position=site_position,
            naked_backbone_sequence=sequence,
            proforma_style_sequence=f"{sequence}[{glycan}]" if glycan else sequence,
            tryptic_background_eligible=False,
            physicochemical_sequence_basis="naked_backbone",
            confidence_tier="Tier C",
        )
        rows.append(row)
    return rows, total_lines
