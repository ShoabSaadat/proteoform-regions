"""Search-engine peptide-table parser registry.

Each adapter is a pure function ``(config, paths) -> (rows, *stats)`` emitting
canonical-schema row dicts. Adding support for a new export format means
adding one module and one registry entry - no study metadata inside the
adapter (that lives in :class:`~proteoform_regions.study.StudyConfig`).

Adapters shipped (the source publication's nine parser functions collapse to
seven format adapters; ``maxquant_peptides`` and ``diann`` each cover two
notebook variants):

======================  =====================================================
adapter key             source format
======================  =====================================================
``maxquant_peptides``   MaxQuant peptides.txt / peptides.csv tables
``maxquant_msms``       MaxQuant msms.txt PSM tables
``diann``               DIA-NN report.tsv (precursor level)
``mzidentml``           mzIdentML (.mzid / .mzid.gz)
``glycopeptide``        pGlyco3/GlyPep-Quant per-run glycopeptide records
``spectronaut``         Spectronaut facility peptide quant TSV
``progenesis``          Progenesis QI peptide CSV
``fragpipe``            FragPipe/MSFragger combined PSM table
======================  =====================================================
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from ..study import StudyConfig
from .diann import parse_diann
from .fragpipe import parse_fragpipe
from .glycopeptide import parse_glycopeptide
from .maxquant import parse_maxquant_msms, parse_maxquant_peptides
from .mzidentml import parse_mzidentml
from .progenesis import parse_progenesis
from .spectronaut import parse_spectronaut

PARSERS: dict[str, Callable] = {
    "maxquant_peptides": parse_maxquant_peptides,
    "maxquant_msms": parse_maxquant_msms,
    "diann": parse_diann,
    "mzidentml": parse_mzidentml,
    "glycopeptide": parse_glycopeptide,
    "spectronaut": parse_spectronaut,
    "progenesis": parse_progenesis,
    "fragpipe": parse_fragpipe,
}


def get_parser(adapter: str) -> Callable:
    try:
        return PARSERS[adapter]
    except KeyError:
        raise ValueError(f"unknown adapter {adapter!r}; registered: {sorted(PARSERS)}") from None


def list_parsers() -> list[str]:
    return sorted(PARSERS)


def parse_study(config: StudyConfig, data_dir: str | Path) -> tuple[list[dict], dict]:
    """Run a study's adapter over its files; returns (canonical rows, stats)."""
    parser = get_parser(config.adapter)
    paths = config.file_paths(data_dir)
    result = parser(config, paths)
    rows, stats = result[0], result[1:]
    stats_out = {
        "dataset_accession": config.dataset_accession,
        "adapter": config.adapter,
        "n_files": len(paths),
        "n_evidence_rows": len(rows),
    }
    # adapters report heterogeneous leading stats (raw rows / SIRs / lines...)
    for i, value in enumerate(stats):
        stats_out[f"stat{i + 1}"] = value
    return rows, stats_out
