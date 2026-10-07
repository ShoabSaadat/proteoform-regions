"""Stage 2: UniProt accession mapping and protein properties.

Logic extracted from the source publication's S2 (``uniprot_mapping.ipynb``):

- ``normalize_accessions``   sp|/tr| prefixes, semicolon groups, isoform
                             suffixes, CON_/REV_ drops - verbatim
- ``ensure_uniprot_entries`` chunked UniProtKB REST stream with per-chunk
                             cache + dated snapshot, retry-once, 0.6 s
                             politeness - verbatim semantics, snapshot path
                             parameterized (no cwd-relative hack)
- ``entry_properties`` / ``contaminant_flag``  protein metadata + the four
                             contaminant rules (TRYPSIN/KERATIN/BSA/IG_CHAIN)
- ``map_peptides``           exact-substring peptide->protein mapping with
                             explicit failure tiers and NO silent drops
                             (join-QC philosophy, review B9 part 1)

Contaminants are flagged, never dropped; degenerate peptides (one sequence
mapping to >1 protein) are flagged per row for downstream ambiguity rules.

**Verbatim-bug pin (divergence #32).** ``entry_properties``'s ``reviewed``
flag reproduces the publication's ``entryType.endswith("Swiss-Prot")``
check, which is always False for real entries ("UniProtKB reviewed
(Swiss-Prot)" ends with a paren). The frozen properties table therefore
carries ``reviewed=False`` on all 5,777 rows; the package preserves this
exactly for parity (pinned by test) - use ``entry_type`` directly if you
need the review status.
"""

from __future__ import annotations

import hashlib
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

from .physchem import biopython_mw, biopython_pi

__all__ = [
    "ISOFORM_SUFFIX",
    "CONTAMINANT_RULES",
    "normalize_accessions",
    "fetch_uniprot_chunk",
    "ensure_uniprot_entries",
    "load_snapshot",
    "entry_properties",
    "contaminant_flag",
    "protein_properties_table",
    "map_peptides",
    "map_to_uniprot",
]

ISOFORM_SUFFIX = re.compile(r"-\d+$")

#: (label, regex) contaminant rules searched against "accession protein_name".
CONTAMINANT_RULES = [
    ("TRYPSIN", re.compile(r"trypsin", re.I)),
    ("KERATIN", re.compile(r"keratin|krt[0-9]", re.I)),
    ("BSA", re.compile(r"^serum albumin$|bovine serum albumin", re.I)),
    ("IG_CHAIN", re.compile(r"immunoglobulin|^ig[hklm]_|^igh|^igk|^igl", re.I)),
]

_UNIPROT_STREAM_URL = "https://rest.uniprot.org/uniprotkb/stream"
_USER_AGENT = "proteoform-regions/0.1 (academic; cached batched UniProtKB REST)"


def normalize_accessions(raw_value) -> list[tuple[str, str | None]]:
    """Split a raw accession field into [(base_accession, isoform_suffix), ...].

    Handles fasta ``sp|P02768|ALBU_HUMAN`` prefixes, ``;``/``,``
    separated groups, ``-2`` isoform suffixes, and drops ``CON_``/``REV_``
    tokens. Verbatim from the publication's S2 (cell 4).
    """
    if raw_value is None or (isinstance(raw_value, float) and np.isnan(raw_value)):
        return []
    tokens: list[tuple[str, str | None]] = []
    for token in str(raw_value).replace(",", ";").split(";"):
        token = token.strip()
        if not token:
            continue
        if "|" in token:
            parts = token.split("|")
            token = parts[1] if len(parts) >= 2 else parts[-1]
        if token.upper().startswith(("CON_", "REV_")):
            continue
        isoform = None
        match = ISOFORM_SUFFIX.search(token)
        if match and token[: match.start()].upper().startswith(("P", "Q", "O", "A")):
            isoform = token[match.start() :]
            token = token[: match.start()]
        tokens.append((token, isoform))
    return tokens


def fetch_uniprot_chunk(accessions: list[str], timeout: float = 120.0, retries: int = 2) -> list[dict] | None:
    """GET one UniProtKB stream chunk (``accession:A OR accession:B ...``).

    Retry-once; a hard failure returns ``None`` (documented, never fabricated).
    """
    query = " OR ".join(f"accession:{a}" for a in accessions)
    url = _UNIPROT_STREAM_URL + "?" + urllib.parse.urlencode({"query": query, "format": "json"})
    headers = {"Accept": "application/json", "User-Agent": _USER_AGENT}
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(
                urllib.request.Request(url, headers=headers), timeout=timeout
            ) as response:
                return json.loads(response.read().decode("utf-8")).get("results", [])
        except Exception:
            if attempt == retries:
                print(f"UNIPROT CHUNK FAILED after retry (documented): {len(accessions)} accessions")
                return None
            time.sleep(2.0)
    return None


def snapshot_path(cache_dir: str | Path, snapshot_date: str) -> Path:
    return Path(cache_dir) / f"uniprot_snapshot_{snapshot_date}.json"


def load_snapshot(cache_dir: str | Path, snapshot_date: str) -> dict[str, dict]:
    """Load a dated snapshot as {primaryAccession: entry}; {} when absent."""
    path = snapshot_path(cache_dir, snapshot_date)
    if not path.exists():
        return {}
    return {e["primaryAccession"]: e for e in json.loads(path.read_text(encoding="utf-8"))}


def ensure_uniprot_entries(
    base_accessions,
    cache_dir: str | Path = "data/raw/uniprot_cache",
    snapshot_date: str = "20260917",
    chunk_size: int = 100,
    politeness_seconds: float = 0.6,
    offline: bool = False,
) -> dict[str, dict]:
    """Resolve base accessions against UniProtKB with cache + snapshot pinning.

    Resolution order per chunk: dated snapshot -> per-chunk sha1-keyed cache ->
    live REST (skipped when ``offline``). The merged snapshot is (re)written so
    later runs are network-free - the paper's reproducible-snapshot semantics.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    entries = load_snapshot(cache_dir, snapshot_date)
    ordered = sorted({a for a in base_accessions if a not in entries})
    chunks = [ordered[i : i + chunk_size] for i in range(0, len(ordered), chunk_size)]
    for chunk in chunks:
        key = hashlib.sha1(";".join(chunk).encode()).hexdigest()[:16]
        chunk_path = cache_dir / f"chunk_{key}.json"
        if chunk_path.exists():
            payload = json.loads(chunk_path.read_text(encoding="utf-8"))
        elif offline:
            continue
        else:
            payload = fetch_uniprot_chunk(chunk)
            if payload is None:
                continue
            chunk_path.write_text(json.dumps(payload), encoding="utf-8")
            time.sleep(politeness_seconds)
        for entry in payload:
            entries[entry["primaryAccession"]] = entry
    snapshot_path(cache_dir, snapshot_date).write_text(json.dumps(list(entries.values())), encoding="utf-8")
    return entries


def entry_properties(entry: dict) -> dict:
    """Canonical sequence, gene/name, entryType, organism, version - verbatim S2."""
    sequence = entry.get("sequence", {}).get("value", "")
    description = entry.get("proteinDescription", {})
    name = (description.get("recommendedName") or {}).get("fullName", {}).get("value")
    if not name:
        names = description.get("submissionNames") or []
        name = names[0].get("fullName", {}).get("value") if names else None
    genes = [g.get("geneName", {}).get("value") for g in entry.get("genes", []) if g.get("geneName")]
    return {
        "sequence": sequence,
        "sequence_length": len(sequence),
        "gene_name": "; ".join(filter(None, genes)) or None,
        "protein_name": name,
        "entry_type": entry.get("entryType"),
        "organism": (entry.get("organism") or {}).get("scientificName"),
        "reviewed": entry.get("entryType", "").endswith("Swiss-Prot"),
        "sequence_version": entry.get("entryAudit", {}).get("sequenceVersion"),
    }


def contaminant_flag(accession: str, protein_name) -> str | None:
    """First matching contaminant rule over "<accession> <protein_name>"."""
    text = f"{accession} {protein_name or ''}"
    for label, rule in CONTAMINANT_RULES:
        if rule.search(text):
            return label
    return None


def protein_properties_table(entries: dict[str, dict]) -> pd.DataFrame:
    """Per-protein properties + contaminant flag + precursor pI/MW (Biopython)."""
    rows = []
    for accession in sorted(entries):
        props = entry_properties(entries[accession])
        rows.append(
            {
                "accession": accession,
                "contaminant_flag": contaminant_flag(accession, props["protein_name"]),
                **props,
            }
        )
    proteins = pd.DataFrame(rows)
    proteins["precursor_pi"] = [biopython_pi(seq) if seq else None for seq in proteins["sequence"]]
    proteins["precursor_mw_da"] = [biopython_mw(seq) if seq else None for seq in proteins["sequence"]]
    return proteins


def map_peptides(evidence: pd.DataFrame, entries: dict[str, dict]) -> pd.DataFrame:
    """Exact-substring peptide->protein mapping with explicit tiers (no silent drops).

    Tiers: ``exact_unique`` / ``exact_multi_position`` on success; failures are
    ``no_accession_in_row`` / ``accession_unmapped_in_uniprot`` /
    ``not_found_in_group_sequences``. Also flags isoform-suffixed rows,
    contaminant hits, and degenerate peptides (sequence mapping to >1 protein).
    """
    prop_by_accession = {}
    for accession, entry in entries.items():
        props = entry_properties(entry)
        prop_by_accession[accession] = {
            "sequence": props["sequence"],
            "contaminant": contaminant_flag(accession, props["protein_name"]),
        }

    map_rows = []
    peptide_proteins: dict[str, set] = {}
    for row_index, record in enumerate(evidence.itertuples(index=False)):
        peptide = record.peptide_sequence
        tokens = normalize_accessions(record.protein_accession_raw)
        mapped_accession, positions, tier, flag = None, [], None, None
        resolved_any = False
        for token, _isoform in tokens:
            props = prop_by_accession.get(token)
            if props is None:
                continue
            resolved_any = True
            flag = flag or props["contaminant"]
            if peptide and peptide in props["sequence"]:
                found = [m.start() + 1 for m in re.finditer(f"(?={re.escape(peptide)})", props["sequence"])]
                if found:
                    mapped_accession, positions = token, found
                    tier = "exact_multi_position" if len(found) > 1 else "exact_unique"
                    break
        if tier is None:
            tier = (
                "no_accession_in_row"
                if not tokens
                else ("accession_unmapped_in_uniprot" if not resolved_any else "not_found_in_group_sequences")
            )
        if peptide and mapped_accession:
            peptide_proteins.setdefault(peptide, set()).add(mapped_accession)
        map_rows.append(
            {
                "peptide_row_id": row_index,
                "dataset_accession": record.dataset_accession,
                "peptide_sequence": peptide,
                "mapped_accession": mapped_accession,
                "positions": ";".join(str(p) for p in positions) if positions else None,
                "mapping_tier": tier,
                "row_had_isoform_suffix": any(i for _, i in tokens),
                "contaminant_flag": flag,
            }
        )
    mapping = pd.DataFrame(map_rows)
    mapping["degenerate_peptide"] = mapping["peptide_sequence"].isin(
        [p for p, accs in peptide_proteins.items() if len(accs) > 1]
    )
    return mapping


def map_to_uniprot(
    evidence: pd.DataFrame,
    cache_dir: str | Path = "data/raw/uniprot_cache",
    snapshot_date: str = "20260917",
    offline: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Full S2 stage: normalize -> resolve -> map. Returns (mapping, proteins).

    ``proteins`` carries per-protein properties (precursor pI/MW via the same
    Biopython engine as peptide features); ``mapping`` carries per-row tiers,
    contaminant and degeneracy flags.
    """
    if "protein_accession_raw" not in evidence.columns:
        raise ValueError("evidence table lacks protein_accession_raw (not a canonical S1 table?)")
    all_bases: set[str] = set()
    for raw in evidence["protein_accession_raw"]:
        all_bases.update(token for token, _ in normalize_accessions(raw))
    entries = ensure_uniprot_entries(all_bases, cache_dir, snapshot_date, offline=offline)
    mapping = map_peptides(evidence, entries)
    proteins = protein_properties_table(entries)
    return mapping, proteins
