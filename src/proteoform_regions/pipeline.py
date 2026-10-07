"""End-to-end orchestration: harmonize -> map -> compute -> stitch -> bias.

Each stage function is a thin, faithful composition of the stage modules
(the publication's S2-S4, S8 notebook flows), so ``run()`` reproduces the
paper's pipeline as one call:

- ``map_to_uniprot``     (stage 2)  accession mapping + protein properties
- ``compute_features``   (stage 3a) peptide physicochemical feature table
- ``build_background``   (stage 3b) digestion-matched tryptic background
- ``stitch_regions``     (stage 4)  single-linkage detected-region inference
- ``quantify_bias``      (stage 8)  RQ1 permutation + RQ2 OR meta-analysis
- ``run``                          stage-selected execution + artifact
                                   writing + sha256 run manifest

Every artifact filename matches the paper's ``data/processed/`` vocabulary;
every stochastic knob (seed, B_perm, B_boot) is a parameter with the paper's
defaults.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from . import stats as st
from .digest import background_peptides
from .guard import write_manifest
from .harmonize import harmonize
from .mapping import map_to_uniprot as _map_stage
from .physchem import biopython_mw, biopython_pi, charge_at_ph, gravy, hh_pi_batch, peptide_features
from .regions import MIN_REGION_PEPTIDES, REGION_GAP_AA, TERMINAL_FRACTION, region_record, stitch_spans
from .schema import classify_pi
from .study import StudyConfig, load_cohort

__all__ = [
    "compute_features",
    "build_background",
    "stitch_regions",
    "quantify_bias",
    "run",
    "RunReport",
    "STAGES",
]

#: Executable pipeline stages, in canonical order (S1 harmonize is run() only).
STAGES = ("harmonize", "map", "compute", "background", "stitch", "bias")

_EXACT_TIERS = ("exact_unique", "exact_multi_position")


# ---------------------------------------------------------------- stage 3a


def compute_features(
    evidence: pd.DataFrame,
    mapping: pd.DataFrame,
    proteins: pd.DataFrame,
    calculator: str = "bjellqvist",
) -> pd.DataFrame:
    """Peptide feature table over exactly-mapped rows (S3 feature engine).

    ``calculator`` selects the peptide pI engine: ``bjellqvist`` (Biopython,
    the paper's primary) or an HH table name (``hh_rodwell_emboss_lineage`` /
    ``hh_lehninger``) for the named sensitivity analysis. Net charge at pH 7.4
    always uses the Rodwell/EMBOSS HH table, as in the source.
    """
    mapped = mapping[mapping["mapping_tier"].isin(_EXACT_TIERS)].copy()
    protein_lookup = proteins.set_index("accession")

    def _pi_map(seqs: list[str]) -> dict[str, float]:
        if calculator == "bjellqvist":
            return {s: biopython_pi(s) for s in seqs}
        if calculator in ("hh_rodwell_emboss_lineage", "hh_lehninger"):
            return dict(zip(seqs, hh_pi_batch(seqs, calculator), strict=True))
        raise ValueError(f"unknown calculator {calculator!r}")

    seqs_needed = pd.unique(mapped["peptide_sequence"]).tolist()
    pi_map = _pi_map([str(s) for s in seqs_needed])
    pi_cache: dict[str, tuple] = {}
    for seq in seqs_needed:
        pi_cache[seq] = (
            pi_map[seq],
            biopython_mw(seq),
            gravy(seq),
            peptide_features(seq),
        )

    feature_rows = []
    for record in mapped.itertuples(index=False):
        seq = record.peptide_sequence
        p_i, mw, gr, feats = pi_cache[seq]
        precursor = (
            protein_lookup.loc[record.mapped_accession]
            if record.mapped_accession in protein_lookup.index
            else None
        )
        row = {
            "peptide_row_id": record.peptide_row_id,
            "dataset_accession": record.dataset_accession,
            "peptide_sequence": seq,
            "mapped_accession": record.mapped_accession,
            "mapping_tier": record.mapping_tier,
            "peptide_pi": p_i,
            "peptide_mw_da": mw,
            "peptide_gravy": gr,
            "contaminant_flag": record.contaminant_flag,
            **feats,
            "peptide_class_anl003": classify_pi(p_i),
        }
        if precursor is not None:
            row.update(
                {
                    "precursor_pi": precursor["precursor_pi"],
                    "precursor_mw_da": precursor["precursor_mw_da"],
                    "precursor_length": precursor["sequence_length"],
                    "gene_name": precursor["gene_name"],
                    "protein_name": precursor["protein_name"],
                    "pi_divergence_peptide_minus_precursor": p_i - precursor["precursor_pi"],
                    "mw_ratio_peptide_over_precursor": (
                        mw / precursor["precursor_mw_da"] if precursor["precursor_mw_da"] else None
                    ),
                }
            )
        feature_rows.append(row)
    features = pd.DataFrame(feature_rows)

    # net charge at pH 7.4 (HH Rodwell lineage; batched as in the source)
    unique_seqs = [str(s) for s in features["peptide_sequence"].unique()]
    features["net_charge_ph74"] = features["peptide_sequence"].map(charge_at_ph(unique_seqs, 7.4))
    features["charge_density"] = features["net_charge_ph74"] / features["peptide_length"]
    return features


# ---------------------------------------------------------------- stage 3b


def build_background(
    protein_sequences: dict[str, str],
    max_missed: int = 2,
    min_len: int = 7,
    max_len: int = 60,
) -> pd.DataFrame:
    """Digestion-matched tryptic background over {accession: sequence} (S3).

    Same digest protocol and same physicochemical columns as the detected
    peptides, so the RQ2 comparison is digestion-matched by construction.
    """
    digested = background_peptides(protein_sequences, max_missed, min_len, max_len)
    pi_cache: dict[str, tuple] = {}

    def props(seq: str) -> tuple:
        if seq not in pi_cache:
            pi_cache[seq] = (
                biopython_pi(seq),
                biopython_mw(seq),
                gravy(seq),
                peptide_features(seq),
            )
        return pi_cache[seq]

    rows = []
    for rec in digested:
        p_i, mw, gr, feats = props(rec["peptide_sequence"])
        rows.append(
            {
                "accession": rec["accession"],
                "peptide_sequence": rec["peptide_sequence"],
                "missed_cleavages": rec["missed_cleavages"],
                "peptide_pi": p_i,
                "peptide_mw_da": mw,
                "peptide_gravy": gr,
                "peptide_class_anl003": classify_pi(p_i),
                **feats,
            }
        )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- stage 4


def stitch_regions(
    features: pd.DataFrame,
    mapping: pd.DataFrame,
    precursor_sequences: dict[str, str] | None = None,
    gap: int = REGION_GAP_AA,
    terminal_fraction: float = TERMINAL_FRACTION,
    min_peptides: int = MIN_REGION_PEPTIDES,
) -> pd.DataFrame:
    """Detected-region inference (S4 walk): position spans -> stitched regions.

    Start positions come from the mapping table (first exact-substring hit);
    rows without exact mapping never reach this stage. ``precursor_sequences``
    (e.g. from the UniProt snapshot) enables region sequence / region pI/MW
    columns; without it they stay null (positions and tiers are unaffected).
    """
    precursor_sequences = precursor_sequences or {}
    mapped = mapping[mapping["mapping_tier"].isin(_EXACT_TIERS)].copy()
    if mapped.empty:
        return pd.DataFrame()
    mapped["start_position"] = mapped["positions"].astype(str).str.split(";").str[0].astype(float)
    mapped = mapped[mapped["start_position"].notna()]
    merged = mapped.merge(
        features[
            [
                "peptide_row_id",
                "peptide_pi",
                "peptide_length",
                "precursor_pi",
                "precursor_mw_da",
                "precursor_length",
                "gene_name",
            ]
        ],
        on="peptide_row_id",
        how="left",
        suffixes=("", "_feat"),
    )
    merged["end_position"] = merged["start_position"] + merged["peptide_length"] - 1

    region_rows = []
    for (study, accession), group in merged.groupby(["dataset_accession", "mapped_accession"]):
        precursor_length = group["precursor_length"].iloc[0]
        if not precursor_length or precursor_length != precursor_length:
            continue
        spans = list(
            zip(
                group["start_position"],
                group["end_position"],
                group["peptide_sequence"],
                group["peptide_pi"],
                group["peptide_row_id"],
                strict=True,
            )
        )
        for cluster in stitch_spans(spans, gap):
            region_rows.append(
                region_record(
                    study,
                    accession,
                    cluster,
                    precursor_length=int(precursor_length),
                    precursor_pi=group["precursor_pi"].iloc[0],
                    precursor_mw=group["precursor_mw_da"].iloc[0],
                    gene_name=group["gene_name"].iloc[0],
                    precursor_sequence=precursor_sequences.get(accession, ""),
                    gap=gap,
                    terminal_fraction=terminal_fraction,
                    min_peptides=min_peptides,
                )
            )
    return pd.DataFrame(region_rows)


# ---------------------------------------------------------------- stage 8


def quantify_bias(
    features: pd.DataFrame,
    background: pd.DataFrame,
    seed: int = st.SEED,
    b_perm: int = st.B_PERM,
    b_boot: int = st.B_BOOT,
    run_rq1: bool = True,
) -> dict:
    """Detection-bias quantification on feature + background tables (S8).

    RQ2: per-study basic-vs-rest OR (Woolf, unconditional 0.5) against each
    study's digestion-matched parent-protein background; inverse-variance
    fixed effect; DerSimonian-Laird random effects; leave-one-study-out
    envelope. RQ1: per-protein divergence permutation null (seed-aware).

    Set ``run_rq1=False`` to skip the permutation null (it dominates runtime
    on paper-scale tables: ~5,400 proteins x B=2,000 draws).
    """
    det = features.drop_duplicates(["dataset_accession", "peptide_sequence"])
    study_parents = (
        features.groupby("dataset_accession")["mapped_accession"].agg(lambda s: set(s.dropna())).to_dict()
    )

    study_rows = []
    for study in sorted(det["dataset_accession"].unique()):
        d_cls = det[det["dataset_accession"] == study]["peptide_class_anl003"].dropna()
        bg = background[background["accession"].isin(study_parents[study])].drop_duplicates(
            "peptide_sequence"
        )
        b_cls = bg["peptide_class_anl003"].dropna()
        r = st.study_or(d_cls.to_numpy(), b_cls.to_numpy(), "basic")
        r["dataset_accession"] = study
        study_rows.append(r)
    or_table = pd.DataFrame(study_rows)

    mu_fe, se_fe, p_fe, Q, I2 = st.fe_summary(or_table)
    dl = st.dl_random_effects(or_table)
    loo_lo, loo_hi = st.loo_envelope(or_table)

    result = {
        "k": int(len(or_table)),
        "or_table": or_table,
        "fe": {
            "log_or": mu_fe,
            "se": se_fe,
            "or": float(np.exp(mu_fe)),
            "ci_low": float(np.exp(mu_fe - 1.96 * se_fe)),
            "ci_high": float(np.exp(mu_fe + 1.96 * se_fe)),
            "p": p_fe,
            "Q": Q,
            "I2": I2,
        },
        "random_effects": dl,
        "loo_envelope": [loo_lo, loo_hi],
    }
    if run_rq1:
        result["rq1"] = st.rq1_permutation(features, background, seed=seed, b_perm=b_perm, b_boot=b_boot)
    return result


# ---------------------------------------------------------------- orchestration


@dataclass
class RunReport:
    """Run provenance: executed stages, artifacts, row counts, parameters."""

    out_dir: Path
    stages: list[str]
    artifacts: dict[str, str] = field(default_factory=dict)
    counts: dict[str, int] = field(default_factory=dict)
    parameters: dict = field(default_factory=dict)
    generated: str = ""
    proteoform_regions_version: str = ""

    def to_dict(self) -> dict:
        return {
            "out_dir": str(self.out_dir),
            "stages": self.stages,
            "artifacts": self.artifacts,
            "counts": self.counts,
            "parameters": self.parameters,
            "generated": self.generated,
            "proteoform_regions_version": self.proteoform_regions_version,
        }


def _selected_stages(stages: str | list[str]) -> list[str]:
    if stages in ("all", "*", ["all"]):
        return list(STAGES)
    selected = stages if isinstance(stages, list) else [s.strip() for s in stages.split(",") if s.strip()]
    unknown = [s for s in selected if s not in STAGES]
    if unknown:
        raise ValueError(f"unknown stages {unknown}; available: {list(STAGES)}")
    return selected


def _sequences_from_snapshot(
    parents: set[str], uniprot_cache: str | Path, snapshot_date: str
) -> dict[str, str]:
    """Parent-protein sequences from the pinned UniProt snapshot (offline).

    The paper's S3/S4 read canonical sequences from the snapshot JSON - the
    properties CSV deliberately drops the sequence column. Offline only: the
    map stage (or a prior fetch) is what populates the snapshot.
    """
    from .mapping import ensure_uniprot_entries

    entries = ensure_uniprot_entries(parents, uniprot_cache, snapshot_date, offline=True)
    sequences = {
        acc: entry.get("sequence", {}).get("value", "") for acc, entry in entries.items() if acc in parents
    }
    missing = len(parents) - len(sequences)
    if missing:
        print(
            f"WARNING: {missing} parent accessions missing from snapshot {snapshot_date} "
            "(background/region sequences partial)"
        )
    return sequences


def run(
    studies: str | Path | list[StudyConfig],
    out_dir: str | Path,
    stages: str | list[str] = "all",
    data_dir: str | Path = "data/raw/peptide_tables",
    download: bool = True,
    uniprot_cache: str | Path = "data/raw/uniprot_cache",
    snapshot_date: str = "20260917",
    offline: bool = False,
    calculator: str = "bjellqvist",
    max_missed: int = 2,
    gap: int = REGION_GAP_AA,
    terminal_fraction: float = TERMINAL_FRACTION,
    min_peptides: int = MIN_REGION_PEPTIDES,
    seed: int = st.SEED,
    b_perm: int = st.B_PERM,
    b_boot: int = st.B_BOOT,
    run_rq1: bool = True,
) -> RunReport:
    """Execute pipeline stages over a cohort, writing paper-vocabulary artifacts.

    ``studies`` is a cohort YAML path or a list of ``StudyConfig``. Artifacts
    land in ``out_dir`` under the paper's ``data/processed/`` filenames; every
    run closes by writing a sha256 manifest (``manifest.sha256``) over its own
    outputs plus a machine-readable ``run_report.json``.
    """
    import proteoform_regions

    if isinstance(studies, (str, Path)):
        studies = load_cohort(studies)
    selected = _selected_stages(stages)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    report = RunReport(
        out_dir=out_dir,
        stages=selected,
        parameters={
            "calculator": calculator,
            "max_missed": max_missed,
            "gap": gap,
            "terminal_fraction": terminal_fraction,
            "min_peptides": min_peptides,
            "seed": seed,
            "b_perm": b_perm,
            "b_boot": b_boot,
            "snapshot_date": snapshot_date,
            "n_studies": len(studies),
        },
        generated=datetime.now().isoformat(timespec="minutes"),
        proteoform_regions_version=proteoform_regions.__version__,
    )
    evidence = mapping_df = features = background = regions = None
    proteins = None

    def _write(frame: pd.DataFrame, name: str, stage: str) -> Path:
        path = out_dir / name
        frame.to_csv(path, index=False)
        report.artifacts[stage] = name
        report.counts[stage] = int(len(frame))
        return path

    if "harmonize" in selected:
        evidence, _stats = harmonize(studies, data_dir, download=download)
        _write(evidence, "peptide_evidence_table.csv", "harmonize")
    elif "map" in selected or "compute" in selected or "stitch" in selected:
        evidence = pd.read_csv(out_dir / "peptide_evidence_table.csv", low_memory=False)

    if "map" in selected:
        mapping_df, proteins = _map_stage(evidence, uniprot_cache, snapshot_date, offline=offline)
        _write(mapping_df, "processed_peptide_protein_map.csv", "map")
        proteins_out = proteins.drop(columns=["sequence"], errors="ignore")
        _write(proteins_out, "processed_protein_properties.csv", "map_properties")
    elif "compute" in selected or "stitch" in selected:
        mapping_df = pd.read_csv(out_dir / "processed_peptide_protein_map.csv")

    if "compute" in selected:
        if proteins is None:
            proteins = pd.read_csv(out_dir / "processed_protein_properties.csv", low_memory=False)
        features = compute_features(evidence, mapping_df, proteins, calculator=calculator)
        _write(features, "processed_feature_table.csv", "compute")
    elif "stitch" in selected or "bias" in selected:
        features = pd.read_csv(out_dir / "processed_feature_table.csv", low_memory=False)

    if "background" in selected:
        if features is None:
            features = pd.read_csv(out_dir / "processed_feature_table.csv", low_memory=False)
        parents = set(features["mapped_accession"].dropna())
        sequences = _sequences_from_snapshot(parents, uniprot_cache, snapshot_date)
        background = build_background(sequences, max_missed=max_missed)
        _write(background, "processed_tryptic_background.csv", "background")
    elif "bias" in selected:
        background = pd.read_csv(out_dir / "processed_tryptic_background.csv")

    if "stitch" in selected:
        parents = set(mapping_df["mapped_accession"].dropna())
        precursor_sequences = _sequences_from_snapshot(parents, uniprot_cache, snapshot_date)
        regions = stitch_regions(
            features,
            mapping_df,
            precursor_sequences,
            gap=gap,
            terminal_fraction=terminal_fraction,
            min_peptides=min_peptides,
        )
        _write(regions, "processed_region_inference.csv", "stitch")

    if "bias" in selected:
        bias = quantify_bias(features, background, seed=seed, b_perm=b_perm, b_boot=b_boot, run_rq1=run_rq1)
        bias_payload = {
            "k": bias["k"],
            "fe": bias["fe"],
            "random_effects": bias["random_effects"],
            "loo_envelope": bias["loo_envelope"],
            "or_table": bias["or_table"].to_dict(orient="records"),
            "rq1": {k: v for k, v in bias.get("rq1", {}).items() if k != "per_protein"},
            "parameters": {"seed": seed, "b_perm": b_perm, "b_boot": b_boot},
        }
        (out_dir / "processed_inference_results.json").write_text(
            json.dumps(bias_payload, indent=2, default=str), encoding="utf-8"
        )
        report.artifacts["bias"] = "processed_inference_results.json"
        report.counts["bias"] = bias["k"]

    (out_dir / "run_report.json").write_text(json.dumps(report.to_dict(), indent=2), encoding="utf-8")
    outputs = [out_dir / name for name in report.artifacts.values()]
    if outputs:
        write_manifest(outputs, out_dir / "manifest.sha256")
    report.artifacts["run_report"] = "run_report.json"
    report.artifacts["manifest"] = "manifest.sha256"
    return report
