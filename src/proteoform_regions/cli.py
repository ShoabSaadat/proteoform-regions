"""Command-line interface (Typer).

Subcommands mirror the pipeline stages; every Supp-Methods constant is a CLI
default, so a full run reproduces the paper's analysis by construction::

    proteoform-regions run --study cohort.yaml --out-dir run/

Individual stages compose on intermediate CSVs (harmonize -> map -> compute ->
stitch / bias), and ``benchmark`` re-runs the frozen ExPASy/UniProt engine
benchmark offline (``--live`` re-queries ExPASy; network).
"""

from __future__ import annotations

import csv
import json
from importlib.resources import files
from pathlib import Path

import pandas as pd
import typer

from . import __version__, pipeline
from .harmonize import harmonize
from .mapping import map_to_uniprot
from .physchem import biopython_mw, biopython_pi, fetch_expasy
from .study import load_cohort

app = typer.Typer(
    name="proteoform-regions",
    help="Proteoform-aware region inference and detection-bias quantification.",
    no_args_is_help=True,
    add_completion=False,
)

_EVIDENCE_HELP = "Canonical evidence table (peptide_evidence_table.csv)"
_FEATURES_HELP = "Feature table (processed_feature_table.csv)"


def _version(value: bool) -> None:
    if value:
        typer.echo(f"proteoform-regions {__version__}")
        raise typer.Exit()


@app.callback()
def root(
    version: bool = typer.Option(
        False, "--version", callback=_version, is_eager=True, help="Print version and exit."
    ),
) -> None:
    """Proteoform-aware region inference and detection-bias quantification."""


@app.command("harmonize")
def harmonize_cmd(
    study: Path = typer.Option(..., "--study", help="Cohort YAML (list of StudyConfig dicts)."),
    data_dir: Path = typer.Option(..., "--data-dir", help="Raw-file cache root (<dir>/<accession>/...)."),
    out: Path = typer.Option(..., "--out", help="Output evidence CSV path."),
    stats_out: Path | None = typer.Option(None, "--stats-out", help="Per-study parse stats CSV path."),
    no_download: bool = typer.Option(False, "--no-download", help="Never fetch; use cached files only."),
    extraction_date: str | None = typer.Option(None, help="Override extraction_date column value."),
) -> None:
    """Stage 1: harmonize heterogeneous peptide tables into the canonical schema."""
    studies = load_cohort(study)
    evidence, stats = harmonize(studies, data_dir, download=not no_download, extraction_date=extraction_date)
    out.parent.mkdir(parents=True, exist_ok=True)
    evidence.to_csv(out, index=False)
    if stats_out is not None:
        stats.to_csv(stats_out, index=False)
    typer.echo(f"evidence: {len(evidence)} rows x {len(evidence.columns)} cols -> {out}")


@app.command("map")
def map_cmd(
    evidence: Path = typer.Option(..., "--evidence", help=_EVIDENCE_HELP),
    uniprot_cache: Path = typer.Option(
        Path("data/raw/uniprot_cache"), "--uniprot-cache", help="UniProt snapshot/chunk cache dir."
    ),
    snapshot: str = typer.Option("20260917", "--snapshot", help="Snapshot date pin (YYYYMMDD)."),
    offline: bool = typer.Option(False, "--offline", help="Snapshot/cache only; no network."),
    out_map: Path = typer.Option(..., "--out-map", help="Output peptide->protein map CSV."),
    out_properties: Path = typer.Option(..., "--out-properties", help="Output protein properties CSV."),
) -> None:
    """Stage 2: map peptides to UniProt accessions (exact substring, tiered)."""
    ev = pd.read_csv(evidence, low_memory=False)
    mapping, proteins = map_to_uniprot(ev, uniprot_cache, snapshot, offline=offline)
    out_map.parent.mkdir(parents=True, exist_ok=True)
    mapping.to_csv(out_map, index=False)
    proteins.drop(columns=["sequence"], errors="ignore").to_csv(out_properties, index=False)
    exact = mapping["mapping_tier"].isin(["exact_unique", "exact_multi_position"]).mean()
    typer.echo(
        f"mapping: {len(mapping)} rows | exact rate {exact:.1%} | "
        f"{proteins['accession'].nunique()} proteins -> {out_map.name}, {out_properties.name}"
    )


@app.command()
def compute(
    evidence: Path = typer.Option(..., "--evidence", help=_EVIDENCE_HELP),
    peptide_map: Path = typer.Option(..., "--map", help="Peptide->protein map (stage 2 output)."),
    proteins: Path = typer.Option(..., "--proteins", help="Protein properties (stage 2 output)."),
    calculator: str = typer.Option(
        "bjellqvist", "--calculator", help="pI engine: bjellqvist | hh_rodwell_emboss_lineage | hh_lehninger."
    ),
    out: Path = typer.Option(..., "--out", help="Output feature table CSV."),
) -> None:
    """Stage 3a: peptide physicochemical features over exactly-mapped rows."""
    ev = pd.read_csv(evidence, low_memory=False)
    mp = pd.read_csv(peptide_map)
    pr = pd.read_csv(proteins, low_memory=False)
    features = pipeline.compute_features(ev, mp, pr, calculator=calculator)
    out.parent.mkdir(parents=True, exist_ok=True)
    features.to_csv(out, index=False)
    typer.echo(f"features: {len(features)} rows -> {out}")


@app.command()
def background(
    proteins: Path | None = typer.Option(
        None, "--proteins", help="Protein properties incl. sequence column (alternative to --uniprot-cache)."
    ),
    uniprot_cache: Path | None = typer.Option(
        None, "--uniprot-cache", help="UniProt snapshot cache dir (sequences from the pinned snapshot)."
    ),
    snapshot: str = typer.Option("20260917", "--snapshot", help="Snapshot date pin (YYYYMMDD)."),
    features: Path | None = typer.Option(None, "--features", help=_FEATURES_HELP + " - restricts parents."),
    max_missed: int = typer.Option(2, "--max-missed", help="Missed-cleavage enumeration bound."),
    out: Path = typer.Option(..., "--out", help="Output tryptic background CSV."),
) -> None:
    """Stage 3b: digestion-matched tryptic background over parent proteins."""
    if proteins is not None:
        pr = pd.read_csv(proteins, low_memory=False)
        if "sequence" not in pr.columns:
            typer.echo(
                "FATAL: protein properties lack a `sequence` column; use --uniprot-cache/--snapshot "
                "or a properties file that keeps sequences.",
                err=True,
            )
            raise typer.Exit(code=2)
        sequences = {a: s for a, s in zip(pr["accession"], pr["sequence"], strict=True) if isinstance(s, str)}
    elif uniprot_cache is not None:
        from .mapping import ensure_uniprot_entries

        if features is None:
            typer.echo("FATAL: --uniprot-cache requires --features (to know the parent proteins).", err=True)
            raise typer.Exit(code=2)
        parents = set(pd.read_csv(features, low_memory=False)["mapped_accession"].dropna())
        entries = ensure_uniprot_entries(parents, uniprot_cache, snapshot, offline=True)
        sequences = {a: e.get("sequence", {}).get("value", "") for a, e in entries.items() if a in parents}
    else:
        typer.echo("FATAL: pass --proteins or --uniprot-cache.", err=True)
        raise typer.Exit(code=2)
    if features is not None and proteins is not None:
        parents = set(pd.read_csv(features, low_memory=False)["mapped_accession"].dropna())
        sequences = {a: s for a, s in sequences.items() if a in parents and isinstance(s, str) and s}
    bg = pipeline.build_background(sequences, max_missed=max_missed)
    out.parent.mkdir(parents=True, exist_ok=True)
    bg.to_csv(out, index=False)
    typer.echo(f"background: {len(bg)} rows over {bg['accession'].nunique()} proteins -> {out}")


@app.command()
def stitch(
    features: Path = typer.Option(..., "--features", help=_FEATURES_HELP),
    peptide_map: Path = typer.Option(..., "--map", help="Peptide->protein map (positions)."),
    proteins: Path | None = typer.Option(None, "--proteins", help="Protein properties incl. sequences."),
    uniprot_cache: Path | None = typer.Option(
        None,
        "--uniprot-cache",
        help="UniProt snapshot cache dir (alternative to --proteins).",
    ),
    snapshot: str = typer.Option("20260917", "--snapshot", help="Snapshot date pin (YYYYMMDD)."),
    gap: int = typer.Option(25, "--gap", help="Max residue gap joining spans (paper: 25)."),
    terminal_fraction: float = typer.Option(
        0.10, "--terminal-fraction", help="N/C-terminal rule (paper: 0.10)."
    ),
    min_peptides: int = typer.Option(2, "--min-peptides", help="region_supported threshold (paper: 2)."),
    out: Path = typer.Option(..., "--out", help="Output region table CSV."),
) -> None:
    """Stage 4: single-linkage detected-region inference."""
    ft = pd.read_csv(features, low_memory=False)
    mp = pd.read_csv(peptide_map)
    precursor_sequences: dict[str, str] = {}
    if proteins is not None:
        pr = pd.read_csv(proteins, low_memory=False)
        if "sequence" in pr.columns:
            precursor_sequences = {
                a: s for a, s in zip(pr["accession"], pr["sequence"], strict=True) if isinstance(s, str)
            }
    elif uniprot_cache is not None:
        from .mapping import ensure_uniprot_entries

        parents = set(mp["mapped_accession"].dropna())
        entries = ensure_uniprot_entries(parents, uniprot_cache, snapshot, offline=True)
        precursor_sequences = {
            a: e.get("sequence", {}).get("value", "") for a, e in entries.items() if a in parents
        }
    regions = pipeline.stitch_regions(
        ft,
        mp,
        precursor_sequences,
        gap=gap,
        terminal_fraction=terminal_fraction,
        min_peptides=min_peptides,
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    regions.to_csv(out, index=False)
    n_supported = int((regions["confidence_label"] == "region_supported").sum()) if len(regions) else 0
    typer.echo(f"regions: {len(regions)} rows ({n_supported} region_supported) -> {out}")


@app.command()
def bias(
    features: Path = typer.Option(..., "--features", help=_FEATURES_HELP),
    background: Path = typer.Option(..., "--background", help="Tryptic background (stage 3b output)."),
    seed: int = typer.Option(42, "--seed", help="RNG seed (paper: 42)."),
    b_perm: int = typer.Option(2000, "--b-perm", help="Per-protein permutation draws (paper: 2000)."),
    b_boot: int = typer.Option(2000, "--b-boot", help="Bootstrap draws (paper: 2000)."),
    no_rq1: bool = typer.Option(
        False, "--no-rq1", help="Skip the RQ1 permutation null (slow on full tables)."
    ),
    out: Path = typer.Option(..., "--out", help="Output inference results JSON."),
) -> None:
    """Stage 8 (S8): RQ1 permutation null + RQ2 basic-class OR meta-analysis."""
    ft = pd.read_csv(features, low_memory=False)
    bg = pd.read_csv(background)
    result = pipeline.quantify_bias(ft, bg, seed=seed, b_perm=b_perm, b_boot=b_boot, run_rq1=not no_rq1)
    payload = {
        "k": result["k"],
        "fe": result["fe"],
        "random_effects": result["random_effects"],
        "loo_envelope": result["loo_envelope"],
        "or_table": result["or_table"].to_dict(orient="records"),
        "rq1": {k: v for k, v in result.get("rq1", {}).items() if k != "per_protein"},
        "parameters": {"seed": seed, "b_perm": b_perm, "b_boot": b_boot},
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    fe = result["fe"]
    typer.echo(
        f"FE OR={fe['or']:.3f} [{fe['ci_low']:.3f}, {fe['ci_high']:.3f}] "
        f"p={fe['p']:.3g} | k={result['k']} | LOO [{result['loo_envelope'][0]:.3f}, "
        f"{result['loo_envelope'][1]:.3f}] -> {out}"
    )


@app.command()
def run(
    study: Path = typer.Option(..., "--study", help="Cohort YAML."),
    out_dir: Path = typer.Option(..., "--out-dir", help="Run directory (all artifacts + manifest)."),
    data_dir: Path = typer.Option(Path("data/raw/peptide_tables"), "--data-dir"),
    stages: str = typer.Option(
        "all", "--stages", help="Comma list or 'all': harmonize,map,compute,background,stitch,bias."
    ),
    uniprot_cache: Path = typer.Option(Path("data/raw/uniprot_cache"), "--uniprot-cache"),
    snapshot: str = typer.Option("20260917", "--snapshot"),
    no_download: bool = typer.Option(False, "--no-download"),
    offline: bool = typer.Option(False, "--offline", help="UniProt snapshot/cache only."),
    calculator: str = typer.Option("bjellqvist", "--calculator"),
    gap: int = typer.Option(25, "--gap"),
    terminal_fraction: float = typer.Option(0.10, "--terminal-fraction"),
    min_peptides: int = typer.Option(2, "--min-peptides"),
    seed: int = typer.Option(42, "--seed"),
    b_perm: int = typer.Option(2000, "--b-perm"),
    b_boot: int = typer.Option(2000, "--b-boot"),
    no_rq1: bool = typer.Option(False, "--no-rq1"),
) -> None:
    """End-to-end run: selected stages + paper-vocabulary artifacts + sha256 manifest."""
    report = pipeline.run(
        study,
        out_dir,
        stages=stages,
        data_dir=data_dir,
        download=not no_download,
        uniprot_cache=uniprot_cache,
        snapshot_date=snapshot,
        offline=offline,
        calculator=calculator,
        gap=gap,
        terminal_fraction=terminal_fraction,
        min_peptides=min_peptides,
        seed=seed,
        b_perm=b_perm,
        b_boot=b_boot,
        run_rq1=not no_rq1,
    )
    typer.echo(f"run complete: stages {report.stages} -> {report.out_dir}")
    for stage, artifact in report.artifacts.items():
        typer.echo(f"  {stage:18s} {artifact}")


@app.command()
def benchmark(
    live: bool = typer.Option(False, "--live", help="Re-query ExPASy online (slow; network)."),
    golden: Path | None = typer.Option(
        None, "--golden", help="Benchmark CSV (default: packaged frozen table)."
    ),
) -> None:
    """Engine benchmark: package pI/MW vs the frozen ExPASy/UniProt reference."""
    source = golden if golden is not None else files("proteoform_regions") / "data" / "physchem_golden.csv"
    rows = list(csv.DictReader(source.open()))
    pi_devs, mw_ppms = [], []
    live_rows = []
    for row in rows:
        seq = row["sequence"]
        pi = biopython_pi(seq)
        mw = biopython_mw(seq)
        if row.get("expasy_pi"):
            pi_devs.append(abs(pi - float(row["expasy_pi"])))
        if row.get("uniprot_mw") and row["kind"] == "protein":
            mw_ppms.append(abs(mw - float(row["uniprot_mw"])) / float(row["uniprot_mw"]) * 1e6)
        if live:
            expasy_pi, _expasy_mw = fetch_expasy(seq)
            live_rows.append((row["name"], pi, expasy_pi, pi - expasy_pi))
    typer.echo(
        f"benchmark rows: {len(rows)} | pI vs ExPASy: n={len(pi_devs)} "
        f"max|d|={max(pi_devs):.4f} median|d|={sorted(pi_devs)[len(pi_devs) // 2]:.4f}"
    )
    typer.echo(f"MW vs UniProt (full proteins): n={len(mw_ppms)} max|ppm|={max(mw_ppms):.1f}")
    if live:
        for name, pi, expasy_pi, dev in live_rows:
            typer.echo(f"  {name}: package {pi:.3f} | ExPASy {expasy_pi:.3f} | d={dev:+.3f}")


def main() -> None:  # console_scripts entry point
    app()


if __name__ == "__main__":
    main()
