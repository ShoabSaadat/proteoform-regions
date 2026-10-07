# Changelog

All notable changes to this project are documented in this file.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and the project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] — 2026-10-07

First release candidate: the full paper arc as a package + docs site.

### Added
- Package scaffold (hatchling src layout, ruff, pytest, CI) — slice P1.
- Core science modules: `guard`, `schema`, `physchem`, `digest`, `regions`
  (verbatim-first extraction, ExPASy 16-seq golden benchmark, D1
  last-member-end stitching semantics pinned) — slice P2.
- Harmonize module: parser registry (7 format adapters = 9 notebook parsers),
  `StudyConfig`/cohort YAML, `ensure_downloads`, freeze-parity harness — slice P3.
- `mapping` (S2): accession normalization, UniProtKB snapshot-pinned fetch,
  exact-substring tiered mapping, contaminant/degeneracy flags; `stats` (S8):
  BH q, Woolf OR (unconditional 0.5), FE + DerSimonian-Laird, LOO envelope,
  seed-aware RQ1 permutation null, bootstrap CIs; `pipeline`: `compute_features`,
  `build_background`, `stitch_regions`, `quantify_bias`, stage-selected
  `run()` with paper-vocabulary artifacts + sha256 run manifest; Typer CLI
  (`harmonize`/`map`/`compute`/`background`/`stitch`/`bias`/`run`/`benchmark`)
  — slice P4.
- Test-suite consolidation: committed freeze fixtures + `parity`-marked
  full-parity layer (harmonize ≡ frozen 70,020-row evidence table on all 10
  studies; bias headlines exact); MaxQuant S6 numeric-quirk flags
  (`drop_zero_intensities`, `psm_count_zero_as_none`) — slice P5.
- `fragpipe` adapter (FragPipe/MSFragger combined PSM tables, PSM-aggregated
  with protein coordinates) — slice P6.
- Docs site (MkDocs Material + mkdocstrings + mkdocs-jupyter): install,
  quickstart, parser registry, study-config reference, CLI reference, API
  reference, methods↔paper-modules mapping, changelog — slice P6.
- Walkthrough notebook (`docs/sadeghi_walkthrough.ipynb`, Colab-badged): a
  new FragPipe dataset end-to-end, fully offline on committed sample data
  (50-PSM Sheet12 sample + dated UniProt snapshot under
  `examples/walkthrough/`); executed headless in CI via nbmake — slice P6.
- Release workflows: trusted-publishing `publish.yml` (TestPyPI via
  workflow_dispatch, PyPI on `v*` tags, artifact attestations) and Pages
  `deploy-pages.yml` (workflow_dispatch until Pages is enabled) — slice P7
  prep.

### Fixed
- `harmonize.ensure_downloads`: `_fetch` referenced a `requests` import that
  lived in the caller's scope (NameError on the first uncached download).
- Ruff/format pass over the P3-committed modules (branch was never pushed, so
  CI had not seen the drift).
- CLI help tests now strip ANSI before asserting — GitHub runners force
  color, and typer/rich then style `--` and the option name as separate
  spans, so contiguous-substring assertions failed on styled bytes.


## [0.1.0] - unreleased

Target: initial release accompanying the manuscript revision
("Basic-peptide depletion and precursor masking in bottom-up proteomics").
