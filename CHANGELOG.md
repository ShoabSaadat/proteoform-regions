# Changelog

All notable changes to this project are documented in this file.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and the project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

### Fixed
- `harmonize.ensure_downloads`: `_fetch` referenced a `requests` import that
  lived in the caller's scope (NameError on the first uncached download).
- Ruff/format pass over the P3-committed modules (branch was never pushed, so
  CI had not seen the drift).


## [0.1.0] - unreleased

Target: initial release accompanying the manuscript revision
("Basic-peptide depletion and precursor masking in bottom-up proteomics").
