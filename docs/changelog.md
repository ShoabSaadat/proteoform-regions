# Changelog

<!-- mkdocs will keep this page in sync with the root CHANGELOG.md at build
     time only if the include plugin is added; until then the root file is
     canonical and this page embeds the same content. -->

See the canonical [CHANGELOG.md](https://github.com/ShoabSaadat/proteoform-regions/blob/main/CHANGELOG.md).

## v0.1.0 (2026-10-07)

First public release candidate.

- **Harmonize** — eight format adapters (MaxQuant peptides/msms, DIA-NN,
  mzIdentML, FragPipe, pGlyco3/GlyPep-Quant glycopeptide, Spectronaut,
  Progenesis) behind a declarative `StudyConfig`/cohort YAML; 44-column
  canonical evidence schema with confidence tiers A–D.
- **Map** — snapshot-pinned UniProtKB resolution (chunked cache, dated
  snapshot, offline mode), contaminant rules, exact-substring tiered mapping
  with no silent drops.
- **Compute** — Bjellqvist pI (Biopython), alternative HH pK tables
  (Rodwell/EMBOSS lineage, Lehninger), MW, net charge at pH 7.4, GRAVY.
- **Background** — Keil tryptic digestion (7–60 aa, configurable missed
  cleavages) over parent proteins.
- **Stitch** — single-linkage detected-region inference, gap 25 aa,
  terminal-inclusive flags, ANL-020 confidence gating.
- **Bias** — RQ1 digestion-matched per-protein permutation null (B=2000,
  seed-aware) and RQ2 per-study Woolf OR meta-analysis (FE + DerSimonian–Laird
  RE + I² + leave-one-out envelope), BH q-values.
- **CLI** — `harmonize/map/compute/background/stitch/bias/run/benchmark` with
  the paper's parameters as defaults; offline engine benchmark.
- **Parity** — harmonizer and bias models reproduce the paper's frozen
  artifacts exactly (70,020-row evidence table across ten studies; FE OR
  0.4519 [0.4431, 0.4608]; RQ1 T = −0.6712, p = 0.0005).
- **Walkthrough** — new-dataset end-to-end notebook (Sadeghi FragPipe PSMs)
  running fully offline on committed sample data.
