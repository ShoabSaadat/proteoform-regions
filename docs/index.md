# proteoform-regions

Proteoform-aware **region inference** and **detection-bias quantification**
from bottom-up proteomics peptide tables.

This package is the software implementation of the pipeline described in
*“Basic-peptide depletion and precursor masking in bottom-up proteomics”*
(Saadat et al., 2026). It takes heterogeneous search-engine peptide exports
(MaxQuant, DIA-NN, FragPipe, mzIdentML, Spectronaut, Progenesis, pGlyco3 …),
harmonizes them into one canonical evidence schema, maps peptides onto
snapshot-pinned UniProt entries, computes physicochemical features
(pI / MW / net charge), infers detected regions by single-linkage stitching,
and quantifies the acidic/basic **detection bias** against a digestion-matched
tryptic background (permutation null + odds-ratio meta-analysis).

## Why

Bottom-up proteomics does not observe proteoforms uniformly: acidic peptides
are under-detected relative to what tryptic digestion should produce. The
paper quantifies this bias (fixed-effect OR ≈ 0.45 across ten public cohorts)
and shows it propagates to region-level conclusions. This package packages
that analysis so any cohort can be re-scored with the same engine — no
notebook surgery.

## The six stages

```text
harmonize -> map -> compute -> background -> stitch -> bias
```

Each stage is importable (`proteoform_regions.run` composes them), invokable
from the CLI, and writes paper-vocabulary artifacts (`peptide_evidence_table.csv`,
`processed_feature_table.csv`, …) plus a sha256 run manifest.

## Guarantees

- **Deterministic** — every stochastic knob (seed, permutation/boot counts) is
  an explicit parameter with the paper's defaults.
- **Freeze-parity tested** — the harmonizer and bias models reproduce the
  paper's frozen artifacts and headline numbers exactly
  (FE OR 0.4519 [0.4431, 0.4608]; RQ1 T = −0.6712, p = 0.0005).
- **No silent drops** — mapping failures are tiered and flagged per row;
  downloads that fail are documented, never fabricated.
- **Offline-capable** — UniProt resolution pins to a dated snapshot;
  the engine benchmark and the walkthrough run without network.

## Where to start

- [Install](install.md) · [Quickstart](quickstart.md) ·
  [Full walkthrough: a new FragPipe dataset end-to-end](sadeghi_walkthrough.ipynb)
- Reference: [parsers](parsers.md) · [study YAML](study-config.md) ·
  [CLI](cli.md) · [Python API](api.md)
- [Methods vs paper modules](methods-mapping.md) maps each module to the
  paper's Methods and Supplementary Methods.
