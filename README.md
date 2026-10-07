# proteoform-regions

[![CI](https://github.com/ShoabSaadat/proteoform-regions/actions/workflows/ci.yml/badge.svg)](https://github.com/ShoabSaadat/proteoform-regions/actions/workflows/ci.yml)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue)](https://www.python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Code style: typed](https://img.shields.io/badge/style-typed-333333.svg)](https://mypy-lang.org)

**Proteoform-aware region inference and detection-bias quantification from bottom-up proteomics peptide tables.**

`proteoform-regions` is the software implementation of the pipeline described in
*"Basic-peptide depletion and precursor masking in bottom-up proteomics"* (Saadat et al., 2026).
It ingests heterogeneous bottom-up peptide-level result tables (MaxQuant, DIA-NN, Spectronaut,
Progenesis, mzIdentML, glycopeptide exports, FragPipe PSM tables), harmonizes them into a
single tiered evidence schema, computes peptide/region physicochemistry (pI, MW, net charge),
stitches detected regions from peptide evidence, and quantifies acidic/basic detection bias
with digestion-matched tryptic backgrounds.

## What this package does and adds

When scientists hunt for disease biomarkers in blood plasma, they use **bottom-up proteomics**: proteins are chopped into peptides, the peptides are measured, and findings are reported as if the whole proteins had been measured. But the instrument never saw the whole protein — it saw fragments, and which fragments are detected is biased (some ionize poorly, some are missed entirely). The "protein" in the report is therefore an annotated guess built from patchy evidence, and its stated chemistry (charge balance, size) can misrepresent the actual molecular evidence underneath.

`proteoform-regions` performs the honest re-reading. Feed it the peptide tables any major search engine produces (MaxQuant, DIA-NN, Spectronaut, Progenesis, mzIdentML, FragPipe, glycopeptide exports — 10 formats), and it:

1. **maps every detected fragment** to its exact address on the parent protein;
2. **stitches overlapping fragments into detected regions** — the stretches of the protein actually observed;
3. **computes each region's real chemistry** (pI, MW, charge) instead of inheriting the whole protein's;
4. **grades every claim by confidence** — from single-peptide to multi-peptide region to top-down-validated;
5. **quantifies detection bias** against digestion-matched tryptic backgrounds with full statistics (permutation nulls, meta-analysis).

**What it adds to science:** (i) a correction lens for the thousands of already-deposited public biomarker datasets — re-read what was reported against what was actually detected, no new experiments required; (ii) a reproducible benchmark — shipped parameters, checksum-pinned test data, and a freeze-parity suite demonstrating it reproduces the source publication's statistics exactly; (iii) methodological transparency — every previously hidden constant (stitching gap, seeds, thresholds) is an exposed, documented parameter.

**Scope honesty:** this tool detects and grades *regions* of proteins — it does not claim to identify complete proteoforms. Top-down benchmarking shows bottom-up regions cover proteoform sequence space well but do not delimit proteoform boundaries; the tool's vocabulary and confidence taxonomy enforce that distinction.

## Install

```bash
pip install proteoform-regions            # core
pip install "proteoform-regions[io]"      # + xlsx PSM tables (FragPipe exports)
pip install "proteoform-regions[plots]"   # + plotting helpers
```

## Quickstart

```python
import proteoform_regions as pfr

studies = pfr.load_cohort("cohort.yaml")  # declarative study configs
evidence, stats = pfr.harmonize(studies, "data/raw")  # tiered evidence table
mapping, proteins = pfr.map_to_uniprot(evidence)  # UniProt snapshot-pinned
features = pfr.compute_features(evidence, mapping, proteins, calculator="bjellqvist")
background = pfr.build_background(protein_sequences)  # digestion-matched tryptic peptides
regions = pfr.stitch_regions(features, mapping)  # gap-25 single-linkage regions
bias = pfr.quantify_bias(features, background)  # RQ1 permutation + RQ2 OR meta-analysis

report = pfr.run("cohort.yaml", "run/")  # all stages + artifacts + sha256 manifest
```

Or from the command line:

```console
$ proteoform-regions run --study cohort.yaml --out-dir run/
$ proteoform-regions bias --features run/processed_feature_table.csv \
    --background run/processed_tryptic_background.csv --out results.json
$ proteoform-regions benchmark    # frozen ExPASy/UniProt engine benchmark, offline
```

## Documentation

Hosted documentation (installation, quickstart, walkthrough on a new dataset, CLI + API
reference): <https://shoabsaadat.github.io/proteoform-regions>

New-dataset walkthrough (FragPipe PSMs → regions, fully offline on committed sample data):
[`docs/sadeghi_walkthrough.ipynb`](docs/sadeghi_walkthrough.ipynb) — the notebook executes
headless in CI, and its sample assets live in
[`examples/walkthrough/`](examples/walkthrough/) (50-PSM Sheet12 sample + dated UniProt
snapshot; never the 181 MB raw sheet).

## Method provenance

Every numerical default in this package (pK tables, the 7–60 aa digest length filter, the
≤25-residue region-stitching gap, seed-42 permutation settings, unconditional 0.5 zero-cell
correction) mirrors a pre-specified parameter of the source publication and is covered by
unit and freeze-parity tests. The package is the parametrized form of the paper's
Supplementary Methods.

## License

MIT — see [LICENSE](LICENSE).
