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
