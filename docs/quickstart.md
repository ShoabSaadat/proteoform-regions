# Quickstart

## 1. Describe your cohort

One YAML file, one entry per dataset — metadata lives in data, not code:

```yaml title="cohort.yaml"
studies:
  - dataset_accession: PXD055218        # your accession
    adapter: diann                       # format adapter (see parser registry)
    sample_type: serum
    disease_context: colorectal cancer vs controls
    acquisition_mode: DIA
    search_engine: DIA-NN
    files:
      - name: DIA-NN_report.tsv
        url: https://...                 # fetched once, then cached (or null)
```

## 2. Run the arc

```python
import proteoform_regions as pfr

report = pfr.run("cohort.yaml", "run/", data_dir="data/raw", uniprot_cache="data/raw/uniprot_cache")
report.counts  # rows per stage
report.artifacts  # paper-vocabulary file names
```

or from the shell:

```console
proteoform-regions run --study cohort.yaml --out-dir run/ --data-dir data/raw
```

`run()` executes `harmonize → map → compute → background → stitch → bias`
and writes each stage's artifact plus `run_report.json` and a sha256
`manifest.sha256`.

## 3. Read the answers

| Artifact | Question it answers |
|---|---|
| `peptide_evidence_table.csv` | what was detected, harmonized to 44 canonical columns |
| `processed_peptide_protein_map.csv` | which UniProt entries each peptide maps to (exact-substring, tiered) |
| `processed_feature_table.csv` | per-peptide pI / MW / net charge (Bjellqvist or HH pK tables) |
| `processed_tryptic_background.csv` | what tryptic digestion *should* produce per parent protein |
| `processed_region_inference.csv` | stitched detected regions (gap 25 aa, single linkage) |
| `processed_inference_results.json` | RQ1 permutation null + RQ2 OR meta-analysis (FE/RE/LOO) |

## Stage-by-stage instead of one shot

```python
studies = pfr.load_cohort("cohort.yaml")
evidence, stats = pfr.harmonize(studies, "data/raw")
mapping, proteins = pfr.map_to_uniprot(evidence, "data/raw/uniprot_cache")
features = pfr.compute_features(evidence, mapping, proteins, calculator="bjellqvist")
background = pfr.build_background(mapping, proteins, features)
regions = pfr.stitch_regions(features, mapping, proteins)
bias = pfr.quantify_bias(features, background, seed=42)
```

!!! note "Offline runs"
    Pass `offline=True` (CLI: `--offline`) and a populated `uniprot_cache`
    to guarantee no network access — UniProt entries are pinned to a dated
    snapshot, exactly as the paper's analysis freeze did.

Next: the full [walkthrough](sadeghi_walkthrough.ipynb) takes a brand-new
FragPipe PSM table end-to-end, offline.
