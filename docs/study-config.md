# Study configuration

`StudyConfig` decouples *what a dataset is* (metadata) from *how its files are
read* (the adapter). Cohorts are YAML; the same config drives Python and CLI
runs.

## Fields

| field | default | meaning |
|---|---|---|
| `dataset_accession` | — | required; also the on-disk folder name under `data_dir` |
| `adapter` | — | required; a [parser registry](parsers.md) key |
| `study_id` | `<accession>_<year>` | written to `study_id` column |
| `study_year` | `"2026"` | used only to build the default `study_id` |
| `parser_family` | per-adapter | value written to `source_parser_family`; override only to match a frozen vocabulary |
| `sample_type` | `""` | e.g. `serum`, `tissue_ffpe_and_plasma_sev` |
| `disease_context` | `""` | free text |
| `assay_type` | `"LC-MS/MS"` | |
| `acquisition_mode` | `""` | e.g. `DDA`, `DIA`, `DDA + FAIMS` |
| `search_engine` | `""` | e.g. `MaxQuant`, `DIA-NN`, `FragPipe/MSFragger` |
| `representative_file_arm` | `"primary result file"` | provenance label |
| `representative_file_species` | `"Homo sapiens"` | |
| `publication_doi` | `None` | |
| `files` | `[]` | list of `{name, url?, sha256?}` |
| adapter flags | see [parsers](parsers.md) | `build_usi`, `extract_protein_names`, `drop_zero_intensities`, `psm_count_zero_as_none` |

## Files

Each entry names one raw file, resolved as `data_dir/<accession>/<name>`:

```yaml
files:
  - name: peptides.txt                       # local file (already curated)
  - name: report.tsv
    url: ftp://.../report.tsv                # fetched on first run, cached
    sha256: 3f2a...                          # optional content pin
```

`pfr.ensure_downloads(studies, data_dir)` fetches missing files with a
retry-once policy; failures are documented, never fabricated.

## Cohort layout

```yaml
studies:
  - dataset_accession: PXD055218
    adapter: diann
    ...
  - dataset_accession: PXD077545
    adapter: fragpipe
    ...
```

`pfr.load_cohort(path)` accepts the `studies:` dict form or a bare YAML list.
`pfr.dump_studies(studies, path)` round-trips configs losslessly (defaults are
elided).

## The paper's k=10 tranche

[`examples/studies/cohort-paper-tranche.yaml`](https://github.com/ShoabSaadat/proteoform-regions/blob/main/examples/studies/cohort-paper-tranche.yaml)
re-declares the paper's ten cohorts exactly (including the S6 numeric-quirk
flags) — it is what the freeze-parity harness runs against.
