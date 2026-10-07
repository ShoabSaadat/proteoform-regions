# Parser registry

Every input format is one adapter: a pure function
`(StudyConfig, paths) -> (rows, *stats)` emitting canonical-schema row dicts.
Adding a dataset is a YAML edit; adding a *format* is one module + one
registry entry.

```python
from proteoform_regions import list_parsers

list_parsers()
# ['diann', 'fragpipe', 'glycopeptide', 'maxquant_msms', 'maxquant_peptides',
#  'mzidentml', 'progenesis', 'spectronaut']
```

## Adapters

| adapter key | format | entity level | tier |
|---|---|---|---|
| `maxquant_peptides` | MaxQuant `peptides.txt` / `.csv` tables | peptide | B |
| `maxquant_msms` | MaxQuant `msms.txt` PSM tables (PSM-aggregated) | peptide | B |
| `diann` | DIA-NN `report.tsv` (precursor level) | peptide | A |
| `mzidentml` | mzIdentML `.mzid` / `.mzid.gz` (two-pass lxml, passThreshold PSMs) | peptide | B |
| `fragpipe` | FragPipe/MSFragger combined PSM table (PSM-aggregated, protein coordinates) | peptide | B |
| `glycopeptide` | pGlyco3 / GlyPep-Quant per-run glycopeptide records | peptidoform | A |
| `spectronaut` | Spectronaut facility peptide-quant TSV | peptide | C |
| `progenesis` | Progenesis QI peptide CSV | peptide | C |

## Adapter-specific `StudyConfig` flags

| flag | adapter | semantics |
|---|---|---|
| `extract_protein_names` | `maxquant_peptides` | map `Protein names` into `protein_name_raw` (S1 pilots: true; S6 extension: false) |
| `drop_zero_intensities` | `maxquant_peptides` | drop exact-0.0 reporters before summing (S6 extension quirk) |
| `psm_count_zero_as_none` | `maxquant_peptides` | `MS/MS Count` 0 → `None` (S6 extension quirk) |
| `build_usi` | `diann` | build USI examples from `MS2.Scan` when present |

The S6 quirk flags exist for one reason: **byte-for-byte freeze parity** with
the paper's frozen evidence table, whose S6 studies were produced with those
semantics. Defaults keep the S1 pilot semantics (zeros kept, counts kept).

## Canonical output

Every adapter emits rows with exactly the 44 canonical columns
([`schema.py`](https://github.com/ShoabSaadat/proteoform-regions/blob/main/src/proteoform_regions/schema.py)).
Study metadata never appears inside adapter code — it is seeded from
`StudyConfig` by `base_row()`.

## Testing philosophy

- **Synthetic unit fixtures** pin each adapter's structural quirks (decoy
  flags, skiprows headers, record blocks, reference graphs, PSM aggregation).
- **Freeze parity** (`tests/test_parity.py`, `parity`-marked) proves the
  adapters + harmonizer reproduce the paper's frozen 70,020-row evidence
  table on all ten studies when the paper repository is available read-only.
