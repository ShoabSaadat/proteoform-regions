# Methods vs paper modules

Where each piece of the package comes from in
*“Basic-peptide depletion and precursor masking in bottom-up proteomics”*
(Saadat et al., 2026) — analysis notebooks (S1–S8), Methods, and
Supplementary Methods Modules 1–7.

## Module map

| package module | paper source | Supp Methods |
|---|---|---|
| `schema.py`, `guard.py` | S1 harmonization (canonical schema, preview guard, freeze hashing) | — |
| `parsers/*`, `study.py`, `harmonize.py` | S1 `peptide_table_harmonization.ipynb` + S6 `tranche_extension` (nine notebook parsers → eight format adapters) | — (Methods, *Peptide-table harmonization*) |
| `mapping.py` | S2 `uniprot_mapping.ipynb` (snapshot-pinned UniProtKB, contaminant rules, exact-substring tiers) | — (Methods, *Accession mapping*) |
| `physchem.py` | S3 `physicochemical_pipeline_validation.ipynb` (Bjellqvist pI, Biopython MW, HH pK tables, charge-at-pH, GRAVY) | **Module 1** — Physicochemical engine |
| `digest.py` | S3 (Keil tryptic digestion, 7–60 aa, digestion-matched background) | **Module 2** — In silico digestion & background |
| `regions.py` | S4 `proteoform_region_inference.ipynb` (single-linkage stitching, gap 25, D1 last-member-end rule, terminal flags, ANL-020 confidence gating) | **Module 3** — Region-stitching algorithm |
| `stats.rq1_permutation` | S8 `statistical_bias_models.ipynb` | **Module 4** — RQ1 permutation null |
| `stats.or_2x2` / `fe_summary` / `bh_q` | S8 (Woolf OR, inverse-variance FE, DerSimonian–Laird RE, I², LOO envelope) | **Module 5** — Meta-analytic models (RQ2) |
| `stats.theoretical_charge_ph26` + charge-state hooks | J4/S8 charge-state layer | **Module 6** — RQ3 & charge-state statistics |
| `pipeline.py` | S2–S4 + S8 composition; `run()` mirrors the notebook flow with paper-vocabulary artifacts | — |

## Deliberately out of scope

**Supp Methods Module 7 — orthogonal-validation metrics (Barnes round).**
The boundary-concordance validation against top-down proteomics spans
(Sadeghi PXD077545) and the ProteomeTools digestion-free confirmation are
*validation analyses of the paper*, not package runtime code. Packaging them
is tracked separately; this package links to the paper for those results.
This decoupling is intentional: pipeline usage must not block on validation,
and validation must not silently change the pipeline.

## Engine parity notes

- The 16-protein ExPASy benchmark is packaged (`proteoform-regions benchmark`)
  and reproduces the paper's engine table exactly (one documented poly-E
  outlier, 0.420 pI; all others ≤ 0.005).
- `max_missed=2` reproduces the published digest enumeration, which yields
  missed-cleavage bins {0,1,2,3}; the paper's “0, 1, 2” sentence describes
  the downstream analysis strata (mc2+ folds mc=3). `max_missed=1` gives a
  strict 0–2 enumeration.
- `mapping.entry_properties`'s `reviewed` flag reproduces a publication bug
  (always False; the entryType check misses the trailing parenthetical) for
  freeze parity — use `entry_type` directly.
