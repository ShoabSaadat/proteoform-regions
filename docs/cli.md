# CLI reference

Install provides one executable: `proteoform-regions`. Every subcommand maps
to a pipeline stage; defaults equal the paper's parameters.

```console
$ proteoform-regions --help
Commands: harmonize, map, compute, background, stitch, bias, run, benchmark
```

## `run` — end-to-end

```console
proteoform-regions run --study cohort.yaml --out-dir run/ \
    [--data-dir data/raw/peptide_tables] [--uniprot-cache data/raw/uniprot_cache] \
    [--snapshot 20260917] [--offline] [--no-download] \
    [--stages all|harmonize,map,compute,background,stitch,bias] \
    [--calculator bjellqvist|hh_rodwell_emboss_lineage|hh_lehninger] \
    [--gap 25] [--max-missed 2] [--seed 42] [--b-perm 2000] [--b-boot 2000]
```

Writes paper-vocabulary artifacts, `run_report.json`, and `manifest.sha256`
into `--out-dir`.

## Stage commands (compose like the notebooks)

```console
proteoform-regions harmonize --study cohort.yaml --data-dir data/raw --out evidence.csv
proteoform-regions map      --evidence evidence.csv --uniprot-cache cache/ --offline \
                            --out-map map.csv --out-properties props.csv
proteoform-regions compute  --evidence evidence.csv --map map.csv --proteins props.csv --out features.csv
proteoform-regions background --uniprot-cache cache/ --features features.csv --out bg.csv
proteoform-regions stitch   --features features.csv --map map.csv --uniprot-cache cache/ --out regions.csv
proteoform-regions bias     --features features.csv --background bg.csv [--no-rq1] --out results.json
```

!!! warning "Sequences come from the snapshot"
    `background` and `stitch` read protein sequences from the pinned UniProt
    snapshot (or cache). Passing a properties CSV produced by `map` fails
    explicitly (it drops sequences by design, as in the paper) — point at
    `--uniprot-cache` instead.

## `benchmark` — offline engine check

```console
proteoform-regions benchmark [--live]
```

Compares the pI/MW engine against the frozen 16-protein ExPASy reference
(packaged with the wheel; `--live` re-queries ExPASy, network). Expected:
`max|d|=0.4200` — one documented poly-E outlier, everything else within
0.005 pI.

## Exit codes and logging

Stage commands print one summary line per artifact (`rows: N`) and exit 0;
input-shape problems exit 2 with an explicit message (no silent partial
output).
