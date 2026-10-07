# Install

```console
pip install proteoform-regions
```

Requires Python 3.10–3.12. Core dependencies: `pandas`, `numpy`, `scipy`,
`biopython`, `pyteomics`, `lxml`, `pyyaml`, `requests`, `typer`.

Optional extras:

```console
pip install "proteoform-regions[plots]"   # plotting helpers (planned)
pip install "proteoform-regions[dev]"     # pytest, ruff (contributors)
pip install "proteoform-regions[docs]"    # mkdocs-material, mkdocstrings
```

## From source (development)

```console
git clone https://github.com/ShoabSaadat/proteoform-regions
cd proteoform-regions
uv sync --extra dev
uv run pytest -m "not network"
```

## Verify

```console
proteoform-regions --version
proteoform-regions benchmark   # offline ExPASy/UniProt engine check
```

`benchmark` compares the package's pI/MW engine against the frozen 16-protein
ExPASy reference table packaged with the wheel — it runs fully offline and
should print `max|d|=0.4200` (the single documented poly-E outlier; all other
proteins agree within 0.005 pI).
