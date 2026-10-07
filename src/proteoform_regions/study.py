"""Study configuration: decouple study metadata from format parsers.

In the source notebooks, study metadata (sample type, disease context, search
engine) was hardcoded INSIDE parser function bodies (``STUDY_META`` /
``SELECTION`` dicts). This module lifts it into a declarative ``StudyConfig``
so the format adapters become pure functions of (config, files) and adding a
new dataset is a YAML edit, not a code change.
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields
from pathlib import Path

#: adapter key -> default value written to the ``source_parser_family`` column.
#: These labels are the source publication's family vocabulary (frozen table).
DEFAULT_PARSER_FAMILY = {
    "maxquant_peptides": "maxquant_or_text_peptide_table",
    "maxquant_msms": "maxquant_or_text_peptide_table",
    "diann": "diann_report",
    "mzidentml": "mzidentml",
    "glycopeptide": "specialized_glycopeptide_csv",
    "spectronaut": "spectronaut_peptide_quant",
    "progenesis": "progenesis_peptide_csv",
}


@dataclass
class FileSpec:
    """One raw result file belonging to a study.

    ``name`` is the local filename under ``<data_dir>/<accession>/``; ``url``
    is the retrieval source (omitted for locally-curated files); ``sha256``
    pins the expected content when known.
    """

    name: str
    url: str | None = None
    sha256: str | None = None

    @classmethod
    def from_dict(cls, d: dict | str) -> FileSpec:
        if isinstance(d, str):
            return cls(name=d)
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


@dataclass
class StudyConfig:
    """Declarative description of one dataset's peptide-evidence extraction."""

    dataset_accession: str
    adapter: str
    study_id: str | None = None  # default: f"{accession}_{year}"
    study_year: str = "2026"
    parser_family: str | None = None  # default: DEFAULT_PARSER_FAMILY[adapter]
    sample_type: str = ""
    disease_context: str = ""
    assay_type: str = "LC-MS/MS"
    acquisition_mode: str = ""
    search_engine: str = ""
    representative_file_arm: str = "primary result file"
    representative_file_species: str = "Homo sapiens"
    publication_doi: str | None = None
    files: list[FileSpec] = field(default_factory=list)
    #: DIA-NN adapter: build USI examples from MS2.Scan when present. The paper's
    #: S1 PXD055218 extraction built USIs; its S6 PXD069732 extraction did not
    #: request MS2.Scan (USI stayed None). ``build_usi=False`` replicates that.
    build_usi: bool | None = None
    #: MaxQuant peptides adapter: map the ``Protein names`` column into
    #: ``protein_name_raw``. The S1 pilot parsers did; the S6 extension parsers
    #: did not (frozen rows carry empty names there). Default True (S1).
    extract_protein_names: bool = True

    def __post_init__(self) -> None:
        self.files = [FileSpec.from_dict(f) for f in self.files]
        if self.parser_family is None:
            if self.adapter not in DEFAULT_PARSER_FAMILY:
                raise ValueError(f"unknown adapter {self.adapter!r}; known: {sorted(DEFAULT_PARSER_FAMILY)}")
            self.parser_family = DEFAULT_PARSER_FAMILY[self.adapter]
        if self.study_id is None:
            self.study_id = f"{self.dataset_accession}_{self.study_year}"

    # -- serialization ---------------------------------------------------------

    @classmethod
    def from_dict(cls, d: dict) -> StudyConfig:
        known = {f.name for f in fields(cls)}
        return cls(**{k: v for k, v in d.items() if k in known})

    @classmethod
    def from_yaml(cls, path: str | Path) -> StudyConfig:
        import yaml

        return cls.from_dict(yaml.safe_load(Path(path).read_text(encoding="utf-8")))

    def to_dict(self) -> dict:
        from dataclasses import asdict

        d = asdict(self)
        if self.parser_family == DEFAULT_PARSER_FAMILY.get(self.adapter):
            d.pop("parser_family")
        if self.build_usi is None:
            d.pop("build_usi")
        return d

    def file_paths(self, data_dir: str | Path) -> list[Path]:
        base = Path(data_dir) / self.dataset_accession
        return [base / f.name for f in self.files]


def load_cohort(path: str | Path) -> list[StudyConfig]:
    """Load a cohort file: a YAML list of study dicts, or a dict with 'studies'."""
    import yaml

    payload = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    if isinstance(payload, dict):
        payload = payload["studies"]
    return [StudyConfig.from_dict(d) for d in payload]


def dump_studies(studies: list[StudyConfig], path: str | Path) -> Path:
    """Write studies as a single-cohort YAML."""
    import yaml

    Path(path).write_text(
        yaml.safe_dump({"studies": [s.to_dict() for s in studies]}, sort_keys=False, allow_unicode=True),
        encoding="utf-8",
    )
    return Path(path)
