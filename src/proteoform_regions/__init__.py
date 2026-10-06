"""proteoform-regions.

Proteoform-aware region inference and detection-bias quantification from
bottom-up proteomics peptide tables.

This package is the software implementation of the pipeline described in
"Basic-peptide depletion and precursor masking in bottom-up proteomics"
(Saadat et al., 2026). Modules:

- ``schema``    -- canonical evidence schema, confidence tiers, pI classification
- ``guard``     -- preview-input guards + sha256 manifest verification
- ``physchem``  -- pI / MW / net-charge engine (multiple pK models) + benchmarks
- ``digest``    -- in-silico tryptic digestion + digestion-matched backgrounds
- ``regions``   -- single-linkage detected-region stitching (gap rule)
- ``parsers``   -- search-engine peptide-table adapters + registry
- ``mapping``   -- UniProt accession mapping + protein properties
- ``stats``     -- detection-bias models (permutation null, OR meta-analysis)
- ``pipeline``  -- end-to-end orchestration
- ``cli``       -- command-line interface
"""

__version__ = "0.1.0"

__all__ = [
    "__version__",
    "StudyConfig",
    "FileSpec",
    "load_cohort",
    "harmonize",
    "ensure_downloads",
    "parse_study",
    "list_parsers",
]

from .study import FileSpec, StudyConfig, load_cohort  # noqa: E402
from .harmonize import ensure_downloads, harmonize  # noqa: E402
from .parsers import list_parsers, parse_study  # noqa: E402
