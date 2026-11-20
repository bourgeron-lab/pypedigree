"""pypedigree - Generate pedigree (family tree) diagrams from pedigree files."""

import importlib
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _metadata_version
from typing import TYPE_CHECKING, Any

from pypedigree.api import PedigreeBuilder
from pypedigree.batch import family_file_names, render_families
from pypedigree.config import (
    DEFAULT_PHENOTYPE_PALETTE,
    PedigreeStyle,
    apply_pgs_spacing,
    is_pgs_column,
    is_ptdt_column,
    resolve_style,
)
from pypedigree.errors import (
    PedigreeBatchError,
    PedigreeDataWarning,
    PedigreeError,
    PedigreeLayoutError,
    PedigreeParseError,
    PedigreeStyleWarning,
    PedigreeValidationError,
    PedigreeWarning,
)
from pypedigree.layout import compute_layout, compute_layout_circular
from pypedigree.models import (
    Adoption,
    Childlessness,
    DuplicateArc,
    EdgeLayout,
    EdgeType,
    Family,
    FamilyLayout,
    Gender,
    Individual,
    LineStyle,
    NodeLayout,
    Notation,
    Pedigree,
    Phenotype,
    PhenotypeDisplayMode,
    Pregnancy,
    ReproductiveRole,
    Sex,
    TwinConnector,
    TwinType,
    Union,
)
from pypedigree.parser import (
    apply_annotations_to_pedigree,
    apply_phenotypes_to_pedigree,
    build_notation,
    export_fam,
    parse_annotation_file,
    parse_fam_file,
    parse_pedigree,
    parse_pedigree_from_dataframe,
    parse_pedigree_from_string,
    parse_phenotype_file,
    parse_plink_ped,
)
from pypedigree.svg_renderer import render_interactive_svg, render_svg

if TYPE_CHECKING:
    from pypedigree.renderer import render_family, render_pedigree

# The matplotlib renderer is imported on first use (PEP 562), so that
# `import pypedigree` for parsing, layout or SVG output does not import
# matplotlib.pyplot: slow, and global state in a threaded server (#8).
_LAZY_ATTRIBUTES = {
    "render_pedigree": "pypedigree.renderer",
    "render_family": "pypedigree.renderer",
}


def __getattr__(name: str) -> Any:
    module = _LAZY_ATTRIBUTES.get(name)
    if module is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(importlib.import_module(module), name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted({*globals(), *_LAZY_ATTRIBUTES})


try:
    __version__ = _metadata_version("pypedigree")
except PackageNotFoundError:  # pragma: no cover - running from a source tree
    __version__ = "0.0.0.dev0"

__all__ = [
    # Models
    "Individual",
    "Union",
    "Family",
    "Pedigree",
    "Sex",
    "Phenotype",
    "PhenotypeDisplayMode",
    "TwinType",
    "Notation",
    "Pregnancy",
    "Adoption",
    "Childlessness",
    "ReproductiveRole",
    "Gender",
    "NodeLayout",
    "EdgeLayout",
    "EdgeType",
    "LineStyle",
    "TwinConnector",
    "DuplicateArc",
    "FamilyLayout",
    # Parser
    "parse_pedigree",
    "parse_pedigree_from_dataframe",
    "parse_pedigree_from_string",
    "parse_phenotype_file",
    "apply_phenotypes_to_pedigree",
    "parse_annotation_file",
    "apply_annotations_to_pedigree",
    "parse_fam_file",
    "parse_plink_ped",
    "export_fam",
    "build_notation",
    # Layout
    "compute_layout",
    "compute_layout_circular",
    # Renderer
    "render_pedigree",
    "render_family",
    "render_interactive_svg",
    "render_svg",
    "render_families",
    "family_file_names",
    # Config
    "PedigreeStyle",
    "DEFAULT_PHENOTYPE_PALETTE",
    "resolve_style",
    "apply_pgs_spacing",
    "is_pgs_column",
    "is_ptdt_column",
    # API
    "PedigreeBuilder",
    # Errors
    "PedigreeError",
    "PedigreeParseError",
    "PedigreeValidationError",
    "PedigreeLayoutError",
    "PedigreeBatchError",
    # Warnings
    "PedigreeWarning",
    "PedigreeStyleWarning",
    "PedigreeDataWarning",
    # Version
    "__version__",
]
