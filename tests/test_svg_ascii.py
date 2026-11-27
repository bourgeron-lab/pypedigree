"""The SVG writer's output stays ASCII.

A Python string takes as many bytes per character as its widest character
needs. One "−" (U+2212) in the legend's PGS icon (1.0.0rc2) made the SVG of a
20,000-person pedigree, and every piece it was built from, take two bytes per
character: peak memory went from 515 to 744 MB. Text, tooltips and `data-*`
values now write a non-ASCII character as a reference (``&#8722;``), which an
SVG reader turns back into the same character. All pedigrees here are synthetic.
"""

from __future__ import annotations

import json
import sys
import xml.etree.ElementTree as ET

import pytest
from conftest import PHENOTYPE_DIR, TEST_DATA_DIR

from pypedigree import PedigreeStyle, parse_pedigree, render_interactive_svg, render_svg

SVG = "{http://www.w3.org/2000/svg}"

CASES = [
    # IDs with accents, in labels, tooltips and data-* attributes.
    (TEST_DATA_DIR / "format" / "unicode_ids.ped", {}),
    # The PGS legend icon ("−2 SD") and values with a minus sign.
    (PHENOTYPE_DIR / "multi_phenotype_pgs.ped", {"show_pgs_values": True}),
    # pTDT captions with no-break spaces; steps "60–65", "≥ 76".
    (PHENOTYPE_DIR / "ptdt.ped", json.loads((PHENOTYPE_DIR / "ptdt.json").read_text())),
    (
        PHENOTYPE_DIR / "gradient_bins.ped",
        json.loads((PHENOTYPE_DIR / "gradient_bins.json").read_text()),
    ),
]


@pytest.mark.parametrize(("path", "style"), CASES, ids=lambda c: getattr(c, "name", ""))
@pytest.mark.parametrize("writer", [render_interactive_svg, render_svg])
def test_output_is_ascii(path, style, writer):
    svg = writer(parse_pedigree(path), style=PedigreeStyle(**style))
    wide = sorted({c for c in svg if ord(c) > 127})
    assert not wide, f"non-ASCII characters: {wide}"
    assert sys.getsizeof(svg) < len(svg) + 100  # one byte per character


def test_text_reads_back_the_same():
    """An XML reader gives back the IDs as written in the file."""
    svg = render_interactive_svg(parse_pedigree(TEST_DATA_DIR / "format" / "unicode_ids.ped"))
    root = ET.fromstring(svg)
    ids = {g.get("data-id") for g in root.iter(f"{SVG}g") if g.get("data-id")}
    assert {"père", "mère", "enfant"} <= ids
    texts = {t.text for t in root.iter(f"{SVG}text")}
    assert "père" in texts
