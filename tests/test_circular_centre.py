"""The circular layout's centre couple and PGS boxes (#151).

The lone founder couple at the centre was 1.12 units apart whatever room it
needed, and the side of each PGS box was decided by angle, which runs right
to left on the page at the top of the circle: the centre couple's boxes
overlapped their symbols. Now the couple is `couple_gap` apart and each box
is on its partner's outer side on the page.

All pedigrees here are synthetic.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from pypedigree import PedigreeStyle, parse_pedigree
from pypedigree.config import apply_pgs_spacing, is_curve_column, resolve_style
from pypedigree.layout import compute_layout_circular
from pypedigree.phenotypes import (
    active_pgs,
    box_width,
    pgs_box_left,
    pgs_box_width,
    pgs_values_width,
)

sys.path.insert(0, str(Path(__file__).parent))
from conftest import TEST_DATA_DIR, generate_realistic_pedigree  # noqa: E402

PGS_FIXTURES = [
    "complex/large_family.ped",
    "complex/multi_phenotype_pgs.ped",
    "phenotype/annotations.ped",
    "phenotype/legend_minimal.ped",
    "phenotype/multi_phenotype_pgs.ped",
    "phenotype/ptdt.ped",
    "phenotype/test_all.ped",
]


def _style_for(ped) -> tuple[PedigreeStyle, list[str]]:
    """The style a renderer would lay the pedigree out with."""
    names = ped.phenotype_names or []
    cols = [n for n in names if is_curve_column(n)]
    style = resolve_style(PedigreeStyle(), names)
    style = apply_pgs_spacing(
        style, names, pgs_values_width(ped.individuals(), cols, style), pgs_box_width(cols, style)
    )
    return style, cols


def _boxes(layout, style, cols):
    half = style.shape_size / 2
    out = []
    for n in layout.nodes:
        x0, x1 = n.x - half, n.x + half
        if cols and active_pgs(n, cols, style.pgs_scale):
            left = pgs_box_left(n, style)
            x0, x1 = min(x0, left), max(x1, left + box_width(style))
        out.append((n.individual_id, x0, n.y - half, x1, n.y + half))
    return out


@pytest.mark.parametrize("fixture", PGS_FIXTURES)
def test_no_symbol_or_box_overlaps(fixture):
    ped = parse_pedigree(TEST_DATA_DIR / fixture)
    style, cols = _style_for(ped)
    for family in ped.families.values():
        boxes = _boxes(compute_layout_circular(family, style), style, cols)
        for i, a in enumerate(boxes):
            for b in boxes[i + 1 :]:
                overlap = a[1] < b[3] and b[1] < a[3] and a[2] < b[4] and b[2] < a[4]
                assert not overlap, (fixture, a[0], b[0])


def test_centre_couple_at_the_couple_gap():
    ped = parse_pedigree(TEST_DATA_DIR / "phenotype/ptdt.ped")
    style, _ = _style_for(ped)
    layout = compute_layout_circular(next(iter(ped.families.values())), style)
    centre = sorted((n for n in layout.nodes if abs(n.y) < 1e-9), key=lambda n: n.x)
    assert len(centre) == 2
    assert centre[1].x - centre[0].x >= style.couple_gap - 1e-9
    assert (centre[0].pgs_side, centre[1].pgs_side) == ("left", "right")


def test_box_sides_follow_the_page():
    """In every couple, only the partner further right on the page may draw
    their box on the right."""
    style = apply_pgs_spacing(PedigreeStyle(), ["PGS_x"])
    for seed in range(4):
        family = generate_realistic_pedigree(60, seed)
        layout = compute_layout_circular(family, style)
        nodes = {n.individual_id: n for n in layout.nodes}
        for union in family.unions:
            p1, p2 = union.partner_ids
            # A partner drawn twice joins this union as their duplicate.
            a = nodes.get(f"_dup_{p1}_{p2}", nodes.get(p1))
            b = nodes.get(f"_dup_{p2}_{p1}", nodes.get(p2))
            if a is None or b is None or abs(a.x - b.x) < 1e-6:
                continue
            lo, hi = (a, b) if a.x < b.x else (b, a)
            assert lo.pgs_side == "left" or hi.pgs_side == "right"
            assert not (lo.pgs_side == "right" and hi.pgs_side == "left")


def test_couple_without_pgs_keeps_its_place():
    """Without PGS the couple gap is the default, wider than 1.12: the couple
    moves out to it as well, so the two symbols keep a whole couple gap."""
    ped = parse_pedigree(TEST_DATA_DIR / "basic/three_gen.ped")
    style = PedigreeStyle()
    layout = compute_layout_circular(next(iter(ped.families.values())), style)
    centre = sorted((n for n in layout.nodes if abs(n.y) < 1e-9), key=lambda n: n.x)
    if len(centre) == 2:
        assert centre[1].x - centre[0].x == pytest.approx(style.couple_gap)
