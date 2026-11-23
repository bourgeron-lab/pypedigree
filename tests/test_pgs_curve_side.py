"""The PGS / pTDT curve in its box: baseline on the right, peak to the left (#126).

Turned a quarter turn clockwise, it reads as the usual bell with positive values
on the right, as the legend's curve. All pedigrees here are synthetic.
"""

from __future__ import annotations

import re

import matplotlib.pyplot as plt
import pytest

from pypedigree import parse_pedigree_from_string, render_interactive_svg, render_pedigree
from pypedigree.phenotypes import PGS_CURVE_MARGIN, PGS_CURVE_REACH, pgs_curve_base, pgs_curve_x

# A couple and a child: the father's box on his left, the mother's on her right (#78).
TRIO = (
    "FID\tIID\tPAT\tMAT\tSEX\tPHENO\tPGS_x\tpTDT_x\n"
    "F\tdad\t0\t0\t1\t1\t-0.5\tNA\n"
    "F\tmum\t0\t0\t2\t1\t1.2\tNA\n"
    "F\tkid\tdad\tmum\t1\t2\tNA\t0.8\n"
)


def test_geometry():
    assert pgs_curve_base(0.0, 1.0) == pytest.approx(1 - PGS_CURVE_MARGIN)
    assert pgs_curve_x(0.0, 1.0, 0.0) == pytest.approx(1 - PGS_CURVE_MARGIN)
    assert pgs_curve_x(0.0, 1.0, 1.0) == pytest.approx(1 - PGS_CURVE_MARGIN - PGS_CURVE_REACH)


def test_svg_every_box_whatever_its_side():
    svg = render_interactive_svg(parse_pedigree_from_string(TRIO))
    boxes = re.findall(
        r'<g class="pgs">\s*<rect x="([-\d.]+)" y="[-\d.]+" width="([-\d.]+)".*?'
        r'<polyline points="([^"]+)"',
        svg,
        flags=re.S,
    )
    assert len(boxes) == 3  # two PGS curves and one pTDT curve
    for left, width, points in boxes:
        left, width = float(left), float(width)
        xs = [float(p.split(",")[0]) for p in points.split()]
        base = pgs_curve_base(left, width)
        # The tails end just left of the baseline (the bell is not 0 at ±2.8 SD).
        assert xs[0] == pytest.approx(base, abs=0.05 * width) and xs[0] <= base
        assert xs[-1] == pytest.approx(xs[0])
        assert min(xs) == pytest.approx(pgs_curve_x(left, width, 1.0), abs=1e-3)
        assert min(xs) < left + width / 2 < base, "peak left of the middle, baseline right"


def test_png_curve_peaks_left():
    fig = render_pedigree(parse_pedigree_from_string(TRIO))
    try:
        curves = [
            line
            for line in fig.axes[0].lines
            if line.get_color() == "#aaaaaa" and len(line.get_xdata()) == 50
        ]
        assert len(curves) == 3
        for line in curves:
            xs = list(line.get_xdata())
            assert min(xs) < xs[0] == pytest.approx(xs[-1]), (
                "ends on the baseline, peak to its left"
            )
    finally:
        plt.close(fig)
