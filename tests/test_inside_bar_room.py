"""A taller "inside" bar keeps its distance from the row above (#94).

With names inside the steps (#90) a bar can be taller than a swatch; the row
made room for the extra below it only, so the bar reached up towards the
captions of the row above. All pedigrees here are synthetic.
"""

from __future__ import annotations

import re

import matplotlib.pyplot as plt
import pytest
from conftest import decode_char_refs

from pypedigree import (
    PedigreeStyle,
    parse_pedigree_from_string,
    render_interactive_svg,
    render_pedigree,
)
from pypedigree.phenotypes import LEG_SWATCH_EDGE

PED = (
    "FID\tIID\tFather\tMother\tSex\tPhenotype\tb_autism\tg_SRS\tgi_IQ\tADHD_RS\n"
    "FAM1\tF\t0\t0\t1\t1\t1\t52\t104\t12\n"
    "FAM1\tM\t0\t0\t2\t1\t1\t63\t96\t20\n"
    "FAM1\tC1\tF\tM\t1\t2\t2\t81\t64\t38\n"
    "FAM1\tC2\tF\tM\t2\t1\t1\t70\t80\t25\n"
)
BASE = dict(
    legend="minimal",
    font_size=6.0,
    min_font_size=6.0,
    gradient_bins={"SRS": [60, 66, 76], "IQ": [70, 85]},
    gradient_bin_labels={
        "SRS": [None, "mild", "moderate", "severe"],
        "IQ": ["ID", "borderline", None],
    },
    gradient_diverging={"ADHD_RS": {"min": 0, "mid": 27, "max": 54}},
)
# Captions under each bar, top to bottom: SRS, IQ, ADHD-RS.
CAPTIONS = [
    ("< 60", "60–65", "66–75", "≥ 76"),
    ("< 70", "70–84", "≥ 85"),
    ("0", "54"),
]
# Both binned rows inside; or SRS's names below its limits, IQ's inside.
POSITIONS = ["inside", {"IQ": "inside"}]


def _png_gaps(position) -> list[float]:
    fig = render_pedigree(
        parse_pedigree_from_string(PED),
        style=PedigreeStyle(**BASE, gradient_bin_label_position=position),
    )
    try:
        fig.canvas.draw()
        ax = fig.axes[1]
        bars = sorted(
            (
                p.get_window_extent()
                for p in ax.patches
                if p.get_facecolor()[3] == 0 and p.get_width() > 1.0
            ),
            key=lambda b: -b.y1,
        )
        assert len(bars) == 3
        texts = {t.get_text(): t.get_window_extent() for t in ax.texts}
        names = {"mild", "moderate", "severe", "ID", "borderline"}
        lowest = []
        for row in CAPTIONS:
            below = [texts[c] for c in row]
            # SRS's names under its limits count as its captions there.
            if row is CAPTIONS[0] and position != "inside":
                below += [texts[n] for n in names & {"mild", "moderate", "severe"}]
            lowest.append(min(b.y0 for b in below))
        return [lowest[0] - bars[1].y1, lowest[1] - bars[2].y1]
    finally:
        plt.close(fig)


@pytest.mark.parametrize("position", POSITIONS)
def test_png_same_gap_above_each_bar(position):
    gap_inside, gap_plain = _png_gaps(position)
    assert gap_inside > 0
    assert gap_inside == pytest.approx(gap_plain, abs=1.0)  # a pixel


@pytest.mark.parametrize("position", POSITIONS)
def test_svg_same_gap_above_each_bar(position):
    svg = render_interactive_svg(
        parse_pedigree_from_string(PED),
        style=PedigreeStyle(**BASE, gradient_bin_label_position=position),
    )
    svg = decode_char_refs(svg)
    legend = svg[svg.index('<g class="legend">') :]
    legend = legend[: legend.index("\n</g>")]
    tops = sorted(
        float(y)
        for y, w in re.findall(
            rf'<rect x="[-0-9.]+" y="([-0-9.]+)" width="([0-9.]+)" height="[0-9.]+" '
            rf'fill="none" stroke="{LEG_SWATCH_EDGE}"',
            legend,
        )
        if float(w) > 0.3
    )
    assert len(tops) == 3
    text_y = {
        t: float(y) for y, t in re.findall(r'<text [^>]*y="([-0-9.]+)"[^>]*>([^<]+)</text>', legend)
    }
    names_below = position != "inside"
    lowest = []
    for row in CAPTIONS:
        ys = [text_y[c.replace("<", "&lt;")] for c in row]
        if row is CAPTIONS[0] and names_below:
            ys += [text_y[n] for n in ("mild", "moderate", "severe")]
        lowest.append(max(ys))  # SVG y points down: the lowest caption
    # Captions are the same size: compare where they start.
    assert tops[1] - lowest[0] == pytest.approx(tops[2] - lowest[1], abs=2e-3)
