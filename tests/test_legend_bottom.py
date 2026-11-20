"""The last legend row's captions stay inside the legend's frame (#77).

A gradient's ends, a binned row's limits and category names are written
under the swatch; below the last row nothing absorbed them, so they ran past
the frame and, at small sizes, the figure's edge. All pedigrees here are
synthetic.
"""

from __future__ import annotations

import json
import re

import matplotlib.pyplot as plt
import pytest
from conftest import PHENOTYPE_DIR

from pypedigree import (
    PedigreeStyle,
    parse_pedigree_from_string,
    render_interactive_svg,
    render_pedigree,
)
from pypedigree.phenotypes import LEG_UNIT_Y_EM, legend_metrics

HEADER = "FID\tIID\tFather\tMother\tSex\tPhenotype"
# The last row is a gradient with a range, a binned column, or categories.
GRADIENT = (
    f"{HEADER}\tb_autism\tg_ADHD_RS\nFAM1\ta\t0\t0\t1\t2\t2\t30\nFAM1\tb\t0\t0\t2\t1\t1\t12\n"
)
BINNED = f"{HEADER}\tb_autism\tg_SRS\nFAM1\ta\t0\t0\t1\t2\t2\t70\nFAM1\tb\t0\t0\t2\t1\t1\t52\n"
CATEGORIES = (
    f"{HEADER}\tb_autism\tc_site\nFAM1\ta\t0\t0\t1\t2\t2\tParis\nFAM1\tb\t0\t0\t2\t1\t1\tLyon\n"
)
# No captions under the last row.
BINARY = f"{HEADER}\tb_autism\tb_ADHD\nFAM1\ta\t0\t0\t1\t2\t2\t1\nFAM1\tb\t0\t0\t2\t1\t1\t2\n"
EXTRA = {"gradient_range": {"ADHD_RS": [0, 54]}, "gradient_bins": {"SRS": [60, 66, 76]}}
STYLES = [{}, {"legend": "minimal", "font_size": 6.0, "min_font_size": 6.0}]


@pytest.mark.parametrize("text", [GRADIENT, BINNED, CATEGORIES])
@pytest.mark.parametrize("style", STYLES)
class TestInsideTheFrame:
    def test_png(self, text, style):
        fig = render_pedigree(
            parse_pedigree_from_string(text), style=PedigreeStyle(**style, **EXTRA)
        )
        try:
            fig.canvas.draw()
            legend = fig.axes[1]
            frame = max(legend.patches, key=lambda p: p.get_window_extent().width)
            bottom = frame.get_window_extent().y0
            for t in legend.texts:
                assert t.get_window_extent().y0 > bottom, t.get_text()
        finally:
            plt.close(fig)

    def test_svg(self, text, style):
        svg = render_interactive_svg(
            parse_pedigree_from_string(text), style=PedigreeStyle(**style, **EXTRA)
        )
        legend = svg[svg.index('<g class="legend">') :]
        legend = legend[: legend.index("\n</g>")]
        frame = re.search(
            r'<rect x="[-0-9.]+" y="([-0-9.]+)" width="[0-9.]+" height="([0-9.]+)"', legend
        )
        assert frame is not None
        frame_bottom = float(frame.group(1)) + float(frame.group(2))
        # Captions are hanging text: top at y (moved by the baseline shift),
        # a line below it at most.
        for y, size in re.findall(r'<text [^>]*y="([-0-9.]+)" font-size="([0-9.]+)"', legend):
            assert float(y) + float(size) * 0.3 < frame_bottom


class TestHeight:
    def test_unchanged_without_captions_under_the_last_row(self):
        for style in STYLES:
            st = PedigreeStyle(**style)
            assert legend_metrics(3, st).total_h == legend_metrics(3, st, bottom_captions=0).total_h

    def test_grows_by_what_the_captions_need(self):
        st = PedigreeStyle(legend="minimal", font_size=6.0, min_font_size=6.0)
        plain, captioned = legend_metrics(3, st), legend_metrics(3, st, bottom_captions=1)
        caption = 1.17 * captioned.caption_fs() / (LEG_UNIT_Y_EM * captioned.font)
        hang = captioned.sw_h / 2 + 0.05 + caption - captioned.row_h / 2
        assert captioned.total_h - plain.total_h == pytest.approx(max(0.0, hang + 0.2))
        assert captioned.total_h > plain.total_h

    def test_binary_last_row_keeps_the_figure_height(self):
        """A legend ending in a binary row is as tall as before #77."""
        st = PedigreeStyle()
        fig = render_pedigree(parse_pedigree_from_string(BINARY), style=st)
        try:
            height = fig.axes[1].get_position().height * fig.get_size_inches()[1]
            # The two-row quadrant key, then the shared binary row (#56).
            m = legend_metrics(3, st, binary_headers=True)
            assert height == pytest.approx(m.total_h * m.unit_y)
        finally:
            plt.close(fig)

    def test_gallery_fixture(self):
        """The issue's measurement: gradient_bins.ped, last row ADHD_RS with a range."""
        from pypedigree import parse_pedigree

        style = json.loads((PHENOTYPE_DIR / "gradient_bins.json").read_text())
        fig = render_pedigree(
            parse_pedigree(PHENOTYPE_DIR / "gradient_bins.ped"),
            style=PedigreeStyle(**style, legend="minimal", font_size=6.0, min_font_size=6.0),
        )
        try:
            fig.canvas.draw()
            legend = fig.axes[1]
            low = next(t for t in legend.texts if t.get_text() == "0").get_window_extent()
            frame = max(legend.patches, key=lambda p: p.get_window_extent().width)
            assert low.y0 > frame.get_window_extent().y0 >= 0
        finally:
            plt.close(fig)
