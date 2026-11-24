"""A PGS score's line: solid and thin, and the legend's key drawn the same way (#158).

Dashed lines 1 pt wide read as an interval rather than one value, and their
dashes were too short to help at a figure's scale. The dot gives the color.
"""

from __future__ import annotations

import re

import matplotlib.pyplot as plt
import pytest
from matplotlib.colors import to_hex

import pypedigree as pp
from pypedigree.phenotypes import PGS_DOT_PT, PGS_MARKER_LINE_PT
from pypedigree.symbols import units_per_point

SCORES = (
    "FID\tIID\tPAT\tMAT\tSEX\tPHENO\tPGS_autism\tPGS_school_years\n"
    "F1\tfa\t0\t0\t1\t1\t0.4\t-1.0\n"
    "F1\tmo\t0\t0\t2\t1\t-0.2\t0.6\n"
    "F1\tkid\tfa\tmo\t1\t2\t2.1\t0.3\n"
)


@pytest.fixture
def scores():
    ped = pp.parse_pedigree_from_string(SCORES)
    return ped, pp.resolve_style(None, ped.phenotype_names)


def test_marker_line_is_thinner_than_before():
    assert PGS_MARKER_LINE_PT < 1.0


class TestPng:
    def test_markers_solid_and_thin(self, scores):
        ped, style = scores
        fig = pp.render_pedigree(ped)
        lines = [ln for ln in fig.axes[0].lines if ln.get_gid() == "pgs-marker"]
        colors = {to_hex(ln.get_color()) for ln in lines}
        plt.close(fig)
        assert len(lines) == 6  # two scores for each of three people
        assert {ln.get_linestyle() for ln in lines} == {"-"}
        assert {ln.get_linewidth() for ln in lines} == {PGS_MARKER_LINE_PT}
        assert colors == {c.lower() for c in style.pgs_colors.values()}

    def test_key_line_with_a_dot(self, scores):
        ped, style = scores
        fig = pp.render_pedigree(ped)
        legend = fig.axes[1]
        keys = [ln for ln in legend.lines if ln.get_gid() == "pgs-key-line"]
        dots = [ln for ln in legend.lines if ln.get_gid() == "pgs-key-dot"]
        plt.close(fig)
        assert len(keys) == len(dots) == 2
        for key, dot in zip(keys, dots, strict=True):
            assert key.get_linestyle() == "-"
            assert key.get_linewidth() == PGS_MARKER_LINE_PT
            assert to_hex(dot.get_markerfacecolor()) == to_hex(key.get_color())
            # The dot in the middle of its line.
            assert dot.get_xdata()[0] == pytest.approx(sum(key.get_xdata()) / 2)
            assert dot.get_ydata()[0] == pytest.approx(key.get_ydata()[0])


class TestSvg:
    def test_markers_solid_and_thin(self, scores):
        ped, _style = scores
        svg = pp.render_svg(ped)
        markers = re.findall(r'<line class="pgs-marker"[^>]*/>', svg)
        assert len(markers) == 6
        width = PGS_MARKER_LINE_PT * units_per_point(1.0)
        for line in markers:
            assert "dasharray" not in line
            assert float(re.search(r'stroke-width="([0-9.]+)"', line).group(1)) == pytest.approx(
                width, abs=1e-3
            )

    def test_no_dashed_line_of_a_score_color(self, scores):
        ped, style = scores
        svg = pp.render_svg(ped)
        for color in style.pgs_colors.values():
            dashed = re.findall(rf'<line [^>]*stroke="{color}"[^>]*dasharray[^>]*/>', svg)
            assert dashed == []

    def test_key_line_with_a_dot(self, scores):
        ped, style = scores
        svg = pp.render_svg(ped)
        legend = svg[svg.index('class="legend"') :]
        for color in style.pgs_colors.values():
            line = re.search(
                rf'<line class="pgs-key-line" x1="([-0-9.]+)" y1="([-0-9.]+)" '
                rf'x2="([-0-9.]+)"[^>]*stroke="{color}"',
                legend,
            )
            dot = re.search(
                rf'<circle class="pgs-key-dot" cx="([-0-9.]+)" cy="([-0-9.]+)" '
                rf'r="([0-9.]+)" fill="{color}"',
                legend,
            )
            assert line and dot
            x1, y1, x2 = (float(v) for v in line.groups())
            cx, cy, r = (float(v) for v in dot.groups())
            assert cx == pytest.approx((x1 + x2) / 2, abs=2e-3)
            assert cy == pytest.approx(y1, abs=2e-3)
            assert 0 < 2 * r < x2 - x1
            assert 2 * r <= PGS_DOT_PT * units_per_point(1.0) * 1.5
