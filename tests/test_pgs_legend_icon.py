"""The legend's PGS / pTDT icon: the box as drawn in the figure, upright (#159).

The legend drew a small horizontal bell, while every box in the figure stands
upright (high values at the top, the peak on the left, #126): the reader had
to turn it in their head, and nothing said which way was high. The icon is
now the figure's box, with "high" above it, "low" below and the ±2 SD lines
named on its right.
"""

from __future__ import annotations

import itertools
import re

import matplotlib.pyplot as plt
import pytest
from conftest import TEST_DATA_DIR, decode_char_refs

import pypedigree as pp
from pypedigree.carriers import CarrierScheme
from pypedigree.config import PedigreeStyle, resolve_style
from pypedigree.errors import PedigreeError
from pypedigree.phenotypes import (
    LEG_CAPTION_GAP,
    LEG_FRAME_INSET,
    LEG_LINE_EM,
    LEG_UNIT_Y_EM,
    LEGEND_TEXT_DEFAULTS,
    MINUS,
    NBSP,
    check_legend,
    legend_layout,
    legend_row_centers,
    legend_rows,
    pgs_curve_x,
    pgs_icon_geometry,
    pgs_row_layout,
    pgs_tail_positions,
    resolve_phenotype_render,
    with_quadrants,
)

SCORES = (
    "FID\tIID\tPAT\tMAT\tSEX\tautism\tPGS_autism\tPGS_ADHD\n"
    "F1\tfa\t0\t0\t1\t1\t0.4\t-1.0\n"
    "F1\tmo\t0\t0\t2\t1\t-0.2\t0.6\n"
    "F1\tkid\tfa\tmo\t1\t2\t2.1\t0.3\n"
)
GRADIENT = (
    "FID\tIID\tPAT\tMAT\tSEX\tautism\tg_SRS\tPGS_autism\n"
    "F1\tfa\t0\t0\t1\t1\t50\t0.4\n"
    "F1\tmo\t0\t0\t2\t1\t60\t-0.2\n"
    "F1\tkid\tfa\tmo\t1\t2\t80\t2.1\n"
)
PTDT = TEST_DATA_DIR / "phenotype" / "ptdt.ped"


def _rows(text: str, **style):
    pedigree = pp.parse_pedigree_from_string(text)
    names = pedigree.phenotype_names
    st = resolve_style(PedigreeStyle(**style), names)
    render = resolve_phenotype_render(pedigree, names, pedigree.phenotype_display_modes, st)
    render = with_quadrants(render, render.band_names(names), CarrierScheme(), st)
    rows = legend_rows(names, st, render, pedigree.phenotype_display_modes, CarrierScheme())
    return legend_layout(rows, st)


def _ptdt_text() -> str:
    return PTDT.read_text(encoding="utf-8")


def _pgs(rows):
    return [r for r in rows if r.kind == "pgs"]


class TestGeometry:
    def test_upright_as_in_the_figure(self):
        """Peak on the left, baseline on the right, the same curve as beside a symbol."""
        rows, m = _rows(SCORES)
        g = pgs_icon_geometry(_pgs(rows)[0], m, 5.0)
        left, bottom, w, h = g.box
        xs = [x for x, _y in g.curve]
        ys = [y for _x, y in g.curve]
        assert (
            ys == sorted(ys)
            and ys[0] == pytest.approx(bottom)
            and ys[-1] == pytest.approx(bottom + h)
        )
        peak = xs.index(min(xs))
        assert ys[peak] == pytest.approx(bottom + h / 2, abs=h / 20)
        assert min(xs) < g.base_x and max(xs) == pytest.approx(g.base_x, abs=w / 20)
        assert left < min(xs) and g.base_x < left + w
        # The figure's curve, scaled to the icon's box.
        assert min(xs) == pytest.approx(pgs_curve_x(left, w, 1.0))

    def test_high_above_low_below_tails_on_the_right(self):
        rows, m = _rows(SCORES, pgs_scale="zscore")
        g = pgs_icon_geometry(_pgs(rows)[0], m, 5.0)
        left, bottom, w, h = g.box
        labels = {text: (x, y, anchor) for x, y, text, anchor in g.labels}
        assert labels["high"][1] > bottom + h
        assert labels["low"][1] < bottom
        assert labels["high"][2] == labels["low"][2] == "middle"
        plus, minus = f"+2{NBSP}SD", f"{MINUS}2{NBSP}SD"
        assert labels[plus][0] > left + w and labels[minus][0] > left + w
        expected = [bottom + p * h for p in pgs_tail_positions("zscore")]
        assert list(g.tails) == pytest.approx(expected)
        assert labels[minus][1] == pytest.approx(expected[0])
        assert labels[plus][1] == pytest.approx(expected[1])
        assert g.right >= labels[plus][0]

    @pytest.mark.parametrize("scale", ["zscore", "percentile"])
    def test_tails_on_fixed_scales(self, scale):
        rows, _m = _rows(SCORES, pgs_scale=scale)
        assert len(_pgs(rows)[0].tails) == 2

    def test_no_tails_on_minmax(self):
        rows, m = _rows(SCORES, pgs_scale="minmax")
        row = _pgs(rows)[0]
        assert row.tails == ()
        g = pgs_icon_geometry(row, m, 5.0)
        assert g.tails == ()
        assert sorted(text for _x, _y, text, _a in g.labels) == ["high", "low"]

    def test_ptdt_row_on_the_z_scale_whatever_pgs_scale(self):
        rows, _m = _rows(_ptdt_text(), pgs_scale="minmax")
        pgs, ptdt = _pgs(rows)
        assert pgs.tails == ()
        assert [p for p, _t in ptdt.tails] == pytest.approx(pgs_tail_positions("zscore"))
        assert ptdt.high == "high" and ptdt.low == "low"


class TestRowHeight:
    """The row holds its icon: the rows around it keep their distance."""

    @pytest.mark.parametrize("spacing", [0.6, 1.0, 3.0])
    @pytest.mark.parametrize("max_width", [None, "25mm"])
    def test_icon_inside_its_row(self, spacing, max_width):
        rows, m = _rows(_ptdt_text(), legend_spacing=spacing, legend_max_width=max_width)
        top = 20.0
        centers = legend_row_centers(rows, m, top)
        y = top
        for row, center in zip(rows, centers, strict=True):
            height = row.slots * m.row_h + (row.lines - 1) * m.line_h + row.extra_h
            if row.kind == "pgs":
                g = pgs_icon_geometry(row, m, center)
                _left, bottom, _w, h = g.box
                middle = bottom + h / 2
                assert middle + g.height / 2 <= y + 1e-9
                assert middle - g.height / 2 >= y - height - 1e-9
                # Centred on the row's text.
                last_line = center - (row.lines - 1) * m.line_h
                assert middle == pytest.approx((center + last_line) / 2)
            y -= height

    @pytest.mark.parametrize("style", [{}, {"legend": "minimal"}, {"legend_spacing": 0.6}])
    def test_clear_of_the_captions_above(self, style):
        """A gradient's Low / High hang under its bar: "high" keeps from them
        the distance a swatch keeps from the top of its row."""
        rows, m = _rows(GRADIENT, **style)
        rows = [r for r in rows if r.kind not in ("quadrants", "key_note")]
        k = next(i for i, r in enumerate(rows) if r.kind == "pgs")
        assert rows[k - 1].kind == "gradient"
        centers = legend_row_centers(rows, m, 20.0)
        caption_h = LEG_LINE_EM * m.label_fs / (LEG_UNIT_Y_EM * m.font)
        captions_bottom = centers[k - 1] - m.sw_h / 2 - LEG_CAPTION_GAP - caption_h
        g = pgs_icon_geometry(rows[k], m, centers[k])
        _left, bottom, _w, h = g.box
        icon_top = bottom + h / 2 + g.height / 2
        assert captions_bottom - icon_top >= (m.row_h - m.sw_h) / 2 - 1e-9

    def test_last_row_inside_the_frame(self):
        rows, m = _rows(SCORES)
        assert rows[-1].kind == "pgs"
        y_top = m.content_top - m.binary_header_h
        center = legend_row_centers(rows, m, y_top)[-1]
        g = pgs_icon_geometry(rows[-1], m, center)
        _left, bottom, _w, h = g.box
        assert bottom + h / 2 - g.height / 2 > LEG_FRAME_INSET


class TestWording:
    def test_roles_have_defaults(self):
        for role in ("pgs_high", "pgs_low", "pgs_plus_2sd", "pgs_minus_2sd"):
            assert LEGEND_TEXT_DEFAULTS[role]
        check_legend(PedigreeStyle(legend_text={"pgs_high": "élevé", "pgs_low": "bas"}))

    def test_misspelt_role_rejected(self):
        with pytest.raises(PedigreeError, match="pgs_high"):
            check_legend(PedigreeStyle(legend_text={"pgs_hihg": "x"}))

    def test_legend_text_renames_the_labels(self):
        rows, m = _rows(
            SCORES,
            pgs_scale="zscore",
            legend_text={
                "pgs_high": "élevé",
                "pgs_low": "bas",
                "pgs_plus_2sd": "+2 écarts-types",
                "pgs_minus_2sd": "−2 écarts-types",
            },
        )
        row = _pgs(rows)[0]
        g = pgs_icon_geometry(row, m, 5.0)
        assert {t for _x, _y, t, _a in g.labels} == {
            "élevé",
            "bas",
            "+2 écarts-types",
            "−2 écarts-types",
        }
        # Longer labels push the row's text right, past them.
        layout = pgs_row_layout(row, m.icon_text_x, m)
        assert layout.text_x > m.icon_text_x
        assert layout.text_x > g.right


def _legend_texts(fig):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    legend = fig.axes[-1]
    return [(t, t.get_window_extent(renderer)) for t in legend.texts if t.get_text().strip()]


class TestPng:
    @pytest.mark.parametrize(
        "style",
        [
            {},
            {"font_size": 6.0, "min_font_size": 6.0, "pgs_scale": "zscore"},
            {"font_size": 14.0, "pgs_scale": "percentile"},
            {"legend": "minimal"},
            {"legend_text": {"pgs_plus_2sd": "+2 standard deviations"}},
        ],
    )
    def test_labels_drawn_and_no_text_overlaps(self, style):
        st = PedigreeStyle(**style)
        fig = pp.render_pedigree(pp.parse_pedigree_from_string(SCORES), style=st)
        try:
            texts = _legend_texts(fig)
            labels = [t for t, _b in texts if t.get_gid() == "pgs-legend-label"]
            assert {"high", "low"} <= {t.get_text() for t in labels}
            for t in labels:
                assert t.get_fontsize() >= (st.min_font_size or 0) - 1e-9
            for (a, ba), (b, bb) in itertools.combinations(texts, 2):
                if "pgs-legend-label" in (a.get_gid(), b.get_gid()):
                    assert not ba.overlaps(bb), (a.get_text(), b.get_text())
        finally:
            plt.close(fig)

    def test_curve_upright(self):
        fig = pp.render_pedigree(pp.parse_pedigree_from_string(SCORES))
        try:
            (curve,) = [ln for ln in fig.axes[-1].lines if ln.get_gid() == "pgs-legend-curve"]
            xs, ys = list(curve.get_xdata()), list(curve.get_ydata())
            peak = xs.index(min(xs))
            assert ys[0] < ys[peak] < ys[-1]  # bottom to top, peak in the middle
            assert xs[0] == pytest.approx(xs[-1]) and min(xs) < xs[0]  # peak on the left
        finally:
            plt.close(fig)

    def test_ptdt_figure(self):
        fig = pp.render_pedigree(pp.parse_pedigree(PTDT))
        try:
            texts = [t.get_text() for t, _b in _legend_texts(fig)]
            assert texts.count("high") == 2 and texts.count("low") == 2
            assert texts.count(f"+2{NBSP}SD") == 2
        finally:
            plt.close(fig)


class TestSvg:
    def test_labels_and_upright_curve(self):
        svg = pp.render_svg(
            pp.parse_pedigree_from_string(SCORES), style=PedigreeStyle(pgs_scale="zscore")
        )
        svg = decode_char_refs(svg)
        legend = svg[svg.index('<g class="legend">') :]
        for label in ("high", "low", f"+2{NBSP}SD", f"{MINUS}2{NBSP}SD"):
            assert re.search(rf">{re.escape(label)}</text>", legend), label
        assert legend.count('class="pgs-legend-tail"') == 2
        points = re.search(r'<polyline class="pgs-legend-curve" points="([^"]+)"', legend)
        assert points
        xy = [tuple(map(float, p.split(","))) for p in points.group(1).split()]
        xs = [x for x, _y in xy]
        ys = [y for _x, y in xy]
        peak = xs.index(min(xs))
        # y grows downward in the SVG: the curve runs from the bottom up.
        assert ys[0] > ys[peak] > ys[-1]
        assert xs[0] == pytest.approx(xs[-1], abs=1e-3) and min(xs) < xs[0]

    def test_figure_tails_not_counted_with_the_legend(self):
        svg = pp.render_svg(
            pp.parse_pedigree_from_string(SCORES), style=PedigreeStyle(pgs_scale="zscore")
        )
        assert svg.count('class="pgs-tail"') == 2 * 3  # two per person's box

    def test_text_starts_past_the_labels(self):
        style = PedigreeStyle(legend_text={"pgs_plus_2sd": "+2 standard deviations"})
        svg = pp.render_svg(pp.parse_pedigree_from_string(SCORES), style=style)
        legend = svg[svg.index('<g class="legend">') :]
        label = re.search(r'<text x="([-0-9.]+)"[^>]*>\+2 standard deviations</text>', legend)
        caption = re.search(r'<text x="([-0-9.]+)"[^>]*>PGS[ <]', legend)
        assert label and caption
        width = len("+2 standard deviations") * 0.5 * 8 * 0.02 * 0.75  # rough, below the true
        assert float(caption.group(1)) > float(label.group(1)) + width
