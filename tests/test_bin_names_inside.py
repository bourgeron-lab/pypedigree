"""Step names inside the legend's colored steps (#90).

All pedigrees here are synthetic.
"""

from __future__ import annotations

import re

import matplotlib.pyplot as plt
import pytest

from pypedigree import (
    PedigreeStyle,
    parse_pedigree_from_string,
    render_interactive_svg,
    render_pedigree,
)
from pypedigree.carriers import CarrierScheme
from pypedigree.config import resolve_style
from pypedigree.errors import PedigreeError
from pypedigree.phenotypes import (
    LEG_X0,
    binned_row_geometry,
    check_gradient_bins,
    legend_layout,
    legend_rows,
    resolve_phenotype_render,
    text_color_on,
    to_hex,
    with_quadrants,
)

PED = (
    "FID\tIID\tFather\tMother\tSex\tPhenotype\tb_autism\tg_SRS\tgi_IQ\tADHD_RS\n"
    "FAM1\tF\t0\t0\t1\t1\t1\t52\t104\t12\n"
    "FAM1\tM\t0\t0\t2\t1\t1\t63\t96\t20\n"
    "FAM1\tC1\tF\tM\t1\t2\t2\t81\t64\t38\n"
    "FAM1\tC2\tF\tM\t2\t1\t1\t70\t80\t25\n"
)
BASE = dict(
    gradient_bins={"SRS": [60, 66, 76], "IQ": [70, 85]},
    gradient_bin_labels={
        "SRS": [None, "mild", "moderate", "severe"],
        "IQ": ["ID", "borderline", None],
    },
    gradient_diverging={"ADHD_RS": {"min": 0, "mid": 27, "max": 54}},
)
SMALL = {"legend": "minimal", "font_size": 6.0, "min_font_size": 6.0}
NAMES = ("mild", "moderate", "severe", "ID", "borderline")


def _style(**extra) -> PedigreeStyle:
    return PedigreeStyle(**{**BASE, **extra})


def _layout(style: PedigreeStyle):
    pedigree = parse_pedigree_from_string(PED)
    names = pedigree.phenotype_names
    st = resolve_style(style, names)
    render = resolve_phenotype_render(pedigree, names, pedigree.phenotype_display_modes, st)
    render = with_quadrants(render, render.band_names(names), CarrierScheme(), st)
    rows = legend_rows(names, st, render, pedigree.phenotype_display_modes, CarrierScheme())
    return legend_layout(rows, st)


@pytest.mark.parametrize("extra", [{}, SMALL])
class TestInside:
    def test_png_names_inside_their_steps_in_a_reading_color(self, extra):
        style = _style(gradient_bin_label_position="inside", **extra)
        rows, m = _layout(style)
        fig = render_pedigree(parse_pedigree_from_string(PED), style=style)
        try:
            fig.canvas.draw()
            ax = fig.axes[1]
            texts = {t.get_text(): t for t in ax.texts}
            for row in (r for r in rows if r.kind == "binned"):
                g = binned_row_geometry(row, m)
                for i, (rgb, name) in enumerate(zip(row.colors, row.bin_names, strict=True)):
                    if not name:
                        continue
                    t = texts[name]
                    x0 = LEG_X0 + i * g.step_w
                    to_data = ax.transData.inverted()
                    box = t.get_window_extent()
                    (bx0, by0), (bx1, by1) = to_data.transform([(box.x0, box.y0), (box.x1, box.y1)])
                    cx, cy = t.get_position()
                    assert x0 <= bx0 and bx1 <= x0 + g.step_w, name
                    assert cy - g.bar_h / 2 <= by0 and by1 <= cy + g.bar_h / 2 + 1e-9, name
                    assert t.get_color() == text_color_on(to_hex(rgb)), name
        finally:
            plt.close(fig)

    def test_svg_names_inside_their_steps(self, extra):
        style = _style(gradient_bin_label_position="inside", **extra)
        svg = render_interactive_svg(parse_pedigree_from_string(PED), style=style)
        legend = svg[svg.index('<g class="legend">') :]
        for name in NAMES:
            assert re.search(
                rf'<text [^>]*fill="(#ffffff|white|#[0-9a-f]{{6}})"[^>]*>{name}</text>', legend
            ), name

    def test_shorter_than_below(self, extra):
        _rows, inside = _layout(_style(gradient_bin_label_position="inside", **extra))
        _rows, below = _layout(_style(**extra))
        assert inside.total_h < below.total_h

    def test_one_caption_line(self, extra):
        rows, m = _layout(_style(gradient_bin_label_position="inside", **extra))
        for row in (r for r in rows if r.kind == "binned"):
            assert row.lines == 1 and row.names_inside


def test_contrasting_colors_on_light_and_dark_steps():
    """With a dark column color: white text on the darkest step, dark text on
    the white first step, in the PNG and the SVG."""
    style = _style(
        gradient_bin_label_position="inside",
        gradient_bin_labels={"SRS": ["normal", "mild", "moderate", "severe"]},
        phenotype_colors={"SRS": "#1f4e79"},
    )
    fig = render_pedigree(parse_pedigree_from_string(PED), style=style)
    try:
        colors = {t.get_text(): t.get_color() for t in fig.axes[1].texts}
    finally:
        plt.close(fig)
    assert colors["severe"] == "white" and colors["normal"] != "white"
    svg = render_interactive_svg(parse_pedigree_from_string(PED), style=style)
    assert re.search(r'fill="white"[^>]*>severe</text>', svg)
    assert not re.search(r'fill="white"[^>]*>normal</text>', svg)


def test_bar_grows_to_hold_the_names():
    rows, m = _layout(_style(gradient_bin_label_position="inside", **SMALL))
    srs = next(r for r in rows if r.label == "SRS")
    g = binned_row_geometry(srs, m)
    assert g.bar_h > m.sw_h and srs.extra_h == pytest.approx(
        g.bar_h - m.sw_h
    )  # all of it, half above (#94)
    # At the default size a swatch is nearly a line of names tall: it grows
    # far less.
    rows, d = _layout(_style(gradient_bin_label_position="inside"))
    srs_d = next(r for r in rows if r.label == "SRS")
    assert binned_row_geometry(srs_d, d).bar_h - d.sw_h < g.bar_h - m.sw_h


def test_captions_clear_of_the_next_row_and_the_frame():
    style = _style(gradient_bin_label_position="inside", **SMALL)
    fig = render_pedigree(parse_pedigree_from_string(PED), style=style)
    try:
        fig.canvas.draw()
        ax = fig.axes[1]
        box = {t.get_text(): t.get_window_extent() for t in ax.texts}
        # The SRS limits end above the IQ bar's names, IQ's above ADHD-RS's label.
        assert box["66–75"].y0 > box["borderline"].y1
        assert box["70–84"].y0 > box["ADHD RS"].y1
        frame = max(ax.patches, key=lambda p: p.get_window_extent().width).get_window_extent()
        assert all(b.y0 > frame.y0 for b in box.values())
    finally:
        plt.close(fig)


def test_by_column():
    rows, _m = _layout(_style(gradient_bin_label_position={"SRS": "inside"}))
    by_label = {r.label: r for r in rows}
    assert by_label["SRS"].names_inside and not by_label["IQ"].names_inside
    assert by_label["IQ"].lines == 2


@pytest.mark.parametrize("bad", ["above", {"SRS": "over"}, 3])
def test_invalid(bad):
    with pytest.raises(PedigreeError, match="gradient_bin_label_position"):
        check_gradient_bins(_style(gradient_bin_label_position=bad))
