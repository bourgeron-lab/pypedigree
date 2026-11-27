"""Names for the steps of ``gradient_bins`` (#81).

All pedigrees here are synthetic.
"""

from __future__ import annotations

import json
import re

import matplotlib.pyplot as plt
import pytest
from conftest import PHENOTYPE_DIR, decode_char_refs

from pypedigree import PedigreeStyle, parse_pedigree, render_interactive_svg, render_pedigree
from pypedigree.carriers import CarrierScheme
from pypedigree.config import resolve_style
from pypedigree.errors import PedigreeError
from pypedigree.phenotypes import (
    bin_names,
    bin_text,
    binned_row_geometry,
    check_gradient_bins,
    legend_layout,
    legend_rows,
    resolve_phenotype_render,
    with_quadrants,
)
from pypedigree.symbols import UNITS_PER_POINT, text_width

FIXTURE = PHENOTYPE_DIR / "gradient_bins.ped"
# The gallery's style, without its step names: each test sets its own.
BASE = {
    k: v
    for k, v in json.loads((PHENOTYPE_DIR / "gradient_bins.json").read_text()).items()
    if k != "gradient_bin_labels"
}
NAMES = {"SRS": [None, "mild", "moderate", "severe"], "IQ": ["intellectual disability", None]}
SMALL = {"legend": "minimal", "font_size": 6.0, "min_font_size": 6.0}


def _style(**extra) -> PedigreeStyle:
    return PedigreeStyle(**{**BASE, "gradient_bin_labels": NAMES, **extra})


def _layout(style: PedigreeStyle):
    pedigree = parse_pedigree(FIXTURE)
    names = pedigree.phenotype_names
    st = resolve_style(style, names)
    render = resolve_phenotype_render(pedigree, names, pedigree.phenotype_display_modes, st)
    render = with_quadrants(render, render.band_names(names), CarrierScheme(), st)
    rows = legend_rows(names, st, render, pedigree.phenotype_display_modes, CarrierScheme())
    return legend_layout(rows, st)


class TestNames:
    def test_names_per_step(self):
        assert bin_names("SRS", _style()) == ("", "mild", "moderate", "severe")
        assert bin_names("ADHD_RS", _style()) == ()

    def test_all_unnamed_is_no_names(self):
        style = _style(gradient_bin_labels={"SRS": [None, "", None, None]})
        assert bin_names("SRS", style) == ()

    def test_tooltip_text(self):
        style = _style()
        assert bin_text("SRS", 81, style) == "≥ 76, severe"
        assert bin_text("SRS", 52, style) == "< 60"  # unnamed step: its limits only
        # A gi_ column names its steps from the lowest up, like its cut-offs.
        assert bin_text("IQ", 64, style) == "< 70, intellectual disability"

    def test_svg_tooltip(self):
        svg = decode_char_refs(render_interactive_svg(parse_pedigree(FIXTURE), style=_style()))
        assert "SRS: 81 (≥ 76, severe)" in svg
        assert "IQ: 64 (&lt; 70, intellectual disability)" in svg


class TestLegend:
    def test_row_takes_a_second_caption_line(self):
        rows, m = _layout(_style())
        srs = next(r for r in rows if r.label == "SRS")
        assert srs.bin_names == ("", "mild", "moderate", "severe") and srs.lines == 2
        plain_rows, plain = _layout(PedigreeStyle(**BASE))
        assert next(r for r in plain_rows if r.label == "SRS").lines == 1
        assert m.total_h > plain.total_h

    def test_steps_as_wide_as_their_names(self):
        rows, m = _layout(_style())
        for row in (r for r in rows if r.kind == "binned"):
            g = binned_row_geometry(row, m)
            for name in row.bin_names:
                width = text_width(name, m.caption_fs()) / UNITS_PER_POINT / 72 / m.unit_x
                assert width < g.step_w

    @pytest.mark.parametrize("extra", [{}, SMALL])
    def test_png_names_under_limits_and_clear_of_the_next_row(self, extra):
        fig = render_pedigree(parse_pedigree(FIXTURE), style=_style(**extra))
        try:
            fig.canvas.draw()
            legend = fig.axes[1]
            box = {t.get_text(): t.get_window_extent() for t in legend.texts}
            for limits, name in (("60–65", "mild"), ("≥ 76", "severe")):
                assert box[name].y1 <= box[limits].y0 + 1  # below, a pixel of rounding
            # The names end above the next row (its label, level with its
            # swatch), and every caption is inside the frame.
            assert box["severe"].y0 >= box["IQ"].y1
            frame = max(legend.patches, key=lambda p: p.get_window_extent().width)
            bottom = frame.get_window_extent().y0
            assert all(t.get_window_extent().y0 > bottom for t in legend.texts)
        finally:
            plt.close(fig)

    def test_svg_writes_the_names(self):
        svg = decode_char_refs(render_interactive_svg(parse_pedigree(FIXTURE), style=_style()))
        legend = svg[svg.index('<g class="legend">') :]
        for name in ("mild", "moderate", "severe", "intellectual disability"):
            assert re.search(rf">{name}</text>", legend)

    def test_last_row_names_inside_the_frame(self):
        """Named steps on the legend's last row: the bottom margin counts both lines."""
        style = PedigreeStyle(
            gradient_bins={"ADHD_RS": [20, 30]},
            gradient_bin_labels={"ADHD_RS": [None, "elevated", "high"]},
        )
        rows, m = _layout(style)
        assert rows[-1].label == "ADHD RS" and rows[-1].lines == 2
        fig = render_pedigree(parse_pedigree(FIXTURE), style=style)
        try:
            fig.canvas.draw()
            legend = fig.axes[1]
            frame = max(legend.patches, key=lambda p: p.get_window_extent().width)
            high = next(t for t in legend.texts if t.get_text() == "high")
            assert high.get_window_extent().y0 > frame.get_window_extent().y0
        finally:
            plt.close(fig)


class TestValidation:
    def test_column_without_bins(self):
        with pytest.raises(PedigreeError, match="has no gradient_bins"):
            check_gradient_bins(PedigreeStyle(gradient_bin_labels={"SRS": ["a", "b"]}))

    @pytest.mark.parametrize("bad", [["mild"], [None, "a", "b"], "mild", None])
    def test_one_name_per_step(self, bad):
        style = PedigreeStyle(gradient_bins={"IQ": [70]}, gradient_bin_labels={"IQ": bad})
        with pytest.raises(PedigreeError, match="one name per step"):
            check_gradient_bins(style)

    def test_text_or_none(self):
        style = PedigreeStyle(gradient_bins={"IQ": [70]}, gradient_bin_labels={"IQ": [1, None]})
        with pytest.raises(PedigreeError, match="text or None"):
            check_gradient_bins(style)

    def test_rendering_checks_them(self):
        style = PedigreeStyle(**BASE, gradient_bin_labels={"SRS": ["a"]})
        with pytest.raises(PedigreeError, match="gradient_bin_labels"):
            render_pedigree(parse_pedigree(FIXTURE), style=style)

    def test_from_json_and_resolve_copy(self):
        style = PedigreeStyle.from_dict({**BASE, "gradient_bin_labels": NAMES})
        assert style.gradient_bin_labels == NAMES
        resolved = resolve_style(style, ["SRS"])
        resolved.gradient_bin_labels["SRS"][0] = "x"
        assert style.gradient_bin_labels["SRS"][0] is None
