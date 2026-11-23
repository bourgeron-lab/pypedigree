"""A more compact legend for journal figures (#67).

The title, the sex symbols and the spacing are options; legend "minimal"
sets them and draws a one-symbol quadrant key whose binary columns are read
from a caption under it. Every option defaults to the legend drawn before.
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
    LEG_CAPTION_GAP,
    LEG_KEY_NOTE_LABEL,
    LEG_LINE_EM,
    LEG_MINIMAL_SPACING,
    LEG_PAD,
    LEG_ROW_H,
    LEG_SEX_BELOW,
    LEG_TITLE_H,
    LEG_TOP_BARE,
    LEG_UNIT_Y_EM,
    check_legend,
    legend_metrics,
    legend_rows,
    quadrant_key_geometry,
    resolve_phenotype_render,
    with_quadrants,
)

# A quartet with one binary column: quadrants by default (#56).
ONE_BINARY = (
    "family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\tb_autism\n"
    "F1\tF\t0\t0\t1\t1\t1\n"
    "F1\tM\t0\t0\t2\t1\t1\n"
    "F1\tC1\tF\tM\t1\t2\t2\n"
    "F1\tC2\tF\tM\t2\t1\t1\n"
)
# Two binary columns, a gradient and two PGS: every kind of row the key
# can sit above.
MANY = (
    "family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\tb_autism\tb_ADHD\t"
    "g_SRS\tPGS_autism\tPGS_ADHD\n"
    "F1\tF\t0\t0\t1\t1\t1\t2\t55\t-0.3\t0.2\n"
    "F1\tM\t0\t0\t2\t1\t1\t1\t48\t0.1\t-0.5\n"
    "F1\tC1\tF\tM\t1\t2\t2\t2\t80\t1.6\t1.1\n"
    "F1\tC2\tF\tM\t2\t1\t1\t0\t50\t-0.2\t0.4\n"
)
BARE = {"legend_title": None, "legend_sex_symbols": False}


def _legend_texts(text: str, **style) -> list[str]:
    fig = render_pedigree(parse_pedigree_from_string(text), style=PedigreeStyle(**style))
    try:
        return [t.get_text() for t in fig.axes[1].texts]
    finally:
        plt.close(fig)


def _height(text: str, **style) -> float:
    fig = render_pedigree(parse_pedigree_from_string(text), style=PedigreeStyle(**style))
    try:
        return float(fig.get_size_inches()[1])
    finally:
        plt.close(fig)


def _svg_legend(text: str, **style) -> str:
    svg = render_interactive_svg(parse_pedigree_from_string(text), style=PedigreeStyle(**style))
    start = svg.index('<g class="legend">')
    return svg[start : svg.index("\n</g>", start)]


def _rows(text: str, **style):
    pedigree = parse_pedigree_from_string(text)
    names = pedigree.phenotype_names
    st = resolve_style(PedigreeStyle(**style), names)
    render = resolve_phenotype_render(pedigree, names, pedigree.phenotype_display_modes, st)
    render = with_quadrants(render, render.band_names(names), CarrierScheme(), st)
    return legend_rows(names, st, render, pedigree.phenotype_display_modes, CarrierScheme())


class TestDefaultsUnchanged:
    def test_metrics_as_before(self):
        """Title and sex row: 0.8 margin + 1.0 + 0.8 above the rows, as before #67."""
        m = legend_metrics(3, PedigreeStyle(), binary_headers=True)
        assert m.total_h == pytest.approx(1.8 + 0.5 + 3 * LEG_ROW_H + LEG_PAD)
        assert m.row_h == LEG_ROW_H
        assert m.title == "Legend"
        assert m.title_y == pytest.approx(m.total_h - LEG_PAD)
        assert m.sex_y == pytest.approx(m.title_y - 1.0)
        assert m.content_top == pytest.approx(m.sex_y - 0.8)

    def test_default_legend_has_title_and_sex_symbols(self):
        texts = _legend_texts(ONE_BINARY)
        assert "Legend" in texts and "= Male" in texts and "= Female" in texts


class TestHeaderOptions:
    def test_no_title_no_sex_symbols(self):
        texts = _legend_texts(ONE_BINARY, **BARE)
        assert "Legend" not in texts
        assert not any("Male" in t or "Female" in t for t in texts)
        svg = _svg_legend(ONE_BINARY, **BARE)
        assert ">Legend<" not in svg and "Male" not in svg and "Female" not in svg

    def test_legend_shorter_by_the_header(self):
        """The figure loses the header's height, not just its text (#67)."""
        full = legend_metrics(3, PedigreeStyle())
        bare = legend_metrics(3, PedigreeStyle(**BARE))
        header = (LEG_PAD + LEG_TITLE_H + LEG_SEX_BELOW) - LEG_TOP_BARE
        assert full.total_h - bare.total_h == pytest.approx(header)
        drop = (_height(ONE_BINARY) - _height(ONE_BINARY, **BARE)) * 2.54
        assert drop == pytest.approx(header * full.unit_y * 2.54, abs=0.01)
        assert drop > 0.5  # cm

    @pytest.mark.parametrize("title", [None, ""])
    def test_empty_title_is_no_title(self, title):
        assert legend_metrics(1, PedigreeStyle(legend_title=title)).title == ""

    def test_own_title(self):
        assert "Key" in _legend_texts(ONE_BINARY, legend_title="Key")
        assert ">Key<" in _svg_legend(ONE_BINARY, legend_title="Key")

    def test_title_without_sex_symbols(self):
        m = legend_metrics(1, PedigreeStyle(legend_sex_symbols=False))
        assert m.sex_y is None
        assert m.content_top == pytest.approx(m.title_y - LEG_TITLE_H)
        assert "Legend" in _legend_texts(ONE_BINARY, legend_sex_symbols=False)

    def test_sex_symbols_without_title(self):
        m = legend_metrics(1, PedigreeStyle(legend_title=None))
        assert m.sex_y is not None and m.sex_y < m.total_h
        assert m.content_top < m.sex_y
        texts = _legend_texts(ONE_BINARY, legend_title=None)
        assert "= Male" in texts and "Legend" not in texts


class TestSpacing:
    def test_spacing_scales_rows_and_gaps(self):
        base = legend_metrics(3, PedigreeStyle(), binary_headers=True)
        tight = legend_metrics(3, PedigreeStyle(legend_spacing=0.8), binary_headers=True)
        assert tight.row_h == pytest.approx(0.8 * base.row_h)
        assert tight.binary_header_h == pytest.approx(0.8 * base.binary_header_h)
        assert tight.total_h == pytest.approx(0.8 * base.total_h)
        assert _height(MANY, legend_spacing=0.7) < _height(MANY)

    def test_tight_rows_stop_at_their_captions(self):
        """At 0.7 the rows would leave the captions no room: they stop shorter."""
        base = legend_metrics(3, PedigreeStyle())
        tight = legend_metrics(3, PedigreeStyle(legend_spacing=0.7))
        assert 0.7 * base.row_h < tight.row_h < base.row_h

    def test_larger_spacing(self):
        assert legend_metrics(3, PedigreeStyle(legend_spacing=1.5)).row_h == pytest.approx(
            1.5 * LEG_ROW_H
        )

    @pytest.mark.parametrize("spacing", [0.6, 0.7, 1.0])
    @pytest.mark.parametrize("minimum", [None, 6.0])
    def test_rows_keep_room_for_their_captions(self, spacing, minimum):
        """A gradient's Low / High stay clear of the next row, even when tight."""
        style = PedigreeStyle(font_size=6.0, min_font_size=minimum, legend_spacing=spacing)
        m = legend_metrics(3, style)
        caption = LEG_LINE_EM * m.caption_fs() / (LEG_UNIT_Y_EM * m.font)
        # From a swatch's centre: half of it, the gap, the caption; then
        # half of the next row's swatch.
        assert m.row_h >= m.sw_h + LEG_CAPTION_GAP + caption

    def test_png_gradient_captions_clear_the_next_row(self):
        style = PedigreeStyle(font_size=6.0, min_font_size=6.0, legend_spacing=0.6)
        fig = render_pedigree(parse_pedigree_from_string(MANY), style=style)
        try:
            legend = fig.axes[1]
            fig.canvas.draw()
            low = next(t for t in legend.texts if t.get_text() == "Low").get_window_extent()
            label = next(t for t in legend.texts if t.get_text().startswith("PGS distribution"))
            assert low.y0 >= label.get_window_extent().y1 - 1  # a pixel of rounding
        finally:
            plt.close(fig)

    @pytest.mark.parametrize("bad", [0.5, 3.5, 0, -1, "tight", True])
    def test_invalid_spacing(self, bad):
        with pytest.raises(PedigreeError, match="legend_spacing"):
            check_legend(PedigreeStyle(legend_spacing=bad))

    def test_invalid_title_and_sex_symbols(self):
        with pytest.raises(PedigreeError, match="legend_title"):
            check_legend(PedigreeStyle(legend_title=3))  # type: ignore[arg-type]
        with pytest.raises(PedigreeError, match="legend_sex_symbols"):
            check_legend(PedigreeStyle(legend_sex_symbols="no"))  # type: ignore[arg-type]

    def test_rendering_checks_the_options(self):
        with pytest.raises(PedigreeError, match="legend_spacing"):
            render_pedigree(
                parse_pedigree_from_string(ONE_BINARY), style=PedigreeStyle(legend_spacing=9)
            )

    def test_from_json_style(self):
        style = PedigreeStyle.from_dict(
            {"legend_title": None, "legend_sex_symbols": False, "legend_spacing": 0.8}
        )
        assert (style.legend_title, style.legend_sex_symbols, style.legend_spacing) == (
            None,
            False,
            0.8,
        )


class TestMinimal:
    def test_sets_the_options(self):
        m = legend_metrics(1, PedigreeStyle(legend="minimal"))
        assert m.title == "" and m.sex_y is None
        assert m.spacing == LEG_MINIMAL_SPACING

    def test_options_set_by_the_user_win(self):
        style = PedigreeStyle(legend="minimal", legend_title="Key", legend_spacing=0.9)
        m = legend_metrics(1, style)
        assert m.title == "Key" and m.spacing == 0.9

    def test_one_binary_column_is_read_from_the_key(self):
        """The key alone, and a caption under it: no Affected / Unaffected row."""
        # One column fills the whole symbol under "auto" (#108): ask for quadrants.
        rows = _rows(ONE_BINARY, legend="minimal", phenotype_layout="quadrants")
        assert [r.kind for r in rows] == ["quadrants", "key_note"]
        assert rows[0].single_symbol
        texts = _legend_texts(ONE_BINARY, legend="minimal", phenotype_layout="quadrants")
        assert LEG_KEY_NOTE_LABEL in texts
        assert "Affected" not in texts and "autism" in texts
        svg = _svg_legend(ONE_BINARY, legend="minimal", phenotype_layout="quadrants")
        assert "Affected" not in svg and "white: unaffected" in svg

    def test_compact_keeps_its_rows(self):
        one = _rows(ONE_BINARY, legend="compact", phenotype_layout="quadrants")
        assert [r.kind for r in one] == ["quadrants", "binary"]
        assert [r.kind for r in _rows(ONE_BINARY, legend="compact")] == ["binary"]  # auto: whole
        kinds = [r.kind for r in _rows(MANY, legend="compact")]
        assert kinds[:2] == ["quadrants", "binary_shared"]

    def test_several_binary_columns_and_other_rows(self):
        """Binary columns get no row; a gradient and the PGS keep theirs."""
        kinds = [r.kind for r in _rows(MANY, legend="minimal")]
        assert kinds == ["quadrants", "key_note", "gradient", "pgs"]

    def test_key_has_one_symbol(self):
        rows = _rows(MANY, legend="minimal")
        m = legend_metrics(4, PedigreeStyle(legend="minimal"))
        g = quadrant_key_geometry(rows[0], m, 0.0)
        assert g.circle is None
        sx, _sy, sw, _sh = g.square
        right = [x for x, _y, _t, anchor in g.labels if anchor == "start"]
        assert right and all(x > sx + sw for x in right)
        svg = _svg_legend(MANY, legend="minimal")
        assert "legend-key-square" in svg and "legend-key-circle" not in svg
        fig = render_pedigree(
            parse_pedigree_from_string(MANY), style=PedigreeStyle(legend="minimal")
        )
        try:
            from matplotlib.patches import Ellipse

            assert not [p for p in fig.axes[1].patches if isinstance(p, Ellipse)]
        finally:
            plt.close(fig)

    def test_shorter_than_compact(self):
        for text in (ONE_BINARY, MANY):
            assert _height(text, legend="minimal") < _height(text, legend="compact")

    def test_svg_and_png_legends_same_height(self):
        """Both renderers read the same metrics: their panels are as tall."""
        for style in ({"legend": "minimal"}, BARE, {"legend_spacing": 0.7}):
            st = PedigreeStyle(**style)
            fig = render_pedigree(parse_pedigree_from_string(MANY), style=st)
            try:
                png_h = fig.axes[1].get_position().height * fig.get_size_inches()[1] * 72
            finally:
                plt.close(fig)
            svg = _svg_legend(MANY, **style)
            border = re.search(r'<rect [^>]*height="([0-9.]+)"', svg)
            assert border is not None
            m_font = legend_metrics(1, st).font
            svg_pt = float(border.group(1)) / 0.02  # layout units -> points
            # The border is the panel less 0.1 unit above and below.
            assert svg_pt + 0.2 * LEG_UNIT_Y_EM * m_font == pytest.approx(png_h, rel=0.01)
