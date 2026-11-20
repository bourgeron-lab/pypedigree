"""The PGS legend row wraps instead of widening the legend (#68).

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
    symbols,
)
from pypedigree.carriers import CarrierScheme
from pypedigree.config import resolve_style
from pypedigree.errors import PedigreeError
from pypedigree.phenotypes import (
    check_legend,
    legend_layout,
    legend_max_width_inches,
    legend_metrics,
    legend_row_right,
    legend_rows,
    legend_slot_count,
    pgs_row_layout,
    resolve_phenotype_render,
    with_quadrants,
)

HEADER = "FID\tIID\tFather\tMother\tSex\tPhenotype"

# Four binary columns in the quadrant key, and four scores with long names:
# the issue's example, whose PGS row set the legend's width (12.2 cm at 6 pt).
FOUR = (
    f"{HEADER}\tb_autism\tb_ADHD\tb_ID\tb_epilepsy\t"
    "PGS_autism\tPGS_ADHD\tPGS_intelligence\tPGS_epilepsy\n"
    "F1\tF\t0\t0\t1\t1\t1\t2\t1\t1\t-0.3\t0.2\t0.5\t-1.0\n"
    "F1\tM\t0\t0\t2\t1\t1\t1\t1\t2\t0.1\t-0.5\t-0.2\t0.3\n"
    "F1\tC1\tF\tM\t1\t2\t2\t2\t1\t1\t1.6\t1.1\t-0.8\t0.4\n"
    "F1\tC2\tF\tM\t2\t1\t1\t0\t2\t1\t-0.2\t0.4\t0.9\t-0.1\n"
)
# The same key, one short score.
ONE_SCORE = (
    f"{HEADER}\tb_autism\tb_ADHD\tb_ID\tb_epilepsy\tPGS_ASD\n"
    "F1\tF\t0\t0\t1\t1\t1\t2\t1\t1\t-0.3\n"
    "F1\tM\t0\t0\t2\t1\t1\t1\t1\t2\t0.1\n"
    "F1\tC1\tF\tM\t1\t2\t2\t2\t1\t1\t1.6\n"
)
# Scores and nothing else: no other row to wrap at.
SCORES_ONLY = (
    f"{HEADER}\tPGS_autism\tPGS_ADHD\tPGS_intelligence\n"
    "F1\tF\t0\t0\t1\t1\t-0.3\t0.2\t0.5\n"
    "F1\tM\t0\t0\t2\t1\t0.1\t-0.5\t-0.2\n"
    "F1\tC1\tF\tM\t1\t2\t1.6\t1.1\t-0.8\n"
)
SIX_PT = {"font_size": 6.0, "min_font_size": 6.0, "pgs_scale": "zscore"}


def _layout(text: str, **style):
    pedigree = parse_pedigree_from_string(text)
    names = pedigree.phenotype_names
    st = resolve_style(PedigreeStyle(**style), names)
    render = resolve_phenotype_render(pedigree, names, pedigree.phenotype_display_modes, st)
    render = with_quadrants(render, render.band_names(names), CarrierScheme(), st)
    rows = legend_rows(names, st, render, pedigree.phenotype_display_modes, CarrierScheme())
    return legend_layout(rows, st)


def _pgs(rows):
    return next(r for r in rows if r.kind == "pgs")


class TestWrapping:
    def test_long_key_wraps_below_the_caption(self):
        rows, m = _layout(FOUR, **SIX_PT)
        layout = pgs_row_layout(_pgs(rows), m.icon_text_x, m)
        assert layout.captions == ("PGS distribution, z-score, ±2 SD marked:",)
        assert layout.lines >= 2
        # The key starts on the line below the caption, under its text.
        assert {i.line for i in layout.items} >= {1}
        first = layout.items[0]
        assert first.line == 1 and first.dash_x0 == pytest.approx(m.icon_text_x)
        assert [i.name for i in layout.items] == ["autism", "ADHD", "intelligence", "epilepsy"]

    def test_no_wider_than_the_other_rows(self):
        rows, m = _layout(FOUR, **SIX_PT)
        others = max(legend_row_right(r, m) for r in rows if r.kind != "pgs")
        assert legend_row_right(_pgs(rows), m) <= others + 1e-9

    def test_breaks_between_scores_only(self):
        """Each name follows its dash on the same line; a line ends before the width."""
        rows, m = _layout(SCORES_ONLY, legend_max_width="4cm", **SIX_PT)
        row = _pgs(rows)
        layout = pgs_row_layout(row, m.icon_text_x, m)
        assert layout.lines > 2
        assert row.wrap_right is not None
        for item in layout.items:
            assert item.dash_x0 < item.dash_x1 < item.text_x
            if item.dash_x0 > m.icon_text_x:  # not the first on its line
                width = symbols.text_width(item.name, m.fs) / symbols.UNITS_PER_POINT
                assert item.text_x + width / 72 / m.unit_x <= row.wrap_right + 1e-9

    def test_single_short_score_stays_on_one_line(self):
        rows, m = _layout(ONE_SCORE, **SIX_PT)
        row = _pgs(rows)
        assert row.lines == 1
        layout = pgs_row_layout(row, m.icon_text_x, m)
        assert layout.lines == 1 and {i.line for i in layout.items} == {0}

    def test_no_other_rows_no_wrapping(self):
        rows, _m = _layout(SCORES_ONLY, **SIX_PT)
        assert _pgs(rows).wrap_right is None and _pgs(rows).lines == 1

    def test_caption_wraps_when_too_long(self):
        rows, m = _layout(SCORES_ONLY, legend_max_width="25mm", **SIX_PT)
        row = _pgs(rows)
        layout = pgs_row_layout(row, m.icon_text_x, m)
        assert len(layout.captions) > 1
        assert " ".join(layout.captions) == row.label
        assert layout.items[0].line == len(layout.captions)

    def test_max_width_lets_the_row_be_wider(self):
        """With room, the row stays on one line however narrow the others."""
        rows, _m = _layout(FOUR, legend_max_width="30cm", **SIX_PT)
        assert _pgs(rows).lines == 1

    def test_max_width_never_wraps_before_the_other_rows(self):
        narrow, _m = _layout(FOUR, legend_max_width="1cm", **SIX_PT)
        default, _m = _layout(FOUR, **SIX_PT)
        assert _pgs(narrow).wrap_right == pytest.approx(_pgs(default).wrap_right)


class TestHeight:
    def test_extra_lines_are_counted(self):
        rows, m = _layout(SCORES_ONLY, legend_max_width="4cm", **SIX_PT)
        extra = _pgs(rows).lines - 1
        assert extra >= 2
        flat = legend_metrics(legend_slot_count(rows), PedigreeStyle(**SIX_PT))
        assert m.total_h == pytest.approx(flat.total_h + extra * m.line_h)

    def test_png_panel_is_as_tall_and_no_wider(self):
        """The 12 cm row of the issue: the panel now takes the other rows' width."""
        style = PedigreeStyle(**SIX_PT)
        fig = render_pedigree(parse_pedigree_from_string(FOUR), style=style)
        try:
            fig.canvas.draw()
            legend = fig.axes[1]
            rows, m = _layout(FOUR, **SIX_PT)
            height = legend.get_position().height * fig.get_size_inches()[1]
            assert height == pytest.approx(m.total_h * m.unit_y)
            texts = {t.get_text(): t.get_window_extent() for t in legend.texts}
            widest_other = texts["Phenotypes in the key (their colors)"].x1
            for name in ("autism", "ADHD", "intelligence", "epilepsy"):
                assert max(b.x1 for t, b in texts.items() if t == name) <= widest_other + 2
            border = max(legend.patches, key=lambda p: p.get_window_extent().width)
            assert border.get_window_extent().width / fig.dpi * 2.54 < 8.5  # cm: one column
        finally:
            plt.close(fig)

    def test_svg_wraps_at_the_same_places(self):
        rows, m = _layout(FOUR, legend_max_width=5.0, **SIX_PT)
        layout = pgs_row_layout(_pgs(rows), m.icon_text_x, m)
        svg = render_interactive_svg(
            parse_pedigree_from_string(FOUR), style=PedigreeStyle(legend_max_width=5.0, **SIX_PT)
        )
        legend = svg[svg.index('<g class="legend">') :]
        ys = {}
        for name in ("intelligence", "epilepsy", "autism", "ADHD"):
            found = re.findall(rf'<text x="([-0-9.]+)" y="([-0-9.]+)"[^>]*>{name}</text>', legend)
            ys[name] = float(found[-1][1])  # the key's, after the quadrant labels
        uy = 2.0 * m.font * 0.02  # LEG_UNIT_Y_EM * font * _PT
        for a in layout.items:
            for b in layout.items:
                expected = (b.line - a.line) * m.line_h * uy
                assert ys[b.name] - ys[a.name] == pytest.approx(expected, abs=1e-3)

    def test_png_lines_where_the_layout_says(self):
        style = PedigreeStyle(legend_max_width=5.0, **SIX_PT)
        rows, m = _layout(FOUR, legend_max_width=5.0, **SIX_PT)
        layout = pgs_row_layout(_pgs(rows), m.icon_text_x, m)
        fig = render_pedigree(parse_pedigree_from_string(FOUR), style=style)
        try:
            texts = [t for t in fig.axes[1].texts if t.get_position()[0] > m.icon_text_x - 1e-6]
            y = {t.get_text(): t.get_position()[1] for t in texts}
            caption_y = y[layout.captions[0]]
            for item in layout.items:
                assert caption_y - y[item.name] == pytest.approx(item.line * m.line_h)
        finally:
            plt.close(fig)


class TestMaxWidthOption:
    @pytest.mark.parametrize(
        ("value", "inches"),
        [
            (8.5, 8.5 / 2.54),
            ("8.5cm", 8.5 / 2.54),
            ("85mm", 8.5 / 2.54),
            ("3.3in", 3.3),
            (None, None),
        ],
    )
    def test_units(self, value, inches):
        got = legend_max_width_inches(PedigreeStyle(legend_max_width=value))
        assert got == (None if inches is None else pytest.approx(inches))

    @pytest.mark.parametrize("bad", [0, -2, "wide", "8.5 pt", True, "0cm"])
    def test_invalid(self, bad):
        with pytest.raises(PedigreeError, match="legend_max_width"):
            check_legend(PedigreeStyle(legend_max_width=bad))

    def test_from_json_style(self):
        assert PedigreeStyle.from_dict({"legend_max_width": "85mm"}).legend_max_width == "85mm"
