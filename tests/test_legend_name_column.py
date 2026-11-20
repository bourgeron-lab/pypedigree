"""The names of the legend's swatch and bar rows start in one column (#88).

Since #81 a binned row with named steps can be longer than the usual name
column, and its name followed its own bar, so the names were staggered.
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
from pypedigree.phenotypes import (
    binned_row_geometry,
    legend_layout,
    legend_row_right,
    legend_rows,
    resolve_phenotype_render,
    with_quadrants,
)

PED = (
    "FID\tIID\tFather\tMother\tSex\tPhenotype\tb_autism\tg_SRS\tgi_IQ\tADHD_RS\n"
    "FAM1\tF\t0\t0\t1\t1\t1\t52\t104\t12\n"
    "FAM1\tM\t0\t0\t2\t1\t1\t63\t96\t20\n"
    "FAM1\tC1\tF\tM\t1\t2\t2\t81\t64\t38\n"
    "FAM1\tC2\tF\tM\t2\t1\t1\t70\t80\t25\n"
)
NAMES = ("SRS T-score", "IQ", "ADHD-RS")
STYLE = dict(
    gradient_bins={"SRS": [60, 66, 76], "IQ": [70, 85]},
    gradient_bin_labels={
        "SRS": [None, "mild", "moderate", "severe"],
        "IQ": ["ID", "borderline", None],
    },
    gradient_diverging={"ADHD_RS": {"min": 0, "mid": 27, "max": 54}},
    phenotype_labels={"SRS": "SRS T-score", "ADHD_RS": "ADHD-RS"},
    # Bars at their own widths, so they end at different places (#95 makes
    # them equal by default).
    legend_bar_width="own",
)
SIZES = [{}, {"legend": "minimal", "font_size": 6.0, "min_font_size": 6.0}]


def _layout(**extra):
    pedigree = parse_pedigree_from_string(PED)
    names = pedigree.phenotype_names
    st = resolve_style(PedigreeStyle(**{**STYLE, **extra}), names)
    render = resolve_phenotype_render(pedigree, names, pedigree.phenotype_display_modes, st)
    render = with_quadrants(render, render.band_names(names), CarrierScheme(), st)
    rows = legend_rows(names, st, render, pedigree.phenotype_display_modes, CarrierScheme())
    return legend_layout(rows, st)


@pytest.mark.parametrize("extra", SIZES)
class TestOneColumn:
    def test_column_past_the_widest_bar(self, extra):
        rows, m = _layout(**extra)
        binned = [r for r in rows if r.kind == "binned"]
        assert len(binned) == 2
        ends = [binned_row_geometry(r, m).label_x for r in binned]
        assert ends[0] != pytest.approx(ends[1])  # bars of different lengths
        assert m.name_x == pytest.approx(max([m.text_x, *ends]))
        assert m.name_x > m.text_x

    def test_png_names_start_together(self, extra):
        fig = render_pedigree(
            parse_pedigree_from_string(PED), style=PedigreeStyle(**STYLE, **extra)
        )
        try:
            xs = {
                t.get_text(): t.get_position()[0]
                for t in fig.axes[1].texts
                if t.get_text() in NAMES and t.get_ha() == "left"
            }
        finally:
            plt.close(fig)
        assert set(xs) == set(NAMES)
        assert max(xs.values()) - min(xs.values()) == pytest.approx(0.0)

    def test_svg_names_start_together(self, extra):
        svg = render_interactive_svg(
            parse_pedigree_from_string(PED), style=PedigreeStyle(**STYLE, **extra)
        )
        legend = svg[svg.index('<g class="legend">') :]
        xs = {}
        for name in NAMES:
            found = re.findall(rf'<text x="([-0-9.]+)"[^>]*>{name}</text>', legend)
            xs[name] = float(found[-1])  # the row's name, after the quadrant key's
        assert max(xs.values()) - min(xs.values()) == pytest.approx(0.0, abs=1e-3)

    def test_widths_count_the_column(self, extra):
        rows, m = _layout(**extra)
        gradient = next(r for r in rows if r.label == "ADHD-RS")
        assert legend_row_right(gradient, m) > m.name_x


def test_unchanged_without_a_wide_binned_row():
    """No named steps: every bar fits before the usual column."""
    rows, m = _layout(gradient_bin_labels={})
    assert m.name_x == pytest.approx(m.text_x)
