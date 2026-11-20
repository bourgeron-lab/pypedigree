"""One width for all the bars of a legend (#95).

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
    LEG_GRAD_W,
    LEG_SWATCH_EDGE,
    LEG_X0,
    binned_row_geometry,
    check_legend,
    legend_layout,
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
BASE = dict(
    gradient_bins={"SRS": [60, 66, 76], "IQ": [70, 85]},
    gradient_bin_labels={
        "SRS": [None, "mild", "moderate", "severe"],
        "IQ": ["ID", "borderline", None],
    },
    gradient_diverging={"ADHD_RS": {"min": 0, "mid": 27, "max": 54}},
)
SIZES = [{}, {"legend": "minimal", "font_size": 6.0, "min_font_size": 6.0}]


def _layout(**extra):
    pedigree = parse_pedigree_from_string(PED)
    names = pedigree.phenotype_names
    st = resolve_style(PedigreeStyle(**{**BASE, **extra}), names)
    render = resolve_phenotype_render(pedigree, names, pedigree.phenotype_display_modes, st)
    render = with_quadrants(render, render.band_names(names), CarrierScheme(), st)
    rows = legend_rows(names, st, render, pedigree.phenotype_display_modes, CarrierScheme())
    return legend_layout(rows, st)


def _png_bars(**extra):
    """Each bar's outline (the edged, unfilled rectangle), top to bottom, in
    legend units; and the legend's width in pixels."""
    fig = render_pedigree(parse_pedigree_from_string(PED), style=PedigreeStyle(**{**BASE, **extra}))
    try:
        fig.canvas.draw()
        ax = fig.axes[1]
        bars = sorted(
            (
                (p.get_x(), p.get_x() + p.get_width(), p.get_y())
                for p in ax.patches
                if p.get_facecolor()[3] == 0 and p.get_width() > 1.0
            ),
            key=lambda b: -b[2],
        )
        frame = max(ax.patches, key=lambda p: p.get_window_extent().width)
        return bars, frame.get_window_extent().width, ax
    finally:
        plt.close(fig)


@pytest.mark.parametrize("extra", SIZES)
class TestEqual:
    def test_png_bars_start_and_end_together(self, extra):
        bars, _width, _ax = _png_bars(**extra)
        assert len(bars) == 3
        assert {round(b[0], 6) for b in bars} == {round(LEG_X0, 6)}
        assert max(b[1] for b in bars) - min(b[1] for b in bars) == pytest.approx(0.0, abs=1e-9)

    def test_svg_bars_start_and_end_together(self, extra):
        svg = render_interactive_svg(
            parse_pedigree_from_string(PED), style=PedigreeStyle(**{**BASE, **extra})
        )
        legend = svg[svg.index('<g class="legend">') :]
        legend = legend[: legend.index("\n</g>")]
        bars = re.findall(
            rf'<rect x="([-0-9.]+)" y="[-0-9.]+" width="([0-9.]+)" height="[0-9.]+" '
            rf'fill="none" stroke="{LEG_SWATCH_EDGE}"',
            legend,
        )
        bars = [(float(x), float(w)) for x, w in bars if float(w) > 0.3]
        assert len(bars) == 3
        assert len({round(x, 3) for x, _ in bars}) == 1
        assert len({round(w, 3) for _, w in bars}) == 1

    def test_steps_equal_and_never_narrower(self, extra):
        rows, m = _layout(**extra)
        _, own = _layout(legend_bar_width="own", **extra)
        assert m.bar_w > LEG_GRAD_W
        for row in (r for r in rows if r.kind == "binned"):
            g = binned_row_geometry(row, m)
            assert g.step_w * len(row.colors) == pytest.approx(m.bar_w)
            assert g.step_w >= binned_row_geometry(row, own).step_w - 1e-12

    def test_captions_within_their_steps(self, extra):
        rows, m = _layout(**extra)
        bars, _w, _ax = _png_bars(**extra)
        fig = render_pedigree(
            parse_pedigree_from_string(PED), style=PedigreeStyle(**{**BASE, **extra})
        )
        try:
            fig.canvas.draw()
            ax = fig.axes[1]
            inv = ax.transData.inverted()
            boxes = {}
            for t in ax.texts:
                e = t.get_window_extent()
                (x0, _), (x1, _) = inv.transform([(e.x0, e.y0), (e.x1, e.y1)])
                boxes[t.get_text()] = (x0, x1)
        finally:
            plt.close(fig)
        for row in (r for r in rows if r.kind == "binned"):
            g = binned_row_geometry(row, m)
            for i, (limits, name) in enumerate(zip(row.bins, row.bin_names, strict=True)):
                lo, hi = LEG_X0 + i * g.step_w, LEG_X0 + (i + 1) * g.step_w
                for caption in (limits, name):
                    if caption:
                        x0, x1 = boxes[caption]
                        assert lo - 1e-6 <= x0 and x1 <= hi + 1e-6, caption
        # The gradient's end captions move with its ends.
        x0, _ = boxes["0"]
        _, x1 = boxes["54"]
        assert x0 == pytest.approx(LEG_X0, abs=0.05)
        assert x1 == pytest.approx(LEG_X0 + m.bar_w, abs=0.05)

    def test_legend_keeps_its_width(self, extra):
        _, equal_w, _ = _png_bars(**extra)
        _, own_w, _ = _png_bars(legend_bar_width="own", **extra)
        assert equal_w == pytest.approx(own_w, abs=1.0)


def test_own_keeps_each_bar_its_width():
    bars, _w, _ax = _png_bars(legend_bar_width="own")
    ends = [b[1] for b in bars]
    assert len(set(round(e, 6) for e in ends)) == 3
    assert min(ends) == pytest.approx(LEG_X0 + LEG_GRAD_W)  # the gradient


def test_a_single_bar_keeps_its_width():
    """A legend with one bar: nothing to match, the bar keeps its width."""
    style_one = dict(gradient_bins={"SRS": [60, 66, 76]}, gradient_bin_labels={})
    pedigree = parse_pedigree_from_string(
        "FID\tIID\tFather\tMother\tSex\tPhenotype\tg_SRS\nFAM1\ta\t0\t0\t1\t2\t70\n"
    )
    names = pedigree.phenotype_names
    st = resolve_style(PedigreeStyle(**style_one), names)
    render = resolve_phenotype_render(pedigree, names, pedigree.phenotype_display_modes, st)
    rows = legend_rows(names, st, render, pedigree.phenotype_display_modes, CarrierScheme())
    _, m = legend_layout(rows, st)
    assert m.bar_w == 0.0


def test_invalid():
    with pytest.raises(PedigreeError, match="legend_bar_width"):
        check_legend(PedigreeStyle(legend_bar_width="wide"))
