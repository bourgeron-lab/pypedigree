"""No blank band above or below a wide pedigree, nor between it and its legend.

The pedigree's axes keep an equal aspect, so a pedigree wider than the figure's
proportions (at least `fig_width` x `fig_height`) was drawn shorter than its
slot, centred, and the legend went under the lower blank band: up to a third of
a three-generation figure was empty. `renderer._trim_figure_height` takes that
height off before the legend is added.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import pytest
from conftest import TEST_DATA_DIR

from pypedigree import PedigreeStyle, parse_pedigree, render_pedigree

QUADRANTS = TEST_DATA_DIR / "phenotype" / "quadrants.ped"
PGS = TEST_DATA_DIR / "phenotype" / "multi_phenotype_pgs.ped"
LARGE = TEST_DATA_DIR / "complex" / "multigenerational.ped"
TRIO = TEST_DATA_DIR / "basic" / "simple_trio.ped"
MULTI = TEST_DATA_DIR / "basic" / "multi_family.ped"

# tight_layout's pad between the pedigree's axes and the legend's, in inches.
PAD_IN = 0.2


@pytest.fixture
def drawn():
    figures = []

    def draw(path, **kwargs):
        fig = render_pedigree(parse_pedigree(path), **kwargs)
        fig.canvas.draw()
        figures.append(fig)
        return fig

    yield draw
    for fig in figures:
        plt.close(fig)


def _heights(fig, ax):
    """(slot, drawn) height of an axes, in inches."""
    height = fig.get_size_inches()[1]
    return ax.get_position(original=True).height * height, ax.get_position().height * height


@pytest.mark.parametrize("path", [QUADRANTS, PGS, LARGE], ids=lambda p: p.stem)
@pytest.mark.parametrize("circular", [False, True], ids=["hierarchical", "circular"])
@pytest.mark.parametrize("legend", [True, False], ids=["legend", "no-legend"])
def test_no_blank_band_above_or_below(drawn, path, circular, legend):
    """Between the figure's top and the pedigree, and between the pedigree and
    the legend (or the figure's bottom), only tight_layout's pad."""
    fig = drawn(path, circular=circular, show_legend=legend)
    height = fig.get_size_inches()[1]
    pedigree = fig.axes[0].get_position()
    below = fig.axes[1].get_position().y1 if len(fig.axes) > 1 else 0.0
    assert (1 - pedigree.y1) * height < PAD_IN
    assert 0 <= (pedigree.y0 - below) * height < PAD_IN


def test_wide_pedigree_shorter_than_fig_height(drawn):
    """A three-generation family is drawn at its own proportions: the figure
    is less tall than `fig_height`, the pedigree keeps the full width."""
    style = PedigreeStyle()
    fig = drawn(QUADRANTS, show_legend=False, style=style)
    width, height = fig.get_size_inches()
    assert height < style.fig_height
    ax = fig.axes[0]
    assert ax.get_position().width == pytest.approx(ax.get_position(original=True).width)


def test_tall_pedigree_keeps_fig_height(drawn):
    """A pedigree narrower than the figure already fills its height."""
    style = PedigreeStyle()
    fig = drawn(TRIO, style=style)
    assert fig.get_size_inches()[1] == pytest.approx(style.fig_height)


def test_several_families_keep_their_titles(drawn):
    """The tallest family fills its height; the titles stay inside the figure."""
    fig = drawn(MULTI)
    fills = [abs(slot - actual) < 0.01 for slot, actual in (_heights(fig, ax) for ax in fig.axes)]
    assert any(fills)
    renderer = fig.canvas.get_renderer()
    top = fig.bbox.y1
    for ax in fig.axes:
        assert ax.title.get_window_extent(renderer).y1 <= top


def test_callers_axes_left_alone():
    """Drawn into a caller's axes, the caller's figure keeps its size."""
    fig, ax = plt.subplots(figsize=(6, 9))
    try:
        render_pedigree(parse_pedigree(QUADRANTS), ax=ax)
        assert tuple(fig.get_size_inches()) == pytest.approx((6, 9))
    finally:
        plt.close(fig)
