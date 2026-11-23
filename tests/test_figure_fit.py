"""A figure as large as its pedigree and legend, with the same small margin all round.

The pedigree's axes keep an equal aspect inside a figure of at least
`fig_width` x `fig_height`, and the view kept 15% of the pedigree's size on each
side: a wide pedigree left blank bands above and below it (the legend went under
the lower one), a large one centimetres of white beside it, while labels or
annotations at the bottom could reach the edge, as text is not in matplotlib's
data limits. Now (`renderer._fit_view_to_content`, `_fit_figure`,
`_fit_legend_width`) each pedigree's axes hug what is drawn plus
`_VIEW_PAD_PT`, at the same scale, families side by side are `_FAMILY_GAP_PT`
apart, the legend's axes hug its frame, and the figure hugs them all.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pytest
from conftest import TEST_DATA_DIR
from matplotlib.transforms import Bbox

from pypedigree import PedigreeStyle, parse_pedigree, parse_pedigree_from_string, render_pedigree
from pypedigree import renderer as renderer_module

QUADRANTS = TEST_DATA_DIR / "phenotype" / "quadrants.ped"
PGS = TEST_DATA_DIR / "phenotype" / "multi_phenotype_pgs.ped"
ANNOTATIONS = TEST_DATA_DIR / "phenotype" / "annotations.ped"
LARGE = TEST_DATA_DIR / "complex" / "multigenerational.ped"
TRIO = TEST_DATA_DIR / "basic" / "simple_trio.ped"
MULTI = TEST_DATA_DIR / "basic" / "multi_family.ped"
# A narrow pedigree under a wider legend: four generations, one child each.
NARROW_WITH_LEGEND = (
    "FID\tIID\tPAT\tMAT\tSEX\tPHENO\tb_autism_spectrum_disorder\tb_ADHD\n"
    "F\tA\t0\t0\t1\t1\t1\t2\nF\tB\t0\t0\t2\t1\t1\t1\n"
    "F\tC\tA\tB\t2\t1\t1\t1\nF\tD\t0\t0\t1\t1\t1\t1\n"
    "F\tE\tD\tC\t1\t1\t1\t2\nF\tG\t0\t0\t2\t1\t1\t1\n"
    "F\tH\tE\tG\t2\t2\t2\t1"
)
FIXTURES = [QUADRANTS, PGS, ANNOTATIONS, LARGE, TRIO]
PT = 1 / 72  # inches


@pytest.fixture
def drawn():
    figures = []

    def draw(ped, **kwargs):
        if not hasattr(ped, "families"):
            ped = parse_pedigree(ped)
        fig = render_pedigree(ped, **kwargs)
        fig.canvas.draw()
        figures.append(fig)
        return fig

    yield draw
    for fig in figures:
        plt.close(fig)


def _drawn_extent(ax) -> Bbox:
    """Everything drawn in ``ax``, each artist measured on its own (inches)."""
    renderer = ax.figure.canvas.get_renderer()
    boxes = []
    for artist in (*ax.patches, *ax.lines, *ax.collections, *ax.texts, *ax.artists):
        if artist.get_visible():
            box = artist.get_window_extent(renderer)
            if np.isfinite(box.get_points()).all() and box.width + box.height > 0:
                boxes.append(box)
    return Bbox.union(boxes).transformed(ax.figure.dpi_scale_trans.inverted())


def _inches(fig, ax) -> Bbox:
    return ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())


@pytest.mark.parametrize("path", FIXTURES, ids=lambda p: p.stem)
@pytest.mark.parametrize("circular", [False, True], ids=["hierarchical", "circular"])
def test_same_white_space_on_every_side(drawn, path, circular):
    """`_VIEW_PAD_PT` between what is drawn and the axes' edge, on every side,
    labels and annotations included."""
    fig = drawn(path, circular=circular)
    ax = fig.axes[0]
    content, box = _drawn_extent(ax), _inches(fig, ax)
    pad = renderer_module._VIEW_PAD_PT * PT
    for gap in (content.x0 - box.x0, box.x1 - content.x1, content.y0 - box.y0, box.y1 - content.y1):
        assert gap == pytest.approx(pad, abs=1.5 * PT)


@pytest.mark.parametrize("path", [QUADRANTS, PGS, LARGE], ids=lambda p: p.stem)
@pytest.mark.parametrize("circular", [False, True], ids=["hierarchical", "circular"])
@pytest.mark.parametrize("legend", [True, False], ids=["legend", "no-legend"])
def test_no_blank_band_above_below_or_beside(drawn, path, circular, legend):
    """Only `_FIGURE_PAD_PT` between the figure's edges and its axes, and
    nothing between the pedigree and the legend."""
    fig = drawn(path, circular=circular, show_legend=legend)
    width, height = fig.get_size_inches()
    pedigree = _inches(fig, fig.axes[0])
    edge = renderer_module._FIGURE_PAD_PT * PT
    assert height - pedigree.y1 == pytest.approx(edge, abs=0.5 * PT)
    used = Bbox.union([_inches(fig, ax) for ax in fig.axes])
    assert used.x0 == pytest.approx(edge, abs=0.5 * PT)
    assert width - used.x1 == pytest.approx(edge, abs=0.5 * PT)
    below = _inches(fig, fig.axes[1]).y1 if legend else 0.0
    assert pedigree.y0 - below == pytest.approx(edge, abs=0.5 * PT)


@pytest.mark.parametrize("path", FIXTURES, ids=lambda p: p.stem)
@pytest.mark.parametrize("circular", [False, True], ids=["hierarchical", "circular"])
def test_drawing_keeps_its_scale(drawn, path, circular, monkeypatch):
    """Only white space goes: one layout unit is as many inches as before."""

    def inches_per_unit(fig):
        ax = fig.axes[0]
        (x0, y0), (x1, y1) = ax.transData.transform([(0, 0), (1, 1)])
        return (x1 - x0) / fig.dpi, (y1 - y0) / fig.dpi

    fitted = inches_per_unit(drawn(path, circular=circular))
    for name in ("_fit_view_to_content", "_fit_figure", "_fit_legend_width"):
        monkeypatch.setattr(renderer_module, name, lambda *args: None)
    unfitted = inches_per_unit(drawn(path, circular=circular))
    assert fitted == pytest.approx(unfitted)
    assert fitted[0] == pytest.approx(fitted[1])  # equal aspect


def test_small_pedigree_smaller_than_the_default_figure(drawn):
    style = PedigreeStyle()
    width, height = drawn(TRIO, style=style).get_size_inches()
    assert width < style.fig_width
    assert height < style.fig_height


def test_legend_wider_than_the_pedigree(drawn):
    """The figure is as wide as the legend's frame, the pedigree in the middle."""
    fig = drawn(parse_pedigree_from_string(NARROW_WITH_LEGEND))
    width = fig.get_size_inches()[0]
    pedigree, legend = fig.axes
    frame = _drawn_extent(legend)
    edge = renderer_module._FIGURE_PAD_PT * PT
    assert frame.x0 >= edge - PT
    assert frame.x1 <= width - edge + PT
    box = _inches(fig, pedigree)
    assert box.width < frame.width
    assert box.x0 == pytest.approx(width - box.x1, abs=PT)


def test_families_side_by_side(drawn):
    """Each family `_FAMILY_GAP_PT` after the previous one, titles inside the figure."""
    fig = drawn(MULTI)
    renderer = fig.canvas.get_renderer()
    boxes = sorted(
        (
            Bbox.union([ax.get_window_extent(), ax.title.get_window_extent(renderer)])
            for ax in fig.axes
        ),
        key=lambda b: b.x0,
    )
    for left, right in zip(boxes, boxes[1:], strict=False):
        assert (right.x0 - left.x1) / fig.dpi == pytest.approx(
            renderer_module._FAMILY_GAP_PT * PT, abs=PT
        )
    assert all(b.y1 <= fig.bbox.y1 for b in boxes)


def test_callers_axes_left_alone():
    """Drawn into a caller's axes: the caller's figure and axes keep their size."""
    fig, ax = plt.subplots(figsize=(6, 9))
    try:
        position = ax.get_position().bounds
        render_pedigree(parse_pedigree(QUADRANTS), ax=ax)
        assert tuple(fig.get_size_inches()) == pytest.approx((6, 9))
        assert ax.get_position(original=True).bounds == pytest.approx(position)
    finally:
        plt.close(fig)


class TestFamiliesAtOneScale:
    """#123: families of different heights side by side are drawn at one scale."""

    TWO = (
        "FID\tIID\tPAT\tMAT\tSEX\tPHENO\n"
        "A\ta1\t0\t0\t1\t1\nA\ta2\t0\t0\t2\t2\nA\ta3\ta1\ta2\t2\t1\n"
        "A\ta4\t0\t0\t1\t1\nA\ta5\ta4\ta3\t1\t2\n"
        "B\tb1\t0\t0\t1\t1\nB\tb2\t0\t0\t2\t1\nB\tb3\tb1\tb2\t1\t2\n"
    )

    def test_same_points_per_unit_and_first_generation_level(self):
        import matplotlib.pyplot as plt

        from pypedigree import parse_pedigree_from_string, render_pedigree

        fig = render_pedigree(parse_pedigree_from_string(self.TWO))
        try:
            fig.canvas.draw()
            a, b = fig.axes[:2]
            scales = []
            for ax in (a, b):
                x0, x1 = ax.get_xlim()
                scales.append(ax.get_window_extent().width / fig.dpi * 72 / abs(x1 - x0))
            assert scales[0] == pytest.approx(scales[1], rel=1e-3)
            # The founders' row (top of each drawing) on one line.
            tops = [ax.transData.transform((0, ax.dataLim.y1))[1] for ax in (a, b)]
            assert tops[0] == pytest.approx(tops[1], abs=1.0)
        finally:
            plt.close(fig)
