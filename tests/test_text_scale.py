"""The PNG never draws a pedigree smaller than the layout's text scale.

The layout keeps room for labels and annotations measured at
`symbols.UNITS_PER_POINT` layout units per point. A wide pedigree used to be
drawn at about 44 points per unit (the figure's margins and its 80-inch cap),
so text took more room than was kept for it and long annotations overlapped
their neighbour's. `renderer._keep_text_scale` grows the figure to that scale,
as long as the drawing stays within `_MAX_FIG_WIDTH_IN` x `_MAX_FIG_HEIGHT_IN`.
"""

from __future__ import annotations

from collections import defaultdict

import matplotlib.pyplot as plt
import pytest
from matplotlib.transforms import Bbox

from pypedigree import parse_pedigree_from_string, render_pedigree, symbols
from pypedigree import renderer as renderer_module

TARGET = 1 / symbols.UNITS_PER_POINT
HEADER = "FID\tIID\tFather\tMother\tSex\tPhenotype\ta_note\n"


def _rows(n_children: int, family: str) -> str:
    """One couple and ``n_children`` children, each with an invented note."""
    rows = [f"{family}\t{family}f\t0\t0\t1\t1\t\n", f"{family}\t{family}m\t0\t0\t2\t1\t\n"]
    for i in range(n_children):
        note = f"invented annotation number {i:03d}"
        rows.append(f"{family}\t{family}c{i}\t{family}f\t{family}m\t{1 + i % 2}\t1\t{note}\n")
    return "".join(rows)


def _sibship(n_children: int, *more: int) -> str:
    """A file of one such family per count, F1, F2..."""
    counts = (n_children, *more)
    return HEADER + "".join(_rows(n, f"F{i + 1}") for i, n in enumerate(counts))


@pytest.fixture
def drawn():
    figures = []

    def draw(text: str, **kwargs):
        fig = render_pedigree(parse_pedigree_from_string(text), **kwargs)
        fig.canvas.draw()
        figures.append(fig)
        return fig

    yield draw
    for fig in figures:
        plt.close(fig)


def _overlaps(fig) -> list[tuple[str, str]]:
    """Pairs of texts of one row (same baseline and size) that overlap."""
    renderer = fig.canvas.get_renderer()
    rows = defaultdict(list)
    for ax in fig.axes:
        for t in ax.texts:
            if t.get_text().strip():
                box = t.get_window_extent(renderer)
                rows[(round(box.y0), t.get_fontsize())].append((box.x0, box.x1, t.get_text()))
    pairs = []
    for row in rows.values():
        row.sort()
        pairs += [(a[2], b[2]) for a, b in zip(row, row[1:], strict=False) if b[0] < a[1]]
    return pairs


def test_wide_pedigree_drawn_at_the_layout_scale(drawn):
    fig = drawn(_sibship(60))
    assert renderer_module._points_per_unit(fig.axes[0]) >= TARGET - 1e-6
    assert _overlaps(fig) == []


def test_text_overlapped_without_it(drawn, monkeypatch):
    """The case above is one where the old scale made annotations overlap."""
    monkeypatch.setattr(renderer_module, "_keep_text_scale", lambda *args: None)
    fig = drawn(_sibship(60))
    assert renderer_module._points_per_unit(fig.axes[0]) < TARGET
    assert _overlaps(fig)


def test_families_side_by_side(drawn):
    """A row of families is squeezed into `fig_width` per family: grown too."""
    fig = drawn(_sibship(20, 20))
    for ax in fig.axes:
        assert renderer_module._points_per_unit(ax) >= TARGET - 1e-6
    assert _overlaps(fig) == []


def test_small_pedigree_keeps_its_larger_scale(drawn, monkeypatch):
    """A pedigree already drawn larger than the layout's scale is not touched."""
    fig = drawn(_sibship(2))
    grown = renderer_module._points_per_unit(fig.axes[0])
    monkeypatch.setattr(renderer_module, "_keep_text_scale", lambda *args: None)
    assert grown > TARGET
    assert renderer_module._points_per_unit(drawn(_sibship(2)).axes[0]) == pytest.approx(grown)


def test_drawing_never_wider_than_the_largest_figure(drawn):
    """Past the largest figure the pedigree is drawn smaller, as before."""
    fig = drawn(_sibship(150), show_labels=False)
    ax = fig.axes[0]
    assert renderer_module._points_per_unit(ax) < TARGET
    width = Bbox.union([ax.dataLim.transformed(ax.transData)]).width / fig.dpi
    # tight_layout, run again on the grown figure, gives the axes a slightly
    # larger share of it.
    assert width <= renderer_module._MAX_FIG_WIDTH_IN * 1.001


def test_circular_layout_left_alone(drawn, monkeypatch):
    """The circular layout keeps no label room, so it keeps its own scale."""
    grown = renderer_module._points_per_unit(drawn(_sibship(60), circular=True).axes[0])
    monkeypatch.setattr(renderer_module, "_keep_text_scale", lambda *args: None)
    plain = renderer_module._points_per_unit(drawn(_sibship(60), circular=True).axes[0])
    assert grown == pytest.approx(plain)
