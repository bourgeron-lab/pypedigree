"""The proband's "P" stays clear of the label on its left.

The marker sits at the lower left of the symbol, level with the labels. With
tight spacing or large text, 0.4.4 drew the "P" over the label of the
individual on the left; the layout now keeps that neighbour far enough
(`symbols.marker_clearance`). At the default style nothing moves.

All pedigrees here are synthetic.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import pytest

from pypedigree import PedigreeStyle, compute_layout, parse_pedigree_from_string, render_pedigree
from pypedigree.svg_renderer import render_interactive_svg
from pypedigree.symbols import marker_clearance, marker_points

TIGHT = PedigreeStyle(sibling_gap=1.0, couple_gap=1.0, font_size=18)
HEADER = "FID\tID\tFather\tMother\tSex\tProband\tConsultand\n"


def _sibship(proband: str = "13", consultand: str = "") -> str:
    """Parents 10 and 11 with sixteen children, 12 to 27.

    Wide enough for the PNG to be drawn at the scale of a large pedigree
    (`renderer._auto_figsize`), where text takes the most room.
    """
    rows = [("10", "0", "0", "1"), ("11", "0", "0", "2")]
    rows += [(str(c), "10", "11", str(1 + c % 2)) for c in range(12, 28)]
    lines = [
        f"F\t{sid}\t{pat}\t{mat}\t{sex}\t{int(sid == proband)}\t{int(sid == consultand)}"
        for sid, pat, mat, sex in rows
    ]
    return HEADER + "\n".join(lines) + "\n"


def _xs(text: str, style: PedigreeStyle) -> dict[str, float]:
    family = parse_pedigree_from_string(text).families["F"]
    return {n.individual_id: n.x for n in compute_layout(family, style).nodes}


def _need(style: PedigreeStyle, letter: bool = True) -> float:
    return marker_clearance(style.shape_size / 2, style.font_size, 2, letter)


class TestLayout:
    def test_left_sibling_keeps_its_distance(self):
        x = _xs(_sibship(), TIGHT)
        assert _need(TIGHT) > TIGHT.sibling_gap
        assert x["13"] - x["12"] == pytest.approx(_need(TIGHT))
        # The other siblings keep the sibling gap.
        assert x["14"] - x["13"] == pytest.approx(TIGHT.sibling_gap)
        assert x["15"] - x["14"] == pytest.approx(TIGHT.sibling_gap)

    def test_couple_keeps_room_for_a_proband_on_the_right(self):
        x = _xs(_sibship(proband="11"), TIGHT)
        assert x["11"] - x["10"] >= _need(TIGHT) - 1e-9

    def test_consultand_arrow_needs_less_room(self):
        x = _xs(_sibship(proband="", consultand="13"), TIGHT)
        assert _need(TIGHT, letter=False) < _need(TIGHT)
        assert x["13"] - x["12"] == pytest.approx(
            max(TIGHT.sibling_gap, _need(TIGHT, letter=False))
        )

    @pytest.mark.parametrize("proband", ["10", "11", "12", "13", "27"])
    def test_default_style_moves_nothing(self, proband):
        style = PedigreeStyle()
        assert _need(style) < style.sibling_gap
        assert _xs(_sibship(proband), style) == _xs(_sibship(""), style)

    def test_leftmost_proband_adds_no_room(self):
        x = _xs(_sibship(proband="12"), TIGHT)
        assert x == _xs(_sibship(""), TIGHT)

    def test_longer_labels_need_more_room(self):
        short = marker_clearance(0.4, 18, 2, True)
        assert marker_clearance(0.4, 18, 6, True) > short
        assert marker_clearance(0.4, 24, 2, True) > short


class TestDrawing:
    def test_p_overlaps_no_label_in_the_png(self):
        fig = render_pedigree(parse_pedigree_from_string(_sibship()), style=TIGHT)
        try:
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            texts = [t for t in fig.axes[0].texts if t.get_text()]
            (p,) = [t for t in texts if t.get_text() == "P"]
            box = p.get_window_extent(renderer)
            others = [t.get_window_extent(renderer) for t in texts if t is not p]
            assert len(others) == 18
            assert not any(box.overlaps(o) for o in others)
        finally:
            plt.close(fig)

    def test_both_renderers_place_the_p_alike(self):
        """The SVG (y down) mirrors the PNG's marker through `marker_points`."""
        (sx, sy), (px, py) = marker_points(1.0, 2.0, 0.4)
        (sx2, sy2), (px2, py2) = marker_points(1.0, -2.0, 0.4, up=-1.0)
        assert (sx, px) == (sx2, px2)
        assert (sy, py) == (-sy2, -py2)
        svg = render_interactive_svg(parse_pedigree_from_string(_sibship()), style=TIGHT)
        assert svg.count(">P</text>") == 1

    def test_circular_layout_still_draws_the_marker(self):
        fig = render_pedigree(parse_pedigree_from_string(_sibship()), style=TIGHT, circular=True)
        try:
            assert [t.get_text() for t in fig.axes[0].texts].count("P") == 1
        finally:
            plt.close(fig)
