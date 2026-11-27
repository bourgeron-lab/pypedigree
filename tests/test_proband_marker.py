"""The proband's "P" stays clear of the label on its left.

The marker sits at the lower left of the symbol, level with the labels. With
tight spacing or large text, 0.4.4 drew the "P" over the label of the
individual on the left; the layout now keeps that neighbour far enough
(`symbols.marker_clearance`). A left neighbour's evaluation asterisk keeps
clear of the arrow too (#143). At the default style nothing moves.

All pedigrees here are synthetic.
"""

from __future__ import annotations

import re

import matplotlib.pyplot as plt
import pytest

from pypedigree import PedigreeStyle, compute_layout, parse_pedigree_from_string, render_pedigree
from pypedigree.svg_renderer import render_interactive_svg
from pypedigree.symbols import (
    MARKER_ARROW_EM,
    MARKER_LETTER_GAP_EM,
    UNITS_PER_POINT,
    descent_clearance,
    marker_clearance,
    marker_points,
    marker_reach,
)

TIGHT = PedigreeStyle(sibling_gap=1.0, couple_gap=1.0, font_size=18)
HEADER = "FID\tID\tFather\tMother\tSex\tProband\tConsultand\n"
TRIO = HEADER + "F\t1\t0\t0\t1\t0\t0\nF\t2\t0\t0\t2\t0\t0\nF\t3\t1\t2\t1\t1\t0\n"


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

    @pytest.mark.parametrize("proband", ["10", "12", "13", "27"])
    def test_default_style_moves_nothing(self, proband):
        """Except for a proband on the right of a couple with children (below)."""
        style = PedigreeStyle()
        assert _need(style) < style.sibling_gap
        assert _xs(_sibship(proband), style) == _xs(_sibship(""), style)


class TestDescentLine:
    """A proband or consultand on the right of a couple with children: their
    marker, at the lower left of the symbol, keeps clear of the line of
    descent hanging from the middle of the gap (the "P" was drawn on it)."""

    def _gap(self, style, proband="11", consultand=""):
        x = _xs(_sibship(proband, consultand), style)
        return x["11"] - x["10"]

    def test_marker_left_of_the_descent_line(self):
        style = PedigreeStyle()
        gap = self._gap(style)
        assert gap / 2 >= marker_reach(0.4, style.font_size, True) + 1e-9
        assert gap == pytest.approx(descent_clearance(0.4, style.font_size, True))
        assert gap > style.couple_gap

    def test_consultand_arrow_fits_the_default_gap(self):
        style = PedigreeStyle()
        assert descent_clearance(0.4, style.font_size, False) < style.couple_gap
        assert self._gap(style, proband="", consultand="11") == pytest.approx(style.couple_gap)

    def test_left_partner_unchanged(self):
        """A marker on the left partner points away from the line."""
        style = PedigreeStyle()
        assert self._gap(style, proband="10") == pytest.approx(style.couple_gap)

    def test_larger_text_wider_gap(self):
        assert self._gap(PedigreeStyle(font_size=14)) > self._gap(PedigreeStyle())

    def test_svg_p_right_of_the_descent_line(self):
        """In the SVG, the whole "P" (right-aligned at its x, about 0.8 em
        wide) stays right of the couple's line of descent."""
        svg = render_interactive_svg(parse_pedigree_from_string(_sibship("11")))
        family = parse_pedigree_from_string(_sibship("11")).families["F"]
        layout = compute_layout(family, PedigreeStyle())
        x = {n.individual_id: n.x for n in layout.nodes}
        middle = (x["10"] + x["11"]) / 2
        letters = [
            float(m.group(1))
            for m in re.finditer(r'<text[^>]*x="([-0-9.]+)"[^>]*text-anchor="end"[^>]*>P<', svg)
        ]
        assert letters, "the proband's P is drawn right-aligned"
        em = PedigreeStyle().font_size * UNITS_PER_POINT
        assert all(lx - 0.8 * em > middle for lx in letters)

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
        """Same offsets from the tip, in points in the PNG and in layout units
        (`UNITS_PER_POINT` per point) in the SVG, whose y points down."""
        sx, sy = marker_points(1.0, 2.0, 8.0)
        sx2, sy2 = marker_points(1.0, -2.0, 8.0, up=-1.0)
        assert (sx, sy) == (sx2, -sy2) == (1.0 - 8 * MARKER_ARROW_EM, 2.0 - 8 * MARKER_ARROW_EM)

        fig = render_pedigree(parse_pedigree_from_string(_sibship()), style=TIGHT)
        try:
            (p,) = [t for t in fig.axes[0].texts if t.get_text() == "P"]
            (arrow,) = [t for t in fig.axes[0].texts if t.get_text() == "" and t.arrow_patch]
            em = TIGHT.font_size
            start = (-MARKER_ARROW_EM * em, -MARKER_ARROW_EM * em)
            assert arrow.xyann == pytest.approx(start)
            assert p.xyann == pytest.approx((start[0] - MARKER_LETTER_GAP_EM * em, start[1]))
            assert (p.get_ha(), p.get_va()) == ("right", "center")
        finally:
            plt.close(fig)

        svg = render_interactive_svg(parse_pedigree_from_string(_sibship()), style=TIGHT)
        (px, py) = re.findall(
            r'<text x="([-\d.]+)" y="([-\d.]+)"[^>]*text-anchor="end"[^>]*>P</text>', svg
        )[0]
        (line,) = re.findall(
            r'<line x1="([-\d.]+)" y1="([-\d.]+)" x2="([-\d.]+)" y2="([-\d.]+)"[^>]*'
            r'marker-end="url\(#arrowhead\)"',
            svg,
        )
        x1, y1, x2, y2 = map(float, line)
        size = TIGHT.font_size * UNITS_PER_POINT
        assert (x2 - x1, y1 - y2) == pytest.approx((MARKER_ARROW_EM * size,) * 2, abs=2e-3)
        assert float(px) == pytest.approx(x1 - MARKER_LETTER_GAP_EM * size, abs=2e-3)

    @pytest.mark.parametrize("figsize", [(10, 8), (4, 3), (2.5, 2), (1.2, 1.0)])
    def test_p_never_covers_its_arrow(self, figsize):
        """0.4.6 put the "P" a fixed distance in layout units from the arrow: with
        text large against the symbols (a pedigree drawn small, as in a
        journal figure), the letter covered the arrow."""
        fig, ax = plt.subplots(figsize=figsize)
        try:
            render_pedigree(parse_pedigree_from_string(TRIO), style=PedigreeStyle(), ax=ax)
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            (p,) = [t for t in ax.texts if t.get_text() == "P"]
            (arrow,) = [t for t in ax.texts if t.get_text() == "" and t.arrow_patch]
            letter = p.get_window_extent(renderer)
            line = arrow.arrow_patch.get_window_extent(renderer)
            assert letter.x1 < line.x0
        finally:
            plt.close(fig)


class TestPgsCaption:
    """The "PGS" caption under a curve met the proband marker (0.4.6)."""

    PED = (
        "FID\tID\tFather\tMother\tSex\tProband\tConsultand\tPGS\n"
        "F\t1\t0\t0\t1\t0\t0\t0.1\n"
        "F\t2\t0\t0\t2\t0\t1\t-0.3\n"
        "F\t3\t1\t2\t1\t1\t0\t0.8\n"
        "F\t4\t1\t2\t2\t0\t0\t-0.5\n"
    )

    def _nodes(self):
        """Laid out with the spacing the renderers apply for PGS curves."""
        from pypedigree.config import apply_pgs_spacing

        family = parse_pedigree_from_string(self.PED).families["F"]
        style = apply_pgs_spacing(PedigreeStyle(), ["PGS"])
        return {n.individual_id: n for n in compute_layout(family, style).nodes}

    @staticmethod
    def _owner(nodes, x, y):
        """The individual whose curve the caption at (x, y) belongs to: the
        nearest symbol in the same generation on the side of its box (#78:
        a right-hand partner's box is on their right)."""

        def side_ok(n):
            dx = n.x - x
            return -1.0 < dx < 0 if n.pgs_side == "right" else 0 < dx < 1.0

        (sid,) = [s for s, n in nodes.items() if side_ok(n) and abs(n.y - y) < 1.0]
        return sid

    @staticmethod
    def _above(nodes, sid):
        """A marked individual's caption goes above a box on the left, where
        their marker is; the consultand "2" is a right-hand partner, whose box
        and caption are on the right, clear of it (#78)."""
        return sid in ("2", "3") and nodes[sid].pgs_side == "left"

    def test_marked_individuals_on_both_sides(self):
        nodes = self._nodes()
        assert nodes["2"].pgs_side == "right" and nodes["3"].pgs_side == "left"

    def test_png_caption_above_the_curve_of_marked_individuals(self):
        nodes = self._nodes()
        fig = render_pedigree(parse_pedigree_from_string(self.PED), show_legend=False)
        try:
            captions = [t for t in fig.axes[0].texts if t.get_text() == "PGS"]
            assert len(captions) == 4
            for caption in captions:
                x, y = caption.get_position()
                sid = self._owner(nodes, x, y)
                marked = self._above(nodes, sid)
                assert (y > nodes[sid].y) is marked, sid
                assert caption.get_va() == ("bottom" if marked else "top")
        finally:
            plt.close(fig)

    def test_svg_caption_above_the_curve_of_marked_individuals(self):
        nodes = self._nodes()
        svg = render_interactive_svg(parse_pedigree_from_string(self.PED), show_legend=False)
        captions = re.findall(r'<text x="([-\d.]+)" y="([-\d.]+)"[^>]*>PGS</text>', svg)
        assert len(captions) == 4
        for x, y in captions:
            # The SVG's y is the layout's, negated (pointing down): above the
            # node means a smaller y.
            sid = self._owner(nodes, float(x), -float(y))
            assert (float(y) < -nodes[sid].y) is self._above(nodes, sid), sid

    def test_circular_layout_still_draws_the_marker(self):
        fig = render_pedigree(parse_pedigree_from_string(_sibship()), style=TIGHT, circular=True)
        try:
            assert [t.get_text() for t in fig.axes[0].texts].count("P") == 1
        finally:
            plt.close(fig)


class TestNeighbourAsterisk:
    """The arrow keeps clear of a left neighbour's evaluation asterisk (#143):
    with tight spacing or large text, the asterisk ran into the arrow's tail."""

    HEAD = "FID\tID\tFather\tMother\tSex\tConsultand\tProband\tevaluated\n"

    def _ped(self, sexes=("1", "1"), proband: bool = False) -> str:
        rows = ["F\tf\t0\t0\t1\t0\t0\t0", "F\tm\t0\t0\t2\t0\t0\t0"]
        rows.append(f"F\ta\tf\tm\t{sexes[0]}\t0\t0\t1")  # evaluated, on the left
        marked = "0\t1" if proband else "1\t0"
        rows.append(f"F\tb\tf\tm\t{sexes[1]}\t{marked}\t0")
        return self.HEAD + "\n".join(rows) + "\n"

    @staticmethod
    def _mark_reach(style: PedigreeStyle, sex) -> float:
        from pypedigree import symbols

        size = style.text_size(symbols.EVALUATION_MARK_FONT_SCALE)
        reach = symbols.evaluation_mark_reach(
            style.shape_size / 2, sex, size, style.shape_line_width
        )
        return reach + symbols.EVALUATION_MARK_GAP_EM * size * UNITS_PER_POINT

    def test_clearance_counts_the_asterisk(self):
        half, size = 0.4, 18
        reach = 0.5
        need = marker_clearance(half, size, 2, False, neighbour_mark=reach)
        assert need >= reach + half + size * UNITS_PER_POINT * MARKER_ARROW_EM - 1e-9
        assert marker_clearance(half, size, 2, False) < need

    @pytest.mark.parametrize("sexes", [("1", "1"), ("2", "2"), ("0", "1")])
    @pytest.mark.parametrize("proband", [False, True])
    def test_layout_keeps_the_arrow_off_the_asterisk(self, sexes, proband):
        from pypedigree.models import Sex

        style = PedigreeStyle(sibling_gap=1.0, couple_gap=1.0, font_size=14)
        x = _xs(self._ped(sexes, proband), style)
        sex = {"1": Sex.MALE, "2": Sex.FEMALE, "0": Sex.UNKNOWN}[sexes[0]]
        tail = style.shape_size / 2 + style.font_size * UNITS_PER_POINT * MARKER_ARROW_EM
        assert x["b"] - x["a"] >= self._mark_reach(style, sex) + tail - 1e-9

    def test_default_style_moves_nothing(self):
        style = PedigreeStyle()
        with_mark = _xs(self._ped(), style)
        without = _xs(self._ped().replace("\t0\t0\t1\n", "\t0\t0\t0\n"), style)
        assert with_mark == without

    def test_png_asterisk_clear_of_the_arrow(self):
        """In the PNG at the text scale (render_pedigree's own figure), the
        asterisk's box ends left of the arrow's tail."""
        style = PedigreeStyle(sibling_gap=1.0, couple_gap=1.0, font_size=14)
        fig = render_pedigree(parse_pedigree_from_string(self._ped()), style=style)
        try:
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            ax = fig.axes[0]
            (star,) = [t for t in ax.texts if t.get_text() == "*"]
            (arrow,) = [t for t in ax.texts if t.get_text() == "" and t.arrow_patch is not None]
            star_box = star.get_window_extent(renderer)
            tail_x = arrow.arrow_patch.get_window_extent(renderer).x0
            assert star_box.x1 < tail_x
        finally:
            plt.close(fig)
