"""The separation marks of a broken relationship (#141).

- Two parallel slashes across the relationship line, which stops on either
  side of them, the same in the PNG and the SVG (the SVG drew a dashed line).
- In the room between the left partner and the line of descent when the
  couple has children, never on the line of descent; at the middle of the
  line when it has none.
- Sized with the symbols and shrunk to their room, so they fit at any scale,
  and stay two marks at 6 pt in a small figure (they merged into one on the
  line of descent).

All pedigrees here are synthetic.
"""

from __future__ import annotations

import math
import xml.etree.ElementTree as ET

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pypedigree import (
    PedigreeBuilder,
    PedigreeStyle,
    render_interactive_svg,
    render_pedigree,
    symbols,
)
from pypedigree.layout import compute_layout, compute_layout_circular
from pypedigree.models import EdgeType, LineStyle

SVG = "{http://www.w3.org/2000/svg}"
STYLE = PedigreeStyle()
HALF = STYLE.shape_size / 2
LINE = ((0.0, 0.0), (0.7, 0.0))  # a couple at the default gaps, outline to outline


def _couple(children: int = 1):
    b = PedigreeBuilder("F")
    b.add_individual("dad", sex="male", phenotype="unaffected")
    b.add_individual("mum", sex="female", phenotype="unaffected")
    for i in range(children):
        b.add_individual(f"c{i}", sex="female", phenotype="unaffected", father="dad", mother="mum")
    b.add_union("dad", "mum", broken=True)
    return b.build()


def _relationship(layout):
    (edge,) = [
        e
        for e in layout.edges
        if e.edge_type == EdgeType.RELATIONSHIP and e.line_style == LineStyle.BROKEN
    ]
    return edge


def _along(p, start=LINE[0], end=LINE[1]) -> float:
    """Position of p along the line, from its start."""
    dx, dy = end[0] - start[0], end[1] - start[1]
    return ((p[0] - start[0]) * dx + (p[1] - start[1]) * dy) / math.hypot(dx, dy)


class TestGeometry:
    def test_two_parallel_slashes_leaning_right(self):
        _, slashes = symbols.separation_marks(*LINE, HALF, descent=False)
        assert len(slashes) == 2
        for (x0, y0), (x1, y1) in slashes:
            assert y1 > y0 and x1 > x0  # "/" with y up
        (a0, a1), (b0, b1) = slashes
        assert (a1[0] - a0[0], a1[1] - a0[1]) == pytest.approx((b1[0] - b0[0], b1[1] - b0[1]))

    def test_line_stops_round_the_slashes(self):
        line, slashes = symbols.separation_marks(*LINE, HALF, descent=True)
        (s0, s1), (e0, e1) = line
        assert s0 == pytest.approx(LINE[0]) and e1 == pytest.approx(LINE[1])
        ends = [_along(p) for seg in slashes for p in seg]
        assert _along(s1) < min(ends) and max(ends) < _along(e0)

    def test_centred_without_descent(self):
        line, _ = symbols.separation_marks(*LINE, HALF, descent=False)
        gap = (_along(line[0][1]) + _along(line[1][0])) / 2
        assert gap == pytest.approx(0.35)

    def test_beside_the_line_of_descent(self):
        """With children, in the middle of the room between the partner and the
        line of descent, which hangs from a whole piece of the line."""
        line, slashes = symbols.separation_marks(*LINE, HALF, descent=True)
        middle = 0.35
        gap = (_along(line[0][1]) + _along(line[1][0])) / 2
        assert gap == pytest.approx(middle / 2)
        assert _along(line[1][0]) < middle
        assert max(_along(p) for seg in slashes for p in seg) < middle

    def test_two_distinct_marks_at_the_text_scale(self):
        """White between the slashes, measured across them, at the SVG's scale."""
        _, ((a0, a1), (b0, _)) = symbols.separation_marks(*LINE, HALF, descent=True)
        ux, uy = a1[0] - a0[0], a1[1] - a0[1]
        across = abs((b0[0] - a0[0]) * uy - (b0[1] - a0[1]) * ux) / math.hypot(ux, uy)
        stroke = STYLE.line_width * symbols.UNITS_PER_POINT
        assert across - stroke >= 2 * stroke

    def test_shrinks_on_a_short_line(self):
        short = ((0.0, 0.0), (0.12, 0.0))
        line, slashes = symbols.separation_marks(*short, HALF, descent=True)
        middle = 0.06
        assert line[0][1][0] > 0.0  # some line is left by the partner
        assert min(p[0] for s in slashes for p in s) > 0.0
        assert max(p[0] for s in slashes for p in s) < line[1][0][0] < middle

    def test_svg_direction_mirrors(self):
        _, up = symbols.separation_marks(*LINE, HALF, descent=True)
        _, down = symbols.separation_marks(*LINE, HALF, descent=True, up=-1.0)
        for (p0, p1), (q0, q1) in zip(up, down, strict=True):
            assert q0 == pytest.approx((p0[0], -p0[1]))
            assert q1 == pytest.approx((p1[0], -p1[1]))

    def test_turned_with_the_line(self):
        """Along an oblique relationship line (circular layout), the marks turn with it."""
        _, flat = symbols.separation_marks(*LINE, HALF, descent=True)
        angle = math.radians(30)
        c, s = math.cos(angle), math.sin(angle)
        _, turned = symbols.separation_marks((0.0, 0.0), (0.7 * c, 0.7 * s), HALF, descent=True)
        for (p0, p1), (q0, q1) in zip(flat, turned, strict=True):
            for p, q in ((p0, q0), (p1, q1)):
                assert q == pytest.approx((p[0] * c - p[1] * s, p[0] * s + p[1] * c))


class TestLayout:
    def test_descent_flag(self):
        assert _relationship(compute_layout(_couple(1).families["F"], STYLE)).has_descent
        assert not _relationship(compute_layout(_couple(0).families["F"], STYLE)).has_descent

    def test_circular_line_from_outline_to_outline(self):
        """The circular layout's line ran centre to centre: the marks fell inside
        the symbols. A separated couple's goes from outline to outline."""
        layout = compute_layout_circular(_couple(2).families["F"], STYLE)
        edge = _relationship(layout)
        assert edge.has_descent
        centres = {n.individual_id: (n.x, n.y) for n in layout.nodes}
        for point in (edge.start, edge.end):
            nearest = min(math.dist(point, c) for c in centres.values())
            assert nearest >= HALF - 1e-9


class TestSvg:
    def _lines(self, ped, cls):
        root = ET.fromstring(render_interactive_svg(ped))
        return root, [e for e in root.iter(f"{SVG}line") if e.get("class") == cls]

    def test_two_slashes_no_dashes(self):
        root, marks = self._lines(_couple(1), "separation-mark")
        assert len(marks) == 2
        assert [e for e in root.iter(f"{SVG}line") if e.get("stroke-dasharray")] == []

    def test_left_of_the_descent_line(self):
        ped = _couple(1)
        _, marks = self._lines(ped, "separation-mark")
        _, pieces = self._lines(ped, "relationship")
        assert len(pieces) == 2
        edge = _relationship(compute_layout(ped.families["F"], STYLE))
        middle = (edge.start[0] + edge.end[0]) / 2
        right = max(max(float(m.get("x1")), float(m.get("x2"))) for m in marks)
        assert right < float(pieces[1].get("x1")) < middle


class TestPng:
    def test_two_marks_clear_at_6pt_in_a_small_figure(self):
        """In a 4 cm figure at 6 pt (as the paper's notation figure), at 600 dpi:
        a row just above the relationship line crosses two separate slashes
        between the father's square and the line of descent, and white on
        either side of them."""
        dpi = 600
        ped = _couple(1)
        style = PedigreeStyle(font_size=6.0)
        fig = plt.figure(figsize=(4 / 2.54, 3 / 2.54), dpi=dpi)
        try:
            ax = fig.add_axes((0, 0, 1, 1))
            render_pedigree(ped, style=style, ax=ax, show_legend=False)
            assert len([a for a in ax.lines if a.get_gid() == "separation-mark"]) == 2
            fig.canvas.draw()
            image = np.asarray(fig.canvas.buffer_rgba())[:, :, :3].mean(axis=2)
            height = image.shape[0]
            layout = compute_layout(ped.families["F"], style)
            edge = _relationship(layout)
            dad = next(n for n in layout.nodes if n.individual_id == "dad")
            sx, _ = ax.transData.transform((dad.x + HALF, 0.0))
            mx, _ = ax.transData.transform(((edge.start[0] + edge.end[0]) / 2, edge.start[1]))
            rise = symbols.SEPARATION_MARK_HEIGHT * HALF / 2
            _, ry = ax.transData.transform((0.0, edge.start[1] + rise / 2))
            stroke = style.line_width * dpi / 72
            row = int(round(height - ry))
            dark = image[row, int(sx + stroke) : int(mx - stroke)] < 128
            runs = int(dark[0]) + int(np.count_nonzero(dark[1:] & ~dark[:-1]))
            assert runs == 2
            assert not dark[0] and not dark[-1]
        finally:
            plt.close(fig)
