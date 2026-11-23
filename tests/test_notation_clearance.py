"""Notation that keeps clear of the symbol: the evaluation asterisk (#110, #134), the
inner letter (#105) and the PGS box beside a diamond (#111).

- The asterisk's gap is in ems of its font, from the outside of the outline's
  stroke to the glyph's ink, on every shape and at any scale; it touched
  circles before. Round adoption brackets it goes outside the right one, with
  the same gap from the bracket.
- The pregnancy "P" (and "D", "S") is written at the centre whatever the
  phenotype columns; band and quadrant text keeps clear of it.
- A diamond reaches `DIAMOND_SCALE` half-sizes, and its PGS box keeps the same
  gap from its tip as from a square's side.

All pedigrees here are synthetic.
"""

from __future__ import annotations

import math
import re
import xml.etree.ElementTree as ET

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.patches import FancyBboxPatch
from matplotlib.text import Annotation, Text

from pypedigree import (
    PedigreeStyle,
    parse_pedigree_from_string,
    render_interactive_svg,
    render_pedigree,
    symbols,
)
from pypedigree.config import apply_pgs_spacing
from pypedigree.layout import compute_layout
from pypedigree.models import Sex
from pypedigree.phenotypes import clear_of_letter, pgs_box_left

SVG_NS = {"svg": "http://www.w3.org/2000/svg"}
HALF = PedigreeStyle().shape_size / 2
MARK = PedigreeStyle().text_size(symbols.EVALUATION_MARK_FONT_SCALE)
LINE = PedigreeStyle().shape_line_width


def _wide(n: int = 12) -> str:
    """The issue's pedigree: two parents and n children, everyone evaluated."""
    rows = ["FID\tID\tFather\tMother\tSex\tPhenotype\tevaluated", "F1\tgf\t0\t0\t1\t1\t1"]
    rows.append("F1\tgm\t0\t0\t2\t2\t1")
    for i in range(n):
        rows.append(f"F1\tc{i}\tgf\tgm\t{i % 2 + 1}\t{1 + (i % 3 == 0)}\t1")
    return "\n".join(rows) + "\n"


def _ink_gap(sex: Sex, upp: float, triangle: bool = False) -> float:
    """Gap from the outside of the outline's stroke to the asterisk's ink."""
    cx, cy = symbols.evaluation_mark_position(0.0, 0.0, HALF, sex, MARK, LINE, upp, triangle)
    hw = symbols.EVALUATION_MARK_INK_HALF_W_EM * MARK * upp
    hh = symbols.EVALUATION_MARK_INK_HALF_H_EM * MARK * upp
    x0, x1, y0, y1 = cx - hw, cx + hw, cy - hh, cy + hh
    stroke = LINE / 2 * upp
    if triangle:
        return x0 - HALF * symbols.TRIANGLE_HALF_WIDTH - stroke
    if sex == Sex.MALE:
        return x0 - HALF - stroke
    if sex == Sex.FEMALE:
        # Nearest point of the ink's box to the centre.
        nx, ny = min(max(0.0, x0), x1), min(max(0.0, y0), y1)
        return math.hypot(nx, ny) - HALF - stroke
    reach = HALF * symbols.DIAMOND_SCALE
    corners = [(x, y) for x in (x0, x1) for y in (y0, y1)]
    return min((abs(x) + abs(y) - reach) / math.sqrt(2) for x, y in corners) - stroke


class TestAsteriskGap:
    @pytest.mark.parametrize("upp", [symbols.UNITS_PER_POINT, 0.06])
    @pytest.mark.parametrize("sex", [Sex.MALE, Sex.FEMALE, Sex.UNKNOWN])
    def test_gap_in_ems_at_any_scale(self, sex, upp):
        """At the SVG's scale and drawn three times smaller, the same gap in points."""
        gap = _ink_gap(sex, upp)
        want = symbols.EVALUATION_MARK_GAP_EM * MARK * upp
        assert gap >= want - 1e-9
        assert gap <= want * 1.15, "not pushed far away either"

    @pytest.mark.parametrize("upp", [symbols.UNITS_PER_POINT, 0.06])
    def test_triangle(self, upp):
        want = symbols.EVALUATION_MARK_GAP_EM * MARK * upp
        assert _ink_gap(Sex.UNKNOWN, upp, triangle=True) == pytest.approx(want)

    @pytest.mark.parametrize("sex", [Sex.MALE, Sex.FEMALE, Sex.UNKNOWN])
    def test_lower_right_and_svg_mirror(self, sex):
        x, y = symbols.evaluation_mark_position(0.0, 0.0, HALF, sex, MARK, LINE, 0.02)
        assert x > 0 and y < 0
        xs, ys = symbols.evaluation_mark_position(0.0, 0.0, HALF, sex, MARK, LINE, 0.02, up=-1.0)
        assert (xs, ys) == pytest.approx((x, -y))

    def test_svg_places_the_ink_middle(self):
        svg = render_interactive_svg(parse_pedigree_from_string(_wide(2)))
        root = ET.fromstring(svg)
        style = PedigreeStyle()
        for group in root.iter("{http://www.w3.org/2000/svg}g"):
            if group.get("class") != "node":
                continue
            mark = group.find("svg:text[@class='documented-evaluation']", SVG_NS)
            assert mark is not None
            sex = {"gf": Sex.MALE, "gm": Sex.FEMALE}.get(group.get("data-id"))
            if sex is None:
                continue
            nx, ny = float(group.get("data-x")), float(group.get("data-y"))
            ex, ey = symbols.evaluation_mark_position(
                nx, ny, HALF, sex, MARK, style.shape_line_width, 0.02, up=-1.0
            )
            size = MARK * 0.02
            assert float(mark.get("x")) == pytest.approx(ex, abs=1e-3)
            baseline = ey + symbols.EVALUATION_MARK_INK_MIDDLE_EM * size
            assert float(mark.get("y")) == pytest.approx(baseline, abs=1e-3)

    def test_png_gap_in_pixels(self):
        """The issue's case at 300 dpi: white between every outline and its asterisk."""
        dpi = 300
        fig = render_pedigree(parse_pedigree_from_string(_wide()), show_legend=False)
        try:
            fig.set_dpi(dpi)
            fig.canvas.draw()
            image = np.asarray(fig.canvas.buffer_rgba())[:, :, :3].mean(axis=2)
            height = image.shape[0]
            ax = fig.axes[0]
            marks = [t for t in ax.texts if isinstance(t, Annotation)]
            assert len(marks) == 14
            px_per_pt = dpi / 72
            want = symbols.EVALUATION_MARK_GAP_EM * MARK * px_per_pt
            for mark in marks:
                sx, sy = ax.transData.transform(mark.xy)
                ox, oy = mark.xyann
                oy += symbols.EVALUATION_MARK_INK_MIDDLE_EM * MARK  # back to the ink's middle
                length = math.hypot(ox, oy)
                ux, uy = ox / length, oy / length
                # From inside the outline outwards: stroke, white, then ink.
                run = []
                for step in range(-12, int(length * px_per_pt) + 4):
                    x = sx + ux * step
                    y = sy + uy * step
                    run.append(image[int(round(height - y)), int(round(x))] < 128)
                end_stroke = next(i for i in range(len(run) - 1) if run[i] and not run[i + 1])
                ink = next(i for i in range(end_stroke + 1, len(run)) if run[i])
                gap = ink - end_stroke - 1
                assert gap >= 0.7 * want, (mark.xy, gap, want)
        finally:
            plt.close(fig)


class TestLetterInside:
    PGS_ONLY = (
        "FID\tIID\tPAT\tMAT\tSEX\tPHENO\tpregnancy\tPGS_x\n"
        "1\tA\t0\t0\t1\t1\t0\t0.1\n"
        "1\tB\t0\t0\t2\t1\t0\t0.5\n"
        "1\tC\tA\tB\t0\t2\tP\t0.9\n"
    )
    BANDS = (
        "FID\tIID\tPAT\tMAT\tSEX\tPHENO\tpregnancy\tart\tg_IQ\tg_SRS\tg_BMI\n"
        "1\tA\t0\t0\t1\t1\t0\t0\t95\t60\t22\n"
        "1\tB\t0\t0\t2\t2\t0\tD\t110\t70\t27\n"
        "1\tC\tA\tB\t0\t1\tP\t0\t102\t55\t24\n"
        "1\tE\tA\tB\t1\t2\tP\t0\t88\t80\t31\n"
    )
    QUADRANTS = (
        "FID\tIID\tPAT\tMAT\tSEX\tPHENO\tpregnancy\tautism\n"
        "1\tA\t0\t0\t1\t1\t0\t2\n"
        "1\tB\t0\t0\t2\t1\t0\t1\n"
        "1\tC\tA\tB\t0\t1\tP\t1\n"
    )
    ALL_NA = (
        "FID\tIID\tPAT\tMAT\tSEX\tPHENO\tpregnancy\tg_IQ\tg_SRS\n"
        "1\tA\t0\t0\t1\t1\t0\t95\t60\n"
        "1\tB\t0\t0\t2\t1\t0\t110\t70\n"
        "1\tC\tA\tB\t0\t1\tP\tNA\tNA\n"
    )

    @staticmethod
    def _png(text: str, **style):
        ped = parse_pedigree_from_string(text)
        fig = render_pedigree(ped, style=PedigreeStyle(**style))
        fig.canvas.draw()
        return fig

    @staticmethod
    def _letters(fig) -> list[Text]:
        return [t for t in fig.axes[0].texts if t.get_text() in ("P", "D", "S")]

    @pytest.mark.parametrize(
        ("text", "style"),
        [
            (PGS_ONLY, {}),
            (QUADRANTS, {}),
            (BANDS, {"phenotype_layout": "bands"}),
            (BANDS, {"phenotype_layout": "quadrants"}),
        ],
    )
    def test_png_letter_at_the_centre(self, text, style):
        spaced = PedigreeStyle()
        if "PGS" in text:
            spaced = apply_pgs_spacing(spaced, ["PGS_x"])
        layout = compute_layout(parse_pedigree_from_string(text).families["1"], spaced)
        centres = {(round(n.x, 6), round(n.y, 6)) for n in layout.nodes}
        fig = self._png(text, **style)
        try:
            letters = self._letters(fig)
            assert letters
            for letter in letters:
                x, y = letter.get_position()
                assert (round(x, 6), round(y, 6)) in centres
                assert letter.get_horizontalalignment() == "center"
        finally:
            plt.close(fig)

    def test_png_backing_only_on_a_divided_symbol(self):
        fig = self._png(self.PGS_ONLY)
        try:
            (letter,) = self._letters(fig)
            assert letter.get_bbox_patch() is None
            # C is affected: a filled diamond, white letter.
            assert letter.get_color() in ("white", "#ffffff", (1.0, 1.0, 1.0))
        finally:
            plt.close(fig)
        fig = self._png(self.BANDS, phenotype_layout="bands")
        try:
            assert all(t.get_bbox_patch() is not None for t in self._letters(fig))
        finally:
            plt.close(fig)

    @pytest.mark.parametrize(
        ("text", "style"),
        [
            (BANDS, {"phenotype_layout": "bands"}),
            (BANDS, {"phenotype_layout": "quadrants"}),
            (ALL_NA, {"phenotype_layout": "bands"}),
            (BANDS, {"phenotype_layout": "bands", "band_text": "name_value"}),
        ],
    )
    def test_png_band_text_clear_of_the_letter(self, text, style):
        fig = self._png(text, **style)
        try:
            renderer = fig.canvas.get_renderer()
            texts = [t for t in fig.findobj(Text) if t.get_text() and t.axes is fig.axes[0]]
            for letter in self._letters(fig):
                box = letter.get_window_extent(renderer)
                for other in texts:
                    if other is letter or other.get_text() in ("P", "D", "S"):
                        continue
                    assert not box.overlaps(other.get_window_extent(renderer)), other.get_text()
        finally:
            plt.close(fig)

    @pytest.mark.parametrize("layout", ["bands", "quadrants"])
    def test_svg_letter_at_the_centre_text_beside_it(self, layout):
        svg = render_interactive_svg(
            parse_pedigree_from_string(self.BANDS), style=PedigreeStyle(phenotype_layout=layout)
        )
        root = ET.fromstring(svg)
        style = PedigreeStyle()
        letter_size = style.text_size() * 0.02
        for group in root.iter("{http://www.w3.org/2000/svg}g"):
            if group.get("class") != "node":
                continue
            letter = group.find("svg:text[@class='symbol-letter']", SVG_NS)
            if letter is None:
                continue
            cx, cy = float(group.get("data-x")), float(group.get("data-y"))
            assert float(letter.get("x")) == pytest.approx(cx, abs=1e-3)
            assert letter.get("text-anchor", "middle") == "middle"
            reach = symbols.LETTER_CLEAR_HALF_W_EM * letter_size
            for band in group.findall("svg:text", SVG_NS):
                if band.get("font-weight") != "bold" or band is letter:
                    continue
                size = float(band.get("font-size"))
                # Centred text: its baseline is 0.35 em below its middle.
                y = float(band.get("y")) - 0.35 * size
                reach_y = symbols.LETTER_CLEAR_HALF_H_EM * letter_size + 0.4 * size
                if abs(y - cy) >= reach_y:
                    continue  # not at the letter's height
                width = symbols.text_width(band.text, size / 0.02, bold=True)
                x = float(band.get("x"))
                assert min(abs(x - width / 2 - cx), abs(x + width / 2 - cx)) >= reach - 1e-6
                assert x - width / 2 >= cx + reach - 1e-6 or x + width / 2 <= cx - reach + 1e-6

    def test_clear_of_letter_geometry(self):
        pt = 0.02
        # A band centred on the letter's height: right of it, in what is left.
        room, dx = clear_of_letter(0.8, 0.1, (0.0, 0.0), 8.0, pt)
        half_w = symbols.LETTER_CLEAR_HALF_W_EM * 8 * pt
        gap = symbols.LETTER_CLEAR_GAP_EM * 8 * pt
        assert room == pytest.approx(0.4 - half_w - gap)
        assert dx - room / 2 == pytest.approx(half_w + gap)
        # A band far above: untouched.
        assert clear_of_letter(0.8, 0.1, (0.0, 0.3), 8.0, pt) == (0.8, 0.0)
        # A quadrant's text: same place, less room.
        room, dx = clear_of_letter(0.4, 0.4, (0.2, 0.2), 8.0, pt)
        assert dx == 0.0 and room == pytest.approx(2 * (0.2 - half_w - gap))
        # No letter: untouched.
        assert clear_of_letter(0.8, 0.1, (0.0, 0.0), None, pt) == (0.8, 0.0)


class TestDiamondPgsBox:
    PED = (
        "FID\tIID\tPAT\tMAT\tSEX\tPHENO\tPGS_x\n"
        "1\tA\t0\t0\t1\t1\t0.1\n"
        "1\tB\t0\t0\t2\t1\t0.5\n"
        "1\tS\tA\tB\t1\t1\t0.3\n"
        "1\tC\tA\tB\t0\t1\t0.9\n"
    )

    def test_box_measured_from_the_tip(self):
        near = symbols.pgs_box_near(HALF, Sex.UNKNOWN, False, MARK, "left")
        assert near == pytest.approx(HALF * (symbols.DIAMOND_SCALE + symbols.PGS_BOX_GAP))
        square = symbols.pgs_box_near(HALF, Sex.MALE, False, MARK, "left")
        assert square == pytest.approx(HALF * (1 + symbols.PGS_BOX_GAP))

    def test_ring_grows_the_diamond(self):
        ring = 0.05
        near = symbols.pgs_box_near(HALF, Sex.UNKNOWN, False, MARK, "left", ring=ring)
        tip = symbols.grown_half(Sex.UNKNOWN, HALF, ring) * symbols.DIAMOND_SCALE
        assert near == pytest.approx(tip + HALF * symbols.PGS_BOX_GAP)

    def test_svg_box_clear_of_the_tip(self):
        svg = render_interactive_svg(parse_pedigree_from_string(self.PED))
        group = next(g for g in svg.split('<g class="node"')[1:] if 'data-id="C"' in g)
        box = re.search(
            r'<g class="pgs">\s*<rect x="([-0-9.]+)" y="[-0-9.]+" width="([0-9.]+)"', group
        )
        right = float(box.group(1)) + float(box.group(2))
        points = re.search(r'<polygon class="shape" points="([^"]+)"', group).group(1)
        tip = min(float(p.split(",")[0]) for p in points.split())
        assert tip - right == pytest.approx(HALF * symbols.PGS_BOX_GAP, abs=2e-3)

    def test_layout_counts_the_extra(self):
        """The box sits as far from the sibling on its left as beside a square."""
        style = apply_pgs_spacing(PedigreeStyle(), ["PGS_x"])

        def room(text: str) -> float:
            nodes = {
                n.individual_id: n
                for n in compute_layout(parse_pedigree_from_string(text).families["1"], style).nodes
            }
            s, c = nodes["S"], nodes["C"]
            if s.x > c.x:
                pytest.skip("the sibling is not on the left")
            return pgs_box_left(c, style) - (s.x + HALF)

        square = self.PED.replace("1\tC\tA\tB\t0", "1\tC\tA\tB\t1")
        assert room(self.PED) >= room(square) - 1e-9


def _adopted(kind: str = "in") -> str:
    """Two parents and three adopted, evaluated children: a man, a woman, a diamond."""
    rows = ["FID\tID\tFather\tMother\tSex\tPhenotype\tadopted\tevaluated"]
    rows += ["F1\tdad\t0\t0\t1\t1\t0\t0", "F1\tmum\t0\t0\t2\t1\t0\t0"]
    for i, sex in enumerate((1, 2, 0)):
        rows.append(f"F1\tc{i}\tdad\tmum\t{sex}\t{1 + (i == 1)}\t{kind}\t1")
    return "\n".join(rows) + "\n"


def _bracket_right(sex: Sex, triangle: bool = False) -> float:
    """X of the right adoption bracket's upright, round a symbol at the origin."""
    lines = symbols.bracket_lines(0.0, 0.0, HALF, sex, triangle)
    return max(x for line in lines for x, _ in line)


class TestAsteriskBesideBrackets:
    """An adopted person's asterisk goes outside the right bracket (#134), with the
    same gap in ems as from an outline (#110); it touched the bracket before."""

    BRACKET = PedigreeStyle().line_width

    @pytest.mark.parametrize("upp", [symbols.UNITS_PER_POINT, 0.06])
    @pytest.mark.parametrize("sex", [Sex.MALE, Sex.FEMALE, Sex.UNKNOWN])
    @pytest.mark.parametrize("triangle", [False, True])
    def test_ink_clear_of_the_bracket(self, sex, upp, triangle):
        cx, _ = symbols.evaluation_mark_position(
            0.0, 0.0, HALF, sex, MARK, LINE, upp, triangle, bracket_width=self.BRACKET
        )
        ink_left = cx - symbols.EVALUATION_MARK_INK_HALF_W_EM * MARK * upp
        outside = _bracket_right(sex, triangle) + self.BRACKET / 2 * upp
        want = symbols.EVALUATION_MARK_GAP_EM * MARK * upp
        assert ink_left - outside == pytest.approx(want)

    @pytest.mark.parametrize("sex", [Sex.MALE, Sex.FEMALE, Sex.UNKNOWN])
    def test_same_height_as_without_brackets(self, sex):
        """Pushed sideways only: still at the symbol's lower right."""
        _, y = symbols.evaluation_mark_position(0.0, 0.0, HALF, sex, MARK, LINE, 0.02)
        _, yb = symbols.evaluation_mark_position(
            0.0, 0.0, HALF, sex, MARK, LINE, 0.02, bracket_width=self.BRACKET
        )
        assert yb == pytest.approx(y)
        _, ys = symbols.evaluation_mark_position(
            0.0, 0.0, HALF, sex, MARK, LINE, 0.02, up=-1.0, bracket_width=self.BRACKET
        )
        assert ys == pytest.approx(-y)

    @pytest.mark.parametrize("kind", ["in", "out"])
    def test_svg_asterisk_outside_the_bracket(self, kind):
        root = ET.fromstring(render_interactive_svg(parse_pedigree_from_string(_adopted(kind))))
        checked = 0
        for group in root.iter("{http://www.w3.org/2000/svg}g"):
            if group.get("class") != "node":
                continue
            mark = group.find("svg:text[@class='documented-evaluation']", SVG_NS)
            brackets = group.findall("svg:polyline[@class='adoption-bracket']", SVG_NS)
            if mark is None:
                continue
            assert len(brackets) == 2
            right = max(float(p.split(",")[0]) for b in brackets for p in b.get("points").split())
            stroke = float(brackets[0].get("stroke-width")) / 2
            size = float(mark.get("font-size"))
            ink_left = float(mark.get("x")) - symbols.EVALUATION_MARK_INK_HALF_W_EM * size
            assert ink_left - (right + stroke) >= symbols.EVALUATION_MARK_GAP_EM * size - 1e-3
            checked += 1
        assert checked == 3

    def test_png_white_between_bracket_and_asterisk(self):
        """At 300 dpi, along a row through the asterisk: bracket, white, then ink."""
        dpi = 300
        fig = render_pedigree(parse_pedigree_from_string(_adopted()), show_legend=False)
        try:
            fig.set_dpi(dpi)
            fig.canvas.draw()
            image = np.asarray(fig.canvas.buffer_rgba())[:, :, :3].mean(axis=2)
            height = image.shape[0]
            ax = fig.axes[0]
            marks = [t for t in ax.texts if isinstance(t, Annotation)]
            assert len(marks) == 3
            px_per_pt = dpi / 72
            want = symbols.EVALUATION_MARK_GAP_EM * MARK * px_per_pt
            for mark in marks:
                sx, sy = ax.transData.transform(mark.xy)
                ox, oy = mark.xyann
                oy += symbols.EVALUATION_MARK_INK_MIDDLE_EM * MARK  # the ink's middle
                row = int(round(height - (sy + oy * px_per_pt)))
                x_ink = sx + ox * px_per_pt
                # Walk left from the ink's middle to the end of its ink.
                xs = range(int(x_ink), int(x_ink - (ox + 6) * px_per_pt), -1)
                dark = [image[row, x] < 128 for x in xs]
                ink_end = next(i for i in range(len(dark) - 1) if dark[i] and not dark[i + 1])
                ink_left = int(x_ink) - ink_end
                # The right bracket as drawn: its upright at x[1] == x[2], ticks to its left.
                uprights = [
                    ln for ln in ax.lines
                    if len(ln.get_xdata()) == 4 and ln.get_xdata()[1] == ln.get_xdata()[2]
                    and ln.get_xdata()[0] < ln.get_xdata()[1]
                    and min(ln.get_ydata()) <= mark.xy[1] <= max(ln.get_ydata())
                    and abs(ln.get_xdata()[1] - mark.xy[0]) < 2 * HALF
                ]  # fmt: skip
                assert len(uprights) == 1
                bx, _ = ax.transData.transform((uprights[0].get_xdata()[1], 0.0))
                outside = bx + uprights[0].get_linewidth() / 2 * px_per_pt
                assert ink_left - outside >= 0.7 * want, (mark.xy, ink_left - outside, want)
        finally:
            plt.close(fig)

    def test_pgs_box_clear_of_the_moved_asterisk(self):
        """A box on the right clears the asterisk, now outside the bracket."""
        near = symbols.pgs_box_near(
            HALF, Sex.FEMALE, True, MARK, "right", line_width=LINE, bracket_width=self.BRACKET
        )
        cx, _ = symbols.evaluation_mark_position(
            0.0, 0.0, HALF, Sex.FEMALE, MARK, LINE, symbols.UNITS_PER_POINT,
            bracket_width=self.BRACKET,
        )  # fmt: skip
        ink_right = cx + symbols.EVALUATION_MARK_INK_HALF_W_EM * MARK * symbols.UNITS_PER_POINT
        assert near >= ink_right + HALF * symbols.PGS_BOX_GAP - 1e-9


class TestPgsBoxBesideBrackets:
    """An adopted person's PGS box keeps its gap from the bracket (#136): it is
    measured from the bracket's stroke, not the shape, on either side."""

    BRACKET = PedigreeStyle().line_width
    PED = (
        "FID\tID\tFather\tMother\tSex\tPhenotype\tadopted\tPGS_x\n"
        "1\tF\t0\t0\t1\t1\t0\t0.1\n"
        "1\tM\t0\t0\t2\t1\t0\t-0.2\n"
        "1\tA\tF\tM\t2\t1\tout\t1.2\n"
        "1\tB\tF\tM\t1\t1\tin\t-0.4\n"
        "1\tC\tF\tM\t0\t1\tin\t0.6\n"
        "1\tD\tF\tM\t2\t1\tout\t-1.1\n"
    )

    @pytest.mark.parametrize("side", ["left", "right"])
    @pytest.mark.parametrize("sex", [Sex.MALE, Sex.FEMALE, Sex.UNKNOWN])
    def test_gap_from_the_bracket(self, sex, side):
        near = symbols.pgs_box_near(
            HALF, sex, False, MARK, side, line_width=LINE, bracket_width=self.BRACKET
        )
        outside = symbols.bracket_offset(HALF, sex) + self.BRACKET / 2 * symbols.UNITS_PER_POINT
        assert near == pytest.approx(outside + HALF * symbols.PGS_BOX_GAP)

    def test_no_brackets_unchanged(self):
        near = symbols.pgs_box_near(HALF, Sex.FEMALE, False, MARK, "left", line_width=LINE)
        assert near == pytest.approx(HALF * (1 + symbols.PGS_BOX_GAP))

    def test_ring_grows_the_brackets_and_the_box(self):
        ring = 0.05
        near = symbols.pgs_box_near(
            HALF, Sex.MALE, False, MARK, "left", ring=ring, bracket_width=self.BRACKET
        )
        grown = symbols.grown_half(Sex.MALE, HALF, ring)
        outside = (
            symbols.bracket_offset(grown, Sex.MALE) + self.BRACKET / 2 * symbols.UNITS_PER_POINT
        )
        assert near == pytest.approx(outside + HALF * symbols.PGS_BOX_GAP)

    def test_svg_box_clear_of_the_bracket(self):
        root = ET.fromstring(render_interactive_svg(parse_pedigree_from_string(self.PED)))
        checked = 0
        for group in root.iter("{http://www.w3.org/2000/svg}g"):
            if group.get("class") != "node":
                continue
            brackets = group.findall("svg:polyline[@class='adoption-bracket']", SVG_NS)
            box = group.find("svg:g[@class='pgs']/svg:rect", SVG_NS)
            if not brackets or box is None:
                continue
            box_left = float(box.get("x"))
            box_right = box_left + float(box.get("width"))
            stroke = float(brackets[0].get("stroke-width")) / 2
            xs = [float(p.split(",")[0]) for b in brackets for p in b.get("points").split()]
            cx = float(group.get("data-x"))
            # A box on the left is measured from the left bracket, on the right from the right one.
            left_gap, right_gap = min(xs) - stroke - box_right, box_left - (max(xs) + stroke)
            gap = left_gap if box_right <= cx else right_gap
            assert gap == pytest.approx(HALF * symbols.PGS_BOX_GAP, abs=2e-3)
            checked += 1
        assert checked == 4

    def test_png_box_clear_of_the_bracket(self):
        """The PNG's box is where pgs_box_left puts it, clear of the drawn bracket."""
        fig = render_pedigree(parse_pedigree_from_string(self.PED), show_legend=False)
        try:
            ax = fig.axes[0]
            uprights = [
                ln.get_xdata()[1] for ln in ax.lines
                if len(ln.get_xdata()) == 4 and ln.get_xdata()[1] == ln.get_xdata()[2]
                and ln.get_xdata()[0] > ln.get_xdata()[1]
            ]  # fmt: skip
            assert len(uprights) == 4  # the left bracket of each adopted child
            boxes = [
                p for p in ax.patches
                if isinstance(p, FancyBboxPatch)
                and p.get_width() == pytest.approx(HALF * symbols.PGS_BOX_WIDTH)
                and p.get_height() == pytest.approx(2 * HALF)
            ]  # fmt: skip
            assert len(boxes) == 6
            for left_bracket in uprights:
                box = min(boxes, key=lambda p: abs(p.get_x() + p.get_width() - left_bracket))
                assert left_bracket - (box.get_x() + box.get_width()) > HALF * symbols.PGS_BOX_GAP
        finally:
            plt.close(fig)

    def test_layout_keeps_the_box_off_the_neighbour(self):
        """Each box ends right of its left neighbour's right bracket."""
        style = apply_pgs_spacing(PedigreeStyle(), ["PGS_x"])
        family = parse_pedigree_from_string(self.PED).families["1"]
        nodes = sorted(
            (n for n in compute_layout(family, style).nodes if n.individual_id in "ABCD"),
            key=lambda n: n.x,
        )
        for left, node in zip(nodes, nodes[1:], strict=False):
            assert node.pgs_side == "left"
            box_left = pgs_box_left(node, style)
            assert box_left > left.x + symbols.bracket_offset(HALF, left.sex), node.individual_id
