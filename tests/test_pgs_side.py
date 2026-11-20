"""The right-hand partner's PGS curve goes on their right (#78).

The box was always drawn left of the symbol: for the partner on the right of
a couple, that is on the line between the partners, which crossed the box,
the curve, the dots and the values. All pedigrees here are synthetic.
"""

from __future__ import annotations

import re

import matplotlib.pyplot as plt
import pytest
from matplotlib.text import Text

from pypedigree import (
    PedigreeStyle,
    parse_pedigree_from_string,
    render_interactive_svg,
    render_pedigree,
)
from pypedigree.config import apply_pgs_spacing
from pypedigree.layout import compute_layout, compute_layout_circular
from pypedigree.phenotypes import pgs_box_left
from pypedigree.symbols import PGS_BOX_WIDTH

HEADER = "FID\tIID\tFather\tMother\tSex\tPhenotype\tPGS_autism\tPGS_ADHD"
# The issue's example.
TRIO = (
    f"{HEADER}\n"
    "FAM1\tF\t0\t0\t1\t1\t0.2\t0.5\n"
    "FAM1\tM\t0\t0\t2\t1\t0.6\t0.4\n"
    "FAM1\tC1\tF\tM\t1\t2\t1.9\t-0.3\n"
)
# Two married siblings: A's spouse (right of A, box on the right) has B,
# whose box is on the left, as right-hand neighbour: two boxes in one gap.
SIBLINGS = (
    f"{HEADER}\n"
    "FAM1\tP1\t0\t0\t1\t1\t0.1\t0.2\n"
    "FAM1\tP2\t0\t0\t2\t1\t0.3\t-0.1\n"
    "FAM1\tA\tP1\tP2\t1\t2\t1.2\t0.4\n"
    "FAM1\tSA\t0\t0\t2\t1\t-1.1\t0.9\n"
    "FAM1\tB\tP1\tP2\t1\t1\t0.5\t-0.6\n"
    "FAM1\tSB\t0\t0\t2\t1\t0.0\t1.3\n"
    "FAM1\tCA\tA\tSA\t2\t1\t0.7\t0.2\n"
    "FAM1\tCB\tB\tSB\t1\t1\t-0.4\t0.1\n"
)
STYLE = PedigreeStyle(pgs_scale="zscore", show_pgs_values=True)


def _layout(text: str, style: PedigreeStyle = STYLE):
    pedigree = parse_pedigree_from_string(text)
    spaced = apply_pgs_spacing(style, pedigree.phenotype_names, 0.3)
    return compute_layout(pedigree.families["FAM1"], spaced), spaced


def _svg_boxes(svg: str) -> dict[str, tuple[float, float, float, float]]:
    """Each node's PGS box with its values: (x0, x1, y0, y1), SVG units."""
    out = {}
    for group in svg.split('<g class="node"')[1:]:
        sid = re.search(r'data-id="([^"]+)"', group).group(1)
        box = re.search(
            r'<g class="pgs">\s*<rect x="([-0-9.]+)" y="([-0-9.]+)" width="([0-9.]+)" '
            r'height="([0-9.]+)"',
            group,
        )
        if not box:
            continue
        x, y, w, h = (float(v) for v in box.groups())
        xs = [x, x + w]
        for vx, anchor, text in re.findall(
            r'<text x="([-0-9.]+)"[^>]*text-anchor="(\w+)"[^>]*class="pgs-value"[^>]*>([^<]*)<',
            group,
        ):
            width = 0.012 * len(text)  # generous: ~0.6 em at 6.4 pt, in SVG units
            v = float(vx)
            xs += [v - width, v] if anchor == "end" else [v, v + width]
        out[sid] = (min(xs), max(xs), y, y + h)
    return out


def _horizontal_lines(svg: str) -> list[tuple[float, float, float]]:
    return [
        (float(x1), float(x2), float(y1))
        for x1, y1, x2, y2 in re.findall(
            r'<line class="edge" x1="([-0-9.]+)" y1="([-0-9.]+)" '
            r'x2="([-0-9.]+)" y2="([-0-9.]+)"',
            svg,
        )
        if abs(float(y1) - float(y2)) < 1e-6
    ]


class TestSide:
    def test_right_hand_partner_only(self):
        layout, _ = _layout(TRIO)
        side = {n.individual_id: n.pgs_side for n in layout.nodes}
        assert side == {"F": "left", "M": "right", "C1": "left"}

    def test_box_left_or_right_of_the_symbol(self):
        layout, style = _layout(SIBLINGS)
        half = style.shape_size / 2
        assert {n.pgs_side for n in layout.nodes} == {"left", "right"}
        for n in layout.nodes:
            left = pgs_box_left(n, style)
            if n.pgs_side == "right":
                assert left > n.x + half
            else:
                assert left + half * PGS_BOX_WIDTH < n.x - half

    def test_circular_layout_too(self):
        """One partner of the couple, the one further round, draws on the right."""
        pedigree = parse_pedigree_from_string(TRIO)
        layout = compute_layout_circular(pedigree.families["FAM1"], STYLE)
        side = {n.individual_id: n.pgs_side for n in layout.nodes}
        assert sorted((side["F"], side["M"])) == ["left", "right"]

    def test_room_kept_right_of_a_right_hand_partner(self):
        """Its right-hand neighbour is pushed out by the box's room (#78)."""
        pedigree = parse_pedigree_from_string(SIBLINGS)
        family = pedigree.families["FAM1"]
        base = apply_pgs_spacing(STYLE, pedigree.phenotype_names, 0.3)
        no_room = PedigreeStyle(**{**base.__dict__, "pgs_box_room": 0.0})

        def gap(style):
            x = {n.individual_id: n.x for n in compute_layout(family, style).nodes}
            (left, right) = sorted(("SA", "B"), key=lambda s: x[s])
            return x[right] - x[left]

        assert gap(base) > gap(no_room)


class TestDrawing:
    @pytest.mark.parametrize("text", [TRIO, SIBLINGS])
    def test_no_partner_line_crosses_a_box_or_its_values(self, text):
        """No line at a box's height runs through it: the partner lines, at the
        symbols' centres, are the only ones there."""
        svg = render_interactive_svg(parse_pedigree_from_string(text), style=STYLE)
        boxes = _svg_boxes(svg)
        assert boxes
        lines = _horizontal_lines(svg)
        assert lines
        for x1, x2, y in lines:
            lo, hi = sorted((x1, x2))
            for sid, (bx0, bx1, by0, by1) in boxes.items():
                if by0 < y < by1:
                    assert bx1 <= lo + 1e-3 or bx0 >= hi - 1e-3, sid

    def test_boxes_and_values_do_not_overlap(self):
        """Two boxes in one gap (a right-hand partner's, the next person's)."""
        svg = render_interactive_svg(parse_pedigree_from_string(SIBLINGS), style=STYLE)
        boxes = list(_svg_boxes(svg).items())
        for i, (a, (ax0, ax1, ay0, ay1)) in enumerate(boxes):
            for b, (bx0, bx1, by0, by1) in boxes[i + 1 :]:
                same_row = ay0 < by1 and by0 < ay1
                assert not (same_row and ax0 < bx1 and bx0 < ax1), (a, b)

    def test_png_values_on_the_far_side(self):
        fig = render_pedigree(parse_pedigree_from_string(TRIO), style=STYLE)
        try:
            fig.canvas.draw()
            values = [t for t in fig.axes[0].findobj(Text) if t.get_text()[:1] in ("+", "−")]
            # M's two values: left-aligned, right of everything else on her row.
            right = [t for t in values if t.get_ha() == "left"]
            assert len(right) == 2
            row_y = right[0].get_position()[1]
            left = [
                t for t in values if t.get_ha() == "right" and abs(t.get_position()[1] - row_y) < 1
            ]
            assert left  # F's, on the left of his box
            assert min(t.get_position()[0] for t in right) > max(t.get_position()[0] for t in left)
        finally:
            plt.close(fig)

    def test_evaluated_right_partner_box_clears_the_asterisk(self):
        text = (
            "FID\tIID\tFather\tMother\tSex\tPhenotype\tevaluated\tPGS\n"
            "FAM1\tF\t0\t0\t1\t1\t0\t0.2\n"
            "FAM1\tW\t0\t0\t2\t1\t1\t0.1\n"
            "FAM1\tC\tF\tW\t1\t1\t0\t0.4\n"
        )
        fig = render_pedigree(parse_pedigree_from_string(text))
        try:
            fig.canvas.draw()
            ax = fig.axes[0]
            star = next(t for t in ax.texts if t.get_text() == "*").get_window_extent()
            # The PGS boxes: rounded, narrower than a symbol (squares are too).
            boxes = [
                p.get_window_extent()
                for p in ax.patches
                if type(p).__name__ == "FancyBboxPatch" and p.get_width() < 0.5
            ]
            assert boxes
            for box in boxes:
                overlap_x = box.x0 < star.x1 and star.x0 < box.x1
                overlap_y = box.y0 < star.y1 and star.y0 < box.y1
                assert not (overlap_x and overlap_y)
        finally:
            plt.close(fig)

    def test_caption_below_a_box_on_the_right(self):
        """The proband marker is at the lower left: a box on the right keeps
        its caption below the curve."""
        text = (
            "FID\tIID\tFather\tMother\tSex\tPhenotype\tproband\tPGS\n"
            "FAM1\tF\t0\t0\t1\t1\t0\t0.2\n"
            "FAM1\tM\t0\t0\t2\t2\t1\t0.6\n"
            "FAM1\tC\tF\tM\t1\t1\t0\t0.4\n"
        )
        svg = render_interactive_svg(parse_pedigree_from_string(text))
        group = next(g for g in svg.split('<g class="node"')[1:] if 'data-id="M"' in g)
        box = re.search(
            r'<g class="pgs">\s*<rect x="[-0-9.]+" y="([-0-9.]+)"[^>]*height="([0-9.]+)"', group
        )
        box_bottom = float(box.group(1)) + float(box.group(2))
        caption_y = float(re.search(r'<text [^>]*y="([-0-9.]+)"[^>]*>PGS</text>', group).group(1))
        assert caption_y > box_bottom  # SVG y points down: below the box
