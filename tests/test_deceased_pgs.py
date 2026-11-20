"""The deceased line stays clear of the PGS box (#80).

The line overshoots the symbol by `DECEASED_OVERSHOOT` half-sizes at each
end, while the box sat `PGS_BOX_GAP` beside it: the line's end reached the
box's corner, where a low score's marker sits. All pedigrees here are
synthetic.
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
from pypedigree.models import Sex
from pypedigree.symbols import (
    DECEASED_OVERSHOOT,
    PGS_BOX_GAP,
    PGS_DECEASED_MARGIN,
    pgs_box_near,
)

# F is deceased with the family's lowest score (marker on the box's bottom
# edge); M, the right-hand partner (box on the right, #78), is deceased with
# the highest (marker on the top edge, where the line's upper end points).
PED = (
    "FID\tIID\tFather\tMother\tSex\tPhenotype\tdeceased\tPGS\n"
    "FAM1\tF\t0\t0\t1\t1\t1\t-1.2\n"
    "FAM1\tM\t0\t0\t2\t1\t1\t1.4\n"
    "FAM1\tC\tF\tM\t1\t2\t0\t0.1\n"
)


# Least distance between the line and a box's edge (its rounded pad
# included): the fix keeps 0.12 half-sizes, 0.048 units at the default size;
# without it the line's end came within 0.02 of the box's corner.
MIN_CLEARANCE = 0.04


def _distance(p, q, x0, y0, x1, y1) -> float:
    """Least distance from segment pq to the rectangle (0 if they meet)."""
    best = float("inf")
    for i in range(201):
        t = i / 200
        x, y = p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])
        dx = max(x0 - x, 0.0, x - x1)
        dy = max(y0 - y, 0.0, y - y1)
        best = min(best, (dx * dx + dy * dy) ** 0.5)
    return best


class TestGap:
    def test_struck_symbol_box_clears_the_overshoot(self):
        half = 0.4
        plain = pgs_box_near(half, Sex.MALE, False, 10.0, "left")
        struck = pgs_box_near(half, Sex.MALE, False, 10.0, "left", struck=True)
        assert plain == pytest.approx(half * (1 + PGS_BOX_GAP))
        assert struck == pytest.approx(half * (1 + DECEASED_OVERSHOOT + PGS_DECEASED_MARGIN))
        assert struck > half * (1 + DECEASED_OVERSHOOT)

    def test_other_symbols_unchanged(self):
        assert pgs_box_near(0.4, Sex.FEMALE, False, 10.0, "right") == pytest.approx(0.4 * 1.15)


class TestDrawing:
    def test_svg_line_clear_of_the_boxes(self):
        svg = render_interactive_svg(parse_pedigree_from_string(PED))
        for sid in ("F", "M"):
            group = next(g for g in svg.split('<g class="node"')[1:] if f'data-id="{sid}"' in g)
            box = re.search(
                r'<g class="pgs">\s*<rect x="([-0-9.]+)" y="([-0-9.]+)" width="([0-9.]+)" '
                r'height="([0-9.]+)"',
                group,
            )
            x, y, w, h = (float(v) for v in box.groups())
            pad = 0.02  # the box's rounded pad
            # The deceased line: the one diagonal of the node.
            diagonals = [
                ((float(a), float(b)), (float(c), float(d)))
                for a, b, c, d in re.findall(
                    r'<line x1="([-0-9.]+)" y1="([-0-9.]+)" x2="([-0-9.]+)" y2="([-0-9.]+)"',
                    group,
                )
                if a != c and b != d
            ]
            assert diagonals, sid
            for p, q in diagonals:
                gap = _distance(p, q, x - pad, y - pad, x + w + pad, y + h + pad)
                assert gap >= MIN_CLEARANCE, (sid, gap)

    def test_png_line_clear_of_the_boxes(self):
        fig = render_pedigree(parse_pedigree_from_string(PED))
        try:
            fig.canvas.draw()
            ax = fig.axes[0]
            # The PGS boxes: rounded, narrower than a symbol (squares are too).
            boxes = [
                p
                for p in ax.patches
                if type(p).__name__ == "FancyBboxPatch" and p.get_width() < 0.5
            ]
            diagonals = [
                line
                for line in ax.lines
                if len(line.get_xdata()) == 2
                and line.get_xdata()[0] != line.get_xdata()[1]
                and line.get_ydata()[0] != line.get_ydata()[1]
            ]
            assert len(diagonals) == 2 and boxes
            for line in diagonals:
                (x1, x2), (y1, y2) = line.get_xdata(), line.get_ydata()
                for box in boxes:
                    bx, by = box.get_x(), box.get_y()
                    bw, bh = box.get_width(), box.get_height()
                    pad = 0.02  # the box's rounded pad
                    gap = _distance(
                        (x1, y1), (x2, y2), bx - pad, by - pad, bx + bw + pad, by + bh + pad
                    )
                    assert gap >= MIN_CLEARANCE, gap
        finally:
            plt.close(fig)

    def test_layout_keeps_the_extra_room(self):
        """A deceased person's box is further out: the placer counts it."""
        from pypedigree.config import apply_pgs_spacing
        from pypedigree.layout import compute_layout

        alive = PED.replace("\t1\t-1.2", "\t0\t-1.2").replace("\t1\t1.4", "\t0\t1.4")
        style = apply_pgs_spacing(PedigreeStyle(), ["PGS"])

        def width(text):
            nodes = compute_layout(parse_pedigree_from_string(text).families["FAM1"], style).nodes
            return max(n.x for n in nodes) - min(n.x for n in nodes)

        assert width(PED) >= width(alive)
