"""Spread-out PGS dots stay inside their box (#89).

Dots of markers closer than a dot's diameter are set side by side (#54); a
group of three or more could be wider than the box, so the outer dots went
past it and reached the symbol. All pedigrees here are synthetic.
"""

from __future__ import annotations

import re

import matplotlib.pyplot as plt
import pytest
from matplotlib.lines import Line2D

from pypedigree import (
    PedigreeStyle,
    parse_pedigree_from_string,
    render_interactive_svg,
    render_pedigree,
)
from pypedigree.phenotypes import PGS_DOT_SPREAD, spread_dots


def _ped(n_scores: int) -> str:
    """A trio whose scores all coincide: one group of ``n_scores`` dots per
    curve; the mother is a right-hand partner, her box on her right (#78)."""
    cols = "\t".join(f"PGS_s{k}" for k in range(n_scores))
    same = "\t".join(["0.3"] * n_scores)
    return (
        f"FID\tIID\tFather\tMother\tSex\tPhenotype\t{cols}\n"
        f"FAM1\tF\t0\t0\t1\t1\t{same}\n"
        f"FAM1\tM\t0\t0\t2\t1\t{same}\n"
        f"FAM1\tC\tF\tM\t2\t1\t{same}\n"
    )


class TestSpread:
    def test_unchanged_when_the_group_fits(self):
        offsets = spread_dots([0.5, 0.5], 0.1, width=1.0)
        assert offsets == pytest.approx([-0.075, 0.075])  # 1.5 diameters apart
        assert spread_dots([0.5, 0.5], 0.1) == offsets

    def test_packed_closer_to_fit_the_box(self):
        offsets = spread_dots([0.5] * 4, 0.1, width=0.3)
        step = (0.3 - 0.1) / 3
        assert step < 0.1 * PGS_DOT_SPREAD
        assert sorted(offsets) == pytest.approx([-1.5 * step, -0.5 * step, 0.5 * step, 1.5 * step])
        # Every dot within the box: centre +/- radius inside +/- width / 2.
        assert max(abs(o) for o in offsets) + 0.05 <= 0.15 + 1e-12

    def test_box_narrower_than_a_dot(self):
        assert spread_dots([0.5, 0.5, 0.5], 0.1, width=0.05) == [0.0, 0.0, 0.0]

    def test_lone_dots_unchanged(self):
        assert spread_dots([0.0, 1.0], 0.1, width=0.01) == [0.0, 0.0]


@pytest.mark.parametrize("n_scores", [2, 3, 4])
class TestDrawn:
    def test_png_in_a_small_axes(self, n_scores):
        """The issue's case: font 6, drawn into a 16 cm x 4 cm axes (box ~8.5 pt)."""
        fig, ax = plt.subplots(figsize=(16 / 2.54, 4 / 2.54), dpi=300)
        style = PedigreeStyle(font_size=6, pgs_scale="zscore", show_pgs_values=True)
        try:
            render_pedigree(
                parse_pedigree_from_string(_ped(n_scores)), style=style, ax=ax, show_legend=False
            )
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            boxes = [
                p.get_window_extent(renderer)
                for p in ax.patches
                if type(p).__name__ == "FancyBboxPatch" and p.get_width() < 0.5
            ]
            assert len(boxes) == 3
            dots = [
                line.get_window_extent(renderer)
                for line in ax.findobj(Line2D)
                if line.get_marker() == "o" and len(line.get_xdata())
            ]
            assert len(dots) == 3 * n_scores
            for d in dots:
                box = min(boxes, key=lambda b: abs((b.x0 + b.x1) - (d.x0 + d.x1)))
                assert box.x0 - 0.5 <= d.x0 and d.x1 <= box.x1 + 0.5  # half a pixel
        finally:
            plt.close(fig)

    def test_svg(self, n_scores):
        svg = render_interactive_svg(
            parse_pedigree_from_string(_ped(n_scores + 1)), style=PedigreeStyle(pgs_scale="zscore")
        )
        groups = svg.split('<g class="pgs">')[1:]
        assert len(groups) == 3
        for group in groups:
            x, w = (
                float(v)
                for v in re.search(
                    r'<rect x="([-0-9.]+)" y="[-0-9.]+" width="([0-9.]+)"', group
                ).groups()
            )
            body = group.split("</g>")[0]
            dots = re.findall(r'<circle cx="([-0-9.]+)" cy="[-0-9.]+" r="([0-9.]+)"', body)
            assert len(dots) == n_scores + 1
            for cx, r in dots:
                assert x - 1e-3 <= float(cx) - float(r) and float(cx) + float(r) <= x + w + 1e-3
