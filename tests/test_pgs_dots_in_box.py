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


class TestDotSize:
    """Dots drawn smaller when the box cannot hold all of an individual's dots
    side by side, `PGS_DOT_SPREAD` diameters apart, so they never cover each
    other (#147); in the paper's Fig. 1b four 3 pt dots were 2.5 pt apart."""

    def test_full_size_when_they_fit(self):
        from pypedigree.phenotypes import PGS_DOT_PT, pgs_dot_size

        pt = 0.02
        width = PGS_DOT_PT * pt * (1 + 3 * PGS_DOT_SPREAD) * 1.01
        assert pgs_dot_size(4, width, pt) == pytest.approx(PGS_DOT_PT * pt)
        assert pgs_dot_size(1, 0.001, pt) == pytest.approx(PGS_DOT_PT * pt)

    def test_smaller_to_fit(self):
        from pypedigree.phenotypes import PGS_DOT_PT, pgs_dot_size

        pt = 0.03
        width = 10.6 * pt  # Fig. 1b's box
        dot = pgs_dot_size(4, width, pt)
        assert dot < PGS_DOT_PT * pt
        assert dot * (1 + 3 * PGS_DOT_SPREAD) == pytest.approx(width)

    def test_floor(self):
        from pypedigree.phenotypes import PGS_DOT_MIN_PT, pgs_dot_size

        assert pgs_dot_size(6, 0.01, 0.02) == pytest.approx(PGS_DOT_MIN_PT * 0.02)

    @staticmethod
    def _apart(centres: list[float], diameters: list[float]) -> bool:
        pairs = sorted(zip(centres, diameters, strict=True))
        return all(
            b[0] - a[0] >= (a[1] + b[1]) / 2 - 1e-6 for a, b in zip(pairs, pairs[1:], strict=False)
        )

    def test_png_dots_apart_when_drawn_small(self):
        """Four coinciding scores, font 6, in an 8 cm x 5 cm axes (a box about
        12 pt wide, as in the paper's Fig. 1b): every dot in its box and none
        covering another."""
        fig, ax = plt.subplots(figsize=(8 / 2.54, 5 / 2.54), dpi=300)
        style = PedigreeStyle(font_size=6, pgs_scale="zscore")
        try:
            render_pedigree(
                parse_pedigree_from_string(_ped(4)), style=style, ax=ax, show_legend=False
            )
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            boxes = [
                p.get_window_extent(renderer)
                for p in ax.patches
                if type(p).__name__ == "FancyBboxPatch" and p.get_width() < 0.5
            ]
            dots = [
                line
                for line in ax.findobj(Line2D)
                if line.get_marker() == "o" and len(line.get_xdata())
            ]
            assert len(boxes) == 3 and len(dots) == 12
            per_box: dict[int, list] = {}
            for d in dots:
                ext = d.get_window_extent(renderer)
                k = min(
                    range(3), key=lambda i: abs((boxes[i].x0 + boxes[i].x1) - (ext.x0 + ext.x1))
                )
                assert boxes[k].x0 - 0.5 <= ext.x0 and ext.x1 <= boxes[k].x1 + 0.5
                diameter = d.get_markersize() * fig.dpi / 72
                per_box.setdefault(k, []).append(((ext.x0 + ext.x1) / 2, diameter))
            for group in per_box.values():
                assert len(group) == 4
                assert group[0][1] < 3.0 * fig.dpi / 72  # smaller than 3 pt
                assert self._apart([c for c, _ in group], [w for _, w in group])
        finally:
            plt.close(fig)

    def test_svg_dots_apart_at_a_larger_text_scale(self):
        svg = render_interactive_svg(
            parse_pedigree_from_string(_ped(4)),
            style=PedigreeStyle(pgs_scale="zscore", text_scale=1.5),
        )
        for group in svg.split('<g class="pgs">')[1:]:
            body = group.split("</g>")[0]
            x, w = (
                float(v)
                for v in re.search(
                    r'<rect x="([-0-9.]+)" y="[-0-9.]+" width="([0-9.]+)"', body
                ).groups()
            )
            dots = [
                (float(c), 2 * float(r))
                for c, r in re.findall(r'<circle cx="([-0-9.]+)" cy="[-0-9.]+" r="([0-9.]+)"', body)
            ]
            assert len(dots) == 4
            for cx, d in dots:
                assert x - 1e-3 <= cx - d / 2 and cx + d / 2 <= x + w + 1e-3
            assert self._apart([c for c, _ in dots], [d for _, d in dots])


class TestBoxWidth:
    """The box grows to hold all of an individual's dots side by side at full
    size, at the style's text scale (#147)."""

    def test_usual_width_when_they_fit(self):
        from pypedigree.phenotypes import pgs_box_width

        style = PedigreeStyle()
        usual = style.shape_size / 2 * 0.8
        assert pgs_box_width(["PGS_a", "PGS_b"], style) == (pytest.approx(usual), 0.0)

    def test_wider_for_many_dots_at_a_larger_scale(self):
        from pypedigree import symbols
        from pypedigree.phenotypes import PGS_DOT_PT, pgs_box_width

        style = PedigreeStyle(text_scale=1.5)
        width, wider = pgs_box_width([f"PGS_{k}" for k in range(4)], style)
        dot = PGS_DOT_PT * symbols.units_per_point(1.5)
        assert width == pytest.approx(dot * (1 + 3 * PGS_DOT_SPREAD))
        assert wider == pytest.approx(width - style.shape_size / 2 * 0.8)

    def test_kinds_counted_apart(self):
        from pypedigree.phenotypes import pgs_dot_count

        assert pgs_dot_count(["PGS_a", "PGS_b", "pTDT_a"]) == 2
        assert pgs_dot_count(["pTDT_a", "pTDT_b", "pTDT_c", "PGS_a"]) == 3

    def test_svg_full_size_dots_in_the_wider_box(self):
        from pypedigree.phenotypes import PGS_DOT_PT

        style = PedigreeStyle(pgs_scale="zscore", text_scale=1.5)
        svg = render_interactive_svg(parse_pedigree_from_string(_ped(4)), style=style)
        pt = 0.02 * 1.5
        for group in svg.split('<g class="pgs">')[1:]:
            body = group.split("</g>")[0]
            width = float(re.search(r'<rect x="[-0-9.]+" y="[-0-9.]+" width="([0-9.]+)"', body)[1])
            assert width == pytest.approx(PGS_DOT_PT * pt * (1 + 3 * PGS_DOT_SPREAD), abs=1e-3)
            radii = [
                float(r)
                for r in re.findall(r'<circle cx="[-0-9.]+" cy="[-0-9.]+" r="([0-9.]+)"', body)
            ]
            assert radii == [pytest.approx(PGS_DOT_PT / 2 * pt, abs=1e-3)] * 4

    def test_layout_keeps_room_for_the_wider_box(self):
        from pypedigree import compute_layout
        from pypedigree.config import apply_pgs_spacing
        from pypedigree.phenotypes import pgs_box_width

        names = [f"PGS_s{k}" for k in range(4)]
        family = parse_pedigree_from_string(_ped(4)).families["FAM1"]
        style = PedigreeStyle(text_scale=1.5)
        usual = apply_pgs_spacing(style, names)
        wide = apply_pgs_spacing(style, names, 0.0, pgs_box_width(names, style))
        assert wide.sibling_gap - usual.sibling_gap == pytest.approx(pgs_box_width(names, style)[1])
        xs = {n.individual_id: n.x for n in compute_layout(family, wide).nodes}
        xu = {n.individual_id: n.x for n in compute_layout(family, usual).nodes}
        assert xs["M"] - xs["F"] > xu["M"] - xu["F"]
