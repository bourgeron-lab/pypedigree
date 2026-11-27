"""The dotted line joining a person's two copies (#149).

A loop broken by drawing a partner twice joins the copies with a dotted
line. It ran from centre to centre as a shallow curve, through the symbols
between them. In the hierarchical layout it now leaves the top of each copy,
runs in the band above the row (between the symbols and their sibship
lines) and comes back down: it touches no other symbol, and lines sharing a
row keep apart. The circular layout keeps its curves.

All pedigrees here are synthetic.
"""

from __future__ import annotations

import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib.pyplot as plt
import pytest

from pypedigree import (
    PedigreeStyle,
    compute_layout,
    parse_pedigree,
    parse_pedigree_from_string,
    render_interactive_svg,
    render_pedigree,
)
from pypedigree.layout import compute_layout_circular
from pypedigree.models import EdgeType

sys.path.insert(0, str(Path(__file__).parent))
from conftest import TEST_DATA_DIR, generate_realistic_pedigree  # noqa: E402

STYLE = PedigreeStyle()
HALF = STYLE.shape_size / 2
SVG = "{http://www.w3.org/2000/svg}"


def _segments(points):
    return list(zip(points, points[1:], strict=False))


def _hits_box(p, q, cx, cy, half, samples: int = 40) -> bool:
    for k in range(samples + 1):
        t = k / samples
        x, y = p[0] + (q[0] - p[0]) * t, p[1] + (q[1] - p[1]) * t
        if abs(x - cx) < half and abs(y - cy) < half:
            return True
    return False


def _crossed(layout) -> int:
    """Symbols, other than its two ends, a routed line passes through."""
    total = 0
    for arc in layout.duplicate_arcs:
        assert arc.waypoints, "hierarchical arcs are routed"
        for n in layout.nodes:
            if (n.x, n.y) in (arc.from_pos, arc.to_pos):
                continue
            if any(_hits_box(p, q, n.x, n.y, HALF) for p, q in _segments(arc.waypoints)):
                total += 1
    return total


@pytest.fixture(scope="module")
def generated():
    return [generate_realistic_pedigree(size, seed) for seed in range(4) for size in (40, 80)]


class TestRoute:
    def test_no_symbol_crossed(self, generated):
        arcs = 0
        for family in generated:
            layout = compute_layout(family, STYLE)
            arcs += len(layout.duplicate_arcs)
            assert _crossed(layout) == 0
        assert arcs > 20

    def test_paper_fixture(self):
        family = next(
            iter(
                parse_pedigree(
                    TEST_DATA_DIR / "complex" / "multigenerational.ped"
                ).families.values()
            )
        )
        layout = compute_layout(family, STYLE)
        assert _crossed(layout) == 0

    def test_ends_on_the_copies_tops(self, generated):
        layout = compute_layout(generated[0], STYLE)
        for arc in layout.duplicate_arcs:
            (sx, sy), (ex, ey) = arc.waypoints[0], arc.waypoints[-1]
            top = arc.from_pos[1] + HALF
            assert sy == pytest.approx(top) and ey == pytest.approx(top)
            assert abs(sx - arc.from_pos[0]) < HALF and abs(ex - arc.to_pos[0]) < HALF
            # Off the copy's own vertical line, which hangs from the top centre.
            assert abs(sx - arc.from_pos[0]) > 0.1 * HALF

    def test_in_the_band_above_the_row(self, generated):
        """Between the symbols' tops and their sibship lines, which are
        `generation_gap` / 2 above the row."""
        for family in generated:
            layout = compute_layout(family, STYLE)
            for arc in layout.duplicate_arcs:
                y = arc.from_pos[1]
                for _, py in arc.waypoints:
                    assert y + HALF - 1e-9 <= py < y + STYLE.generation_gap / 2

    def test_lines_of_a_row_apart(self, generated):
        """Two lines whose spans overlap run at different heights."""
        for family in generated:
            layout = compute_layout(family, STYLE)
            runs = []
            for arc in layout.duplicate_arcs:
                lane = max(py for _, py in arc.waypoints)
                xs = sorted((arc.from_pos[0], arc.to_pos[0]))
                runs.append((arc.from_pos[1], lane, xs))
            for i, (row, lane, xs) in enumerate(runs):
                for row2, lane2, xs2 in runs[i + 1 :]:
                    if row == row2 and xs[0] < xs2[1] and xs2[0] < xs[1]:
                        assert lane != pytest.approx(lane2)

    def test_circular_keeps_its_curves(self, generated):
        layout = compute_layout_circular(generated[0], STYLE)
        assert layout.duplicate_arcs
        assert all(not a.waypoints for a in layout.duplicate_arcs)

    def test_sibship_lines_where_expected(self, generated):
        """The band's top: sibship lines sit `generation_gap` / 2 above their row."""
        layout = compute_layout(generated[0], STYLE)
        rows = {round(n.y, 6) for n in layout.nodes}
        for e in layout.edges:
            if e.edge_type == EdgeType.SIBSHIP:
                assert round(e.start[1] - STYLE.generation_gap / 2, 6) in rows


# A man whose second partner is his niece: the parser draws him a second
# time for that partnership, on her row, a generation below him.
ROWS_APART = (
    "FID\tIID\tPAT\tMAT\tSEX\tPHENO\n"
    "FAM1\tg1\t0\t0\t1\t1\n"
    "FAM1\tg2\t0\t0\t2\t1\n"
    "FAM1\ta\tg1\tg2\t1\t2\n"
    "FAM1\tb\tg1\tg2\t2\t1\n"
    "FAM1\tbh\t0\t0\t1\t1\n"
    "FAM1\tc\tbh\tb\t2\t1\n"
    "FAM1\td\tbh\tb\t1\t1\n"
    "FAM1\tw1\t0\t0\t2\t1\n"
    "FAM1\tx\ta\tw1\t1\t1\n"
    "FAM1\ty\ta\tc\t2\t2\n"
)


class TestCopiesInTwoRows:
    """The line between copies a generation apart reaches both (the line
    once ran in the upper copy's band only and stopped a row short)."""

    @pytest.fixture(scope="class")
    def layout(self):
        family = parse_pedigree_from_string(ROWS_APART).get_family("FAM1")
        return compute_layout(family, STYLE)

    def test_copies_in_two_rows(self, layout):
        (arc,) = layout.duplicate_arcs
        assert arc.from_pos[1] - arc.to_pos[1] == pytest.approx(STYLE.generation_gap)

    def test_ends_on_both_tops(self, layout):
        (arc,) = layout.duplicate_arcs
        (sx, sy), (ex, ey) = arc.waypoints[0], arc.waypoints[-1]
        assert sy == pytest.approx(arc.from_pos[1] + HALF)
        assert ey == pytest.approx(arc.to_pos[1] + HALF)
        assert abs(sx - arc.from_pos[0]) < HALF and abs(ex - arc.to_pos[0]) < HALF

    def test_no_symbol_crossed(self, layout):
        assert _crossed(layout) == 0

    def test_not_up_a_couples_line_of_descent(self, layout):
        """It goes up between w1 and bh, not between a and his wife w1."""
        (arc,) = layout.duplicate_arcs
        x = {n.individual_id: n.x for n in layout.nodes if not n.is_duplicate}
        risers = [p[0] for p, q in _segments(arc.waypoints) if abs(p[0] - q[0]) < 1e-9]
        tall = [
            p[0]
            for p, q in _segments(arc.waypoints)
            if abs(p[0] - q[0]) < 1e-9 and abs(p[1] - q[1]) > STYLE.generation_gap / 2
        ]
        assert risers and tall
        assert x["w1"] < tall[0] < x["bh"]

    def test_drawn_in_both_outputs(self, layout):
        family = parse_pedigree_from_string(ROWS_APART)
        svg = render_interactive_svg(family)
        root = ET.fromstring(svg)
        links = [e for e in root.iter(f"{SVG}polyline") if "duplicate-link" in e.get("class", "")]
        assert len(links) == 1
        fig = render_pedigree(family)
        plt.close(fig)


class TestDrawn:
    def test_svg_draws_the_route(self, generated):
        from pypedigree.models import Pedigree

        family = generated[0]
        ped = Pedigree(families={family.family_id: family})
        root = ET.fromstring(render_interactive_svg(ped))
        dotted = [e for e in root.iter(f"{SVG}polyline") if e.get("class") == "duplicate-link"]
        assert len(dotted) == len(compute_layout(family, STYLE).duplicate_arcs)

    def test_png_draws_the_route(self, generated):
        from pypedigree.models import Pedigree

        family = generated[0]
        fig = render_pedigree(Pedigree(families={family.family_id: family}), show_legend=False)
        try:
            dotted = [ln for ln in fig.axes[0].lines if ln.get_gid() == "duplicate-link"]
            assert len(dotted) == len(compute_layout(family, STYLE).duplicate_arcs)
            assert all(ln.get_linestyle() == ":" for ln in dotted)
        finally:
            plt.close(fig)
