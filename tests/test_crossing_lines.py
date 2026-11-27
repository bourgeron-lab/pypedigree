"""Loops drawn with crossing lines instead of duplicates (`loops` "cross").

A loop that would stretch a sibship line over another is drawn by default
with one partner twice, joined by a dotted line. With "cross" the partner
is drawn once: overlapping sibship lines of a row run at different heights,
and a horizontal line hops over each vertical line it crosses with a small
half-circle. All pedigrees here are synthetic.
"""

from __future__ import annotations

import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib.pyplot as plt
import pytest
from click.testing import CliRunner

from pypedigree import (
    PedigreeStyle,
    compute_layout,
    compute_layout_circular,
    parse_pedigree,
    render_interactive_svg,
    render_pedigree,
)
from pypedigree.cli import main
from pypedigree.errors import PedigreeError
from pypedigree.models import EdgeType
from pypedigree.symbols import hop_path

sys.path.insert(0, str(Path(__file__).parent))
from conftest import TEST_DATA_DIR, generate_realistic_pedigree  # noqa: E402

CROSS = PedigreeStyle(loops="cross")
CROSS_BIRTH = PedigreeStyle(loops="cross", sibling_order="birth")
MARRIED_IN = TEST_DATA_DIR / "complex" / "married_in_founders.ped"
MULTIGENERATIONAL = TEST_DATA_DIR / "complex" / "multigenerational.ped"
SVG = "{http://www.w3.org/2000/svg}"
EPS = 1e-6


def _family(path: Path):
    return next(iter(parse_pedigree(path).families.values()))


@pytest.fixture(scope="module")
def generated():
    return [generate_realistic_pedigree(size, seed) for seed in range(4) for size in (40, 80)]


def _horizontal(e) -> tuple[float, float, float]:
    """(y, left, right) of a horizontal line, hops or not."""
    return e.start[1], min(e.start[0], e.end[0]), max(e.start[0], e.end[0])


def _crossings(layout, e) -> list[float]:
    """x of each vertical line crossing (not ending on) horizontal line ``e``."""
    y, lo, hi = _horizontal(e)
    return sorted(
        v.start[0]
        for v in layout.edges
        if v.edge_type in (EdgeType.DESCENT, EdgeType.INDIVIDUAL)
        and lo + EPS < v.start[0] < hi - EPS
        and min(v.start[1], v.end[1]) + EPS < y < max(v.start[1], v.end[1]) - EPS
    )


def _hops(e) -> int:
    """Half-circles drawn on a horizontal line."""
    y = e.start[1]
    above = [p[1] > y + EPS for p in e.waypoints]
    return sum(1 for a, b in zip([False, *above], above, strict=False) if b and not a)


class TestMarriedInFounders:
    def test_no_copy(self):
        layout = compute_layout(_family(MARRIED_IN), CROSS)
        assert not any(n.is_duplicate for n in layout.nodes)
        assert not layout.duplicate_arcs

    def test_a_hop_over_the_line_it_crosses(self):
        layout = compute_layout(_family(MARRIED_IN), CROSS)
        hopped = [e for e in layout.edges if e.waypoints]
        assert hopped
        for e in hopped:
            assert _hops(e) == len(_crossings(layout, e))

    def test_default_unchanged(self):
        layout = compute_layout(_family(MARRIED_IN))
        assert sum(n.is_duplicate for n in layout.nodes) == 1
        assert not any(e.waypoints for e in layout.edges)


class TestGeometry:
    @pytest.mark.parametrize("style", [CROSS, CROSS_BIRTH], ids=["auto", "birth"])
    def test_sibship_lines_never_overlap(self, generated, style):
        """Rule 3 still holds: sibship lines of one height never share a stretch."""
        for family in [*generated, _family(MULTIGENERATIONAL)]:
            layout = compute_layout(family, style)
            lines = [_horizontal(e) for e in layout.edges if e.edge_type is EdgeType.SIBSHIP]
            for i, (y, lo, hi) in enumerate(lines):
                for y2, lo2, hi2 in lines[i + 1 :]:
                    if abs(y - y2) < EPS:
                        assert hi <= lo2 + EPS or hi2 <= lo + EPS

    @pytest.mark.parametrize("style", [CROSS, CROSS_BIRTH], ids=["auto", "birth"])
    def test_every_crossing_hopped(self, generated, style):
        hops = 0
        for family in [*generated, _family(MULTIGENERATIONAL)]:
            layout = compute_layout(family, style)
            for e in layout.edges:
                if e.edge_type in (EdgeType.SIBSHIP, EdgeType.RELATIONSHIP):
                    assert _hops(e) == len(_crossings(layout, e))
                    hops += _hops(e)
        assert hops > 0

    def test_children_hang_from_their_own_line(self, generated):
        """Each individual line starts on a sibship line, at its height."""
        for family in generated:
            layout = compute_layout(family, CROSS)
            sibships = [_horizontal(e) for e in layout.edges if e.edge_type is EdgeType.SIBSHIP]
            for e in layout.edges:
                if e.edge_type is EdgeType.INDIVIDUAL:
                    x, y = e.start
                    assert any(
                        abs(y - sy) < EPS and lo - EPS <= x <= hi + EPS for sy, lo, hi in sibships
                    ) or any(
                        d.edge_type is EdgeType.DESCENT and abs(d.end[1] - y) < EPS
                        for d in layout.edges
                    )

    def test_fewer_duplicates(self, generated):
        before = after = 0
        for family in generated:
            before += sum(n.is_duplicate for n in compute_layout(family).nodes)
            after += sum(n.is_duplicate for n in compute_layout(family, CROSS).nodes)
        assert after < before

    def test_circular_keeps_duplicates(self, generated):
        plain = compute_layout_circular(generated[0])
        cross = compute_layout_circular(generated[0], CROSS)
        assert sum(n.is_duplicate for n in cross.nodes) == sum(n.is_duplicate for n in plain.nodes)

    def test_hop_path(self):
        points = hop_path(0.0, 4.0, 1.0, [1.0, 3.0], 0.2)
        assert points[0] == (0.0, 1.0) and points[-1] == (4.0, 1.0)
        assert max(p[1] for p in points) == pytest.approx(1.2)
        assert all(p[1] >= 1.0 - EPS for p in points)

    def test_close_hops_shrink(self):
        points = hop_path(0.0, 4.0, 0.0, [1.0, 1.1], 0.2)
        assert max(p[1] for p in points) == pytest.approx(0.05)


class TestDrawn:
    def test_svg_draws_the_hops(self):
        svg = render_interactive_svg(parse_pedigree(MARRIED_IN), style=CROSS)
        root = ET.fromstring(svg)
        lines = [e for e in root.iter(f"{SVG}polyline") if e.get("class") == "edge"]
        assert lines

    def test_png(self):
        fig = render_pedigree(parse_pedigree(MARRIED_IN), style=CROSS)
        plt.close(fig)

    @pytest.mark.parametrize("path", sorted(TEST_DATA_DIR.rglob("*.ped")), ids=lambda p: p.name)
    def test_every_fixture(self, path):
        try:
            pedigree = parse_pedigree(path)
        except Exception:
            pytest.skip("not a drawable pedigree")
        render_interactive_svg(pedigree, style=CROSS_BIRTH)


class TestOption:
    def test_unknown_value(self):
        with pytest.raises(PedigreeError, match="loops"):
            compute_layout(_family(MARRIED_IN), PedigreeStyle(loops="jump"))

    def test_negative_radius(self):
        with pytest.raises(PedigreeError, match="hop_radius"):
            compute_layout(_family(MARRIED_IN), PedigreeStyle(loops="cross", hop_radius=-1))

    def test_cli(self, tmp_path):
        out = tmp_path / "out.svg"
        result = CliRunner().invoke(
            main, [str(MARRIED_IN), "--loops", "cross", "--interactive", "-o", str(out)]
        )
        assert result.exit_code == 0, result.output
        assert 'class="edge" points=' in out.read_text(encoding="utf-8")
