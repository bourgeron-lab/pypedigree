"""Every outline (o_) column an individual is affected in gets a border (#46).

All pedigrees here are synthetic.
"""

from __future__ import annotations

import itertools
import math
import re
import warnings

import matplotlib.pyplot as plt
import pytest
from click.testing import CliRunner
from matplotlib import patches as mpatches
from matplotlib.colors import to_rgb

import pypedigree as pp
from pypedigree import symbols
from pypedigree.cli import main
from pypedigree.config import apply_outline_spacing
from pypedigree.errors import PedigreeStyleWarning
from pypedigree.models import Sex
from pypedigree.phenotypes import (
    LEG_OUTLINE_RINGS_LABEL,
    MAX_QUIET_OUTLINE_RINGS,
    OUTLINE_RING_GAP_PT,
    OUTLINE_WIDTH_SCALE,
    PhenotypeRender,
    outline_colors,
    outline_reach,
    outline_ring_offsets,
    resolve_phenotype_render,
    with_outline_rings,
)

PT = symbols.UNITS_PER_POINT

# The issue's example, with two more children: which data each person has.
DATA = (
    "FID\tIID\tFather\tMother\tSex\tPhenotype\tb_ASD\to_WGS\to_WES\to_Array\n"
    "F\ta\t0\t0\t1\t1\t1\t2\t2\t1\n"
    "F\tb\t0\t0\t2\t1\t1\t1\t2\t2\n"
    "F\tc\ta\tb\t1\t2\t2\t1\t1\t2\n"
    "F\td\ta\tb\t2\t1\t1\t2\t2\t2\n"
    "F\te\ta\tb\t0\t1\t1\t2\t1\t2\n"
)
# The same columns, nobody in more than one.
ONE_EACH = (
    "FID\tIID\tFather\tMother\tSex\tPhenotype\tb_ASD\to_WGS\to_WES\to_Array\n"
    "F\ta\t0\t0\t1\t1\t1\t2\t1\t1\n"
    "F\tb\t0\t0\t2\t1\t1\t1\t2\t1\n"
    "F\tc\ta\tb\t1\t2\t2\t1\t1\t2\n"
)


def _ped(text: str = DATA) -> pp.Pedigree:
    return pp.parse_pedigree_from_string(text)


def _style(ped: pp.Pedigree) -> pp.PedigreeStyle:
    return pp.resolve_style(None, ped.phenotype_names)


def _nodes(svg: str) -> dict[str, str]:
    groups = svg.split('<g class="node"')[1:]
    return {re.search(r'data-id="([^"]+)"', g).group(1): g for g in groups}


def _step(style: pp.PedigreeStyle) -> float:
    return style.shape_line_width * OUTLINE_WIDTH_SCALE + OUTLINE_RING_GAP_PT


class TestColumns:
    def test_colors_in_column_order(self):
        ped = _ped()
        style = _style(ped)
        names = ["WGS", "WES", "Array"]
        d = ped.families["F"].individuals["d"]
        assert outline_colors(d, names, style) == [style.phenotype_colors[n] for n in names]
        b = ped.families["F"].individuals["b"]
        assert outline_colors(b, names, style) == [
            style.phenotype_colors["WES"],
            style.phenotype_colors["Array"],
        ]

    def test_offsets_and_reach(self):
        style = pp.PedigreeStyle()
        assert outline_ring_offsets(1, style) == []
        assert outline_ring_offsets(3, style) == pytest.approx([_step(style), 2 * _step(style)])
        assert outline_reach(0, style) == outline_reach(1, style) == 0.0
        assert outline_reach(3, style) == pytest.approx(2 * _step(style) * PT)

    def test_most_rings(self):
        ped = _ped()
        style = _style(ped)
        render = resolve_phenotype_render(
            ped, ped.phenotype_names, ped.phenotype_display_modes, style
        )
        assert with_outline_rings(render, ped.individuals(), ped.phenotype_names).outline_rings == 3
        assert with_outline_rings(PhenotypeRender(), ped.individuals(), []).outline_rings == 0

    def test_warns_above_three(self):
        header = "FID\tIID\tFather\tMother\tSex\tPhenotype\to_A\to_B\to_C\to_D\n"
        ped = _ped(header + "F\tx\t0\t0\t1\t1\t2\t2\t2\t2\n")
        assert MAX_QUIET_OUTLINE_RINGS == 3
        with pytest.warns(PedigreeStyleWarning, match="x has 4 colored borders"):
            pp.render_svg(ped)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            pp.render_svg(_ped())  # three: quiet


class TestGeometry:
    @pytest.mark.parametrize("sex", [Sex.MALE, Sex.FEMALE])
    def test_square_and_circle_grow_by_the_distance(self, sex):
        assert symbols.grown_half(sex, 0.4, 0.1) == pytest.approx(0.5)

    def test_diamond_sides_move_by_the_distance(self):
        half, d = 0.4, 0.1
        grown = symbols.grown_half(Sex.UNKNOWN, half, d)
        # Distance from the centre to a side: the vertex reach over sqrt(2).
        side = lambda h: h * symbols.DIAMOND_SCALE / math.sqrt(2)  # noqa: E731
        assert side(grown) - side(half) == pytest.approx(d)

    def test_spacing(self):
        style = pp.PedigreeStyle()
        assert apply_outline_spacing(style, 0.0) is style
        wide = apply_outline_spacing(style, 0.2)
        assert wide.couple_gap == pytest.approx(style.couple_gap + 0.4)
        assert wide.sibling_gap == pytest.approx(style.sibling_gap + 0.4)
        assert wide.label_offset == pytest.approx(style.label_offset + 0.2)


class TestSvg:
    def test_one_border_per_column(self):
        ped = _ped()
        style = _style(ped)
        nodes = _nodes(pp.render_interactive_svg(ped))
        rings = {
            k: re.findall(r'class="outline-ring"[^>]*stroke="([^"]+)"', g) for k, g in nodes.items()
        }
        assert rings == {
            "a": [style.phenotype_colors["WES"]],
            "b": [style.phenotype_colors["Array"]],
            "c": [],
            "d": [style.phenotype_colors["WES"], style.phenotype_colors["Array"]],
            "e": [style.phenotype_colors["Array"]],
        }
        # The first border stays on the outline, in the first column's color.
        shape = re.search(r'class="shape"[^>]*stroke="([^"]+)"', nodes["d"]).group(1)
        assert shape == style.phenotype_colors["WGS"]

    def test_borders_a_step_apart(self):
        ped = _ped()
        style = _style(ped)
        group = _nodes(pp.render_interactive_svg(ped))["d"]
        shape_r = float(re.search(r'<circle class="shape"[^>]*r="([0-9.]+)"', group).group(1))
        radii = [
            float(r) for r in re.findall(r'<circle class="outline-ring"[^>]*r="([0-9.]+)"', group)
        ]
        step = _step(style) * PT
        assert radii == pytest.approx([shape_r + step, shape_r + 2 * step], abs=2e-3)

    def test_neighbours_and_labels_clear(self):
        ped = _ped()
        style = _style(ped)
        svg = pp.render_interactive_svg(ped)
        nodes = _nodes(svg)
        reach = outline_reach(3, style) + style.shape_line_width * OUTLINE_WIDTH_SCALE * PT / 2
        half = style.shape_size / 2 * symbols.DIAMOND_SCALE  # the widest shape
        xy = {
            k: (
                float(re.search(r'data-x="([-0-9.]+)"', g).group(1)),
                float(re.search(r'data-y="([-0-9.]+)"', g).group(1)),
            )
            for k, g in nodes.items()
        }
        for (_, (x1, y1)), (_, (x2, y2)) in itertools.combinations(xy.items(), 2):
            if y1 == y2:
                assert abs(x1 - x2) > 2 * (half + reach)
        label_y = {}
        for key, group in nodes.items():
            y = float(
                re.search(rf'<text x="[-0-9.]+" y="([-0-9.]+)"[^>]*>{key}</text>', group).group(1)
            )
            label_y[key] = y - xy[key][1]
        # Every label at the same distance under its symbol, below the outer border.
        assert len({round(v, 6) for v in label_y.values()}) == 1
        assert label_y["c"] > style.shape_size / 2 + reach

    def test_tooltip_names_the_borders(self):
        nodes = _nodes(pp.render_interactive_svg(_ped()))
        assert "Borders: WGS, WES, Array" in nodes["d"]
        assert "Borders: WES, Array" in nodes["b"]
        assert "Borders:" not in nodes["c"]

    def test_legend_row(self):
        svg = pp.render_svg(_ped())
        assert LEG_OUTLINE_RINGS_LABEL in svg[svg.index('class="legend"') :]

    def test_pgs_box_clear_of_the_borders(self):
        text = DATA.replace("o_Array\n", "o_Array\tPGS\n")
        rows = text.splitlines()
        text = "\n".join([rows[0]] + [r + "\t0.5" for r in rows[1:]]) + "\n"
        ped = _ped(text)
        style = _style(ped)
        nodes = _nodes(pp.render_interactive_svg(ped))
        group = nodes["d"]
        x = float(re.search(r'data-x="([-0-9.]+)"', group).group(1))
        box = float(re.search(r'<g class="pgs">\s*<rect x="([-0-9.]+)"', group).group(1))
        half = style.shape_size / 2
        assert box == pytest.approx(x - half * 1.95 - outline_reach(3, style), abs=2e-3)

    def test_nothing_changes_with_one_border_each(self):
        ped = _ped(ONE_EACH)
        svg = pp.render_interactive_svg(ped)
        assert "outline-ring" not in svg
        assert "Borders:" not in svg
        assert LEG_OUTLINE_RINGS_LABEL not in svg
        # Same layout as without the o_ columns' borders at all.
        plain = ONE_EACH.replace("\t2\t1\t1\n", "\t1\t1\t1\n").replace("\t1\t2\t1\n", "\t1\t1\t1\n")
        xs = lambda s: re.findall(r'data-x="([-0-9.]+)"', s)  # noqa: E731
        assert xs(svg) == xs(pp.render_interactive_svg(_ped(plain)))


class TestPng:
    def _rings(self, fig) -> list[mpatches.Patch]:
        from pypedigree.renderer import _OutlineRings

        artists = [a for a in fig.axes[0].get_children() if isinstance(a, _OutlineRings)]
        return [p for a in artists for p in a.get_children()]

    def test_borders_drawn(self):
        ped = _ped()
        style = _style(ped)
        fig = pp.render_pedigree(ped)
        fig.canvas.draw()
        colors = sorted(tuple(p.get_edgecolor()[:3]) for p in self._rings(fig))
        plt.close(fig)
        expected = sorted(
            to_rgb(style.phenotype_colors[n]) for n in ("WES", "Array", "WES", "Array", "Array")
        )
        assert colors == pytest.approx(expected)

    def test_borders_a_step_apart_in_points(self):
        ped = _ped()
        style = _style(ped)
        fig = pp.render_pedigree(ped)
        fig.canvas.draw()
        ax = fig.axes[0]
        circles = sorted(
            (p for p in self._rings(fig) if isinstance(p, mpatches.Circle)), key=lambda p: p.radius
        )
        assert len(circles) == 3  # b's one, d's two
        # d's two share a centre; b's one is alone.
        d_rings = sorted(
            (c for c in circles if sum(1 for o in circles if o.center == c.center) == 2),
            key=lambda c: c.radius,
        )
        origin, unit = ax.transData.transform([(0.0, 0.0), (1.0, 0.0)])
        pt_per_unit = abs(unit[0] - origin[0]) / (fig.dpi / 72)
        half = style.shape_size / 2
        gaps = [(c.radius - half) * pt_per_unit for c in d_rings]
        plt.close(fig)
        assert gaps == pytest.approx([_step(style), 2 * _step(style)], rel=1e-3)

    def test_legend_row(self):
        fig = pp.render_pedigree(_ped())
        texts = {t.get_text() for t in fig.axes[1].texts}
        plt.close(fig)
        assert LEG_OUTLINE_RINGS_LABEL in texts

    def test_circular_turns_with_the_node(self):
        fig = pp.render_pedigree(_ped(), circular=True)
        fig.canvas.draw()
        assert self._rings(fig)
        plt.close(fig)


def test_cli_warning(tmp_path):
    header = "FID\tIID\tFather\tMother\tSex\tPhenotype\to_A\to_B\to_C\to_D\n"
    ped = tmp_path / "four.ped"
    ped.write_text(header + "F\tx\t0\t0\t1\t1\t2\t2\t2\t2\n", encoding="utf-8")
    out = tmp_path / "f.svg"
    result = CliRunner().invoke(main, [str(ped), "--interactive", "-o", str(out)])
    assert result.exit_code == 0, result.output
    assert "Warning: x has 4 colored borders" in result.output
