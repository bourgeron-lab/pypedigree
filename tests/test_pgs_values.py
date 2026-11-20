"""PGS values in the figure, and markers that no longer hide each other (#54).

All pedigrees here are synthetic.
"""

from __future__ import annotations

import itertools
import math
import re

import matplotlib.pyplot as plt
import pytest
from click.testing import CliRunner
from matplotlib.lines import Line2D

import pypedigree as pp
from pypedigree import symbols
from pypedigree.cli import main
from pypedigree.config import apply_pgs_spacing
from pypedigree.layout import compute_layout
from pypedigree.phenotypes import (
    MINUS,
    PGS_DOT_PT,
    PGS_DOT_SPREAD,
    pgs_marker_text,
    pgs_values_width,
    spread_dots,
    spread_labels,
)

PT = symbols.UNITS_PER_POINT

# Three scores; the siblings' are within a fifth of a standard deviation of
# each other, as in the issue.
CLOSE = (
    "FID\tIID\tPAT\tMAT\tSEX\tPHENO\tb_autism\tPGS_autism\tPGS_ADHD\tPGS_intelligence\n"
    "F1\tfa\t0\t0\t1\t1\t1\t0.4\t-1.2\t0.8\n"
    "F1\tmo\t0\t0\t2\t1\t1\t-0.6\t0.3\t-0.5\n"
    "F1\ts1\tfa\tmo\t1\t2\t2\t-0.14\t0.14\t-0.21\n"
    "F1\ts2\tfa\tmo\t2\t1\t1\t0.21\t-0.14\t2.35\n"
)
ZSCORE = pp.PedigreeStyle(pgs_scale="zscore", show_pgs_values=True)


def _ped() -> pp.Pedigree:
    return pp.parse_pedigree_from_string(CLOSE)


def _nodes(svg: str) -> dict[str, str]:
    groups = svg.split('<g class="node"')[1:]
    return {re.search(r'data-id="([^"]+)"', g).group(1): g for g in groups}


class TestText:
    @pytest.mark.parametrize(
        ("raw", "scale", "text"),
        [
            (2.12, "zscore", "+2.1"),
            (-0.4, "zscore", f"{MINUS}0.4"),
            (0.0, "zscore", "+0.0"),
            (97.6, "percentile", "98"),
            (0.71234, "minmax", "0.712"),
            (-1.5, "minmax", f"{MINUS}1.5"),
            (1234.0, "minmax", "1.23e+03"),
        ],
    )
    def test_format(self, raw, scale, text):
        assert pgs_marker_text(raw, scale) == text

    def test_minus_is_measured(self):
        """The minus has its DejaVu Sans width, not the fallback's."""
        assert symbols._CHAR_EM[MINUS] == pytest.approx(0.838)


class TestSpreadDots:
    def test_apart_dots_stay(self):
        assert spread_dots([0.0, 1.0, 2.0], 0.1) == [0.0, 0.0, 0.0]

    def test_close_dots_side_by_side(self):
        offsets = spread_dots([0.5, 0.52, 0.0], 0.1)
        step = 0.1 * PGS_DOT_SPREAD
        assert offsets == pytest.approx([-step / 2, step / 2, 0.0])

    def test_three_close_dots_in_input_order(self):
        offsets = spread_dots([0.53, 0.5, 0.56], 0.1)
        step = 0.1 * PGS_DOT_SPREAD
        assert offsets == pytest.approx([0.0, -step, step])


class TestSpreadLabels:
    def test_apart_labels_stay(self):
        assert spread_labels([0.1, 0.5, 0.9], 0.1, 0.0, 1.0) == pytest.approx([0.1, 0.5, 0.9])

    def test_close_labels_spread_around_their_mean(self):
        ys = spread_labels([0.5, 0.52, 0.48], 0.1, 0.0, 1.0)
        assert sorted(ys) == pytest.approx([0.4, 0.5, 0.6])
        assert ys[2] < ys[0] < ys[1]  # order kept

    def test_kept_inside(self):
        ys = spread_labels([0.99, 0.98, 0.97], 0.1, 0.0, 1.0)
        assert max(ys) <= 1.0 + 1e-9
        assert min(b - a for a, b in itertools.pairwise(sorted(ys))) >= 0.1 - 1e-9

    def test_chains_of_groups_merge(self):
        ys = spread_labels([0.0, 0.08, 0.16, 0.24], 0.1, -1.0, 1.0)
        assert min(b - a for a, b in itertools.pairwise(sorted(ys))) >= 0.1 - 1e-9


class TestSpacing:
    def test_no_room_without_the_option(self):
        ped = _ped()
        style = pp.PedigreeStyle(pgs_scale="zscore")
        assert pgs_values_width(ped.individuals(), ["PGS_autism"], style) == 0.0

    def test_room_for_the_widest_value(self):
        ped = _ped()
        names = ["PGS_autism", "PGS_ADHD", "PGS_intelligence"]
        width = pgs_values_width(ped.individuals(), names, ZSCORE)
        size = ZSCORE.text_size(symbols.PGS_CAPTION_SCALE)
        assert width > symbols.text_width(f"{MINUS}1.2", size)

    def test_symbols_move_apart(self):
        ped = _ped()
        family = next(iter(ped.families.values()))
        names = ped.phenotype_names

        def gap(style: pp.PedigreeStyle) -> float:
            width = pgs_values_width(ped.individuals(), names[1:], style)
            spaced = apply_pgs_spacing(style, names, width)
            xs = {n.individual_id: n.x for n in compute_layout(family, spaced).nodes}
            return xs["s2"] - xs["s1"]

        assert gap(ZSCORE) > gap(pp.PedigreeStyle(pgs_scale="zscore"))


class TestSvg:
    def test_values_written_in_their_colors(self):
        svg = pp.render_svg(_ped(), style=ZSCORE)
        style = pp.resolve_style(ZSCORE, _ped().phenotype_names)
        values = re.findall(r'<text [^>]*fill="([^"]+)"[^>]*class="pgs-value"[^>]*>([^<]+)<', svg)
        assert len(values) == 12
        assert (style.pgs_colors["PGS_intelligence"], "+2.4") in values
        assert (style.pgs_colors["PGS_ADHD"], f"{MINUS}1.2") in values

    def test_no_values_by_default(self):
        svg = pp.render_svg(_ped(), style=pp.PedigreeStyle(pgs_scale="zscore"))
        assert 'class="pgs-value"' not in svg

    def test_values_beyond_the_box_and_apart(self):
        """On the far side of the box from the symbol: left of a box on the
        left, right of a right-hand partner's box on the right (#78)."""
        svg = pp.render_svg(_ped(), style=ZSCORE)
        size = ZSCORE.text_size(symbols.PGS_CAPTION_SCALE) * PT
        for group in _nodes(svg).values():
            box = re.search(
                r'<g class="pgs">\s*<rect x="([-0-9.]+)" y="[-0-9.]+" width="([0-9.]+)"', group
            )
            box_left, box_w = float(box.group(1)), float(box.group(2))
            node_x = float(re.search(r'data-x="([-0-9.]+)"', group).group(1))
            values = re.findall(r'<text x="([-0-9.]+)" y="([-0-9.]+)"[^>]*class="pgs-value"', group)
            assert values
            if box_left > node_x:  # a box on the right
                assert all(float(x) > box_left + box_w for x, _y in values)
            else:
                assert all(float(x) < box_left for x, _y in values)
            ys = sorted(float(y) for _x, y in values)
            assert min(b - a for a, b in itertools.pairwise(ys)) >= 1.15 * size - 2e-3

    def test_close_dots_do_not_overlap(self):
        svg = pp.render_svg(_ped(), style=pp.PedigreeStyle(pgs_scale="zscore"))
        group = _nodes(svg)["s1"]
        dots = [
            (float(x), float(y))
            for x, y in re.findall(r'<circle cx="([-0-9.]+)" cy="([-0-9.]+)" r="0.030"', group)
        ]
        assert len(dots) == 3
        for a, b in itertools.combinations(dots, 2):
            assert math.dist(a, b) >= PGS_DOT_PT * PT - 2e-3

    def test_viewbox_holds_the_values(self):
        svg = pp.render_svg(_ped(), style=ZSCORE)
        vb_x = float(re.search(r'viewBox="([-0-9.]+)', svg).group(1))
        xs = [float(x) for x in re.findall(r'<text x="([-0-9.]+)"[^>]*class="pgs-value"', svg)]
        size = ZSCORE.text_size(symbols.PGS_CAPTION_SCALE)
        assert min(xs) - symbols.text_width(f"{MINUS}1.2", size) >= vb_x - 1e-3


class TestPng:
    def _drawn(self, style: pp.PedigreeStyle):
        fig = pp.render_pedigree(_ped(), style=style)
        fig.canvas.draw()
        return fig

    def test_values_written(self):
        fig = self._drawn(ZSCORE)
        texts = {t.get_text() for t in fig.findobj(plt.Text)}
        plt.close(fig)
        assert {"+2.4", f"{MINUS}1.2", "+0.1"} <= texts

    def test_close_dots_do_not_overlap(self):
        fig = self._drawn(pp.PedigreeStyle(pgs_scale="zscore"))
        ax = fig.axes[0]
        dots = [
            line
            for line in ax.findobj(Line2D)  # not the legend's symbols
            if line.get_marker() == "o" and len(line.get_xdata()) == 1
        ]
        points = [ax.transData.transform((d.get_xdata()[0], d.get_ydata()[0])) for d in dots]
        px_per_pt = fig.dpi / 72
        close = [
            (a, b)
            for a, b in itertools.combinations(points, 2)
            if math.dist(a, b) < PGS_DOT_PT * px_per_pt - 0.5
        ]
        plt.close(fig)
        assert len(dots) == 12
        assert close == []

    def test_no_values_by_default(self):
        fig = self._drawn(pp.PedigreeStyle(pgs_scale="zscore"))
        texts = {t.get_text() for t in fig.findobj(plt.Text)}
        plt.close(fig)
        assert "+2.4" not in texts


class TestOptions:
    def test_cli(self, tmp_path):
        ped = tmp_path / "close.ped"
        ped.write_text(CLOSE, encoding="utf-8")
        out = tmp_path / "f.svg"
        args = [str(ped), "--interactive", "--pgs-scale", "zscore", "--show-pgs-values"]
        result = CliRunner().invoke(main, [*args, "-o", str(out)])
        assert result.exit_code == 0, result.output
        assert 'class="pgs-value"' in out.read_text(encoding="utf-8")

    def test_min_font_size(self):
        style = pp.PedigreeStyle(
            pgs_scale="zscore", show_pgs_values=True, font_size=6.0, min_font_size=6.0
        )
        svg = pp.render_svg(_ped(), style=style)
        sizes = re.findall(r'font-size="([0-9.]+)"[^>]*class="pgs-value"', svg)
        assert sizes
        assert min(float(v) for v in sizes) >= 6.0 * PT - 1e-9
