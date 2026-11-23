"""`PedigreeStyle.text_scale`: a pedigree drawn smaller than its text (#145).

The layout keeps room for what is sized in points (labels, markers, the
evaluation asterisk, PGS boxes and values) at `symbols.UNITS_PER_POINT`
layout units per point. A pedigree drawn into one's own smaller axes needs
that room at its real scale: `text_scale` multiplies it, the SVG draws its
text and lines at it, and `render_pedigree` keeps its own figures no
smaller. 1.0 changes nothing.

All pedigrees here are synthetic.
"""

from __future__ import annotations

import math
import xml.etree.ElementTree as ET

import matplotlib.pyplot as plt
import pytest

from pypedigree import (
    PedigreeError,
    PedigreeStyle,
    compute_layout,
    parse_pedigree_from_string,
    render_interactive_svg,
    render_pedigree,
    symbols,
)
from pypedigree.models import Sex

SVG = "{http://www.w3.org/2000/svg}"
PAPER = 2.8  # the paper's Fig. 1a: drawn at about 0.056 units per point
HEAD = "FID\tID\tFather\tMother\tSex\tConsultand\tevaluated\n"
# An evaluated son left of his consultand sister, as IV-1 and IV-2 of Fig. 1a.
SIBS = HEAD + (
    "F\tf\t0\t0\t1\t0\t0\nF\tm\t0\t0\t2\t0\t0\n"
    "F\ta\tf\tm\t1\t0\t1\nF\tb\tf\tm\t2\t1\t0\nF\tc\tf\tm\t1\t0\t0\n"
)


def _family(text: str = SIBS):
    return parse_pedigree_from_string(text).families["F"]


def _xs(style: PedigreeStyle, text: str = SIBS) -> dict[str, float]:
    return {n.individual_id: n.x for n in compute_layout(_family(text), style).nodes}


class TestDefault:
    def test_one_changes_nothing(self):
        ped = parse_pedigree_from_string(SIBS)
        assert _xs(PedigreeStyle(text_scale=1.0)) == _xs(PedigreeStyle())
        assert render_interactive_svg(ped, style=PedigreeStyle(text_scale=1.0)) == (
            render_interactive_svg(ped, style=PedigreeStyle())
        )

    def test_units_per_point(self):
        assert symbols.units_per_point(1.0) == symbols.UNITS_PER_POINT
        assert symbols.units_per_point(2.0) == pytest.approx(2 * symbols.UNITS_PER_POINT)


class TestCheck:
    @pytest.mark.parametrize("value", [0, -1.0, math.inf, math.nan, "2", True, None])
    def test_rejected(self, value):
        style = PedigreeStyle(text_scale=value)
        with pytest.raises(PedigreeError, match="text_scale"):
            compute_layout(_family(), style)
        with pytest.raises(PedigreeError, match="text_scale"):
            render_interactive_svg(parse_pedigree_from_string(SIBS), style=style)

    def test_int_accepted(self):
        compute_layout(_family(), PedigreeStyle(text_scale=2))


class TestLayout:
    def test_asterisk_and_arrow_planned_at_the_scale(self):
        """At 2.8 times the text scale, the consultand's left neighbour keeps the
        asterisk and the whole arrow apart, measured at that scale."""
        style = PedigreeStyle(text_scale=PAPER)
        upp = symbols.units_per_point(PAPER)
        x = _xs(style)
        size = style.text_size(symbols.EVALUATION_MARK_FONT_SCALE)
        mark = symbols.evaluation_mark_reach(
            style.shape_size / 2, Sex.MALE, size, style.shape_line_width, upp=upp
        )
        mark += symbols.EVALUATION_MARK_GAP_EM * size * upp
        tail = style.shape_size / 2 + style.text_size() * upp * symbols.MARKER_ARROW_EM
        assert x["b"] - x["a"] >= mark + tail - 1e-9
        # The default scale needs no more than the sibling gap.
        assert _xs(PedigreeStyle())["b"] - _xs(PedigreeStyle())["a"] == pytest.approx(1.5)

    def test_labels_planned_at_the_scale(self):
        """Long labels keep apart at the scale they will be drawn at."""
        text = "FID\tID\tFather\tMother\tSex\tPhenotype\n" + "".join(
            f"F\t{name}\t0\t0\t1\t1\n" for name in ("long_label_one", "long_label_two")
        )
        rows = text.splitlines()
        rows += ["F\tkid_a\tlong_label_one\tmum\t1\t1", "F\tmum\t0\t0\t2\t1"]
        text = "\n".join(rows) + "\n"
        style = PedigreeStyle(text_scale=PAPER)
        upp = symbols.units_per_point(PAPER)
        nodes = compute_layout(parse_pedigree_from_string(text).families["F"], style).nodes
        xs = sorted(n.x for n in nodes if n.individual_id in ("long_label_one", "long_label_two"))
        width = symbols.label_width("long_label_one", [], style.text_size(), upp=upp)
        assert xs[1] - xs[0] >= width - 1e-9


class TestSvg:
    def test_text_and_lines_at_the_scale(self):
        ped = parse_pedigree_from_string(SIBS)
        root = ET.fromstring(render_interactive_svg(ped, style=PedigreeStyle(text_scale=2.0)))
        upp = symbols.units_per_point(2.0)
        label = next(t for t in root.iter(f"{SVG}text") if t.text == "a" and t.get("font-size"))
        assert float(label.get("font-size")) == pytest.approx(8.0 * upp, abs=1e-3)
        edges = [e for e in root.iter(f"{SVG}line") if e.get("class") == "edge"]
        assert edges
        for e in edges:
            assert float(e.get("stroke-width")) == pytest.approx(1.0 * upp, abs=1e-3)


def _paper_axes(ped, style: PedigreeStyle):
    """Draw into axes sized so the pedigree is at 1 / (0.02 x PAPER) points per
    unit, as the paper's Fig. 1a is."""
    layout = compute_layout(next(iter(ped.families.values())), style)
    xs = [n.x for n in layout.nodes]
    ys = [n.y for n in layout.nodes]
    pad = 1.5
    width_u = max(xs) - min(xs) + 2 * pad
    height_u = max(ys) - min(ys) + 2 * pad
    pt_per_unit = 1 / symbols.units_per_point(PAPER)
    fig = plt.figure(figsize=(width_u * pt_per_unit / 72, height_u * pt_per_unit / 72), dpi=300)
    ax = fig.add_axes((0, 0, 1, 1))
    render_pedigree(ped, style=style, ax=ax, show_legend=False)
    ax.set_xlim(min(xs) - pad, max(xs) + pad)
    ax.set_ylim(min(ys) - pad, max(ys) + pad)
    return fig, ax


class TestPng:
    @staticmethod
    def _star_and_tail(fig, ax):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        (star,) = [t for t in ax.texts if t.get_text() == "*"]
        (arrow,) = [t for t in ax.texts if t.get_text() == "" and t.arrow_patch is not None]
        return (
            star.get_window_extent(renderer).x1,
            arrow.arrow_patch.get_window_extent(renderer).x0,
        )

    def test_paper_scale_asterisk_clear_of_the_arrow(self):
        """Drawn into small axes at the paper's scale: without the option the
        asterisk runs into the arrow; with it, the layout made room."""
        ped = parse_pedigree_from_string(SIBS)
        for scale, clear in ((1.0, False), (PAPER, True)):
            fig, ax = _paper_axes(ped, PedigreeStyle(text_scale=scale))
            try:
                star_right, tail_left = self._star_and_tail(fig, ax)
                assert bool(star_right < tail_left) is clear, scale
            finally:
                plt.close(fig)

    def test_own_figure_kept_at_the_scale(self):
        """`render_pedigree`'s own figure grows to the text scale: less with a
        larger text_scale, which allows a smaller drawing."""
        rows = ["FID\tID\tFather\tMother\tSex\tPhenotype", "F\tf\t0\t0\t1\t1", "F\tm\t0\t0\t2\t1"]
        rows += [f"F\tc{i:02d}\tf\tm\t{1 + i % 2}\t1" for i in range(40)]
        ped = parse_pedigree_from_string("\n".join(rows) + "\n")
        sizes = {}
        for scale in (1.0, 2.0):
            fig = render_pedigree(ped, style=PedigreeStyle(text_scale=scale), show_legend=False)
            sizes[scale] = fig.get_size_inches()[0]
            plt.close(fig)
        assert sizes[2.0] < sizes[1.0]


class TestCli:
    def test_text_scale_option(self, tmp_path):
        from click.testing import CliRunner

        from pypedigree.cli import main

        src = tmp_path / "sibs.ped"
        src.write_text(SIBS, encoding="utf-8")
        out = tmp_path / "out.svg"
        result = CliRunner().invoke(
            main, [str(src), "-o", str(out), "--interactive", "--text-scale", "2"]
        )
        assert result.exit_code == 0, result.output
        root = ET.fromstring(out.read_text(encoding="utf-8"))
        label = next(t for t in root.iter(f"{SVG}text") if t.text == "a" and t.get("font-size"))
        assert float(label.get("font-size")) == pytest.approx(
            8.0 * symbols.units_per_point(2.0), abs=1e-3
        )
        bad = CliRunner().invoke(main, [str(src), "-o", str(out), "--text-scale", "0"])
        assert bad.exit_code != 0
