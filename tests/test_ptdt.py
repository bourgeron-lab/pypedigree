"""pTDT deviation columns (#52).

A ``pTDT_`` column holds each child's polygenic transmission deviation,
(PGS_child - PGS_mid-parent) / SD(PGS_mid-parent): drawn with the PGS curve
on the z scale, captioned "pTDT", with its own legend row and tooltip.
All pedigrees here are synthetic.
"""

from __future__ import annotations

import json
import re
import warnings

import matplotlib.pyplot as plt
import pytest
from conftest import PHENOTYPE_DIR
from matplotlib.text import Text

from pypedigree import (
    PedigreeStyle,
    is_pgs_column,
    is_ptdt_column,
    parse_pedigree,
    parse_pedigree_from_string,
    render_interactive_svg,
    render_pedigree,
)
from pypedigree.carriers import CarrierScheme
from pypedigree.config import apply_pgs_spacing, is_curve_column, resolve_style
from pypedigree.errors import PedigreeDataWarning
from pypedigree.phenotypes import (
    active_pgs,
    column_scale,
    curve_caption,
    curve_value_text,
    legend_rows,
    pgs_marker_values,
    pgs_position,
    pgs_trait_name,
    ptdt_legend_description,
    resolve_phenotype_render,
)

FIXTURE = PHENOTYPE_DIR / "ptdt.ped"
STYLE = json.loads((PHENOTYPE_DIR / "ptdt.json").read_text())


def _individual(pedigree, sample):
    return next(i for i in pedigree.individuals() if i.sample_id == sample)


class TestColumns:
    @pytest.mark.parametrize("name", ["pTDT_autism", "PTDT_IQ", "ptdt", "pTDT"])
    def test_ptdt_columns(self, name):
        assert is_ptdt_column(name) and is_curve_column(name)
        assert not is_pgs_column(name)

    def test_pgs_columns_unchanged(self):
        assert is_pgs_column("PGS_autism") and not is_ptdt_column("PGS_autism")
        assert not is_curve_column("SRS")

    def test_always_on_the_z_scale(self):
        for scale in ("minmax", "zscore", "percentile"):
            assert column_scale("pTDT_autism", scale) == "zscore"
            assert column_scale("PGS_autism", scale) == scale

    def test_trait_name(self):
        style = PedigreeStyle()
        assert pgs_trait_name("pTDT_autism", style) == "autism"
        assert pgs_trait_name("pTDT_school_years", style) == "school years"

    def test_drawn_as_a_curve(self):
        pedigree = parse_pedigree(FIXTURE)
        names = pedigree.phenotype_names
        render = resolve_phenotype_render(
            pedigree, names, pedigree.phenotype_display_modes, PedigreeStyle()
        )
        assert render.effective("pTDT_autism") == "pgs"
        assert render.pgs_names(names) == ["PGS_autism", "pTDT_autism"]

    def test_curves_get_their_room(self):
        text = "FID\tIID\tFather\tMother\tSex\tPhenotype\tpTDT_autism\nFAM1\ta\t0\t0\t1\t2\t1.2\n"
        names = parse_pedigree_from_string(text).phenotype_names
        style = PedigreeStyle()
        assert apply_pgs_spacing(style, names).couple_gap > style.couple_gap


class TestCurves:
    def test_parents_show_their_pgs_children_their_ptdt(self):
        pedigree = parse_pedigree(FIXTURE)
        names = ["PGS_autism", "pTDT_autism"]
        parent = active_pgs(_individual(pedigree, "I-1"), names, "zscore")
        child = active_pgs(_individual(pedigree, "II-1"), names, "zscore")
        assert [n for n, _ in parent] == ["PGS_autism"] and curve_caption(parent) == "PGS"
        assert [n for n, _ in child] == ["pTDT_autism"] and curve_caption(child) == "pTDT"

    @pytest.mark.parametrize("scale", ["minmax", "zscore", "percentile"])
    def test_deviation_placed_as_a_z_score(self, scale):
        child = _individual(parse_pedigree(FIXTURE), "II-1")  # +2.1
        [(_name, position)] = active_pgs(child, ["pTDT_autism"], scale)
        assert position == pytest.approx(pgs_position(2.1))

    def test_marker_value(self):
        child = _individual(parse_pedigree(FIXTURE), "II-4")  # -0.2
        assert pgs_marker_values(child, ["pTDT_autism"], "minmax") == {"pTDT_autism": "−0.2"}

    def test_both_kinds_warn_and_show_the_ptdt(self):
        text = (
            "FID\tIID\tFather\tMother\tSex\tPhenotype\tPGS_autism\tpTDT_autism\n"
            "FAM1\ta\t0\t0\t1\t1\t0.3\tNA\n"
            "FAM1\tb\t0\t0\t2\t1\t0.1\tNA\n"
            "FAM1\tc\ta\tb\t1\t2\t1.5\t2.0\n"
        )
        pedigree = parse_pedigree_from_string(text)
        child = _individual(pedigree, "c")
        assert curve_caption(active_pgs(child, ["PGS_autism", "pTDT_autism"])) == "pTDT"
        with pytest.warns(PedigreeDataWarning, match=r"1 individual\(s\) have both.*\(c\)"):
            fig = render_pedigree(pedigree)
        plt.close(fig)

    def test_no_warning_when_kinds_are_apart(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error", PedigreeDataWarning)
            plt.close(render_pedigree(parse_pedigree(FIXTURE), style=PedigreeStyle(**STYLE)))


class TestDrawing:
    def test_png_captions(self):
        fig = render_pedigree(parse_pedigree(FIXTURE), style=PedigreeStyle(**STYLE))
        try:
            texts = [t.get_text() for t in fig.axes[0].findobj(Text)]
        finally:
            plt.close(fig)
        assert texts.count("PGS") == 2 and texts.count("pTDT") == 4

    def test_svg_captions_and_tooltips(self):
        svg = render_interactive_svg(parse_pedigree(FIXTURE), style=PedigreeStyle(**STYLE))
        assert len(re.findall(r">pTDT</text>", svg)) == 4
        assert len(re.findall(r">PGS</text>", svg)) == 2
        assert "pTDT autism: pTDT = +2.10" in svg
        assert "pTDT autism: pTDT = −0.20" in svg
        assert not re.search(r"pTDT = [^|]*percentile", svg)
        assert "PGS autism: z = +0.40, 66th percentile" in svg

    def test_tooltip_text(self):
        assert curve_value_text("pTDT_autism", 2.123, "minmax") == "pTDT = +2.12"
        assert curve_value_text("PGS_autism", 1.1, "zscore") == "z = +1.10, 86th percentile"


class TestLegend:
    def _rows(self, **style):
        pedigree = parse_pedigree(FIXTURE)
        names = pedigree.phenotype_names
        st = resolve_style(PedigreeStyle(**style), names)
        render = resolve_phenotype_render(pedigree, names, pedigree.phenotype_display_modes, st)
        return legend_rows(names, st, render, pedigree.phenotype_display_modes, CarrierScheme())

    def test_own_row(self):
        rows = [r for r in self._rows(pgs_scale="minmax") if r.kind == "pgs"]
        assert [r.label for r in rows] == ["PGS distribution:", ptdt_legend_description()]
        assert rows[1].key == (("autism", rows[1].pgs["pTDT_autism"]),)
        description = ptdt_legend_description()
        assert description.startswith("pTDT deviation")
        assert "mid-parent" in description and description.endswith("/\u00a0SD:")

    def test_row_on_the_z_scale_whatever_pgs_scale(self):
        for scale in ("minmax", "zscore", "percentile"):
            rows = [r for r in self._rows(pgs_scale=scale) if r.kind == "pgs"]
            assert rows[-1].label == ptdt_legend_description()

    def test_png_and_svg_legends(self):
        fig = render_pedigree(parse_pedigree(FIXTURE), style=PedigreeStyle(**STYLE))
        try:
            texts = " ".join(t.get_text() for t in fig.axes[1].texts)
        finally:
            plt.close(fig)
        svg = render_interactive_svg(parse_pedigree(FIXTURE), style=PedigreeStyle(**STYLE))
        legend = svg[svg.index('<g class="legend">') :]
        for text in (texts, legend):
            assert "pTDT deviation" in text and "PGS distribution" in text


class TestColors:
    def test_same_trait_may_share_its_color(self):
        names = ["b_autism", "PGS_autism", "PGS_ADHD", "pTDT_autism", "pTDT_ADHD"]
        style = resolve_style(PedigreeStyle(), names)
        assert style.pgs_colors["pTDT_autism"] == style.pgs_colors["PGS_autism"]
        assert style.pgs_colors["pTDT_ADHD"] == style.pgs_colors["PGS_ADHD"]
        assert style.pgs_colors["pTDT_autism"] != style.pgs_colors["pTDT_ADHD"]

    def test_several_ptdt_columns_get_distinct_colors(self):
        style = resolve_style(PedigreeStyle(), ["pTDT_a", "pTDT_b", "pTDT_c"])
        assert len({style.pgs_colors[n] for n in ("pTDT_a", "pTDT_b", "pTDT_c")}) == 3
