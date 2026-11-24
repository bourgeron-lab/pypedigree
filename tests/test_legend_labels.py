"""Names shown for columns, a compact legend with quadrants, the PGS key (#56).

All pedigrees here are synthetic.
"""

from __future__ import annotations

import re
from dataclasses import replace

import matplotlib.pyplot as plt
import pytest
from matplotlib.colors import to_hex

import pypedigree as pp
from pypedigree.carriers import CarrierScheme
from pypedigree.errors import PedigreeError
from pypedigree.phenotypes import (
    LEG_SHARED_BINARY_LABEL,
    display_name,
    legend_key_items,
    legend_metrics,
    legend_rows,
    resolve_phenotype_render,
    with_quadrants,
)

# Four binary columns (quadrants), one with underscores; a score; a gene.
FOUR = (
    "FID\tIID\tPAT\tMAT\tSEX\tPHENO\tb_autism\tb_intellectual_disability\tb_ADHD\tb_epilepsy\n"
    "F1\tfa\t0\t0\t1\t1\t1\t1\t2\t1\n"
    "F1\tmo\t0\t0\t2\t1\t1\t2\t1\t1\n"
    "F1\tch\tfa\tmo\t1\t2\t2\t2\t1\t2\n"
)
SCORES = (
    "FID\tIID\tPAT\tMAT\tSEX\tPHENO\tb_autism\tg_SRS\tPGS_autism\tPGS_school_years\tcr_SHANK3\n"
    "F1\tfa\t0\t0\t1\t1\t1\t52\t0.2\t0.5\t1\n"
    "F1\tmo\t0\t0\t2\t1\t1\t58\t-0.4\t-1.0\t0\n"
    "F1\tch\tfa\tmo\t1\t2\t2\t78\t1.8\t0.3\t1\n"
)


def _rows(text: str, style: pp.PedigreeStyle | None = None):
    ped = pp.parse_pedigree_from_string(text)
    names, modes = ped.phenotype_names, ped.phenotype_display_modes
    style = pp.resolve_style(style, names)
    render = resolve_phenotype_render(ped, names, modes, style)
    render = with_quadrants(render, render.band_names(names), CarrierScheme(), style)
    return legend_rows(names, style, render, modes, CarrierScheme()), style


class TestDisplayName:
    def test_underscores_become_spaces(self):
        assert display_name("intellectual_disability", pp.PedigreeStyle()) == (
            "intellectual disability"
        )

    def test_label_wins(self):
        style = pp.PedigreeStyle(phenotype_labels={"SRS": "SRS T-score"})
        assert display_name("SRS", style) == "SRS T-score"
        assert display_name("IQ", style) == "IQ"

    def test_legend_and_key_use_it(self):
        rows, _ = _rows(FOUR, pp.PedigreeStyle(legend="full"))
        assert "intellectual disability" in {c.label for c in rows[0].quadrants}
        assert "intellectual disability" in {r.label for r in rows}

    def test_band_text_and_tooltip(self):
        ped = pp.parse_pedigree_from_string(SCORES)
        style = pp.PedigreeStyle(
            phenotype_layout="bands",
            band_text="name_value",
            phenotype_labels={"SRS": "SRS T-score", "PGS_school_years": "EA"},
            pgs_scale="minmax",  # the raw value in the tooltip
        )
        svg = pp.render_interactive_svg(ped, style=style)
        assert ">SRS T-score: 78</text>" in svg
        assert "SRS T-score: 78 |" in svg or "SRS T-score: 78<" in svg  # tooltip
        assert "EA: 0.3" in svg
        assert "PGS school years" not in svg

    def test_carrier_label(self):
        ped = pp.parse_pedigree_from_string(SCORES)
        style = pp.PedigreeStyle(phenotype_labels={"SHANK3": "SHANK3 (paternal)"})
        svg = pp.render_interactive_svg(ped, style=style)
        assert "Carrier: SHANK3 (paternal)" in svg

    def test_invalid(self):
        ped = pp.parse_pedigree_from_string(FOUR)
        for style in (
            pp.PedigreeStyle(legend="short"),
            pp.PedigreeStyle(phenotype_labels={"autism": 3}),
        ):
            with pytest.raises(PedigreeError):
                pp.render_svg(ped, style=style)


class TestCompactLegend:
    def test_quadrants_share_one_row(self):
        rows, style = _rows(FOUR)
        kinds = [r.kind for r in rows]
        assert kinds == ["quadrants", "binary_shared"]
        shared = rows[1]
        assert shared.label == LEG_SHARED_BINARY_LABEL
        expected = [
            style.phenotype_colors[n]
            for n in ("autism", "intellectual_disability", "ADHD", "epilepsy")
        ]
        assert [to_hex(c) for c in shared.colors] == [c.lower() for c in expected]

    def test_full(self):
        rows, _ = _rows(FOUR, pp.PedigreeStyle(legend="full"))
        assert [r.kind for r in rows] == ["quadrants"] + ["binary"] * 4

    def test_bands_keep_one_row_per_column(self):
        for legend in ("auto", "compact"):
            rows, _ = _rows(FOUR, pp.PedigreeStyle(phenotype_layout="bands", legend=legend))
            assert [r.kind for r in rows] == ["binary"] * 4

    def test_other_rows_kept(self):
        rows, _ = _rows(SCORES)
        assert [r.kind for r in rows] == ["quadrants", "binary", "gradient", "pgs"]

    def test_drawn_in_both_renderers(self):
        ped = pp.parse_pedigree_from_string(FOUR)
        style = pp.resolve_style(None, ped.phenotype_names)
        colors = [style.phenotype_colors[n] for n in ped.phenotype_names]
        svg = pp.render_svg(ped)
        legend = svg[svg.index('class="legend"') :]
        for color in colors:
            assert f'fill="{color.lower()}"' in legend.lower()
        assert LEG_SHARED_BINARY_LABEL in legend
        fig = pp.render_pedigree(ped)
        texts = {t.get_text() for t in fig.axes[1].texts}
        fills = {to_hex(p.get_facecolor()) for p in fig.axes[1].patches}
        plt.close(fig)
        assert LEG_SHARED_BINARY_LABEL in texts
        assert {c.lower() for c in colors} <= fills

    def test_shorter(self):
        ped = pp.parse_pedigree_from_string(FOUR)
        compact = pp.render_pedigree(ped)
        full = pp.render_pedigree(ped, style=pp.PedigreeStyle(legend="full"))
        heights = compact.get_size_inches()[1], full.get_size_inches()[1]
        plt.close(compact)
        plt.close(full)
        assert heights[0] < heights[1]


class TestPgsKey:
    def test_key_rows(self):
        rows, style = _rows(SCORES)
        pgs = next(r for r in rows if r.kind == "pgs")
        assert pgs.label == "PGS distribution, z-score:"
        assert [name for name, _ in pgs.key] == ["autism", "school years"]
        assert [color for _, color in pgs.key] == [
            style.pgs_colors["PGS_autism"],
            style.pgs_colors["PGS_school_years"],
        ]

    def test_items_follow_each_other(self):
        rows, style = _rows(SCORES)
        pgs = next(r for r in rows if r.kind == "pgs")
        m = legend_metrics(len(rows), style)
        items = legend_key_items(pgs, m.icon_text_x, m)
        assert items[0].dash_x0 > m.icon_text_x
        for item in items:
            assert item.dash_x0 < item.dash_x1 < item.text_x
        assert items[1].dash_x0 > items[0].text_x

    def test_svg_dash_in_each_color_before_its_name(self):
        ped = pp.parse_pedigree_from_string(SCORES)
        svg = pp.render_svg(ped)
        style = pp.resolve_style(None, ped.phenotype_names)
        legend = svg[svg.index('class="legend"') :]
        for name, column in (("autism", "PGS_autism"), ("school years", "PGS_school_years")):
            color = style.pgs_colors[column]
            dash = re.search(
                rf'<line class="pgs-key-line" x1="[-0-9.]+" [^>]*x2="([-0-9.]+)"[^>]*'
                rf'stroke="{color}"',
                legend,
            )
            text = re.search(rf'<text x="([-0-9.]+)"[^>]*>{name}</text>', legend)
            assert dash and text
            assert float(dash.group(1)) < float(text.group(1))

    def test_png_key_lines(self):
        ped = pp.parse_pedigree_from_string(SCORES)
        style = pp.resolve_style(None, ped.phenotype_names)
        fig = pp.render_pedigree(ped)
        legend = fig.axes[1]
        keys = {to_hex(ln.get_color()) for ln in legend.lines if ln.get_gid() == "pgs-key-line"}
        texts = {t.get_text() for t in legend.texts}
        plt.close(fig)
        assert {style.pgs_colors["PGS_autism"].lower()} <= keys
        assert "school years" in texts

    def test_single_pgs_column_keyed_by_its_name(self):
        """A lone "PGS" column is keyed like the others (#159): the icon has no
        marker, so the key gives the marker's color."""
        text = (
            "FID\tIID\tPAT\tMAT\tSEX\tPHENO\tPGS\n"
            "F1\tfa\t0\t0\t1\t1\t0.2\n"
            "F1\tmo\t0\t0\t2\t1\t-0.4\n"
        )
        rows, _ = _rows(text)
        pgs = next(r for r in rows if r.kind == "pgs")
        assert [name for name, _ in pgs.key] == ["PGS"]


def test_style_json_round_trip():
    style = pp.PedigreeStyle.from_dict(
        {"legend": "full", "phenotype_labels": {"SRS": "SRS T-score"}}
    )
    resolved = pp.resolve_style(style, ["SRS"])
    assert resolved.legend == "full"
    assert resolved.phenotype_labels == {"SRS": "SRS T-score"}
    resolved.phenotype_labels["x"] = "y"
    assert "x" not in style.phenotype_labels  # resolve_style copies it
    assert replace(style, legend="auto").legend == "auto"
