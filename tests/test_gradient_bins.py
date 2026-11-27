"""Clinical thresholds on quantitative scores (#57).

``gradient_bins`` cuts a gradient column into steps at its cut-offs; the
legend shows the steps with their limits. With ``gradient_range`` the
legend's gradient bar writes the range's ends instead of Low / High.
All pedigrees here are synthetic.
"""

from __future__ import annotations

import json
import re

import matplotlib.pyplot as plt
import pytest
from conftest import PHENOTYPE_DIR, decode_char_refs

from pypedigree import PedigreeStyle, parse_pedigree, render_interactive_svg, render_pedigree
from pypedigree.carriers import CarrierScheme
from pypedigree.config import resolve_style
from pypedigree.errors import PedigreeError
from pypedigree.models import PhenotypeDisplayMode
from pypedigree.phenotypes import (
    LEG_X0,
    band_cell,
    bin_color,
    bin_index,
    bin_labels,
    binned_row_geometry,
    blend_with_white,
    check_gradient_bins,
    legend_metrics,
    legend_rows,
    resolve_phenotype_render,
    to_rgb,
    with_quadrants,
)
from pypedigree.symbols import UNITS_PER_POINT, text_width

FIXTURE = PHENOTYPE_DIR / "gradient_bins.ped"
STYLE = json.loads((PHENOTYPE_DIR / "gradient_bins.json").read_text())
SRS = [60, 66, 76]


def _style(**extra) -> PedigreeStyle:
    return PedigreeStyle(**{**STYLE, **extra})


def _resolved(style: PedigreeStyle):
    pedigree = parse_pedigree(FIXTURE)
    names = pedigree.phenotype_names
    st = resolve_style(style, names)
    render = resolve_phenotype_render(pedigree, names, pedigree.phenotype_display_modes, st)
    render = with_quadrants(render, render.band_names(names), CarrierScheme(), st)
    return pedigree, st, render


def _rows(style: PedigreeStyle):
    pedigree, st, render = _resolved(style)
    names = pedigree.phenotype_names
    return legend_rows(names, st, render, pedigree.phenotype_display_modes, CarrierScheme())


def _cell(sample: str, column: str, style: PedigreeStyle):
    pedigree, st, render = _resolved(style)
    node = next(i for i in pedigree.individuals() if i.sample_id == sample)
    return band_cell(node, column, st, render, pedigree.phenotype_display_modes)


class TestSteps:
    @pytest.mark.parametrize(
        ("value", "step"),
        [(40, 0), (59.9, 0), (60, 1), (65, 1), (66, 2), (75.5, 2), (76, 3), (99, 3)],
    )
    def test_value_at_a_cutoff_goes_up(self, value, step):
        assert bin_index(value, SRS) == step

    def test_labels(self):
        assert bin_labels(SRS) == ("< 60", "60–65", "66–75", "≥ 76")
        assert bin_labels([70]) == ("< 70", "≥ 70")
        assert bin_labels([1, 2, 5]) == ("< 1", "1", "2–4", "≥ 5")
        assert bin_labels([0.5, 1.5]) == ("< 0.5", "0.5–1.5", "≥ 1.5")

    def test_colors_from_white_to_the_column_color(self):
        rgb = to_rgb("#e67e22")
        steps = [bin_color(rgb, i, SRS, inverted=False) for i in range(4)]
        assert steps[0] == pytest.approx((1.0, 1.0, 1.0))
        assert steps[-1] == pytest.approx(rgb)
        assert steps[1] == pytest.approx(blend_with_white(rgb, 1 / 3))

    def test_inverted_column_bins_the_other_way(self):
        rgb = to_rgb("#2ecc71")
        assert bin_color(rgb, 0, [70], inverted=True) == pytest.approx(rgb)
        assert bin_color(rgb, 1, [70], inverted=True) == pytest.approx((1.0, 1.0, 1.0))


class TestBands:
    def test_band_takes_its_step_color(self):
        style = _style()
        _ped, st, _render = _resolved(style)
        color = to_rgb(st.phenotype_colors["SRS"])
        assert to_rgb(_cell("C2", "SRS", style).color) == pytest.approx(
            bin_color(color, 2, SRS, inverted=False)
        )  # 70: 66-75
        assert to_rgb(_cell("F", "SRS", style).color) == pytest.approx((1.0, 1.0, 1.0))  # 52

    def test_value_still_written(self):
        assert _cell("C1", "SRS", _style()).short == "81"

    def test_gi_iq_below_70_takes_the_full_color(self):
        style = _style()
        _ped, st, _render = _resolved(style)
        iq = to_rgb(st.phenotype_colors["IQ"])
        assert to_rgb(_cell("C1", "IQ", style).color) == pytest.approx(iq)  # 64
        assert to_rgb(_cell("C3", "IQ", style).color) == pytest.approx((1.0, 1.0, 1.0))  # 71

    def test_bins_take_precedence_over_the_range(self):
        style = _style(gradient_range={"SRS": [0, 200]})
        full = to_rgb(_cell("C1", "SRS", style).color)  # 81 >= 76: the last step
        _ped, st, _render = _resolved(style)
        assert full == pytest.approx(to_rgb(st.phenotype_colors["SRS"]))

    def test_bins_make_a_column_a_gradient(self):
        pedigree = parse_pedigree(FIXTURE)
        modes = dict(pedigree.phenotype_display_modes)
        modes["SRS"] = PhenotypeDisplayMode.AUTO
        render = resolve_phenotype_render(
            pedigree, ["SRS"], modes, PedigreeStyle(gradient_bins={"SRS": [1, 2]})
        )
        assert render.effective("SRS") == "gradient"


class TestLegend:
    def test_binned_rows(self):
        rows = {r.label: r for r in _rows(_style())}
        assert rows["SRS"].kind == "binned"
        assert rows["SRS"].bins == ("< 60", "60–65", "66–75", "≥ 76")
        assert len(rows["SRS"].colors) == 4
        assert rows["IQ"].kind == "binned" and rows["IQ"].bins == ("< 70", "≥ 70")

    def test_range_ends_instead_of_low_high(self):
        rows = {r.label: r for r in _rows(_style())}
        assert (rows["ADHD RS"].low, rows["ADHD RS"].high) == ("0", "54")
        plain = {r.label: r for r in _rows(PedigreeStyle())}
        assert (plain["ADHD RS"].low, plain["ADHD RS"].high) == ("Low", "High")

    def test_png_and_svg_write_the_limits(self):
        fig = render_pedigree(parse_pedigree(FIXTURE), style=_style())
        try:
            texts = [t.get_text() for t in fig.axes[1].texts]
        finally:
            plt.close(fig)
        svg = decode_char_refs(render_interactive_svg(parse_pedigree(FIXTURE), style=_style()))
        legend = svg[svg.index('<g class="legend">') :]
        for limits in ("< 60", "60–65", "66–75", "≥ 76", "< 70", "≥ 70"):
            assert limits in texts
            assert f">{limits.replace('<', '&lt;')}<" in legend
        assert "0" in texts and "54" in texts and "Low" not in texts

    def test_svg_tooltip_names_the_step(self):
        svg = decode_char_refs(render_interactive_svg(parse_pedigree(FIXTURE), style=_style()))
        # The gallery's style also names the steps (#81).
        assert re.search(r"SRS: 70 \(66–75, moderate\)", svg)
        assert re.search(r"IQ: 64 \(&lt; 70, intellectual disability\)", svg)

    def test_steps_wide_enough_for_their_limits(self):
        """Many steps with long limits widen the bar; the label follows it."""
        style = PedigreeStyle(gradient_bins={"SRS": [10.25, 20.5, 30.75, 40.5, 50.25, 60.5]})
        rows = {r.label: r for r in _rows(style)}
        row = rows["SRS"]
        m = legend_metrics(4, style)
        g = binned_row_geometry(row, m)
        size = m.caption_fs()
        for limits in row.bins:
            width = text_width(limits, size) / UNITS_PER_POINT / 72 / m.unit_x
            assert width < g.step_w
        assert g.label_x > LEG_X0 + g.step_w * len(row.colors)
        assert binned_row_geometry(_rows(_style())[1], m).label_x >= m.text_x


class TestValidation:
    @pytest.mark.parametrize("bad", [[66, 60], [], [60, 60], ["60"], [True], 60, [float("nan")]])
    def test_invalid(self, bad):
        with pytest.raises(PedigreeError, match="gradient_bins"):
            check_gradient_bins(PedigreeStyle(gradient_bins={"SRS": bad}))

    def test_rendering_checks_them(self):
        with pytest.raises(PedigreeError, match="gradient_bins"):
            render_pedigree(
                parse_pedigree(FIXTURE), style=PedigreeStyle(gradient_bins={"SRS": [66, 60]})
            )
        with pytest.raises(PedigreeError, match="gradient_bins"):
            render_interactive_svg(
                parse_pedigree(FIXTURE), style=PedigreeStyle(gradient_bins={"SRS": [66, 60]})
            )

    def test_from_json_style(self):
        assert PedigreeStyle.from_dict(STYLE).gradient_bins == {"SRS": [60, 66, 76], "IQ": [70]}

    def test_resolve_style_copies_them(self):
        style = _style()
        resolved = resolve_style(style, ["SRS"])
        resolved.gradient_bins["SRS"].append(90)
        assert style.gradient_bins["SRS"] == [60, 66, 76]
