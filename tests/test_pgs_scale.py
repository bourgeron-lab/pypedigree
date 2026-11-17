"""PGS on a fixed scale (issue #7).

``pgs_scale="minmax"`` (the default) places a PGS value between the
pedigree's lowest and highest score. ``"zscore"`` and ``"percentile"`` place
it on the drawn standard normal, the same for every pedigree, mark the
+/-2 SD tails and give z and percentile in the SVG tooltip. PGS curves are
drawn even when the pedigree has no other phenotype column.
"""

from __future__ import annotations

import re

import matplotlib.pyplot as plt
import pytest
from click.testing import CliRunner

from pypedigree import PedigreeStyle, parse_pedigree_from_string, render_pedigree
from pypedigree.cli import main
from pypedigree.errors import PedigreeError
from pypedigree.phenotypes import (
    PGS_CURVE_SIGMA,
    active_pgs,
    check_pgs_scale,
    pgs_legend_description,
    pgs_position,
    pgs_tail_positions,
    pgs_value_text,
)
from pypedigree.svg_renderer import render_interactive_svg

# The trio of the issue: z = -0.2 (42nd percentile), +1.1, +2.5.
TRIO = (
    "FID\tIID\tFather\tMother\tSex\tPhenotype\tPGS_autism\n"
    "F\tdad\t0\t0\t1\t1\t-0.2\n"
    "F\tmum\t0\t0\t2\t1\t1.1\n"
    "F\tkid\tdad\tmum\t1\t2\t2.5\n"
)
ZSCORE = PedigreeStyle(pgs_scale="zscore")


def _positions(text, scale):
    family = parse_pedigree_from_string(text).families["F"]
    return {
        sid: active_pgs(ind, ["PGS_autism"], scale)[0][1]
        for sid, ind in family.individuals.items()
        if active_pgs(ind, ["PGS_autism"], scale)
    }


class TestPositions:
    def test_zscore_is_placed_on_the_drawn_normal(self):
        assert pgs_position(0.0) == pytest.approx(0.5)
        assert pgs_position(2.0) == pytest.approx(0.5 + 2 * PGS_CURVE_SIGMA)
        assert pgs_position(10.0) == 1.0 and pgs_position(-10.0) == 0.0

    def test_minmax_spans_the_curve_whatever_the_values(self):
        pos = _positions(TRIO, "minmax")
        assert pos["dad"] == pytest.approx(0.0) and pos["kid"] == pytest.approx(1.0)

    def test_zscore_keeps_the_population_meaning(self):
        pos = _positions(TRIO, "zscore")
        assert pos["dad"] == pytest.approx(0.5 - 0.2 * PGS_CURVE_SIGMA)
        assert pos["mum"] == pytest.approx(0.5 + 1.1 * PGS_CURVE_SIGMA)
        assert pos["kid"] > pgs_tail_positions("zscore")[1], "beyond +2 SD"

    def test_close_z_scores_stay_close(self):
        text = TRIO.replace("-0.2", "0.05").replace("1.1", "-0.05").replace("2.5", "0.0")
        pos = _positions(text, "zscore")
        assert max(pos.values()) - min(pos.values()) < 0.05
        spread = _positions(text, "minmax")
        assert max(spread.values()) - min(spread.values()) == pytest.approx(1.0)

    def test_percentile(self):
        text = TRIO.replace("-0.2", "50").replace("1.1", "97.72").replace("2.5", "100")
        pos = _positions(text, "percentile")
        assert pos["dad"] == pytest.approx(0.5)
        assert pos["mum"] == pytest.approx(0.5 + 2 * PGS_CURVE_SIGMA, abs=1e-3)
        assert pos["kid"] == 1.0  # 100 is clipped, not an error

    def test_tails_only_on_a_fixed_scale(self):
        assert pgs_tail_positions("minmax") == ()
        low, high = pgs_tail_positions("percentile")
        assert (low, high) == pytest.approx((0.5 - 2 * PGS_CURVE_SIGMA, 0.5 + 2 * PGS_CURVE_SIGMA))

    def test_unknown_scale_is_rejected(self):
        with pytest.raises(PedigreeError, match="Unknown pgs_scale"):
            check_pgs_scale(PedigreeStyle(pgs_scale="rank"))


class TestText:
    @pytest.mark.parametrize(
        ("raw", "scale", "expected"),
        [
            (-0.2, "zscore", "z = -0.20, 42nd percentile"),
            (1.1, "zscore", "z = +1.10, 86th percentile"),
            (2.5, "zscore", "z = +2.50, 99th percentile"),
            (0.0, "zscore", "z = +0.00, 50th percentile"),
            (11.0, "percentile", "11th percentile, z = -1.23"),
            (0.72, "minmax", "0.72"),
        ],
    )
    def test_values(self, raw, scale, expected):
        assert pgs_value_text(raw, scale) == expected

    def test_legend_caption_names_the_scale(self):
        assert pgs_legend_description(["PGS_autism"], "zscore").endswith("z-score, ±2 SD marked")
        assert "percentile" in pgs_legend_description(["PGS_autism"], "percentile")
        assert pgs_legend_description(["PGS_autism"]) == "PGS distribution (autism)"


class TestRendering:
    def test_png_curves_without_other_phenotypes_and_tails(self):
        for style, tails in ((PedigreeStyle(), 0), (ZSCORE, 2 * 3)):
            fig = render_pedigree(parse_pedigree_from_string(TRIO), style=style)
            try:
                ax = fig.axes[0]
                assert [t.get_text() for t in ax.texts].count("PGS") == 3, "a curve per person"
                dotted = [ln for ln in ax.lines if ln.get_linestyle() == ":"]
                assert len(dotted) == tails
                assert len(fig.axes) == 2, "the legend is drawn"
            finally:
                plt.close(fig)

    def test_svg_curves_tails_tooltip_and_legend(self):
        svg = render_interactive_svg(parse_pedigree_from_string(TRIO), style=ZSCORE)
        assert svg.count('<g class="pgs">') == 3
        assert svg.count('class="pgs-tail"') == 6
        assert "PGS_autism: z = +2.50, 99th percentile" in svg
        assert "z-score, ±2 SD marked" in svg
        minmax = render_interactive_svg(parse_pedigree_from_string(TRIO))
        assert 'class="pgs-tail"' not in minmax
        assert re.search(r"PGS_autism: 2\.5\b", minmax)

    def test_cli_and_style_json(self, tmp_path):
        ped = tmp_path / "trio.ped"
        ped.write_text(TRIO)
        out = tmp_path / "trio.svg"
        args = [str(ped), "--interactive", "--pgs-scale", "zscore", "-o", str(out)]
        result = CliRunner().invoke(main, args)
        assert result.exit_code == 0, result.output
        assert 'class="pgs-tail"' in out.read_text(encoding="utf-8")
        assert PedigreeStyle.from_dict({"pgs_scale": "percentile"}).pgs_scale == "percentile"
