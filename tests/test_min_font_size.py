"""A minimum size for every text, for journal figures (#51).

All pedigrees here are synthetic.
"""

from __future__ import annotations

import re

import matplotlib.pyplot as plt
import pytest
from click.testing import CliRunner
from conftest import TEST_DATA_DIR
from matplotlib.text import Text

import pypedigree as pp
from pypedigree import symbols
from pypedigree.cli import main
from pypedigree.errors import PedigreeError
from pypedigree.layout import compute_layout
from pypedigree.phenotypes import QUADRANT_FONT_SCALE, legend_metrics

# Every kind of text: quadrants and bands, gradient values, "NA", PGS
# captions, annotations, notation notes and letters, evaluation marks, carrier
# sections, twins of unknown zygosity, a couple with no child ("?"), legend.
FIXTURES = [
    "phenotype/test_all.ped",
    "phenotype/multi_phenotype_inline.ped",
    "phenotype/carrier_genes.ped",
    "phenotype/prefix_modes.ped",
    "nomenclature/nomenclature.ped",
    "twins/twins_unknown.ped",
    "basic/couple_only.ped",
    "extended/extended_full.ped",
]

SMALL = pp.PedigreeStyle(font_size=6.0, min_font_size=6.0)


def _sizes_png(fig) -> list[float]:
    fig.canvas.draw()
    return [t.get_fontsize() for t in fig.findobj(Text) if t.get_text().strip()]


def _sizes_svg(svg: str) -> list[float]:
    return [float(v) for v in re.findall(r'<text [^>]*font-size="([0-9.]+)"', svg)]


class TestTextSize:
    def test_no_minimum(self):
        style = pp.PedigreeStyle(font_size=10.0)
        assert style.text_size() == 10.0
        assert style.text_size(0.6) == pytest.approx(6.0)

    def test_minimum_raises_small_text_only(self):
        style = pp.PedigreeStyle(font_size=6.0, min_font_size=5.0)
        assert style.text_size() == 6.0
        assert style.text_size(0.6) == 5.0
        assert style.text_size(2.0) == 12.0

    def test_minimum_above_font_size(self):
        assert pp.PedigreeStyle(font_size=5.0, min_font_size=6.0).text_size() == 6.0

    def test_from_json(self):
        assert pp.PedigreeStyle.from_dict({"min_font_size": 6}).text_size(0.5) == 6.0

    @pytest.mark.parametrize("bad", [0, -1, "6", True])
    def test_invalid(self, bad):
        ped = pp.parse_pedigree(TEST_DATA_DIR / "basic/simple_trio.ped")
        style = pp.PedigreeStyle(min_font_size=bad)
        with pytest.raises(PedigreeError, match="min_font_size"):
            pp.render_svg(ped, style=style)
        with pytest.raises(PedigreeError, match="min_font_size"):
            pp.render_pedigree(ped, style=style)


class TestEveryTextReachesTheMinimum:
    @pytest.mark.parametrize("fixture", FIXTURES)
    @pytest.mark.parametrize("layout", ["auto", "bands"])
    def test_png(self, fixture, layout):
        ped = pp.parse_pedigree(TEST_DATA_DIR / fixture)
        style = pp.PedigreeStyle(font_size=6.0, min_font_size=6.0, phenotype_layout=layout)
        fig = pp.render_pedigree(ped, style=style)
        sizes = _sizes_png(fig)
        plt.close(fig)
        assert sizes
        assert min(sizes) >= 6.0

    @pytest.mark.parametrize("fixture", FIXTURES)
    @pytest.mark.parametrize("layout", ["auto", "bands"])
    def test_svg(self, fixture, layout):
        ped = pp.parse_pedigree(TEST_DATA_DIR / fixture)
        style = pp.PedigreeStyle(font_size=6.0, min_font_size=6.0, phenotype_layout=layout)
        sizes = _sizes_svg(pp.render_svg(ped, style=style))
        assert sizes
        assert min(sizes) >= 6.0 * symbols.UNITS_PER_POINT - 1e-9

    def test_circular(self):
        ped = pp.parse_pedigree(TEST_DATA_DIR / "phenotype/test_all.ped")
        fig = pp.render_pedigree(ped, style=SMALL, circular=True)
        sizes = _sizes_png(fig)
        plt.close(fig)
        assert min(sizes) >= 6.0

    def test_without_minimum_text_is_smaller(self):
        """The fixtures do have text below 6 pt at 6 pt, so the tests above test something."""
        ped = pp.parse_pedigree(TEST_DATA_DIR / "phenotype/test_all.ped")
        sizes = _sizes_svg(pp.render_svg(ped, style=pp.PedigreeStyle(font_size=6.0)))
        assert min(sizes) == pytest.approx(6.0 * 0.6 * symbols.UNITS_PER_POINT, abs=1e-3)

    def test_secondary_text_keeps_its_size_above_the_minimum(self):
        ped = pp.parse_pedigree(TEST_DATA_DIR / "phenotype/test_all.ped")
        style = pp.PedigreeStyle(font_size=20.0, min_font_size=6.0)
        sizes = _sizes_svg(pp.render_svg(ped, style=style))
        quadrant = 20.0 * QUADRANT_FONT_SCALE * symbols.UNITS_PER_POINT
        assert any(s == pytest.approx(quadrant, abs=1e-3) for s in sizes)


class TestLayoutAndLegend:
    ANNOTATED = (
        "FID\tIID\tPAT\tMAT\tSEX\tPHENO\ta_variant\n"
        "F1\tf\t0\t0\t1\t1\t\n"
        "F1\tm\t0\t0\t2\t1\t\n"
        "F1\tc1\tf\tm\t1\t2\tNM_000000.1:c.1234_1240delinsACGTACG\n"
        "F1\tc2\tf\tm\t2\t1\tNM_000000.1:c.4321_4327delinsTTTTGGG\n"
    )

    def test_larger_notes_spread_the_symbols(self):
        """The layout measures the lines below a label at their drawn size (#41)."""
        family = next(iter(pp.parse_pedigree_from_string(self.ANNOTATED).families.values()))

        def gap(style: pp.PedigreeStyle) -> float:
            xs = {n.individual_id: n.x for n in compute_layout(family, style).nodes}
            return xs["c2"] - xs["c1"]

        assert gap(pp.PedigreeStyle(font_size=6.0, min_font_size=8.0)) > gap(
            pp.PedigreeStyle(font_size=6.0)
        )

    def test_legend_grows_with_its_text(self):
        plain = legend_metrics(3, pp.PedigreeStyle(font_size=6.0))
        raised = legend_metrics(3, SMALL)
        assert raised.fs >= 6.0
        assert raised.caption_fs(0.8) >= 6.0
        assert raised.title_fs > plain.title_fs
        assert raised.unit_y > plain.unit_y
        assert raised.font == pytest.approx(6.0 / 0.85)

    def test_legend_unchanged_without_minimum(self):
        m = legend_metrics(3, pp.PedigreeStyle(font_size=8.0))
        assert (m.font, m.fs, m.label_fs, m.title_fs) == pytest.approx((8.0, 6.8, 5.1, 7.8))


class TestCli:
    def test_option(self, tmp_path):
        out = tmp_path / "f.svg"
        args = [str(TEST_DATA_DIR / "phenotype/test_all.ped"), "--interactive", "-o", str(out)]
        result = CliRunner().invoke(main, [*args, "--min-font-size", "7"])
        assert result.exit_code == 0, result.output
        sizes = _sizes_svg(out.read_text(encoding="utf-8"))
        assert min(sizes) >= 7.0 * symbols.UNITS_PER_POINT - 1e-9

    @pytest.mark.parametrize("bad", ["0", "-2", "big"])
    def test_invalid(self, tmp_path, bad):
        args = [str(TEST_DATA_DIR / "basic/simple_trio.ped"), "-o", str(tmp_path / "f.png")]
        result = CliRunner().invoke(main, [*args, "--min-font-size", bad])
        assert result.exit_code == 2
        assert "--min-font-size" in result.output
