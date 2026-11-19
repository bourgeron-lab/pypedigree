"""Text in phenotype bands and quadrants fits its place (#55).

All pedigrees here are synthetic.
"""

from __future__ import annotations

import math
import re

import matplotlib.pyplot as plt
import pytest
from conftest import TEST_DATA_DIR
from matplotlib.text import Text

import pypedigree as pp
from pypedigree import symbols
from pypedigree.errors import PedigreeError
from pypedigree.models import Sex
from pypedigree.phenotypes import (
    FIT_FLOOR_SCALE,
    BandCell,
    band_room,
    cell_texts,
    fit_text,
    quadrant_room,
)

PT = symbols.UNITS_PER_POINT
HALF = 0.4  # half the default shape size

# Three gradient columns with a long name, as in Freddy's journal figure.
SCORES = (
    "FID\tIID\tPAT\tMAT\tSEX\tPHENO\tg_SRS\tg_IQ\tg_ADHD-RS\n"
    "F1\tfa\t0\t0\t1\t1\t61\t95\t14\n"
    "F1\tmo\t0\t0\t2\t1\t58\t110\t11\n"
    "F1\tc1\tfa\tmo\t1\t2\t78\t72\t31\n"
    "F1\tc2\tfa\tmo\t2\t2\t69\t85\t22\n"
    "F1\tc3\tfa\tmo\t0\t1\t55\t108\t8\n"
)

FIXTURES = [
    "phenotype/test_all.ped",
    "phenotype/multi_phenotype_inline.ped",
    "phenotype/prefix_modes.ped",
    "phenotype/carrier_genes.ped",
    "complex/multi_phenotype_pgs.ped",
]


def _wide(_h: float) -> float:
    return 10.0


class TestBoldWidths:
    def test_table_matches_dejavu_sans_bold(self):
        from matplotlib import font_manager, ft2font
        from matplotlib.font_manager import FontProperties

        path = font_manager.findfont(FontProperties(family="DejaVu Sans", weight="bold"))
        font = ft2font.FT2Font(path)
        if hasattr(ft2font, "LoadFlags"):
            flags = ft2font.LoadFlags.NO_SCALE
        else:
            flags = ft2font.LOAD_NO_SCALE
        for code in range(32, 127):
            glyph = font.load_char(code, flags=flags)
            assert symbols._CHAR_EM_BOLD[chr(code)] == pytest.approx(
                glyph.horiAdvance / font.units_per_EM, abs=5e-4
            ), chr(code)

    def test_bold_is_wider(self):
        assert symbols.text_width("IQ: 103", 8, bold=True) > symbols.text_width("IQ: 103", 8)


class TestGeometry:
    def test_square_keeps_its_width(self):
        assert symbols.shape_half_width_at(Sex.MALE, 1.0, 0.9) == 1.0
        assert symbols.shape_half_width_at(Sex.MALE, 1.0, 1.1) == 0.0

    def test_circle_chord(self):
        assert symbols.shape_half_width_at(Sex.FEMALE, 1.0, 0.6) == pytest.approx(0.8)
        assert symbols.shape_half_width_at(Sex.FEMALE, 1.0, -0.6) == pytest.approx(0.8)

    def test_diamond(self):
        reach = symbols.DIAMOND_SCALE
        assert symbols.shape_half_width_at(Sex.UNKNOWN, 1.0, 0.0) == pytest.approx(reach)
        assert symbols.shape_half_width_at(Sex.UNKNOWN, 1.0, 0.5) == pytest.approx(reach - 0.5)

    def test_band_room(self):
        assert band_room(Sex.FEMALE, 1.0, 0.5)(0.2) == pytest.approx(1.6)
        assert band_room(Sex.MALE, 1.0, 0.5)(0.2) == pytest.approx(2.0)

    def test_quadrant_room_stops_at_the_axis_and_the_edge(self):
        assert quadrant_room(Sex.MALE, 1.0, "top-left")(0.1) == pytest.approx(1.0)
        edge = math.sqrt(1 - 0.55**2)
        assert quadrant_room(Sex.FEMALE, 1.0, "top-right")(0.1) == pytest.approx(2 * (edge - 0.5))


class TestCellTexts:
    GRADIENT = BandCell("#ffffff", "SRS: 78", False, "78")

    @pytest.mark.parametrize(
        ("mode", "texts"),
        [
            ("auto", ("SRS: 78", "78")),
            ("name_value", ("SRS: 78",)),
            ("value", ("78",)),
            ("none", ()),
        ],
    )
    def test_modes(self, mode, texts):
        assert cell_texts(self.GRADIENT, mode) == texts

    def test_na_whatever_the_mode(self):
        na = BandCell("#d9d9d9", "NA", True)
        assert cell_texts(na, "none") == ("NA",)

    def test_text_without_short_form(self):
        assert cell_texts(BandCell("#ffffff", "Control", False), "auto") == ("Control",)

    def test_no_text(self):
        assert cell_texts(BandCell("#ffffff", None, False), "auto") == ()


class TestFitText:
    STYLE = pp.PedigreeStyle(font_size=8.0)

    def test_the_name_when_it_fits(self):
        fit = fit_text(("SRS: 78", "78"), 6.0, self.STYLE, PT, _wide, 1.0)
        assert (fit.text, fit.size) == ("SRS: 78", 6.0)

    def test_the_value_alone_when_the_name_does_not_fit(self):
        width = symbols.text_width("78", 6.0, bold=True) * 1.5
        fit = fit_text(("SRS: 78", "78"), 6.0, self.STYLE, PT, lambda h: width, 1.0)
        assert (fit.text, fit.size) == ("78", 6.0)

    def test_smaller_when_the_value_does_not_fit(self):
        width = symbols.text_width("78", 6.0, bold=True) * 0.9
        fit = fit_text(("SRS: 78", "78"), 6.0, self.STYLE, PT, lambda h: width, 1.0)
        assert fit.text == "78"
        assert 8.0 * FIT_FLOOR_SCALE <= fit.size < 6.0
        used = symbols.text_width("78", fit.size, bold=True) + 0.3 * fit.size * PT
        assert used <= width + 1e-9

    def test_nothing_below_the_floor(self):
        assert fit_text(("78",), 6.0, self.STYLE, PT, lambda h: 0.01, 1.0) is None

    def test_the_floor_is_the_minimum_size_when_set(self):
        style = pp.PedigreeStyle(font_size=8.0, min_font_size=6.0)
        width = symbols.text_width("78", 6.0, bold=True) * 0.9
        assert fit_text(("78",), 6.0, style, PT, lambda h: width, 1.0) is None

    def test_too_low_a_band(self):
        fit = fit_text(("78",), 6.0, self.STYLE, PT, _wide, 4.0 * PT)
        assert fit.size <= 5.0 + 1e-6

    @pytest.mark.parametrize("mode", ["name_value", "value", "none"])
    def test_other_modes_are_not_fitted(self, mode):
        style = pp.PedigreeStyle(band_text=mode)
        fit = fit_text(("SRS: 78",), 6.0, style, PT, lambda h: 0.0, 0.0)
        assert (fit.text, fit.size) == ("SRS: 78", 6.0)

    def test_invalid_mode(self):
        ped = pp.parse_pedigree(TEST_DATA_DIR / "basic/simple_trio.ped")
        style = pp.PedigreeStyle(band_text="short")
        with pytest.raises(PedigreeError, match="band_text"):
            pp.render_svg(ped, style=style)


def _svg_overflows(svg: str) -> list[str]:
    """Bold texts of the SVG that reach outside their symbol or quadrant."""
    bad = []
    for group in svg.split('<g class="node"')[1:]:
        # The band and quadrant texts come before the symbol's outline; the
        # proband's "P" and labels after it.
        group = group.split('class="shape"')[0]
        nx = float(re.search(r'data-x="([-0-9.]+)"', group).group(1))
        ny = float(re.search(r'data-y="([-0-9.]+)"', group).group(1))
        shape = re.search(r"<(rect|circle|polygon)[^>]*$", group)
        if shape is None:
            continue
        sex = {"rect": Sex.MALE, "circle": Sex.FEMALE}.get(shape.group(1), Sex.UNKNOWN)
        for m in re.finditer(
            r'<text x="([-0-9.]+)" y="([-0-9.]+)" font-size="([0-9.]+)"[^>]*'
            r'font-weight="bold">([^<]*)</text>',
            group,
        ):
            x, y, size, text = float(m[1]), float(m[2]), float(m[3]), m[4]
            yc = y - 0.35 * size  # the "central" baseline shift
            width = symbols.text_width(text, size / PT, bold=True)
            dx, h = abs(x - nx), 0.8 * size
            edge = symbols.shape_half_width_at(sex, HALF, abs(yc - ny) + h / 2)
            if dx + width / 2 > edge + 1e-3 or (dx > 1e-6 and dx - width / 2 < -1e-3):
                bad.append(text)
    return bad


class TestEveryTextFits:
    @pytest.mark.parametrize("fixture", FIXTURES)
    @pytest.mark.parametrize("layout", ["auto", "bands"])
    def test_svg(self, fixture, layout):
        ped = pp.parse_pedigree(TEST_DATA_DIR / fixture)
        svg = pp.render_svg(ped, style=pp.PedigreeStyle(phenotype_layout=layout))
        assert _svg_overflows(svg) == []

    def test_svg_journal_size(self):
        ped = pp.parse_pedigree_from_string(SCORES)
        style = pp.PedigreeStyle(font_size=6.0, min_font_size=6.0, phenotype_layout="bands")
        assert _svg_overflows(pp.render_svg(ped, style=style)) == []

    def test_the_check_finds_overflowing_text(self):
        ped = pp.parse_pedigree_from_string(SCORES)
        style = pp.PedigreeStyle(phenotype_layout="bands", band_text="name_value")
        assert "ADHD-RS: 14" in _svg_overflows(pp.render_svg(ped, style=style))

    def test_svg_shortens_the_text_where_needed(self):
        ped = pp.parse_pedigree_from_string(SCORES)
        svg = pp.render_svg(ped, style=pp.PedigreeStyle(phenotype_layout="bands"))
        assert ">ADHD-RS: 14</text>" not in svg  # still in the tooltip title
        assert ">14</text>" in svg
        assert ">IQ: 95</text>" in svg  # the square's middle band has room


class TestPng:
    def _band_texts(self, fig) -> list[Text]:
        fig.canvas.draw()
        return [t for t in fig.axes[0].texts if t.get_fontweight() == "bold"]

    def test_journal_figure(self):
        """A 17 cm figure drawn into the caller's axes: the text is chosen at that scale."""
        ped = pp.parse_pedigree_from_string(SCORES)
        style = pp.PedigreeStyle(font_size=6.0, min_font_size=6.0, phenotype_layout="bands")
        fig, ax = plt.subplots(figsize=(17 / 2.54, 4 / 2.54))
        pp.render_pedigree(ped, style=style, ax=ax, show_legend=False)
        texts = [t.get_text() for t in self._band_texts(fig)]
        plt.close(fig)
        assert "14" in texts
        assert not any(t.startswith("ADHD-RS") for t in texts)

    def test_large_figure_keeps_the_name(self):
        ped = pp.parse_pedigree_from_string(SCORES)
        fig = pp.render_pedigree(ped, style=pp.PedigreeStyle(phenotype_layout="bands"))
        texts = [t.get_text() for t in self._band_texts(fig)]
        plt.close(fig)
        assert "IQ: 95" in texts

    def test_refitted_when_the_figure_is_resized(self):
        ped = pp.parse_pedigree_from_string(SCORES)
        fig, ax = plt.subplots(figsize=(30, 12))
        pp.render_pedigree(ped, style=pp.PedigreeStyle(phenotype_layout="bands"), ax=ax)
        before = {t.get_text() for t in self._band_texts(fig)}
        fig.set_size_inches(5, 2)
        after = {t.get_text() for t in self._band_texts(fig)}
        plt.close(fig)
        assert "ADHD-RS: 14" in before
        assert "ADHD-RS: 14" not in after

    def test_name_value_draws_plain_text(self):
        ped = pp.parse_pedigree_from_string(SCORES)
        style = pp.PedigreeStyle(phenotype_layout="bands", band_text="name_value")
        fig = pp.render_pedigree(ped, style=style)
        texts = {t.get_text() for t in self._band_texts(fig)}
        plt.close(fig)
        assert "ADHD-RS: 14" in texts

    def test_none_writes_no_value(self):
        ped = pp.parse_pedigree_from_string(SCORES)
        style = pp.PedigreeStyle(phenotype_layout="bands", band_text="none")
        fig = pp.render_pedigree(ped, style=style)
        texts = [t.get_text() for t in self._band_texts(fig) if t.get_text()]
        plt.close(fig)
        assert texts == []
