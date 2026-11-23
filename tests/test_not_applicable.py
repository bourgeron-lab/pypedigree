"""Not applicable (N/A) as a value in phenotype columns (#107).

A phenotype that cannot apply to someone (ovarian cancer for someone assigned
male at birth, prostate cancer for someone assigned female at birth) is
written ``N/A`` or ``not_applicable``; ``NA`` keeps meaning no data. It is
drawn white, crossed by a thin grey diagonal, explained by a legend row, and
not counted as data. All pedigrees here are synthetic.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import pytest
from matplotlib.colors import to_hex

from pypedigree import (
    PedigreeBuilder,
    PedigreeStyle,
    parse_pedigree_from_string,
    render_interactive_svg,
    render_pedigree,
)
from pypedigree.carriers import CarrierScheme
from pypedigree.config import resolve_style
from pypedigree.layout import compute_layout
from pypedigree.phenotypes import (
    NOT_APPLICABLE_STROKE,
    all_bands_na,
    band_cell,
    check_legend,
    detect_auto_mode,
    legend_rows,
    not_applicable_stroke,
    quadrant_cell,
    resolve_phenotype_render,
)

CANCERS = (
    "FID\tIID\tPAT\tMAT\tSEX\tPHENO\tb_breast_cancer\tb_ovarian_cancer\tb_prostate_cancer\n"
    "FAM1\tA\t0\t0\t1\t1\t1\tN/A\t2\n"
    "FAM1\tB\t0\t0\t2\t2\t2\t1\tn/a\n"
    "FAM1\tC\tA\tB\t2\t2\t2\t2\tnot_applicable\n"
    "FAM1\tD\tA\tB\t1\t1\t-9\tNot Applicable\t-9\n"
    "FAM1\tE\tA\tB\t1\t1\tNA\tN/A\tN/A\n"
)
NA_CELLS = 6  # A, B, C, D, and E twice


def _people(text: str = CANCERS):
    ped = parse_pedigree_from_string(text)
    return ped, ped.families["FAM1"].individuals


def _without_na() -> str:
    """The same pedigree with every N/A replaced by 1 (unaffected)."""
    text = CANCERS
    for spelling in ("N/A", "n/a", "not_applicable", "Not Applicable"):
        text = text.replace(spelling, "1")
    return text


def _render(ped):
    names = ped.phenotype_names
    style = resolve_style(PedigreeStyle(), names)
    render = resolve_phenotype_render(ped, names, ped.phenotype_display_modes, style)
    return style, render


class TestParsing:
    def test_spellings(self):
        _, people = _people()
        assert people["A"].not_applicable == {"ovarian_cancer"}
        assert people["B"].not_applicable == {"prostate_cancer"}
        assert people["C"].not_applicable == {"prostate_cancer"}
        assert people["D"].not_applicable == {"ovarian_cancer"}
        assert people["E"].not_applicable == {"ovarian_cancer", "prostate_cancer"}
        assert people["A"].extra_phenotypes_raw["ovarian_cancer"] is None

    def test_na_still_means_no_data(self):
        _, people = _people()
        assert "breast_cancer" not in people["E"].not_applicable  # "NA"
        assert "breast_cancer" not in people["D"].not_applicable  # -9

    def test_positional_header(self):
        """Columns 11+ of a standard ten-column header are read too."""
        text = (
            "FID\tIID\tPAT\tMAT\tSEX\tPHENO\tDOB\tTWIN\tDEAD\tPROBAND\tOC\n"
            "FAM1\tA\t0\t0\t1\t1\t0\t0\t0\t0\tN/A\n"
            "FAM1\tB\t0\t0\t2\t1\t0\t0\t0\t0\t2\n"
        )
        _, people = _people(text)
        assert people["A"].not_applicable == {"OC"}

    def test_auto_column_stays_binary(self):
        """An N/A cell is not a category: the column is still binary."""
        text = CANCERS.replace("b_ovarian_cancer", "ovarian_cancer")
        ped, _ = _people(text)
        individuals = list(ped.individuals())
        assert detect_auto_mode("ovarian_cancer", individuals) == "binary"

    def test_gradient_range_ignores_na(self):
        text = (
            "FID\tIID\tPAT\tMAT\tSEX\tPHENO\tg_score\n"
            "FAM1\tA\t0\t0\t1\t1\t10\n"
            "FAM1\tB\t0\t0\t2\t1\t20\n"
            "FAM1\tC\tA\tB\t1\t1\tN/A\n"
        )
        _, people = _people(text)
        assert people["A"].extra_phenotypes["score"] == 0.0
        assert people["B"].extra_phenotypes["score"] == 1.0
        assert people["C"].not_applicable == {"score"}

    def test_layout_and_builder_keep_it(self):
        ped, _ = _people()
        nodes = {n.individual_id: n for n in compute_layout(ped.families["FAM1"]).nodes}
        assert nodes["A"].not_applicable == {"ovarian_cancer"}
        b = PedigreeBuilder("F")
        b.add_individual("m", sex="male", phenotypes={"OC": "N/A"})
        b.add_individual("f", sex="female", phenotypes={"OC": 2})
        built = b.build().families["F"].individuals
        assert built["m"].not_applicable == {"OC"}
        assert not built["f"].not_applicable

    def test_dataframe(self):
        pd = pytest.importorskip("pandas")
        from pypedigree import parse_pedigree_from_dataframe

        df = pd.DataFrame(
            {
                "FID": ["FAM1", "FAM1"],
                "IID": ["A", "B"],
                "PAT": ["0", "0"],
                "MAT": ["0", "0"],
                "SEX": [1, 2],
                "PHENO": [1, 1],
                "b_OC": ["N/A", "2"],
            }
        )
        people = parse_pedigree_from_dataframe(df).families["FAM1"].individuals
        assert people["A"].not_applicable == {"OC"}


class TestDecisions:
    def test_band_and_quadrant_cells(self):
        ped, people = _people()
        style, render = _render(ped)
        modes = ped.phenotype_display_modes
        for cell_of in (band_cell, quadrant_cell):
            cell = cell_of(people["A"], "ovarian_cancer", style, render, modes)
            assert cell.not_applicable and cell.text is None and not cell.is_na
            assert to_hex(cell.color) == "#ffffff"
            other = cell_of(people["A"], "prostate_cancer", style, render, modes)
            assert not other.not_applicable

    def test_one_centred_na_only_with_missing_data(self):
        _, people = _people()
        cols = ["breast_cancer", "ovarian_cancer", "prostate_cancer"]
        assert all_bands_na(people["D"], cols), "no data and N/A: one 'NA'"
        assert all_bands_na(people["E"], cols)
        assert not all_bands_na(people["E"], ["ovarian_cancer", "prostate_cancer"]), "all N/A"
        assert not all_bands_na(people["A"], cols)

    def test_stroke(self):
        assert not_applicable_stroke(0, 0, 2, 1) == ((0, 0), (2, 1))
        assert not_applicable_stroke(0, 0, 2, 1, up=-1.0) == ((0, 1), (2, 0)), "SVG: same /"

    def test_legend_row_only_when_needed(self):
        ped, _ = _people()
        style, render = _render(ped)
        assert render.not_applicable
        rows = legend_rows(
            ped.phenotype_names, style, render, ped.phenotype_display_modes, CarrierScheme()
        )
        assert [r.label for r in rows if r.kind == "not_applicable"] == ["Not applicable"]
        plain = parse_pedigree_from_string(_without_na())
        style, render = _render(plain)
        assert not render.not_applicable

    def test_legend_text(self):
        ped, _ = _people()
        style = resolve_style(
            PedigreeStyle(legend_text={"not_applicable": "Does not apply"}), ped.phenotype_names
        )
        check_legend(style)
        render = resolve_phenotype_render(
            ped, ped.phenotype_names, ped.phenotype_display_modes, style
        )
        rows = legend_rows(
            ped.phenotype_names, style, render, ped.phenotype_display_modes, CarrierScheme()
        )
        assert "Does not apply" in [r.label for r in rows]


class TestDrawing:
    @pytest.mark.parametrize("layout", ["bands", "quadrants"])
    def test_svg(self, layout):
        svg = render_interactive_svg(
            parse_pedigree_from_string(CANCERS), style=PedigreeStyle(phenotype_layout=layout)
        )
        assert svg.count('class="not-applicable"') == NA_CELLS + 1  # and the legend row
        assert "ovarian cancer: not applicable" in svg
        assert "Not applicable" in svg

    @pytest.mark.parametrize("layout", ["bands", "quadrants"])
    def test_png(self, layout):
        fig = render_pedigree(
            parse_pedigree_from_string(CANCERS), style=PedigreeStyle(phenotype_layout=layout)
        )
        try:
            strokes = [
                line
                for ax in fig.axes
                for line in ax.lines
                if to_hex(line.get_color()) == NOT_APPLICABLE_STROKE
            ]
            assert len(strokes) == NA_CELLS + 1
            texts = [t.get_text() for ax in fig.axes for t in ax.texts]
            assert "Not applicable" in texts
            assert texts.count("NA") == 2, "D and E: one centred 'NA' each"
        finally:
            plt.close(fig)

    def test_no_stroke_without_na(self):
        svg = render_interactive_svg(parse_pedigree_from_string(_without_na()))
        assert 'class="not-applicable"' not in svg
        assert "Not applicable" not in svg
