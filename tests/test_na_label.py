"""The centred "NA" only on a symbol with no data, never on "unknown" or a category.

"NA" means no data (-9, NA, or no value at all); a 0 means "unknown", drawn grey
without text. `phenotypes.all_bands_na` used to count both as no data, and every
categorical value too (the parser keeps categories apart, with None as the
numeric value): a symbol whose diagnoses were all 0, and every symbol of a
figure with only categorical columns, got an "NA" over it.
"""

from __future__ import annotations

import re

import matplotlib.pyplot as plt
import pytest
from matplotlib.text import Text

from pypedigree import (
    PedigreeStyle,
    apply_phenotypes_to_pedigree,
    parse_pedigree_from_string,
    render_pedigree,
    render_svg,
)
from pypedigree.errors import PedigreeDataWarning
from pypedigree.phenotypes import all_bands_na, band_has_no_data

HEAD = "FID\tIID\tPAT\tMAT\tSEX\tPHENO\t"


def _family(columns: str, a: str, b: str, c: str) -> str:
    """A trio with the given values of the extra columns for A, B and C."""
    return f"{HEAD}{columns}\nF\tA\t0\t0\t1\t1\t{a}\nF\tB\t0\t0\t2\t1\t{b}\nF\tC\tA\tB\t1\t2\t{c}"


CASES = {
    # name: (file, number of "NA" written)
    "binary_all_unknown": (_family("b_x\tb_y", "0\t0", "2\t1", "2\t2"), 0),
    "binary_all_no_data": (_family("b_x\tb_y", "-9\t-9", "2\t1", "2\t2"), 1),
    "binary_unknown_and_no_data": (_family("b_x\tb_y", "0\t-9", "2\t1", "2\t2"), 1),
    "categorical_only": (_family("c_sub", "x", "y", "x"), 0),
    "categorical_missing": (_family("c_sub", "-9", "y", "x"), 1),
    "gradient_no_data": (_family("g_s\tb_x", "-9\t-9", "4.5\t1", "6\t2"), 1),
}


def _png_na_count(ped, style) -> int:
    fig = render_pedigree(ped, style=style, show_legend=False)
    try:
        fig.canvas.draw()
        return sum(1 for t in fig.findobj(Text) if t.get_text() == "NA" and t.get_visible())
    finally:
        plt.close(fig)


@pytest.mark.parametrize("layout", ["bands", "quadrants"])
@pytest.mark.parametrize("case", list(CASES))
def test_na_written_only_without_data(case, layout):
    text, expected = CASES[case]
    ped = parse_pedigree_from_string(text)
    style = PedigreeStyle(phenotype_layout=layout)
    svg = render_svg(ped, style=style, show_legend=False)
    assert len(re.findall(r">NA<", svg)) == expected
    assert _png_na_count(ped, style) == expected


def test_individual_missing_from_the_phenotype_file_is_no_data():
    """An individual a separate phenotype file does not list has no value."""
    ped = parse_pedigree_from_string(_family("b_x", "1", "1", "2"))
    with pytest.warns(PedigreeDataWarning):  # the extra row matches nobody
        apply_phenotypes_to_pedigree(ped, {"B": {"z": 2.0}, "C": {"z": 1.0}, "Q": {}}, ["z"])
    a = ped.get_family("F").individuals["A"]
    assert band_has_no_data(a, "z")
    assert all_bands_na(a, ["z"])


def test_rule_per_value():
    ped = parse_pedigree_from_string(_family("b_x\tg_y\tc_z", "0\t-9\t-9", "2\t3.5\tq", "1\t5\t0"))
    a, b, c = (ped.get_family("F").individuals[i] for i in "ABC")
    assert not band_has_no_data(a, "x")  # 0: unknown, a value
    assert band_has_no_data(a, "y")  # -9
    assert band_has_no_data(a, "z")  # missing category
    assert not band_has_no_data(b, "z")  # a category
    assert not band_has_no_data(c, "z")  # "0" is a category name in a c_ column
    assert not band_has_no_data(c, "x")
