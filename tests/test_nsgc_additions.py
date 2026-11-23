"""Three NSGC notations added in #104: several individuals in one symbol, the
asymptomatic-carrier line and a couple's childlessness. All pedigrees here are
synthetic.
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import pytest

from pypedigree import (
    PedigreeBuilder,
    PedigreeStyle,
    export_fam,
    parse_pedigree,
    parse_pedigree_from_string,
    render_interactive_svg,
    render_pedigree,
    symbols,
)
from pypedigree.carriers import CarrierScheme
from pypedigree.config import resolve_style
from pypedigree.errors import PedigreeDataWarning, PedigreeParseError, PedigreeStyleWarning
from pypedigree.layout import compute_layout, compute_layout_circular
from pypedigree.models import Childlessness, EdgeType
from pypedigree.phenotypes import resolve_phenotype_render, with_quadrants

FIXTURE = Path(__file__).parent / "data" / "nomenclature" / "group_carrier_couple.ped"
HEADER = "FID\tIID\tFather\tMother\tSex\tPhenotype\tcount\tasymptomatic\tchildless\n"


def _ped(*rows: str):
    return parse_pedigree_from_string(HEADER + "\n".join(rows) + "\n")


def _people():
    ped = parse_pedigree(FIXTURE)
    return ped, ped.families["FAM1"]


class TestSeveralIndividuals:
    @pytest.mark.parametrize(("value", "count"), [("5", "5"), ("n", "n"), ("N", "n"), ("12", "12")])
    def test_values(self, value, count):
        ped = _ped(f"FAM1\ta\t0\t0\t0\t1\t{value}\t0\t0")
        assert ped.families["FAM1"].individuals["a"].notation.count == count

    @pytest.mark.parametrize("value", ["1", "many", "2.5"])
    def test_other_values_make_a_phenotype_column(self, value):
        """#118: a column named count that does not hold counts is a phenotype, as in 0.9.1."""
        with pytest.warns(PedigreeDataWarning, match=f"line 2 has '{re.escape(value)}'"):
            ped = _ped(f"FAM1\ta\t0\t0\t0\t1\t{value}\t0\t0")
        assert "count" in ped.phenotype_names
        assert ped.families["FAM1"].individuals["a"].notation.count is None

    def test_a_parent_makes_a_phenotype_column(self):
        """#118: counts on parents (4 / 2 / 3) are a quantitative trait, not groups."""
        with pytest.warns(PedigreeDataWarning, match="g would stand for several individuals"):
            ped = _ped(
                "FAM1\tg\t0\t0\t1\t1\t3\t0\t0",
                "FAM1\tm\t0\t0\t2\t1\t0\t0\t0",
                "FAM1\tc\tg\tm\t1\t1\t0\t0\t0",
            )
        assert ped.phenotype_names == ["count"]

    def test_number_column_of_0_9_1_draws(self):
        """The file of #118: drawn in 0.9.1, an error in 0.10.0."""
        text = (
            "FID\tIID\tPAT\tMAT\tSEX\tPHENO\tnumber\n"
            "F\tf\t0\t0\t1\t1\t1\nF\tm\t0\t0\t2\t1\t0\nF\tc\tf\tm\t1\t2\t3\n"
        )
        with pytest.warns(PedigreeDataWarning, match="'number' is read as a phenotype"):
            ped = parse_pedigree_from_string(text)
        assert ped.phenotype_names == ["number"]
        assert 'class="node"' in render_interactive_svg(ped)

    def test_builder_stays_strict(self):
        with pytest.raises(PedigreeParseError, match="count must be 2 or more"):
            PedigreeBuilder("F").add_individual("a", count=1).build()

    def test_cannot_be_a_parent(self):
        """Through the builder, a group given as a parent is still an error."""
        b = PedigreeBuilder("F")
        b.add_individual("g", sex="male", count=3)
        b.add_individual("m", sex="female")
        b.add_individual("c", sex="male", father="g", mother="m")
        with pytest.raises(PedigreeParseError, match="cannot be the parent of c"):
            b.build()

    def test_written_inside_tooltip_and_export(self, tmp_path):
        ped, _ = _people()
        assert symbols.inner_letter(ped.families["FAM1"].individuals["II-3"].notation) == "5"
        svg = render_interactive_svg(ped)
        assert 'class="symbol-letter">5<' in svg and "5 individuals" in svg
        assert 'class="symbol-letter">n<' in svg
        with pytest.warns(PedigreeDataWarning, match="II-3, II-4 stand"):
            export_fam(ped, tmp_path / "out.fam")
        assert "II-3" not in (tmp_path / "out.fam").read_text(encoding="utf-8")

    def test_png_number_at_the_centre(self):
        ped, family = _people()
        nodes = {n.individual_id: n for n in compute_layout(family).nodes}
        fig = render_pedigree(ped, style=PedigreeStyle(phenotype_layout="bands"))
        try:
            texts = {t.get_text(): t for t in fig.axes[0].texts}
            assert "5" in texts and "n" in texts
        finally:
            plt.close(fig)
        assert nodes["II-3"].notation.count == "5"


class TestAsymptomatic:
    @pytest.mark.parametrize("value", ["1", "yes", "TRUE"])
    def test_values_and_alias(self, value):
        ped = parse_pedigree_from_string(
            f"FID\tIID\tFather\tMother\tSex\tPhenotype\tpresymptomatic\nFAM1\ta\t0\t0\t2\t1\t{value}\n"
        )
        assert ped.families["FAM1"].individuals["a"].notation.is_asymptomatic

    @pytest.mark.parametrize(("name", "value"), [("asymptomatic", "2"), ("presymptomatic", "3.5")])
    def test_other_values_make_a_phenotype_column(self, name, value):
        """#118: a 1/2-coded or quantitative column of that name stays a phenotype."""
        text = (
            f"FID\tIID\tFather\tMother\tSex\tPhenotype\t{name}\n"
            f"FAM1\ta\t0\t0\t2\t1\t1\nFAM1\tb\t0\t0\t1\t1\t{value}\n"
        )
        with pytest.warns(PedigreeDataWarning, match=f"line 3 has '{value}'"):
            ped = parse_pedigree_from_string(text)
        assert ped.phenotype_names == [name]
        assert not ped.families["FAM1"].individuals["a"].notation.is_asymptomatic

    def test_line_through_and_past_the_symbol(self):
        (x1, y1), (x2, y2) = symbols.asymptomatic_line(0.0, 0.0, 0.4, symbols.Sex.FEMALE)
        assert x1 == x2 == 0.0
        assert -y1 == y2 == pytest.approx(0.4 * (1 + symbols.ASYMPTOMATIC_OVERSHOOT))

    def test_svg_lines_legend_tooltip(self):
        svg = render_interactive_svg(parse_pedigree(FIXTURE))
        assert svg.count('class="asymptomatic"') == 2 + 1  # II-1, II-6 and the legend row
        assert "Asymptomatic or presymptomatic carrier" in svg

    def test_auto_layout_falls_back_to_bands(self):
        ped = parse_pedigree(FIXTURE)
        names = ped.phenotype_names
        for layout, quadrants in (("auto", False), ("quadrants", True)):
            style = resolve_style(PedigreeStyle(phenotype_layout=layout), names)
            render = resolve_phenotype_render(ped, names, ped.phenotype_display_modes, style)
            assert render.asymptomatic
            if quadrants:
                with pytest.warns(PedigreeStyleWarning, match="between quadrants"):
                    out = with_quadrants(render, render.band_names(names), CarrierScheme(), style)
            else:
                out = with_quadrants(render, render.band_names(names), CarrierScheme(), style)
            assert bool(out.quadrants) == quadrants

    def test_warning_with_carrier_sections(self):
        text = (
            "FID\tIID\tFather\tMother\tSex\tPhenotype\tasymptomatic\tcr_CFTR\n"
            "FAM1\ta\t0\t0\t2\t1\t1\t1\n"
        )
        with pytest.warns(PedigreeDataWarning, match="also have carrier sections"):
            render_interactive_svg(parse_pedigree_from_string(text))


class TestChildlessCouple:
    def test_union_without_children(self):
        _, family = _people()
        union = family.get_union_for_parents("II-1", "II-2")
        assert union is not None and union.children_ids == []
        assert union.childlessness == Childlessness.INFERTILITY
        assert family.get_union_for_parents("II-6", "II-5").childlessness == Childlessness.BY_CHOICE

    def test_declared_on_both_rows_one_union(self):
        ped = _ped(
            "FAM1\ta\t0\t0\t1\t1\t0\t0\tchoice:b",
            "FAM1\tb\t0\t0\t2\t1\t0\t0\tchoice:a",
        )
        unions = ped.families["FAM1"].unions
        assert len(unions) == 1 and unions[0].childlessness == Childlessness.BY_CHOICE

    def test_unknown_partner(self):
        with pytest.raises(PedigreeParseError, match="childless partner 'zz'"):
            _ped("FAM1\ta\t0\t0\t1\t1\t0\t0\tchoice:zz")

    def test_self_reference_is_an_error(self):
        """#122: infertility:s on s's own row was silently lost."""
        with pytest.raises(PedigreeParseError, match="'s' is the individual themselves"):
            parse_pedigree_from_string(
                "FID\tIID\tPAT\tMAT\tSEX\tPHENO\tchildless\nF\ts\t0\t0\t2\t1\tinfertility:s\n"
            )

    def test_different_kinds_warn_and_keep_the_first(self):
        """#122: choice on one row, infertility on the other."""
        with pytest.warns(PedigreeDataWarning, match="different kinds of childlessness"):
            ped = _ped(
                "FAM1\ta\t0\t0\t1\t1\t0\t0\tchoice:b",
                "FAM1\tb\t0\t0\t2\t1\t0\t0\tinfertility:a",
            )
        (union,) = ped.families["FAM1"].unions
        assert union.childlessness == Childlessness.BY_CHOICE

    def test_couple_with_children_warns(self):
        with pytest.warns(PedigreeDataWarning, match="declared without children"):
            ped = _ped(
                "FAM1\ta\t0\t0\t1\t1\t0\t0\tchoice:b",
                "FAM1\tb\t0\t0\t2\t1\t0\t0\t0",
                "FAM1\tc\ta\tb\t1\t1\t0\t0\t0",
            )
        assert ped.families["FAM1"].unions[0].childlessness is None

    def test_marker_hangs_from_the_middle_of_the_line(self):
        _, family = _people()
        layout = compute_layout(family)
        nodes = {n.individual_id: n for n in layout.nodes}
        mid = (nodes["II-1"].x + nodes["II-2"].x) / 2
        stems = [
            e
            for e in layout.edges
            if e.edge_type == EdgeType.CHILDLESS and e.start[0] == e.end[0] == pytest.approx(mid)
        ]
        assert stems and stems[0].start[1] == pytest.approx(nodes["II-1"].y)
        bars = [
            e
            for e in layout.edges
            if e.edge_type == EdgeType.CHILDLESS
            and e.start[1] == e.end[1]
            and abs(e.start[0] + e.end[0] - 2 * mid) < 1e-9
        ]
        assert len(bars) == 2, "infertility: two bars"
        # No marker under either partner, and no room kept below them for one.
        assert nodes["II-1"].marker_depth == 0.0 and nodes["II-2"].marker_depth == 0.0
        assert nodes["II-7"].marker_depth > 0, "a bare 'choice' keeps the individual marker"

    def test_circular_layout_draws_it(self):
        _, family = _people()
        layout = compute_layout_circular(family)
        assert len([e for e in layout.edges if e.edge_type == EdgeType.CHILDLESS]) >= 2 + 3

    def test_tooltip(self):
        svg = render_interactive_svg(parse_pedigree(FIXTURE))
        assert "Infertility (with II-2)" in svg


class TestBuilder:
    def test_keywords(self):
        b = PedigreeBuilder("F")
        b.add_individual("a", sex="male", count=4, asymptomatic=True)
        b.add_individual("m", sex="male")
        b.add_individual("w", sex="female", childless="choice:m")
        b.add_individual("x", sex="male")
        b.add_individual("y", sex="female")
        b.add_union("x", "y", childless="infertility")
        family = b.build().families["F"]
        a = family.individuals["a"].notation
        assert a.count == "4" and a.is_asymptomatic
        assert family.get_union_for_parents("m", "w").childlessness == Childlessness.BY_CHOICE
        assert family.get_union_for_parents("x", "y").childlessness == Childlessness.INFERTILITY
