"""The quadrant phenotype layout: one quadrant of the symbol per column.

With ``phenotype_layout="quadrants"`` (``--phenotype-layout quadrants``) a
symbol is cut in four and each phenotype column, or carrier source, fills
one quadrant. Autism, ADHD, intellectual disability and epilepsy take their
conventional quadrant and color; the legend starts with a key showing which
quadrant holds which column. The default layout stays horizontal bands.
"""

from __future__ import annotations

import re
import warnings
from dataclasses import replace

import matplotlib.pyplot as plt
import pytest
from click.testing import CliRunner
from conftest import TEST_DATA_DIR
from matplotlib.colors import to_hex

from pypedigree import PedigreeStyle, parse_pedigree, parse_pedigree_from_string, render_pedigree
from pypedigree.carriers import CarrierScheme
from pypedigree.cli import main
from pypedigree.config import resolve_style
from pypedigree.errors import PedigreeError, PedigreeStyleWarning
from pypedigree.phenotypes import (
    NA_FILL,
    QUADRANT_REACH,
    assign_quadrants,
    check_phenotype_layout,
    legend_rows,
    legend_slot_count,
    quadrant_box,
    quadrant_cell,
    quadrant_center,
    resolve_phenotype_render,
)
from pypedigree.svg_renderer import render_interactive_svg

FIXTURE = TEST_DATA_DIR / "phenotype" / "quadrants.ped"
QUADRANTS = PedigreeStyle(phenotype_layout="quadrants")

TRIO = (
    "FID\tIID\tFather\tMother\tSex\tb_Epilepsy\tb_ID\tb_ADHD\tb_Autism\n"
    "F\tdad\t0\t0\t1\t1\t1\t1\t2\n"  # autism only
    "F\tmum\t0\t0\t2\t2\t1\t2\t1\n"  # ADHD and epilepsy
    "F\tkid\tdad\tmum\t1\t-9\t2\t1\t1\n"  # ID, epilepsy unknown
)
COLORS = {"Autism": "#386eb2", "ID": "#d34436", "ADHD": "#e59635", "Epilepsy": "#5e287e"}


def _quadrant_fills(ax, node_x, node_y, half):
    """Position -> fill of the quadrant rectangles drawn for one node."""
    reach = half * QUADRANT_REACH
    fills = {}
    for patch in ax.patches:
        if not hasattr(patch, "get_width") or abs(patch.get_width() - reach) > 1e-9:
            continue
        x, y = patch.get_xy()
        on_x = abs(x - (node_x - reach)) < 1e-9 or abs(x - node_x) < 1e-9
        on_y = abs(y - node_y) < 1e-9 or abs(y - (node_y - reach)) < 1e-9
        if on_x and on_y:
            vertical = "top" if abs(y - node_y) < 1e-9 else "bottom"
            horizontal = "left" if x < node_x - 1e-9 else "right"
            fills[f"{vertical}-{horizontal}"] = to_hex(patch.get_facecolor())
    return fills


class TestAssignment:
    def test_known_phenotypes_take_their_quadrant_whatever_the_order(self):
        style = resolve_style(QUADRANTS, ["Epilepsy", "ID", "ADHD", "Autism"])
        assert assign_quadrants(["Epilepsy", "ID", "ADHD", "Autism"], style) == {
            "Epilepsy": "bottom-right",
            "ID": "bottom-left",
            "ADHD": "top-right",
            "Autism": "top-left",
        }

    def test_other_columns_fill_free_quadrants_in_reading_order(self):
        assignment = assign_quadrants(["IQ", "ADHD", "SRS"], QUADRANTS)
        assert assignment == {"IQ": "top-left", "ADHD": "top-right", "SRS": "bottom-left"}

    def test_explicit_quadrants_win(self):
        style = PedigreeStyle(
            phenotype_layout="quadrants", phenotype_quadrants={"ADHD": "bottom-right"}
        )
        assignment = assign_quadrants(["Autism", "ADHD", "Epilepsy"], style)
        # Epilepsy's own quadrant is taken, so it gets the first free one.
        assert assignment == {"Autism": "top-left", "ADHD": "bottom-right", "Epilepsy": "top-right"}

    def test_bands_layout_assigns_nothing(self):
        assert assign_quadrants(["Autism", "ADHD"], PedigreeStyle(phenotype_layout="bands")) == {}

    def test_default_is_auto_and_assigns_up_to_four(self):
        assert PedigreeStyle().phenotype_layout == "auto"
        assert assign_quadrants(["Autism", "ADHD"], PedigreeStyle()) == {
            "Autism": "top-left",
            "ADHD": "top-right",
        }

    def test_auto_falls_back_to_bands_beyond_four_without_a_warning(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            assert assign_quadrants(["A", "B", "C", "D", "E"], PedigreeStyle()) == {}

    def test_more_than_four_columns_fall_back_to_bands(self):
        with pytest.warns(PedigreeStyleWarning, match="at most 4 columns"):
            assert assign_quadrants(["A", "B", "C", "D", "E"], QUADRANTS) == {}

    @pytest.mark.parametrize(
        ("style", "message"),
        [
            (PedigreeStyle(phenotype_layout="pie"), "Unknown phenotype_layout"),
            (PedigreeStyle(phenotype_quadrants={"ASD": "top"}), "expected one of"),
            (
                PedigreeStyle(phenotype_quadrants={"ASD": "top-left", "ID": "top-left"}),
                "several columns in top-left",
            ),
        ],
    )
    def test_invalid_settings_are_rejected(self, style, message):
        with pytest.raises(PedigreeError, match=message):
            check_phenotype_layout(style)


class TestGeometry:
    @pytest.mark.parametrize(
        ("position", "up", "expected"),
        [
            ("top-left", 1.0, (-1.5, 0.0)),
            ("top-right", 1.0, (0.0, 0.0)),
            ("bottom-left", 1.0, (-1.5, -1.5)),
            ("bottom-right", 1.0, (0.0, -1.5)),
            # SVG, y pointing down: the top quadrants lie at smaller y.
            ("top-left", -1.0, (-1.5, -1.5)),
            ("bottom-right", -1.0, (0.0, 0.0)),
        ],
    )
    def test_quadrant_box(self, position, up, expected):
        x, y, w, h = quadrant_box(position, 0.0, 0.0, 1.0, up)
        assert (x, y) == pytest.approx(expected)
        assert (w, h) == pytest.approx((QUADRANT_REACH, QUADRANT_REACH))

    def test_text_sits_inside_its_quadrant(self):
        assert quadrant_center("top-right", 0.0, 0.0, 1.0) == pytest.approx((0.5, 0.5))
        assert quadrant_center("top-right", 0.0, 0.0, 1.0, up=-1.0) == pytest.approx((0.5, -0.5))
        assert quadrant_center("bottom-left", 0.0, 0.0, 1.0) == pytest.approx((-0.5, -0.5))

    def test_gradient_value_is_written_alone(self):
        text = "FID\tID\tFather\tMother\tSex\tg_IQ\nF\ta\t0\t0\t1\t103\nF\tb\t0\t0\t2\t87\n"
        pedigree = parse_pedigree_from_string(text)
        style = resolve_style(QUADRANTS, pedigree.phenotype_names)
        render = resolve_phenotype_render(
            pedigree, pedigree.phenotype_names, pedigree.phenotype_display_modes, style
        )
        a = pedigree.families["F"].individuals["a"]
        cell = quadrant_cell(a, "IQ", style, render, pedigree.phenotype_display_modes)
        assert cell.text == "103"


class TestRendering:
    def _layout(self, svg=False):
        pedigree = parse_pedigree_from_string(TRIO)
        if svg:
            return render_interactive_svg(pedigree, style=QUADRANTS)
        return render_pedigree(pedigree, style=QUADRANTS)

    def test_png_quadrants_follow_the_convention(self):
        fig = self._layout()
        try:
            ax = fig.axes[0]
            half = QUADRANTS.shape_size / 2
            nodes = {t.get_text(): t.get_position() for t in ax.texts}
            # Labels sit below their node: find each node by its quadrants.
            fills = []
            for patch in ax.patches:
                if abs(getattr(patch, "get_width", lambda: 0)() - half * QUADRANT_REACH) < 1e-9:
                    fills.append(to_hex(patch.get_facecolor()))
            assert set(nodes) >= {"dad", "mum", "kid"}
            assert len(fills) == 12  # three people, four quadrants each
            assert COLORS["Autism"] in fills and COLORS["ID"] in fills
            assert NA_FILL in fills  # kid's epilepsy is -9
        finally:
            plt.close(fig)

    def test_png_each_color_in_its_quadrant(self):
        from pypedigree.layout import compute_layout

        pedigree = parse_pedigree_from_string(TRIO)
        style = resolve_style(QUADRANTS, pedigree.phenotype_names)
        layout = compute_layout(pedigree.families["F"], style)
        fig = render_pedigree(pedigree, style=QUADRANTS)
        try:
            half = QUADRANTS.shape_size / 2
            fills = {
                n.individual_id: _quadrant_fills(fig.axes[0], n.x, n.y, half) for n in layout.nodes
            }
        finally:
            plt.close(fig)
        white = "#ffffff"
        assert fills["dad"] == {
            "top-left": COLORS["Autism"],
            "top-right": white,
            "bottom-left": white,
            "bottom-right": white,
        }
        assert fills["mum"] == {
            "top-left": white,
            "top-right": COLORS["ADHD"],
            "bottom-left": white,
            "bottom-right": COLORS["Epilepsy"],
        }
        assert fills["kid"] == {
            "top-left": white,
            "top-right": white,
            "bottom-left": COLORS["ID"],
            "bottom-right": NA_FILL,
        }

    def test_svg_draws_the_same_quadrants(self):
        svg = self._layout(svg=True)
        fills = re.findall(r'<rect class="band" [^>]*fill="([^"]+)"', svg)
        assert len(fills) == 12
        assert fills.count(COLORS["Autism"]) == 1  # dad
        assert fills.count(COLORS["ADHD"]) == 1  # mum
        assert fills.count(COLORS["Epilepsy"]) == 1  # mum
        assert fills.count(COLORS["ID"]) == 1  # kid
        assert fills.count(NA_FILL) == 1  # kid's epilepsy

    def test_bands_layout_draws_equal_bands(self):
        pedigree = parse_pedigree_from_string(TRIO)
        svg = render_interactive_svg(pedigree, style=PedigreeStyle(phenotype_layout="bands"))
        heights = set(re.findall(r'<rect class="band" [^>]*height="([^"]+)"', svg))
        assert len(heights) == 1, "equal horizontal bands"
        assert float(heights.pop()) == pytest.approx(PedigreeStyle().shape_size / 4, abs=1e-3)

    def test_carrier_source_takes_a_hatched_quadrant(self):
        text = (
            "FID\tIID\tFather\tMother\tSex\tb_Autism\tcr_HEXA\n"
            "F\ta\t0\t0\t1\t2\t1\n"
            "F\tb\t0\t0\t2\t1\t0\n"
        )
        fig = render_pedigree(parse_pedigree_from_string(text), style=QUADRANTS)
        try:
            half = QUADRANTS.shape_size / 2
            hatched = [
                p
                for p in fig.axes[0].patches
                if p.get_hatch() and abs(p.get_width() - half * QUADRANT_REACH) < 1e-9
            ]
            assert len(hatched) == 1
        finally:
            plt.close(fig)

    def test_circular_layout_draws_quadrants(self):
        """Shapes are larger in the circular layout; four equal quadrants each."""
        from collections import Counter

        fig = render_pedigree(parse_pedigree_from_string(TRIO), style=QUADRANTS, circular=True)
        try:
            widths = Counter(
                round(p.get_width(), 6)
                for p in fig.axes[0].patches
                if type(p).__name__ == "Rectangle"
            )
        finally:
            plt.close(fig)
        assert widths.most_common(1)[0][1] == 12

    def test_too_many_columns_warn_and_draw_bands(self):
        text = (
            "FID\tIID\tFather\tMother\tSex\tb_A\tb_B\tb_C\tb_D\tb_E\nF\ta\t0\t0\t1\t2\t1\t1\t1\t1\n"
        )
        with pytest.warns(PedigreeStyleWarning, match="drawn as bands"):
            svg = render_interactive_svg(parse_pedigree_from_string(text), style=QUADRANTS)
        assert len(re.findall(r'<rect class="band" ', svg)) == 5


class TestLegend:
    def test_key_comes_first_with_every_column(self):
        pedigree = parse_pedigree_from_string(TRIO)
        names = pedigree.phenotype_names
        style = resolve_style(QUADRANTS, names)
        render = resolve_phenotype_render(pedigree, names, pedigree.phenotype_display_modes, style)
        from pypedigree.phenotypes import with_quadrants

        render = with_quadrants(render, render.band_names(names), CarrierScheme(), style)
        rows = legend_rows(names, style, render, pedigree.phenotype_display_modes, CarrierScheme())
        assert rows[0].kind == "quadrants"
        assert {(c.position, c.label, c.color) for c in rows[0].quadrants} == {
            ("top-left", "Autism", COLORS["Autism"]),
            ("top-right", "ADHD", COLORS["ADHD"]),
            ("bottom-left", "ID", COLORS["ID"]),
            ("bottom-right", "Epilepsy", COLORS["Epilepsy"]),
        }
        # The key is two rows tall, then one shared row for the binary
        # columns, whose colors the key gives (compact legend, #56).
        assert legend_slot_count(rows) == 2 + 1
        assert rows[1].kind == "binary_shared"
        full = replace(style, legend="full")
        rows = legend_rows(names, full, render, pedigree.phenotype_display_modes, CarrierScheme())
        assert legend_slot_count(rows) == 2 + 4

    def test_png_and_svg_legends_name_the_quadrants(self):
        fig = self._fig()
        try:
            texts = {t.get_text() for t in fig.axes[1].texts}
        finally:
            plt.close(fig)
        assert {"Autism", "ADHD", "ID", "Epilepsy"} <= texts
        svg = render_interactive_svg(parse_pedigree_from_string(TRIO), style=QUADRANTS)
        assert 'clipPath id="legend-key-square"' in svg
        assert 'clipPath id="legend-key-circle"' in svg

    @pytest.mark.parametrize(
        "text",
        [
            # The structural carrier column (carrier.ped).
            None,
            # Two carrier sources and no phenotype column.
            "FID\tIID\tFather\tMother\tSex\tcr_HEXA\tcr_CFTR\n"
            "F\ta\t0\t0\t1\t1\t0\n"
            "F\tb\t0\t0\t2\t2\t1\n",
        ],
        ids=["carrier-column", "cr-columns"],
    )
    def test_carriers_alone_are_not_cut_in_quadrants(self, text):
        """With no phenotype column, the symbol and the legend stay as with bands.

        0.4.4 hatched the whole symbol but put a quadrant key in the legend.
        """
        if text is None:
            pedigree = parse_pedigree(str(TEST_DATA_DIR / "phenotype" / "carrier.ped"))
        else:
            pedigree = parse_pedigree_from_string(text)
        assert render_interactive_svg(pedigree, style=QUADRANTS) == render_interactive_svg(
            pedigree, style=PedigreeStyle()
        )
        quadrants = render_pedigree(pedigree, style=QUADRANTS)
        bands = render_pedigree(pedigree, style=PedigreeStyle())
        try:
            for fig in (quadrants, bands):
                fig.canvas.draw()
            assert quadrants.get_size_inches().tolist() == bands.get_size_inches().tolist()
            assert len(quadrants.axes[1].patches) == len(bands.axes[1].patches)
            assert [t.get_text() for t in quadrants.axes[1].texts] == [
                t.get_text() for t in bands.axes[1].texts
            ]
        finally:
            plt.close(quadrants)
            plt.close(bands)

    def _fig(self):
        return render_pedigree(parse_pedigree_from_string(TRIO), style=QUADRANTS)


class TestOptions:
    def test_style_json_keys(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            style = PedigreeStyle.from_dict(
                {"phenotype_layout": "quadrants", "phenotype_quadrants": {"IQ": "top-left"}}
            )
        assert style.phenotype_layout == "quadrants"
        assert PedigreeStyle().phenotype_layout == "auto"

    def test_cli_flag(self, tmp_path):
        runner = CliRunner()
        out = tmp_path / "q.svg"
        args = [str(FIXTURE), "--interactive", "--phenotype-layout", "quadrants", "-o", str(out)]
        result = runner.invoke(main, args)
        assert result.exit_code == 0, result.output
        assert 'clipPath id="legend-key-square"' in out.read_text(encoding="utf-8")

    def test_cli_rejects_an_unknown_layout(self, tmp_path):
        result = CliRunner().invoke(
            main, [str(FIXTURE), "--phenotype-layout", "pie", "-o", str(tmp_path / "x.png")]
        )
        assert result.exit_code != 0

    def test_fixture_renders(self):
        pedigree = parse_pedigree(FIXTURE)
        assert pedigree.phenotype_names == ["Autism", "ADHD", "ID", "Epilepsy", "WGS"]
        fig = render_pedigree(pedigree, style=QUADRANTS)
        plt.close(fig)


class TestOneSourceFillsTheSymbol:
    """#108: under "auto", one column fills the whole symbol, as NSGC draws one condition."""

    ONE = "FID\tIID\tPAT\tMAT\tSEX\tPHENO\tb_autism\nF\ta\t0\t0\t1\t1\t2\nF\tb\t0\t0\t2\t1\t1\n"

    def _quadrants(self, text, layout):
        from pypedigree import PedigreeStyle, parse_pedigree_from_string
        from pypedigree.carriers import resolve_carrier_scheme
        from pypedigree.config import resolve_style
        from pypedigree.phenotypes import resolve_phenotype_render, with_quadrants

        ped = parse_pedigree_from_string(text)
        names = ped.phenotype_names
        style = resolve_style(PedigreeStyle(phenotype_layout=layout), names)
        render = resolve_phenotype_render(ped, names, ped.phenotype_display_modes, style)
        modes = ped.phenotype_display_modes
        cr = [n for n in names if modes.get(n) is not None and modes[n].value == "carrier"]
        carriers = resolve_carrier_scheme(ped.individuals(), cr, style.carrier_notation)
        return with_quadrants(render, render.band_names(names), carriers, style).quadrants

    def test_auto_one_column_no_quadrant(self):
        assert self._quadrants(self.ONE, "auto") == {}
        assert self._quadrants(self.ONE, "quadrants") == {"autism": "top-left"}

    def test_a_carrier_column_makes_two_sources(self):
        text = (
            self.ONE.replace("b_autism\n", "b_autism\tcr_CFTR\n")
            .replace("\t2\n", "\t2\t1\n")
            .replace("\t1\n", "\t1\t0\n")
        )
        assert len(self._quadrants(text, "auto")) == 2
