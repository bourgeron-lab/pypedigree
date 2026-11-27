"""Normally distributed traits drawn on a curve, as a PGS is (R6).

A ``z_`` column (IQ...) gets a curve box of its own beside the symbol, with
its values made into z-scores by ``curve_reference`` and written as they are.
A kind of curve is a caption group: PGS, pTDT, IQ... An individual gets one
box per kind it has a value of, side by side outward from the symbol, and the
legend a row per kind.
"""

from __future__ import annotations

import json
import re
import sys
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import pytest
from click.testing import CliRunner
from matplotlib.text import Text

from pypedigree import (
    PedigreeBuilder,
    PedigreeStyle,
    parse_pedigree,
    parse_pedigree_from_string,
    render_interactive_svg,
    render_pedigree,
)
from pypedigree.carriers import CarrierScheme
from pypedigree.cli import main
from pypedigree.config import (
    apply_pgs_spacing,
    curve_kind,
    is_curve_column,
    is_pgs_column,
    is_ptdt_column,
    resolve_style,
)
from pypedigree.errors import PedigreeDataWarning, PedigreeError
from pypedigree.models import PhenotypeDisplayMode
from pypedigree.phenotypes import (
    NBSP,
    active_curves,
    box_tail_positions,
    box_width,
    check_curve_options,
    column_scale,
    curve_spacing,
    curve_value_text,
    curve_zscore,
    legend_rows,
    pgs_box_left,
    pgs_marker_values,
    pgs_position,
    pgs_tail_positions,
    resolve_phenotype_render,
)

sys.path.insert(0, str(Path(__file__).parent))
from conftest import TEST_DATA_DIR  # noqa: E402

FIXTURE = TEST_DATA_DIR / "phenotype" / "curve_traits.ped"
STYLE = json.loads((TEST_DATA_DIR / "phenotype" / "curve_traits.json").read_text("utf-8"))
IQ = {"IQ": {"mean": 100, "sd": 15}}


def _resolved(pedigree, **style):
    return resolve_style(
        PedigreeStyle(**style), pedigree.phenotype_names, pedigree.phenotype_display_modes
    )


def _individual(pedigree, sample):
    return next(i for i in pedigree.individuals() if i.sample_id == sample)


class TestColumns:
    def test_z_prefix_is_a_curve_column(self):
        pedigree = parse_pedigree(FIXTURE)
        assert "IQ" in pedigree.phenotype_names
        assert pedigree.phenotype_display_modes["IQ"] == PhenotypeDisplayMode.CURVE
        style = _resolved(pedigree)
        assert curve_kind("IQ", style) == "IQ"
        render = resolve_phenotype_render(
            pedigree, pedigree.phenotype_names, pedigree.phenotype_display_modes, style
        )
        assert render.curve_names(pedigree.phenotype_names) == ["PGS_autism", "PGS_ADHD", "IQ"]
        assert render.band_names(pedigree.phenotype_names) == ["autism"]

    def test_pgs_and_ptdt_kinds_by_name(self):
        assert curve_kind("PGS_autism") == "PGS" and curve_kind("pTDT_autism") == "pTDT"
        assert curve_kind("IQ") is None and not is_curve_column("IQ")
        # The 1.0 wrappers keep their meaning.
        assert is_pgs_column("PGS_IQ") and not is_pgs_column("IQ")
        assert is_ptdt_column("pTDT_autism")

    def test_a_reference_or_a_caption_makes_a_curve(self):
        assert curve_kind("IQ", PedigreeStyle(curve_reference=IQ)) == "IQ"
        labelled = PedigreeStyle(curve_reference=IQ, phenotype_labels={"IQ": "Full-scale IQ"})
        assert curve_kind("IQ", labelled) == "Full-scale IQ"
        assert curve_kind("SRS", PedigreeStyle(curve_captions={"SRS": "SRS"})) == "SRS"
        assert curve_kind("PGS_a", PedigreeStyle(curve_captions={"PGS_a": "PRS"})) == "PRS"

    def test_resolution_is_pure_and_idempotent(self):
        pedigree = parse_pedigree(FIXTURE)
        style = PedigreeStyle()
        once = _resolved(pedigree)
        assert style.curve_captions == {}
        assert once.curve_captions == {"IQ": "IQ"}
        again = resolve_style(once, pedigree.phenotype_names, pedigree.phenotype_display_modes)
        assert again == once

    def test_builder(self):
        builder = PedigreeBuilder()
        builder.add_individual("a", sex=1, phenotypes={"z_IQ": 112, "PGS_x": 0.5})
        builder.add_individual("b", sex=2, phenotypes={"z_IQ": 95, "PGS_x": -0.5})
        builder.add_union("a", "b")
        svg = builder.render_interactive_svg(style=PedigreeStyle(curve_reference=IQ))
        assert len(re.findall(r">IQ</text>", svg)) == 2
        assert "IQ: 112 (z = +0.80, 79th percentile)" in svg


class TestScale:
    def test_reference_makes_z_scores(self):
        style = PedigreeStyle(curve_reference=IQ)
        assert curve_zscore("IQ", 115, style) == pytest.approx(1.0)
        assert column_scale("IQ", style) == "zscore"

    def test_without_reference_values_are_z_scores_whatever_pgs_scale(self):
        style = PedigreeStyle(pgs_scale="minmax", curve_captions={"IQ": "IQ"})
        assert column_scale("IQ", style) == "zscore"
        assert curve_zscore("IQ", 1.5, style) == 1.5
        assert column_scale("PGS_a", style) == "minmax"

    def test_positions(self):
        pedigree = parse_pedigree(FIXTURE)
        style = _resolved(pedigree, curve_reference=IQ)
        [(kind, [(name, position)])] = active_curves(
            _individual(pedigree, "II-3"), ["IQ"], style
        )  # 133
        assert (kind, name) == ("IQ", "IQ")
        assert position == pytest.approx(pgs_position(33 / 15))

    def test_tails_per_box(self):
        minmax = PedigreeStyle(pgs_scale="minmax", curve_reference=IQ)
        assert box_tail_positions(["PGS_a"], minmax) == ()
        assert box_tail_positions(["IQ"], minmax) == pgs_tail_positions("zscore")
        # A pTDT box has its ±2 SD lines on "minmax" too, as its legend row has.
        assert box_tail_positions(["pTDT_a"], minmax) == pgs_tail_positions("zscore")


class TestBoxes:
    def test_one_box_per_kind_in_column_order(self):
        pedigree = parse_pedigree(FIXTURE)
        style = _resolved(pedigree, curve_reference=IQ)
        names = ["PGS_autism", "PGS_ADHD", "IQ"]
        boxes = active_curves(_individual(pedigree, "II-1"), names, style)
        assert [(k, [n for n, _ in m]) for k, m in boxes] == [
            ("PGS", ["PGS_autism", "PGS_ADHD"]),
            ("IQ", ["IQ"]),
        ]
        assert [k for k, _m in active_curves(_individual(pedigree, "I-1"), names, style)] == ["PGS"]
        assert [k for k, _m in active_curves(_individual(pedigree, "III-2"), names, style)] == [
            "IQ"
        ]

    def test_one_caption_one_box(self):
        text = (
            "FID\tIID\tFather\tMother\tSex\tPhenotype\tz_IQ_verbal\tz_IQ_performance\n"
            "F\ta\t0\t0\t1\t1\t110\t95\n"
        )
        pedigree = parse_pedigree_from_string(text)
        captions = {"IQ_verbal": "IQ", "IQ_performance": "IQ"}
        reference = {n: {"mean": 100, "sd": 15} for n in captions}
        style = _resolved(pedigree, curve_captions=captions, curve_reference=reference)
        [(kind, markers)] = active_curves(next(pedigree.individuals()), list(captions), style)
        assert kind == "IQ" and [n for n, _ in markers] == list(captions)
        svg = render_interactive_svg(
            pedigree, style=PedigreeStyle(curve_captions=captions, curve_reference=reference)
        )
        assert svg.count('<g class="pgs">') == 1
        assert svg.count('class="pgs-marker"') == 2

    def test_spacing_unchanged_when_nobody_mixes_kinds(self):
        text = (
            "FID\tIID\tFather\tMother\tSex\tPhenotype\tPGS_a\tz_IQ\n"
            "F\ta\t0\t0\t1\t1\t0.5\tNA\n"
            "F\tb\t0\t0\t2\t1\tNA\t110\n"
        )
        pedigree = parse_pedigree_from_string(text)
        style = _resolved(pedigree, curve_reference=IQ)
        names = pedigree.phenotype_names
        spacing = curve_spacing(pedigree.individuals(), names, style)
        assert spacing["boxes"] == 1 and spacing["box_gap"] == 0.0
        alone = apply_pgs_spacing(style, ["PGS_a"])
        assert apply_pgs_spacing(style, names, **spacing).sibling_gap == alone.sibling_gap

    def test_room_for_the_most_kinds(self):
        pedigree = parse_pedigree(FIXTURE)
        style = _resolved(pedigree, curve_reference=IQ)
        names = pedigree.phenotype_names
        spacing = curve_spacing(pedigree.individuals(), names, style)
        assert spacing["boxes"] == 2 and spacing["box_gap"] > 0
        one = apply_pgs_spacing(style, names, **{**spacing, "boxes": 1})
        two = apply_pgs_spacing(style, names, **spacing)
        width = spacing["box_width"][0]
        assert two.pgs_box_room - one.pgs_box_room == pytest.approx(width + spacing["box_gap"])
        assert two.pgs_box_gap == spacing["box_gap"] and one.pgs_box_gap == 0.0

    @pytest.mark.parametrize("show_values", [False, True])
    def test_boxes_clear_of_each_other_and_of_neighbours(self, show_values):
        """Boxes in a row outward, apart; none reaching another symbol or box."""
        from pypedigree import compute_layout

        pedigree = parse_pedigree(FIXTURE)
        style = _resolved(pedigree, curve_reference=IQ, show_pgs_values=show_values)
        names = pedigree.phenotype_names
        style = apply_pgs_spacing(
            style, names, **curve_spacing(pedigree.individuals(), names, style)
        )
        layout = compute_layout(pedigree.families["FAM1"], style)
        half = style.shape_size / 2
        w = box_width(style)
        extents = []
        for node in layout.nodes:
            boxes = active_curves(node, names, style)
            lefts = [pgs_box_left(node, style, index=i) for i in range(len(boxes))]
            outward = lefts if node.pgs_side == "right" else [-x for x in lefts]
            assert outward == sorted(outward), "each box further out than the last"
            for a, b in zip(lefts, lefts[1:], strict=False):
                assert abs(a - b) >= w + style.pgs_box_gap - 1e-9
            xs = [node.x - half, node.x + half, *lefts, *(x + w for x in lefts)]
            extents.append((node.y, min(xs), max(xs)))
        for i, (y, x0, x1) in enumerate(extents):
            for y2, u0, u1 in extents[i + 1 :]:
                if abs(y - y2) < 1e-6:
                    assert x1 <= u0 or u1 <= x0, "two individuals' boxes or symbols overlap"


class TestText:
    def test_tooltip(self):
        style = PedigreeStyle(curve_reference=IQ)
        assert curve_value_text("IQ", 115, style) == "115 (z = +1.00, 84th percentile)"
        assert curve_value_text("IQ", 85, style) == "85 (z = -1.00, 16th percentile)"
        z = PedigreeStyle(curve_captions={"IQ": "IQ"})
        assert curve_value_text("IQ", 1.0, z) == "z = +1.00, 84th percentile"

    def test_values_written_raw(self):
        pedigree = parse_pedigree(FIXTURE)
        style = _resolved(pedigree, curve_reference=IQ)
        child = _individual(pedigree, "III-1")
        assert pgs_marker_values(child, ["IQ", "PGS_autism"], style) == {
            "IQ": "96",
            "PGS_autism": "+0.8",
        }

    def test_svg(self):
        svg = render_interactive_svg(parse_pedigree(FIXTURE), style=PedigreeStyle(**STYLE))
        # Five people with IQ, five with PGS, three with both.
        assert len(re.findall(r">IQ</text>", svg)) == 6
        assert len(re.findall(r">PGS</text>", svg)) == 5
        assert svg.count('<g class="pgs">') == 11
        assert "IQ: 133 (z = +2.20, 99th percentile)" in svg
        assert re.search(r'class="pgs-value"[^>]*>111<', svg)

    def test_png(self):
        fig = render_pedigree(parse_pedigree(FIXTURE), style=PedigreeStyle(**STYLE))
        try:
            texts = [t.get_text() for t in fig.findobj(Text)]
        finally:
            plt.close(fig)
        assert texts.count("IQ") == 6 and texts.count("PGS") == 5
        assert "111" in texts and "133" in texts
        # Wrapped before its reference, which stays on one line (#68).
        assert f"(mean{NBSP}100,{NBSP}SD{NBSP}15)" in texts


class TestLegend:
    def _rows(self, pedigree, **style):
        names = pedigree.phenotype_names
        modes = pedigree.phenotype_display_modes
        st = resolve_style(PedigreeStyle(**style), names, modes)
        render = resolve_phenotype_render(pedigree, names, modes, st)
        return [
            r for r in legend_rows(names, st, render, modes, CarrierScheme()) if r.kind == "pgs"
        ]

    def test_a_row_per_kind(self):
        pgs, iq = self._rows(parse_pedigree(FIXTURE), curve_reference=IQ)
        assert pgs.label == "PGS distribution, z-score:"
        assert [name for name, _c in pgs.key] == ["autism", "ADHD"]
        assert iq.label == f"IQ distribution (mean{NBSP}100,{NBSP}SD{NBSP}15)"
        assert iq.key == () and [p for p, _t in iq.tails] == list(pgs_tail_positions("zscore"))

    def test_without_reference(self):
        with pytest.warns(PedigreeDataWarning, match="Curve column"):
            pass_through = render_pedigree(parse_pedigree(FIXTURE))
        plt.close(pass_through)
        _pgs, iq = self._rows(parse_pedigree(FIXTURE))
        assert iq.label == "IQ distribution, z-score"

    def test_grouped_columns_keep_their_key(self):
        text = (
            "FID\tIID\tFather\tMother\tSex\tPhenotype\tz_IQ_verbal\tz_IQ_performance\n"
            "F\ta\t0\t0\t1\t1\t0.5\t-0.5\n"
        )
        captions = {"IQ_verbal": "IQ", "IQ_performance": "IQ"}
        [row] = self._rows(parse_pedigree_from_string(text), curve_captions=captions)
        assert row.label == "IQ distribution, z-score:"
        assert [n for n, _c in row.key] == ["IQ verbal", "IQ performance"]
        assert len({c for _n, c in row.key}) == 2, "distinct colors within a kind"


class TestChecks:
    @pytest.mark.parametrize(
        "reference",
        [
            {"IQ": {"mean": 100}},
            {"IQ": {"mean": 100, "sd": 0}},
            {"IQ": {"mean": "100", "sd": 15}},
            {"IQ": {"mean": 100, "sd": 15, "skew": 1}},
            {"IQ": [100, 15]},
        ],
    )
    def test_invalid_reference(self, reference):
        with pytest.raises(PedigreeError, match="curve_reference for 'IQ'"):
            check_curve_options(PedigreeStyle(curve_reference=reference))

    def test_invalid_caption(self):
        with pytest.raises(PedigreeError, match="curve_captions for 'IQ'"):
            check_curve_options(PedigreeStyle(curve_captions={"IQ": ""}))

    def test_checked_before_drawing(self):
        bad = PedigreeStyle(curve_reference={"IQ": {"mean": 100, "sd": -1}})
        with pytest.raises(PedigreeError):
            render_interactive_svg(parse_pedigree(FIXTURE), style=bad)
        with pytest.raises(PedigreeError):
            render_pedigree(parse_pedigree(FIXTURE), style=bad)

    def test_raw_values_without_reference_warn(self):
        with pytest.warns(PedigreeDataWarning, match=r"Curve column\(s\) IQ \(a value of 133\)"):
            render_interactive_svg(parse_pedigree(FIXTURE))

    def test_no_warning_with_reference(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error", PedigreeDataWarning)
            render_interactive_svg(parse_pedigree(FIXTURE), style=PedigreeStyle(**STYLE))


def test_cli(tmp_path):
    out = tmp_path / "curve.svg"
    style = TEST_DATA_DIR / "phenotype" / "curve_traits.json"
    result = CliRunner().invoke(
        main, [str(FIXTURE), "--style", str(style), "--interactive", "-o", str(out)]
    )
    assert result.exit_code == 0, result.output
    assert "Warning" not in result.output
    assert ">IQ</text>" in out.read_text("utf-8")
