"""Symbols are spaced by their labels' width, so long IDs do not overlap (#41)."""

from __future__ import annotations

import xml.etree.ElementTree as ET
from collections import defaultdict

import pytest

import pypedigree as pp
from pypedigree import symbols
from pypedigree.layout import compute_layout

# Generated IDs as long as genetics-viz's "dummy" members (16-17 characters).
LONG_IDS = (
    "FID\tIID\tFather\tMother\tSex\tPhenotype\n"
    "F2\tGRANDFATHER_0001\t0\t0\t1\t1\n"
    "F2\tGRANDMOTHER_0002\t0\t0\t2\t1\n"
    "F2\tDUMMY_S0001234_MF\tGRANDFATHER_0001\tGRANDMOTHER_0002\t1\t2\n"
    "F2\tSIBLING_S0005678\tGRANDFATHER_0001\tGRANDMOTHER_0002\t2\t1\n"
    "F2\tPARTNER_S0009012\t0\t0\t2\t1\n"
    "F2\tPROBAND_S0003456\tDUMMY_S0001234_MF\tPARTNER_S0009012\t1\t2\n"
)

# Two families side by side under one founder couple: long IDs meet
# across neighbouring sub-trees, not only within a sibship or a couple.
COUSINS = (
    "FID\tIID\tFather\tMother\tSex\tPhenotype\n"
    "F3\tg1\t0\t0\t1\t1\n"
    "F3\tg2\t0\t0\t2\t1\n"
    "F3\ta\tg1\tg2\t1\t1\n"
    "F3\tb\tg1\tg2\t1\t1\n"
    "F3\tpa\t0\t0\t2\t1\n"
    "F3\tpb\t0\t0\t2\t1\n"
    "F3\tCOUSIN_LEFT_000000001\ta\tpa\t1\t2\n"
    "F3\tCOUSIN_RIGHT_00000002\tb\tpb\t2\t1\n"
)


def _family(text: str):
    ped = pp.parse_pedigree_from_string(text)
    return next(iter(ped.families.values()))


def _rows(layout):
    rows = defaultdict(list)
    for node in layout.nodes:
        rows[round(node.y, 6)].append(node)
    return [sorted(row, key=lambda n: n.x) for row in rows.values()]


def _assert_no_overlap(layout, font_size: float) -> None:
    """Neighbours on a row are at least half of each label plus the gap apart."""
    pad = symbols.LABEL_GAP_EM * font_size * symbols.UNITS_PER_POINT
    for row in _rows(layout):
        for a, b in zip(row, row[1:], strict=False):
            wa = symbols.label_width(a.label, a.notes + a.annotations, font_size)
            wb = symbols.label_width(b.label, b.notes + b.annotations, font_size)
            need = (wa + wb + 2 * pad) / 2 if a.label and b.label else 0.0
            assert b.x - a.x >= need - 1e-9, (a.label, b.label, b.x - a.x, need)


class TestTextWidth:
    def test_table_matches_dejavu_sans(self):
        """The widths are those of matplotlib's DejaVu Sans, to the thousandth."""
        from matplotlib import font_manager, ft2font

        font = ft2font.FT2Font(font_manager.findfont("DejaVu Sans"))
        # matplotlib 3.10 moved the flag into an enum; 3.7 to 3.9 have the constant.
        if hasattr(ft2font, "LoadFlags"):
            flags = ft2font.LoadFlags.NO_SCALE
        else:
            flags = ft2font.LOAD_NO_SCALE
        for code in range(32, 127):
            glyph = font.load_char(code, flags=flags)
            assert symbols._CHAR_EM[chr(code)] == pytest.approx(
                glyph.horiAdvance / font.units_per_EM, abs=5e-4
            ), chr(code)

    def test_width_scales_with_font_size(self):
        assert symbols.text_width("ABC", 16) == pytest.approx(2 * symbols.text_width("ABC", 8))

    def test_narrow_and_wide_characters(self):
        assert symbols.text_width("iiii", 8) < symbols.text_width("1111", 8)
        assert symbols.text_width("1111", 8) < symbols.text_width("WWWW", 8)

    def test_accented_letter_counts_as_its_base(self):
        assert symbols.text_width("père", 8) == pytest.approx(symbols.text_width("pere", 8))

    def test_wide_characters_count_one_em(self):
        assert symbols.text_width("家", 10) == pytest.approx(10 * symbols.UNITS_PER_POINT)

    def test_lines_below_use_smaller_type(self):
        """Notes and annotations count at 70 % of the label's size."""
        line = "a much longer annotation line"
        width = symbols.label_width("A", [line], 8)
        assert width == pytest.approx(symbols.text_width(line, 8 * symbols.SMALL_TEXT_SCALE))

    def test_no_label_takes_no_room(self):
        assert symbols.label_width("", ["annotation"], 8) == 0.0


class TestLayout:
    @pytest.mark.parametrize("font_size", [8.0, 12.0])
    def test_long_ids_do_not_overlap(self, font_size):
        style = pp.PedigreeStyle(font_size=font_size)
        _assert_no_overlap(compute_layout(_family(LONG_IDS), style), font_size)

    @pytest.mark.parametrize("font_size", [8.0, 12.0])
    def test_across_neighbouring_subtrees(self, font_size):
        style = pp.PedigreeStyle(font_size=font_size)
        _assert_no_overlap(compute_layout(_family(COUSINS), style), font_size)

    def test_larger_font_spreads_further(self):
        family = _family(LONG_IDS)
        widths = []
        for size in (8.0, 12.0):
            xs = [n.x for n in compute_layout(family, pp.PedigreeStyle(font_size=size)).nodes]
            widths.append(max(xs) - min(xs))
        assert widths[1] > widths[0]

    def test_short_ids_keep_the_usual_gaps(self):
        """Labels that fit change nothing: partners stay `couple_gap` apart."""
        text = (
            "FID\tIID\tFather\tMother\tSex\tPhenotype\n"
            "F1\tf\t0\t0\t1\t1\nF1\tm\t0\t0\t2\t1\nF1\tc\tf\tm\t1\t2\n"
        )
        layout = compute_layout(_family(text), pp.PedigreeStyle())
        xs = {n.individual_id: n.x for n in layout.nodes}
        assert xs["m"] - xs["f"] == pytest.approx(pp.PedigreeStyle().couple_gap)

    def test_hidden_labels_take_no_room(self):
        family = _family(LONG_IDS)
        style = pp.PedigreeStyle(font_size=12.0)
        xs = {n.individual_id: n.x for n in compute_layout(family, style, labels=False).nodes}
        gap = xs["GRANDMOTHER_0002"] - xs["GRANDFATHER_0001"]
        assert gap == pytest.approx(style.couple_gap)

    def test_long_annotation_counts(self):
        text = (
            "FID\tIID\tFather\tMother\tSex\tPhenotype\ta_variant\n"
            "F1\tf\t0\t0\t1\t1\t\n"
            "F1\tm\t0\t0\t2\t1\t\n"
            "F1\tc1\tf\tm\t1\t2\tNM_000000.1:c.1234_1240delinsACGTACG\n"
            "F1\tc2\tf\tm\t2\t1\tNM_000000.1:c.4321_4327delinsTTTTGGG\n"
        )
        layout = compute_layout(_family(text), pp.PedigreeStyle())
        xs = {n.individual_id: n.x for n in layout.nodes}
        assert xs["c2"] - xs["c1"] > pp.PedigreeStyle().sibling_gap
        _assert_no_overlap(layout, 8.0)

    def test_couple_counts_both_partners_labels(self):
        """A long label beside a short one: the couple is as wide as both reaches.

        Only labels reaching past half the sibling gap need room among
        siblings, but a couple's gap adds the partners' reaches, so the
        shorter one must count too (it counted as nothing before).
        """
        text = (
            "FID\tIID\tFather\tMother\tSex\tPhenotype\ta_name\n"
            "F1\tf\t0\t0\t1\t1\tAn invented name here\n"
            "F1\tm\t0\t0\t2\t1\tA much longer invented annotation\n"
            "F1\tc\tf\tm\t1\t2\t\n"
        )
        style = pp.PedigreeStyle()
        size, small = style.text_size(), style.text_size(symbols.SMALL_TEXT_SCALE)
        pad = symbols.LABEL_GAP_EM * size * symbols.UNITS_PER_POINT

        def reach(label, note):
            return (symbols.label_width(label, [note], size, small) + pad) / 2

        short = reach("f", "An invented name here")
        long = reach("m", "A much longer invented annotation")
        # The case: one reach within half the sibling gap, the other beyond,
        # together wider than the couple gap.
        assert short <= style.sibling_gap / 2 < long
        assert short + long > style.couple_gap
        layout = compute_layout(_family(text), style)
        xs = {n.individual_id: n.x for n in layout.nodes}
        # With a child, each text also stays clear of the descent line (#130).
        assert xs["m"] - xs["f"] == pytest.approx(max(short + long, 2 * long))
        _assert_no_overlap(layout, 8.0)

    def test_proband_marker_room_still_kept(self):
        """The marker's room and the label's room combine (the larger wins)."""
        text = (
            "FID\tIID\tFather\tMother\tSex\tPhenotype\tProband\n"
            "F1\tf\t0\t0\t1\t1\t0\n"
            "F1\tm\t0\t0\t2\t1\t0\n"
            "F1\tSIBLING_LONG_ID_01\tf\tm\t1\t1\t0\n"
            "F1\tPROBAND_LONG_ID_02\tf\tm\t2\t2\t1\n"
        )
        style = pp.PedigreeStyle(font_size=12.0)
        layout = compute_layout(_family(text), style)
        _assert_no_overlap(layout, 12.0)
        xs = {n.individual_id: n.x for n in layout.nodes}
        clearance = symbols.marker_clearance(
            style.shape_size / 2, style.font_size, len("SIBLING_LONG_ID_01"), True
        )
        assert xs["PROBAND_LONG_ID_02"] - xs["SIBLING_LONG_ID_01"] >= clearance - 1e-9

    def test_builder_pedigrees_are_spaced_too(self):
        b = pp.PedigreeBuilder(family_id="F1")
        b.add_individual("FATHER_WITH_A_LONG_ID", sex="male")
        b.add_individual("MOTHER_WITH_A_LONG_ID", sex="female")
        b.add_individual(
            "CHILD", sex="male", father="FATHER_WITH_A_LONG_ID", mother="MOTHER_WITH_A_LONG_ID"
        )
        family = b.build().get_family("F1")
        _assert_no_overlap(compute_layout(family, pp.PedigreeStyle()), 8.0)


class TestSvgFrame:
    @pytest.mark.parametrize("font_size", [8.0, 12.0])
    def test_labels_stay_inside_the_viewbox(self, font_size):
        ped = pp.parse_pedigree_from_string(LONG_IDS)
        style = pp.PedigreeStyle(font_size=font_size)
        root = ET.fromstring(pp.render_svg(ped, style=style, show_legend=False))
        vb_x, _vb_y, vb_w, _vb_h = map(float, root.get("viewBox").split())
        layout = compute_layout(next(iter(ped.families.values())), style)
        for node in layout.nodes:
            half = symbols.label_width(node.label, [], font_size) / 2
            assert vb_x <= node.x - half
            assert node.x + half <= vb_x + vb_w


class TestTextClearOfTheDescentLine:
    """#130: a long line under one partner crossed the couple's descent line."""

    TEXT = (
        "FID\tIID\tPAT\tMAT\tSEX\tPHENO\ta_genes\n"
        "F\tfa\t0\t0\t1\t1\t0\n"
        "F\tmo\t0\t0\t2\t2\tBRCA1, BRCA2, PALB2, ATM, CHEK2 - (search)\n"
        "F\td1\tfa\tmo\t2\t1\t0\n"
        "F\td2\tfa\tmo\t2\t1\t0\n"
    )

    def test_each_partners_text_stays_on_its_side(self):
        from pypedigree import PedigreeStyle, compute_layout, parse_pedigree_from_string, symbols

        style = PedigreeStyle()
        family = parse_pedigree_from_string(self.TEXT).families["F"]
        nodes = {n.individual_id: n for n in compute_layout(family, style).nodes}
        mid = (nodes["fa"].x + nodes["mo"].x) / 2  # where the descent line starts
        for sid in ("fa", "mo"):
            node = nodes[sid]
            width = symbols.label_width(
                node.label,
                node.notes + node.annotations,
                style.text_size(),
                style.text_size(symbols.SMALL_TEXT_SCALE),
            )
            assert abs(node.x - mid) > width / 2, sid

    def test_childless_couple_keeps_the_sum(self):
        """With no descent line (a couple without children), the labels only need not meet (#41)."""
        from pypedigree import PedigreeStyle, compute_layout, parse_pedigree_from_string

        rows = self.TEXT.splitlines()[:3]
        text = rows[0] + "\tchildless\n" + rows[1] + "\tchoice:mo\n" + rows[2] + "\t0\n"
        with_children = {
            n.individual_id: n.x
            for n in compute_layout(
                parse_pedigree_from_string(self.TEXT).families["F"], PedigreeStyle()
            ).nodes
        }
        alone = {
            n.individual_id: n.x
            for n in compute_layout(
                parse_pedigree_from_string(text).families["F"], PedigreeStyle()
            ).nodes
        }
        assert alone["mo"] - alone["fa"] < with_children["mo"] - with_children["fa"]
