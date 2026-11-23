"""NSGC numbering of individuals from the drawing: generation and position (#125).

Generations in Roman numerals from the top, individuals in Arabic numerals from
the left within their generation (I-1, II-3). All pedigrees here are synthetic.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from pypedigree import (
    PedigreeStyle,
    compute_layout,
    compute_layout_circular,
    parse_pedigree,
    parse_pedigree_from_string,
    render_interactive_svg,
)
from pypedigree.cli import main
from pypedigree.errors import PedigreeDataWarning, PedigreeError, PedigreeParseError
from pypedigree.symbols import generation_numbers, roman

DATA = Path(__file__).parent / "data"
NOMENCLATURE = DATA / "nomenclature" / "nomenclature.ped"

# A child whose father has no row (a named phantom), and a mother with two
# partners (drawn twice).
FAMILY = (
    "FID\tIID\tPAT\tMAT\tSEX\tPHENO\n"
    "F\tgf\t0\t0\t1\t1\nF\tgm\t0\t0\t2\t1\n"
    "F\tmum\tgf\tgm\t2\t1\n"
    "F\tp1\t0\t0\t1\t1\nF\tp2\t0\t0\t1\t1\n"
    "F\tc1\tp1\tmum\t1\t2\nF\tc2\tp2\tmum\t2\t1\n"
    "F\tlone\tdad\t0\t1\t1\n"
)


def _nodes(text: str, labels: str, **kw):
    with pytest.warns(PedigreeDataWarning, match="drawn as a phantom"):
        family = parse_pedigree_from_string(text).families["F"]
    return compute_layout(family, PedigreeStyle(labels=labels, **kw)).nodes


class TestNumbers:
    @pytest.mark.parametrize(
        ("n", "text"), [(1, "I"), (4, "IV"), (9, "IX"), (12, "XII"), (40, "XL")]
    )
    def test_roman(self, n, text):
        assert roman(n) == text

    def test_rank_by_level_then_x(self):
        levels = {"a": 0, "b": 0, "c": 2, "d": 2, "e": 2}
        xs = {"a": 5.0, "b": 1.0, "c": 3.0, "d": 1.0, "e": 2.0}
        assert generation_numbers(levels, xs, levels) == {
            "b": "I-1", "a": "I-2", "d": "II-1", "e": "II-2", "c": "II-3",
        }  # fmt: skip


class TestLayout:
    def test_default_keeps_ids(self):
        nodes = _nodes(FAMILY, "id")
        assert {n.label for n in nodes} >= {"gf", "mum", "c1"}
        assert all(n.number is None for n in nodes)

    def test_generation_labels_follow_the_drawing(self):
        nodes = _nodes(FAMILY, "generation")
        for gen, row in _rows(nodes).items():
            labels = [n.label for n in row if n.label and not n.is_duplicate]
            assert labels == [f"{roman(gen)}-{i}" for i in range(1, len(labels) + 1)], labels

    def test_duplicate_takes_its_originals_number(self):
        nodes = _nodes(FAMILY, "generation")
        mum = next(n for n in nodes if n.individual_id == "mum")
        dups = [n for n in nodes if n.is_duplicate]
        assert dups and all(d.label == mum.label for d in dups)

    def test_unnamed_phantoms_unnumbered_named_ones_numbered(self):
        nodes = {n.individual_id: n for n in _nodes(FAMILY, "generation")}
        assert nodes["dad"].is_phantom and nodes["dad"].number  # named in the file: labelled
        unnamed = [n for n in nodes.values() if n.individual_id.startswith("_phantom_")]
        assert unnamed and all(not n.label and n.number is None for n in unnamed)

    def test_both_writes_the_id_first_under_the_number(self):
        family = parse_pedigree(NOMENCLATURE).families["FAM1"]
        nodes = {
            n.individual_id: n for n in compute_layout(family, PedigreeStyle(labels="both")).nodes
        }
        plain = {n.individual_id: n for n in compute_layout(family).nodes}
        assert nodes["II-4"].notes[0] == "II-4" and nodes["II-4"].notes[1:] == plain["II-4"].notes
        # The file's II-2 and II-3 are drawn the other way round.
        assert nodes["II-2"].label == "II-3" and nodes["II-3"].label == "II-2"

    def test_both_keeps_room_for_the_extra_line(self):
        def gap(labels):
            ys = sorted({n.y for n in _nodes(FAMILY, labels)}, reverse=True)
            return ys[0] - ys[1]

        assert gap("both") > gap("generation") == pytest.approx(gap("id"))

    def test_circular_layout_has_the_same_numbers(self):
        family = parse_pedigree_from_string(FAMILY).families["F"]
        style = PedigreeStyle(labels="generation")
        flat = {n.individual_id: n.number for n in compute_layout(family, style).nodes}
        ring = {n.individual_id: n.number for n in compute_layout_circular(family, style).nodes}
        assert {k: v for k, v in ring.items() if v} == {k: v for k, v in flat.items() if v}

    def test_unknown_mode(self):
        with pytest.raises(PedigreeError, match="Unknown labels 'roman'"):
            _nodes(FAMILY, "roman")


def _rows(nodes):
    levels = sorted({n.y for n in nodes}, reverse=True)
    return {
        i + 1: sorted(
            (n for n in nodes if n.y == y and not n.is_unknown_offspring), key=lambda n: n.x
        )
        for i, y in enumerate(levels)
    }


class TestOutputs:
    def test_svg_keeps_the_id(self):
        svg = render_interactive_svg(
            parse_pedigree_from_string(FAMILY), style=PedigreeStyle(labels="generation")
        )
        assert 'data-id="c1"' in svg and 'data-number="III-' in svg
        assert "ID: c1" in svg and "Number: III-" in svg

    def test_data_number_column_is_reserved(self):
        with pytest.raises(PedigreeParseError, match="data_number"):
            parse_pedigree_from_string(
                "FID\tIID\tPAT\tMAT\tSEX\tPHENO\tdata_number\nF\ta\t0\t0\t1\t1\t3\n"
            )

    def test_cli(self, tmp_path):
        out = tmp_path / "n.svg"
        result = CliRunner().invoke(
            main, [str(NOMENCLATURE), "--labels", "both", "--interactive", "-o", str(out)]
        )
        assert result.exit_code == 0, result.output
        assert 'data-number="III-7"' in out.read_text(encoding="utf-8")
