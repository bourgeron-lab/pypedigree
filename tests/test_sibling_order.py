"""Strict birth order of siblings (`sibling_order`).

With "auto" (the default) a child whose line leads to a loop is moved to the
end of their sibship facing the other branch. With "birth" every sibship is
drawn in birth order (date of birth, then file order) and such a loop is
broken with a duplicate instead. All pedigrees here are synthetic.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest
from click.testing import CliRunner

from pypedigree import PedigreeStyle, compute_layout, compute_layout_circular, parse_pedigree
from pypedigree.cli import main
from pypedigree.errors import PedigreeError
from pypedigree.layout import _order_children

sys.path.insert(0, str(Path(__file__).parent))
from conftest import TEST_DATA_DIR, generate_realistic_pedigree  # noqa: E402

BIRTH = PedigreeStyle(sibling_order="birth")
MULTIGENERATIONAL = TEST_DATA_DIR / "complex" / "multigenerational.ped"


def _family(path: Path):
    return next(iter(parse_pedigree(path).families.values()))


def _out_of_order(family, style) -> list[tuple[str, ...]]:
    """Sibships whose children are not drawn in `_order_children`'s order.

    A sibship in which two siblings are partners is left out: the couple is
    drawn side by side, whatever their birth order.
    """
    layout = compute_layout(family, style)
    x = {n.individual_id: n.x for n in layout.nodes if not n.is_duplicate}
    couples = {frozenset(u.partner_ids) for u in family.unions}
    wrong = []
    for union in family.unions:
        kids = set(union.children_ids)
        if any(c <= kids for c in couples):
            continue  # two siblings partners: drawn side by side
        expected = [c for c in _order_children(family, union.children_ids) if c in x]
        drawn = sorted(expected, key=lambda c: x[c])
        if drawn != expected:
            wrong.append(tuple(expected))
    return wrong


class TestBirthOrder:
    def test_auto_moves_siblings_of_a_loop(self):
        """The gallery's five-generation pedigree has loops: by default some
        sibships leave birth order to close them."""
        assert _out_of_order(_family(MULTIGENERATIONAL), PedigreeStyle())

    def test_birth_keeps_every_sibship_in_order(self):
        assert _out_of_order(_family(MULTIGENERATIONAL), BIRTH) == []

    @pytest.mark.parametrize("seed", range(4))
    def test_generated_pedigrees_in_order(self, seed):
        assert _out_of_order(generate_realistic_pedigree(60, seed), BIRTH) == []

    def test_loops_broken_with_duplicates_instead(self):
        family = _family(MULTIGENERATIONAL)
        auto = sum(n.is_duplicate for n in compute_layout(family).nodes)
        birth = sum(n.is_duplicate for n in compute_layout(family, BIRTH).nodes)
        assert birth > auto

    def test_no_loop_no_change(self):
        """Without loops nothing is reordered, so "birth" draws the same."""
        family = _family(TEST_DATA_DIR / "extended" / "dob_ordering.ped")
        auto = [(n.individual_id, n.x, n.y) for n in compute_layout(family).nodes]
        birth = [(n.individual_id, n.x, n.y) for n in compute_layout(family, BIRTH).nodes]
        assert auto == birth

    def test_circular_layout_too(self):
        layout = compute_layout_circular(_family(MULTIGENERATIONAL), BIRTH)
        assert layout.nodes


class TestOption:
    def test_unknown_value(self):
        with pytest.raises(PedigreeError, match="sibling_order"):
            compute_layout(_family(MULTIGENERATIONAL), PedigreeStyle(sibling_order="dob"))

    def test_from_a_style_dict(self):
        assert PedigreeStyle.from_dict({"sibling_order": "birth"}).sibling_order == "birth"

    def test_cli(self, tmp_path):
        out = tmp_path / "out.svg"
        result = CliRunner().invoke(
            main,
            [str(MULTIGENERATIONAL), "--sibling-order", "birth", "--interactive", "-o", str(out)],
        )
        assert result.exit_code == 0, result.output
        auto = tmp_path / "auto.svg"
        args = [str(MULTIGENERATIONAL), "--interactive", "-o", str(auto)]
        result = CliRunner().invoke(main, args)
        assert result.exit_code == 0, result.output
        assert out.read_text(encoding="utf-8") != auto.read_text(encoding="utf-8")
