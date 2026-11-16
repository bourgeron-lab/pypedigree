"""A couple alone in its family is drawn with a "?" child.

A family holding only a man and a woman, both founders, with no child in
the file (a trio whose child is missing, say) used to be drawn as two
unrelated symbols. The layout now draws their couple, father on the left,
with a descent line to a "?". The pedigree itself is not changed.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import pytest

from pypedigree import PedigreeStyle, parse_pedigree_from_string, render_pedigree
from pypedigree.layout import (
    UNKNOWN_OFFSPRING_ID,
    compute_layout,
    compute_layout_circular,
    with_unknown_offspring,
)
from pypedigree.models import EdgeType
from pypedigree.svg_renderer import render_interactive_svg

# Mother first in the file: the father must still be drawn on the left.
COUPLE = "FAM1\tmum\t0\t0\t2\t1\nFAM1\tdad\t0\t0\t1\t2\n"


def _family(text: str, fid: str = "FAM1"):
    return parse_pedigree_from_string(text).families[fid]


class TestWhenItApplies:
    def test_couple_alone_gets_a_placeholder_child(self):
        family = _family(COUPLE)
        drawn = with_unknown_offspring(family)
        child = drawn.individuals[UNKNOWN_OFFSPRING_ID]
        assert child.is_unknown_offspring
        assert (child.father_id, child.mother_id) == ("dad", "mum")
        assert [(u.partner_ids, u.children_ids) for u in drawn.unions] == [
            (("dad", "mum"), [UNKNOWN_OFFSPRING_ID])
        ]

    def test_the_pedigree_itself_is_not_changed(self):
        family = _family(COUPLE)
        with_unknown_offspring(family)
        assert set(family.individuals) == {"mum", "dad"}
        assert family.unions == []

    @pytest.mark.parametrize(
        "text",
        [
            # A child in the file: an ordinary trio.
            COUPLE + "FAM1\tkid\tdad\tmum\t1\t1\n",
            # Two women, or a sex unknown: no father and mother to pair.
            "FAM1\ta\t0\t0\t2\t1\nFAM1\tb\t0\t0\t2\t1\n",
            "FAM1\ta\t0\t0\t1\t1\nFAM1\tb\t0\t0\t0\t1\n",
            # Three unrelated people.
            COUPLE + "FAM1\tc\t0\t0\t1\t1\n",
            # A single person.
            "FAM1\ta\t0\t0\t1\t1\n",
        ],
    )
    def test_other_families_are_left_alone(self, text):
        family = _family(text)
        assert with_unknown_offspring(family) is family

    def test_placeholder_id_never_collides(self):
        family = _family(f"FAM1\t{UNKNOWN_OFFSPRING_ID}\t0\t0\t1\t1\nFAM1\tmum\t0\t0\t2\t1\n")
        drawn = with_unknown_offspring(family)
        placeholders = [i for i in drawn.individuals.values() if i.is_unknown_offspring]
        assert len(placeholders) == 1
        assert placeholders[0].sample_id != UNKNOWN_OFFSPRING_ID


class TestLayout:
    def test_father_left_mother_right_question_mark_below(self):
        layout = compute_layout(_family(COUPLE), PedigreeStyle())
        nodes = {n.individual_id: n for n in layout.nodes}
        dad, mum, child = nodes["dad"], nodes["mum"], nodes[UNKNOWN_OFFSPRING_ID]
        assert dad.x < mum.x
        assert child.is_unknown_offspring and child.label == ""
        assert dad.x < child.x < mum.x
        assert child.y < dad.y  # one generation below
        kinds = {e.edge_type for e in layout.edges}
        assert EdgeType.RELATIONSHIP in kinds

    def test_circular_layout_too(self):
        layout = compute_layout_circular(_family(COUPLE), PedigreeStyle())
        assert any(n.is_unknown_offspring for n in layout.nodes)


class TestRendering:
    def test_png_writes_a_question_mark_and_no_symbol_for_it(self):
        pedigree = parse_pedigree_from_string(COUPLE)
        fig = render_pedigree(pedigree)
        try:
            ax = fig.axes[0]
            texts = [t.get_text() for t in ax.texts]
            assert texts.count("?") == 1
            assert "dad" in texts and "mum" in texts
            # Two symbols only, the square and the circle (a trio has three).
            assert len(ax.patches) == 2
        finally:
            plt.close(fig)

    def test_svg_writes_a_question_mark_without_a_node(self):
        svg = render_interactive_svg(parse_pedigree_from_string(COUPLE))
        assert svg.count('class="unknown-offspring"') == 1
        assert svg.count('<g class="node"') == 2
        assert UNKNOWN_OFFSPRING_ID not in svg

    def test_multi_family_file(self):
        text = COUPLE + "FAM2\tf\t0\t0\t1\t1\nFAM2\tm\t0\t0\t2\t1\nFAM2\tk\tf\tm\t1\t2\n"
        svg = render_interactive_svg(parse_pedigree_from_string(text))
        assert svg.count('class="unknown-offspring"') == 1  # FAM1 only
        assert svg.count('<g class="node"') == 5
