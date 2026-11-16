"""Individuals with a row but never sampled (issue #6).

"Dummy" individuals link real members (the shared parents of MZ twins, the
grandparents joining an aunt to her niece) without having been sampled. A
``not_sampled`` / ``inferred`` column, or the builder's ``inferred=True``,
draws them with a dashed outline like a phantom parent, while keeping their
row: label, parents, phenotypes, and their place in a ``.fam`` export.
"""

from __future__ import annotations

import re

import matplotlib.pyplot as plt
import pytest

from pypedigree import parse_pedigree_from_string, render_pedigree
from pypedigree.api import PedigreeBuilder
from pypedigree.layout import compute_layout
from pypedigree.parser import export_fam
from pypedigree.svg_renderer import render_interactive_svg

TEXT = (
    "FID\tIID\tFather\tMother\tSex\tPhenotype\tNot_sampled\tb_Autism\n"
    "FAM1\tgf\t0\t0\t1\t1\t1\t-9\n"
    "FAM1\tgm\t0\t0\t2\t1\tyes\t-9\n"
    "FAM1\taunt\tgf\tgm\t2\t2\t0\t2\n"
    "FAM1\tdad\tgf\tgm\t1\t1\t1\t1\n"
    "FAM1\tmum\t0\t0\t2\t1\t\t1\n"
    "FAM1\tniece\tdad\tmum\t2\t2\t0\t2\n"
)
INFERRED = {"gf", "gm", "dad"}


def _family():
    return parse_pedigree_from_string(TEXT).families["FAM1"]


class TestParsing:
    def test_column_marks_the_individuals(self):
        flags = {sid: ind.is_inferred for sid, ind in _family().individuals.items()}
        assert {sid for sid, flag in flags.items() if flag} == INFERRED

    @pytest.mark.parametrize("header", ["inferred", "Not_sampled", "unsampled", "dummy"])
    def test_column_names(self, header):
        text = f"FID\tIID\tFather\tMother\tSex\tPhenotype\t{header}\nF\ta\t0\t0\t1\t1\t1\n"
        assert parse_pedigree_from_string(text).families["F"].individuals["a"].is_inferred

    def test_they_are_not_phantoms(self):
        family = _family()
        assert not any(ind.is_phantom for ind in family.individuals.values())
        assert family.individuals["dad"].father_id == "gf"

    def test_builder_keyword(self):
        b = PedigreeBuilder("F")
        b.add_individual("dummy", sex="male", inferred=True)
        b.add_individual("real", sex="female")
        family = b.build().families["F"]
        assert family.individuals["dummy"].is_inferred
        assert not family.individuals["real"].is_inferred


class TestDrawing:
    def test_layout_keeps_labels(self):
        nodes = {n.individual_id: n for n in compute_layout(_family()).nodes}
        for sid in INFERRED:
            assert nodes[sid].is_inferred and nodes[sid].label == sid
        assert not nodes["niece"].is_inferred

    def test_png_dashed_outline_and_bands_kept(self):
        fig = render_pedigree(parse_pedigree_from_string(TEXT))
        try:
            ax = fig.axes[0]
            dashed = [p for p in ax.patches if p.get_linestyle() not in ("-", "solid")]
            # gf, gm, dad: one dashed outline each.
            assert len(dashed) == 3
            texts = [t.get_text() for t in ax.texts]
            assert {"gf", "gm", "dad"} <= set(texts)
        finally:
            plt.close(fig)

    def test_svg_dashed_outline_and_tooltip(self):
        svg = render_interactive_svg(parse_pedigree_from_string(TEXT))
        groups = re.findall(r'<g class="node[^"]*" data-id="([^"]+)".*?</g>', svg, re.S)
        bodies = dict(
            zip(groups, re.split(r'<g class="node[^"]*" data-id="[^"]+"', svg)[1:], strict=True)
        )
        for sid, body in bodies.items():
            assert ("stroke-dasharray" in body.split("</g>")[0]) == (sid in INFERRED), sid
        assert svg.count("Not sampled (inferred)") == 3

    def test_fam_export_keeps_them(self, tmp_path):
        out = export_fam(parse_pedigree_from_string(TEXT), tmp_path / "out.fam")
        ids = {line.split()[1] for line in out.read_text().splitlines()}
        assert ids >= INFERRED
