"""The legend's wording can be changed, and a minimal legend's frame hugs it (#79).

All pedigrees here are synthetic.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import pytest
from conftest import PHENOTYPE_DIR

from pypedigree import (
    PedigreeStyle,
    parse_pedigree,
    parse_pedigree_from_string,
    render_interactive_svg,
    render_pedigree,
)
from pypedigree.errors import PedigreeError
from pypedigree.phenotypes import (
    LEGEND_TEXT_DEFAULTS,
    check_legend,
    legend_metrics,
    legend_sex_items,
    legend_text,
)

FIXTURE = PHENOTYPE_DIR / "legend_minimal.ped"  # quadrants, a gradient, two PGS
PTDT = PHENOTYPE_DIR / "ptdt.ped"
CARRIER = (
    "FID\tIID\tFather\tMother\tSex\tPhenotype\tcr_GJB2\n"
    "FAM1\ta\t0\t0\t1\t1\t1\nFAM1\tb\t0\t0\t2\t1\t1\nFAM1\tc\ta\tb\t1\t2\t2\n"
)
FRENCH = {
    "shared_binary": "Phénotypes de la clé",
    "affected": "Atteint",
    "unaffected": "Non atteint",
    "unknown": "Inconnu",
    "male": "Homme",
    "female": "Femme",
    "low": "Bas",
    "high": "Haut",
    "pgs": "Score polygénique :",
}


def _png_texts(path, **style) -> list[str]:
    fig = render_pedigree(parse_pedigree(path), style=PedigreeStyle(**style))
    try:
        return [t.get_text() for t in fig.axes[1].texts]
    finally:
        plt.close(fig)


def _svg_legend(path, **style) -> str:
    svg = render_interactive_svg(parse_pedigree(path), style=PedigreeStyle(**style))
    start = svg.index('<g class="legend">')
    return svg[start : svg.index("\n</g>", start)]


class TestWording:
    def test_defaults_unchanged(self):
        texts = _png_texts(FIXTURE)
        for word in ("Affected", "Unaffected", "Unknown", "= Male", "= Female", "Low", "High"):
            assert word in texts
        assert "Phenotypes in the key (their colors)" in texts
        assert legend_text(PedigreeStyle(), "key_note") == LEGEND_TEXT_DEFAULTS["key_note"]

    def test_every_role_in_png_and_svg(self):
        texts = _png_texts(FIXTURE, legend_text=FRENCH)
        legend = _svg_legend(FIXTURE, legend_text=FRENCH)
        for role, word in FRENCH.items():
            shown = f"= {word}" if role in ("male", "female") else word
            assert shown in texts, role
            assert f">{shown}<" in legend, role
        for old in ("Affected", "= Male", "Low", "PGS distribution:"):
            assert old not in texts and f">{old}<" not in legend

    def test_key_note_and_ptdt(self):
        note = "Coloured: affected · white: unaffected · grey: unknown"
        assert note in _png_texts(FIXTURE, legend="minimal", legend_text={"key_note": note})
        ptdt = "Écart pTDT :"
        assert ptdt in _png_texts(PTDT, legend_text={"ptdt": ptdt})

    def test_scores_key_still_follows_the_caption(self):
        texts = _png_texts(FIXTURE, legend_text={"pgs": "Polygenic scores:"})
        assert "Polygenic scores:" in texts and "autism" in texts and "ADHD" in texts

    def test_wrapping_measures_the_replaced_text(self):
        """A longer caption wraps where the default fitted on one line (#68)."""
        long = "Polygenic score of each individual, standardised in the cohort:"
        fig = render_pedigree(
            parse_pedigree(FIXTURE), style=PedigreeStyle(legend_text={"pgs": long})
        )
        try:
            legend = fig.axes[1]
            y = {t.get_text(): t.get_position()[1] for t in legend.texts}
            # The caption, wrapped itself, then the key on the line below it.
            caption = [v for t, v in y.items() if t.startswith("Polygenic")]
            assert caption and y["autism"] < min(caption)
        finally:
            plt.close(fig)


class TestSexRow:
    def test_long_male_text_moves_the_female_symbol(self):
        m = legend_metrics(1, PedigreeStyle())
        (_, _, short), (female_x, _, _) = legend_sex_items(PedigreeStyle(), m)
        assert short == "= Male" and female_x == 2.5
        style = PedigreeStyle(legend_text={"male": "Homme (sexe assigné à la naissance)"})
        (_, _, _), (moved_x, _, female) = legend_sex_items(style, m)
        assert moved_x > 2.5 and female == "= Female"

    def test_no_overlap_in_png(self):
        style = PedigreeStyle(legend_text={"male": "Homme (sexe assigné à la naissance)"})
        fig = render_pedigree(parse_pedigree(FIXTURE), style=style)
        try:
            fig.canvas.draw()
            legend = fig.axes[1]
            male = next(t for t in legend.texts if t.get_text().startswith("= Homme"))
            circle = next(
                line for line in legend.lines if line.get_marker() == "o"
            ).get_window_extent()
            assert male.get_window_extent().x1 < circle.x0
        finally:
            plt.close(fig)


class TestValidation:
    def test_unknown_role_suggests_a_spelling(self):
        with pytest.raises(PedigreeError, match=r"'afected' \(did you mean 'affected'\?\)"):
            check_legend(PedigreeStyle(legend_text={"afected": "x"}))

    def test_text_only(self):
        with pytest.raises(PedigreeError, match="legend_text must map roles to text"):
            check_legend(PedigreeStyle(legend_text={"male": 1}))  # type: ignore[dict-item]

    def test_rendering_checks_it(self):
        with pytest.raises(PedigreeError, match="legend_text"):
            render_pedigree(parse_pedigree(FIXTURE), style=PedigreeStyle(legend_text={"x": "y"}))

    def test_from_json(self):
        assert PedigreeStyle.from_dict({"legend_text": FRENCH}).legend_text == FRENCH


class TestFrame:
    def _frame_and_content(self, legend: str) -> tuple[float, float]:
        fig = render_pedigree(
            parse_pedigree_from_string(CARRIER), style=PedigreeStyle(legend=legend)
        )
        try:
            fig.canvas.draw()
            ax = fig.axes[1]
            frame = max(ax.patches, key=lambda p: p.get_window_extent().width).get_window_extent()
            content = max(t.get_window_extent().x1 for t in ax.texts)
            return frame.width, frame.x1 - content
        finally:
            plt.close(fig)

    def test_minimal_frame_hugs_its_content(self):
        compact_w, compact_margin = self._frame_and_content("compact")
        minimal_w, minimal_margin = self._frame_and_content("minimal")
        assert minimal_w < compact_w
        assert minimal_margin < compact_margin
        assert 0 < minimal_margin < 40  # pixels: the left padding, mirrored

    def test_svg_frame_hugs_too(self):
        import re

        def width(legend: str) -> float:
            svg = render_interactive_svg(
                parse_pedigree_from_string(CARRIER), style=PedigreeStyle(legend=legend)
            )
            part = svg[svg.index('<g class="legend">') :]
            frame = re.search(r'<rect [^>]*width="([0-9.]+)"[^>]*rx=', part)
            assert frame is not None
            return float(frame.group(1))

        assert width("minimal") < width("compact")
