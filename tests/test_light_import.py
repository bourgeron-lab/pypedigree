"""Parsing, layout and SVG output without matplotlib (#8).

All pedigrees here are synthetic.
"""

from __future__ import annotations

import subprocess
import sys
import textwrap
from pathlib import Path

import matplotlib.colors as mcolors
import pytest

import pypedigree
from pypedigree.colors import _BASE_COLORS, _NAMED_COLORS, parse_rgb, rgb_to_hex
from pypedigree.errors import PedigreeError
from pypedigree.phenotypes import to_hex, to_rgb

PED = (
    "FID\tID\tFather\tMother\tSex\tb_ASD\tPGS_x\n"
    "F\t1\t0\t0\t1\t2\t0.3\n"
    "F\t2\t0\t0\t2\t1\t-1\n"
    "F\t3\t1\t2\t1\t2\t1\n"
)

# Run in a fresh interpreter: matplotlib is already imported in this one.
_HEAVY = ("matplotlib", "matplotlib.pyplot", "numpy", "pandas")


def _run(code: str, *, block_matplotlib: bool = False) -> str:
    prelude = ""
    if block_matplotlib:
        # As if matplotlib and numpy were not installed.
        prelude = textwrap.dedent(
            """
            import sys
            class _Block:
                def find_spec(self, name, path=None, target=None):
                    if name.split(".")[0] in ("matplotlib", "numpy"):
                        raise ModuleNotFoundError(f"No module named {name!r}")
                    return None
            sys.meta_path.insert(0, _Block())
            """
        )
    result = subprocess.run(
        [sys.executable, "-c", prelude + textwrap.dedent(code)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    return result.stdout


class TestImport:
    def test_import_loads_neither_matplotlib_nor_numpy(self) -> None:
        out = _run(
            f"""
            import sys
            import pypedigree
            print([m for m in {_HEAVY!r} if m in sys.modules])
            """
        )
        assert out.strip() == "[]"

    def test_parse_layout_and_svg_load_no_matplotlib(self) -> None:
        out = _run(
            f"""
            import sys
            from pypedigree import (
                PedigreeStyle, compute_layout, parse_pedigree_from_string,
                render_interactive_svg,
            )
            ped = parse_pedigree_from_string({PED!r})
            compute_layout(ped.families["F"], PedigreeStyle())
            svg = render_interactive_svg(ped, style=PedigreeStyle(fill_color="tab:red"))
            assert "<svg" in svg
            print([m for m in {_HEAVY!r} if m in sys.modules])
            """
        )
        assert out.strip() == "[]"

    def test_cli_interactive_loads_no_matplotlib(self, tmp_path: Path) -> None:
        ped = tmp_path / "f.ped"
        ped.write_text(PED, encoding="utf-8")
        out = _run(
            f"""
            import sys
            from pypedigree.cli import main
            try:
                main([{str(ped)!r}, "--interactive", "-o", {str(tmp_path / "f.svg")!r}])
            except SystemExit as e:
                assert e.code in (0, None), e.code
            print([m for m in {_HEAVY!r} if m in sys.modules])
            """
        )
        assert out.strip().splitlines()[-1] == "[]"
        assert (tmp_path / "f.svg").exists()

    def test_svg_works_without_matplotlib_installed(self) -> None:
        out = _run(
            f"""
            from pypedigree import PedigreeStyle, parse_pedigree_from_string
            from pypedigree import render_interactive_svg
            ped = parse_pedigree_from_string({PED!r})
            svg = render_interactive_svg(ped, style=PedigreeStyle(line_color="0.3"))
            print("ok" if "<svg" in svg else "no svg")
            """,
            block_matplotlib=True,
        )
        assert out.strip() == "ok"

    def test_without_matplotlib_its_own_colors_are_invalid(self) -> None:
        out = _run(
            """
            from pypedigree.errors import PedigreeError
            from pypedigree.phenotypes import to_rgb
            try:
                to_rgb("xkcd:sky blue")
            except PedigreeError:
                print("refused")
            """,
            block_matplotlib=True,
        )
        assert out.strip() == "refused"

    def test_without_matplotlib_render_pedigree_says_so(self) -> None:
        out = _run(
            """
            import pypedigree
            try:
                pypedigree.render_pedigree
            except ImportError as e:
                print(type(e).__name__, "matplotlib" in str(e))
            """,
            block_matplotlib=True,
        )
        assert out.strip() == "ModuleNotFoundError True"


class TestLazyRenderer:
    def test_render_functions_from_the_package_root(self) -> None:
        from pypedigree import render_family, render_pedigree
        from pypedigree.renderer import render_family as rf
        from pypedigree.renderer import render_pedigree as rp

        assert render_pedigree is rp
        assert render_family is rf

    def test_listed_by_dir_and_all(self) -> None:
        assert {"render_pedigree", "render_family"} <= set(dir(pypedigree))
        assert {"render_pedigree", "render_family"} <= set(pypedigree.__all__)

    def test_unknown_attribute(self) -> None:
        with pytest.raises(AttributeError, match="no_such_thing"):
            pypedigree.no_such_thing  # noqa: B018


class TestColorsMatchMatplotlib:
    """`pypedigree.colors` reads every color as matplotlib.colors does."""

    @pytest.mark.parametrize("name", sorted({*_BASE_COLORS, *_NAMED_COLORS}))
    def test_every_name(self, name: str) -> None:
        assert parse_rgb(name) == mcolors.to_rgb(name)

    def test_tables_are_complete(self) -> None:
        assert set(_NAMED_COLORS) == {*mcolors.CSS4_COLORS, *mcolors.TABLEAU_COLORS}
        assert set(_BASE_COLORS) == set(mcolors.BASE_COLORS)

    @pytest.mark.parametrize(
        "value",
        [
            "#e74c3c",
            "#E74C3C",
            "#abc",
            "#abcd",
            "#11223344",
            "Red",
            "GREY",
            "DarkSlateGray",
            "Tab:Blue",
            "none",
            "None",
            "0",
            "1",
            "0.35",
            " 0.5",
            (0.1, 0.2, 0.3),
            (0.1, 0.2, 0.3, 0.5),
            [1, 0, 0],
            "xkcd:sky blue",
            "C1",
        ],
    )
    def test_same_rgb(self, value: object) -> None:
        assert parse_rgb(value) == mcolors.to_rgb(value)  # type: ignore[arg-type]

    @pytest.mark.parametrize(
        "value",
        [
            "R",
            "notacolor",
            "#12",
            "#12345",
            "#gggggg",
            "1.5",
            "-0.1",
            "nan",
            "",
            (1, 2, 3),
            (0.1, 0.2),
            3,
            None,
        ],
    )
    def test_same_refusals(self, value: object) -> None:
        with pytest.raises(ValueError):
            mcolors.to_rgb(value)  # type: ignore[arg-type]
        with pytest.raises(ValueError):
            parse_rgb(value)

    @pytest.mark.parametrize(
        "rgb",
        [
            (0.0, 0.0, 0.0),
            (1.0, 1.0, 1.0),
            (0.5, 0.25, 0.75),
            (0.002, 0.998, 0.501),
            (1 / 3, 2 / 3, 0.1),
        ],
    )
    def test_same_hex(self, rgb: tuple[float, float, float]) -> None:
        assert rgb_to_hex(rgb) == mcolors.to_hex(rgb)

    def test_phenotypes_helpers(self) -> None:
        assert to_hex("tab:red") == mcolors.to_hex("tab:red")
        assert to_hex((0.5, 0.25, 0.75, 1.0)) == mcolors.to_hex((0.5, 0.25, 0.75))
        assert to_rgb("grey") == mcolors.to_rgb("grey")
        with pytest.raises(PedigreeError, match="Invalid color"):
            to_rgb("notacolor")
