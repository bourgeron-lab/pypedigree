"""The SVG can be inlined in a page next to other SVGs, or printed (#4).

All pedigrees here are synthetic.
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest
from click.testing import CliRunner

from pypedigree import render_interactive_svg, render_svg
from pypedigree.api import PedigreeBuilder
from pypedigree.cli import main
from pypedigree.config import PedigreeStyle
from pypedigree.models import Pedigree
from pypedigree.parser import parse_pedigree_from_string
from pypedigree.svg_renderer import ROOT_CLASS, check_id_prefix

# Bands (clip paths), a proband arrow (marker), carriers (patterns) and a
# legend: every kind of id the renderer writes.
PED = (
    "FID\tID\tFather\tMother\tSex\tb_ASD\tb_ID\to_Proband\tPGS_x\n"
    "F\t1\t0\t0\t1\t2\t1\t1\t0.4\n"
    "F\t2\t0\t0\t2\t1\t2\t1\t-0.2\n"
    "F\t3\t1\t2\t1\t2\t2\t2\t0.9\n"
    "F\t4\t1\t2\t2\t1\t1\t1\t0.1\n"
)

ID_RE = re.compile(r'\sid="([^"]+)"')
REF_RE = re.compile(r'="url\(#([^)]+)\)"|\shref="#([^"]+)"')


def _pedigree() -> Pedigree:
    pedigree = parse_pedigree_from_string(PED)
    fam = pedigree.families["F"]
    fam.individuals["2"].is_carrier = True
    return pedigree


def _ids(svg: str) -> list[str]:
    return ID_RE.findall(svg)


def _refs(svg: str) -> set[str]:
    return {a or b for a, b in REF_RE.findall(svg)}


def _parse(svg: str) -> ET.Element:
    return ET.fromstring(svg.replace(' xmlns="http://www.w3.org/2000/svg"', ""))


class TestRoot:
    def test_root_has_the_pypedigree_class(self) -> None:
        for render in (render_interactive_svg, render_svg):
            root = _parse(render(_pedigree()))
            assert root.get("class") == ROOT_CLASS

    def test_css_is_scoped_under_the_root_class(self) -> None:
        svg = render_interactive_svg(_pedigree())
        css = svg[svg.index("<style>") + len("<style>") : svg.index("</style>")]
        rules = [r.split("{")[0].strip() for r in css.split("}") if "{" in r]
        assert rules
        for selectors in rules:
            for selector in selectors.split(","):
                assert selector.strip().startswith(f".{ROOT_CLASS} "), selector


class TestIdPrefix:
    def test_default_ids_are_unchanged(self) -> None:
        ids = _ids(render_interactive_svg(_pedigree()))
        assert {"main-group", "tooltip", "tooltip-bg", "tooltip-text", "arrowhead"} <= set(ids)

    def test_every_id_and_reference_is_prefixed(self) -> None:
        plain = render_interactive_svg(_pedigree())
        svg = render_interactive_svg(_pedigree(), id_prefix="fam-a-")
        ids = _ids(svg)
        assert ids == ["fam-a-" + i for i in _ids(plain)]
        refs = _refs(svg)
        assert refs, "the fixture should draw clip paths, markers or patterns"
        assert all(r.startswith("fam-a-") for r in refs)
        assert refs <= set(ids), "every reference points to an id in the same SVG"

    def test_data_attributes_are_not_prefixed(self) -> None:
        svg = render_interactive_svg(_pedigree(), id_prefix="p_")
        assert 'data-id="1"' in svg
        assert 'data-id="p_1"' not in svg

    def test_two_figures_share_no_id(self) -> None:
        a = render_interactive_svg(_pedigree(), id_prefix="a-")
        b = render_interactive_svg(_pedigree(), id_prefix="b-")
        assert not set(_ids(a)) & set(_ids(b))
        assert len(_ids(a)) == len(set(_ids(a))), "ids are unique within one figure"

    def test_script_looks_up_the_prefixed_ids(self) -> None:
        svg = render_interactive_svg(_pedigree(), id_prefix="fam-a-")
        assert 'const P = "fam-a-";' in svg

    def test_prefix_is_kept_when_the_svg_is_parsed(self) -> None:
        root = _parse(render_svg(_pedigree(), id_prefix="x."))
        assert root.find(".//g[@id='x.main-group']") is not None

    @pytest.mark.parametrize("bad", ["1a", "-a", ".a", "a b", 'a"', "a#", "é"])
    def test_invalid_prefix_is_refused(self, bad: str) -> None:
        with pytest.raises(ValueError, match="id_prefix"):
            render_interactive_svg(_pedigree(), id_prefix=bad)

    @pytest.mark.parametrize("good", ["", "a", "_a", "fam1-", "a.b_c-"])
    def test_valid_prefix_is_accepted(self, good: str) -> None:
        check_id_prefix(good)


class TestScript:
    def test_script_is_scoped_to_its_own_root(self) -> None:
        svg = render_interactive_svg(_pedigree())
        script = svg[svg.index("<script>") : svg.index("</script>")]
        assert "document.currentScript" in script
        assert "closest('svg')" in script
        assert "document.querySelector('svg')" not in script
        assert "document.getElementById" not in script
        assert "document.querySelectorAll('.node')" not in script

    def test_no_listener_on_the_window(self) -> None:
        svg = render_interactive_svg(_pedigree())
        assert "window.addEventListener" not in svg
        assert "document.addEventListener" not in svg

    def test_wheel_zooms_only_with_ctrl_when_inlined(self) -> None:
        svg = render_interactive_svg(_pedigree())
        assert "svg === document.documentElement" in svg
        assert "if (!standalone && !e.ctrlKey && !e.metaKey) return;" in svg


class TestStatic:
    def test_no_script_and_no_tooltip(self) -> None:
        svg = render_svg(_pedigree())
        assert "<script" not in svg
        assert 'id="tooltip' not in svg
        assert "cursor: pointer" not in svg

    def test_same_drawing_as_the_interactive_svg(self) -> None:
        def drawing(svg: str) -> str:
            return svg[svg.index('<g id="main-group">') : svg.index("</g>\n<!-- Tooltip")]

        interactive = render_interactive_svg(_pedigree())
        static = render_svg(_pedigree())
        assert static.index('<g id="main-group">')
        assert drawing(interactive) in static

    def test_titles_and_data_attributes_are_kept(self) -> None:
        svg = render_svg(_pedigree())
        assert svg.count("<title>") == 4
        assert 'data-parents="1,2"' in svg
        assert ".node.dimmed" in svg

    def test_builder_render_svg(self, tmp_path: Path) -> None:
        builder = PedigreeBuilder("F")
        builder.add_individual("a", sex="male")
        out = tmp_path / "static.svg"
        svg = builder.render_svg(output=str(out), id_prefix="k-", responsive=True)
        assert out.read_text(encoding="utf-8") == svg
        assert "<script" not in svg
        assert 'id="k-main-group"' in svg
        assert 'width="100%"' in svg


class TestResponsive:
    def test_fixed_size_by_default(self) -> None:
        root = _parse(render_interactive_svg(_pedigree()))
        assert root.get("width", "").isdigit()
        assert root.get("height", "").isdigit()

    def test_responsive_fills_the_container_width(self) -> None:
        for render in (render_interactive_svg, render_svg):
            root = _parse(render(_pedigree(), responsive=True))
            assert root.get("width") == "100%"
            assert root.get("height") is None
            assert len(root.get("viewBox", "").split()) == 4

    def test_empty_pedigree(self) -> None:
        svg = render_svg(Pedigree(), id_prefix="e-", responsive=True)
        assert 'id="e-main-group"' in svg

    def test_style_is_honoured(self) -> None:
        svg = render_svg(_pedigree(), style=PedigreeStyle(line_color="#123456"))
        assert "#123456" in svg


class TestCli:
    def _ped(self, tmp_path: Path) -> Path:
        path = tmp_path / "f.ped"
        path.write_text(PED, encoding="utf-8")
        return path

    def test_no_script_writes_a_static_svg(self, tmp_path: Path) -> None:
        out = tmp_path / "s.svg"
        args = [str(self._ped(tmp_path)), "--interactive", "--no-script", "-o", str(out)]
        result = CliRunner().invoke(main, args)
        assert result.exit_code == 0, result.output
        assert "Saved static SVG" in result.output
        assert "<script" not in out.read_text(encoding="utf-8")

    def test_id_prefix_and_responsive(self, tmp_path: Path) -> None:
        out = tmp_path / "e.svg"
        args = [str(self._ped(tmp_path)), "--interactive", "--id-prefix", "fam-a-", "--responsive"]
        result = CliRunner().invoke(main, [*args, "-o", str(out)])
        assert result.exit_code == 0, result.output
        svg = out.read_text(encoding="utf-8")
        assert 'id="fam-a-main-group"' in svg
        assert 'width="100%"' in svg

    def test_invalid_prefix_is_a_usage_error(self, tmp_path: Path) -> None:
        args = [str(self._ped(tmp_path)), "--interactive", "--id-prefix", "1a"]
        result = CliRunner().invoke(main, [*args, "-o", str(tmp_path / "x.svg")])
        assert result.exit_code == 2
        assert "--id-prefix" in result.output

    @pytest.mark.parametrize("flag", [["--no-script"], ["--responsive"], ["--id-prefix", "a-"]])
    def test_options_need_interactive(self, tmp_path: Path, flag: list[str]) -> None:
        args = [str(self._ped(tmp_path)), *flag, "-o", str(tmp_path / "x.png")]
        result = CliRunner().invoke(main, args)
        assert result.exit_code == 2
        assert "need --interactive" in result.output
