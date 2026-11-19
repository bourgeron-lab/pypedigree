"""The web app (``web/``): its Python glue, its examples and the site build.

The page runs ``web/pypedigree_web.py`` in Pyodide, with pypedigree installed
without its dependencies; here it runs under CPython. The JavaScript is
checked in a browser, not here; only its references to the page are.

The web app is not shipped in the sdist: these tests are skipped there.
All pedigrees here are synthetic.
"""

from __future__ import annotations

import importlib.util
import json
import re
import subprocess
import sys
import textwrap
from pathlib import Path
from types import ModuleType

import pytest

WEB = Path(__file__).resolve().parent.parent / "web"
if not (WEB / "pypedigree_web.py").exists():
    pytest.skip("web/ is not in the sdist", allow_module_level=True)


def _load(name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, WEB / f"{name}.py")
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


web = _load("pypedigree_web")
build = _load("build")

PED = (
    "FID\tID\tFather\tMother\tSex\tProband\tb_ASD\tPGS_x\n"
    "F\t1\t0\t0\t1\t0\t1\t0.3\n"
    "F\t2\t0\t0\t2\t0\t1\t-1\n"
    "F\t3\t1\t2\t1\t1\t2\t1\n"
    "G\t4\t0\t0\t1\t0\t1\t0\n"
)


def _render(text: str = PED, **options: object) -> dict:
    return json.loads(web.render(text, json.dumps(options)))


class TestRender:
    def test_draws_every_family(self):
        out = _render()
        assert out["ok"] is True
        assert out["families"] == ["F", "G"]
        assert out["family"] == ""
        assert out["individuals"] == 4
        assert out["warnings"] == []

    def test_interactive_svg_is_prefixed_and_responsive(self):
        svg = _render()["svg"]
        assert 'class="pypedigree"' in svg
        assert 'width="100%"' in svg
        assert 'id="pv-main-group"' in svg
        assert "<script>" in svg

    def test_static_svg_has_a_size_and_no_script(self):
        svg = _render()["static_svg"]
        assert re.search(r'width="\d+" height="\d+"', svg)
        assert "<script" not in svg

    def test_one_family(self):
        out = _render(family="G")
        assert out["family"] == "G"
        assert out["families"] == ["F", "G"]
        assert out["individuals"] == 1

    def test_min_font_size(self):
        """#51: the form's "Smallest text" raises every text of the drawing."""
        svg = _render(min_font_size="6")["static_svg"]
        sizes = [float(v) for v in re.findall(r'<text [^>]*font-size="([0-9.]+)"', svg)]
        assert min(sizes) >= 6 * 0.02 - 1e-9
        assert _render(min_font_size="")["ok"] is True

    @pytest.mark.parametrize("bad", ["six", "0"])
    def test_bad_min_font_size(self, bad):
        out = _render(min_font_size=bad)
        assert out["ok"] is False
        assert "smallest text" in out["error"].lower() or "min_font_size" in out["error"]

    def test_one_family_keeps_the_columns(self):
        """#50: the family is drawn with its bands, PGS curves and legend."""
        for key in ("svg", "static_svg"):
            svg = _render(family="F")[key]
            assert "PGS distribution (x)" in svg
            assert svg.count('class="pgs"') >= 3

    def test_options_reach_the_drawing(self):
        default = _render()["static_svg"]
        assert ">3</text>" in default and "Legend" in default
        assert _render(phenotype_layout="quadrants")["static_svg"] == default  # auto, one column
        assert _render(phenotype_layout="bands")["static_svg"] != default
        assert ">3</text>" not in _render(show_labels=False)["static_svg"]
        assert "Legend" not in _render(show_legend=False)["static_svg"]

    def test_style_json_as_text_or_object(self):
        as_text = _render(style='{"font_size": 14}')["static_svg"]
        as_object = _render(style={"font_size": 14})["static_svg"]
        assert as_text == as_object != _render()["static_svg"]

    def test_phantom_parent_warns_and_is_not_counted(self):
        out = _render("F 1 0 0 1 1\nF 2 0 0 2 1\nF 3 1 2 1 2\nF 4 9 2 1 1\n")
        assert out["ok"] is True
        assert out["individuals"] == 4
        assert any("'9'" in w for w in out["warnings"])

    def test_unknown_style_key_warns(self):
        out = _render(style={"font_size": 12, "fontsize": 3})
        assert out["ok"] is True
        assert any("fontsize" in w for w in out["warnings"])

    def test_version(self):
        import pypedigree

        assert web.version() == pypedigree.__version__


class TestErrors:
    @pytest.mark.parametrize(
        ("text", "options", "message"),
        [
            ("   \n", {}, "empty"),
            ("not a pedigree\n", {}, "separator"),
            (PED, {"family": "Z"}, "Family 'Z' not found. Available: F, G"),
            (PED, {"phenotype_layout": "rings"}, "Unknown phenotype_layout 'rings'"),
            (PED, {"style": "{font_size: 12}"}, "not valid JSON"),
            (PED, {"style": "[1, 2]"}, "must be a JSON object"),
            (PED, {"style": {"phenotype_colors": {"ASD": "nocolor"}}}, "Invalid color"),
        ],
    )
    def test_error_is_reported_not_raised(self, text, options, message):
        out = _render(text, **options)
        assert out["ok"] is False
        assert message in out["error"]
        assert "svg" not in out


class TestWithoutDependencies:
    def test_renders_with_matplotlib_numpy_and_click_missing(self):
        """As in Pyodide, where pypedigree is installed with deps=False."""
        code = f"""
            import sys
            class _Block:
                def find_spec(self, name, path=None, target=None):
                    if name.split(".")[0] in ("matplotlib", "numpy", "click", "pandas"):
                        raise ModuleNotFoundError(f"No module named {{name!r}}")
                    return None
            sys.meta_path.insert(0, _Block())
            sys.path.insert(0, {str(WEB)!r})
            import json, pypedigree_web
            out = json.loads(pypedigree_web.render({PED!r}, '{{"phenotype_layout": "quadrants"}}'))
            assert out["ok"], out
            print(len(out["svg"]) > 1000, len(out["static_svg"]) > 1000)
        """
        result = subprocess.run(
            [sys.executable, "-c", textwrap.dedent(code)],
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == 0, result.stderr
        assert result.stdout.split() == ["True", "True"]


class TestExamples:
    @pytest.mark.parametrize("example", build.EXAMPLES, ids=lambda e: str(e["file"]))
    def test_example_draws_without_warnings(self, example):
        text = (build.DATA / str(example["file"])).read_text(encoding="utf-8")
        options = dict(example.get("options", {}))  # type: ignore[call-overload]
        if "style" in example:
            options["style"] = json.loads(
                (build.DATA / str(example["style"])).read_text(encoding="utf-8")
            )
        out = _render(text, **options)
        assert out["ok"] is True, out.get("error")
        assert out["warnings"] == []


class TestSite:
    def test_build_with_a_given_wheel(self, tmp_path):
        wheel = tmp_path / "pypedigree-9.9.9-py3-none-any.whl"
        wheel.write_bytes(b"not a real wheel")
        out = build.build(tmp_path / "site", wheel)
        config = json.loads((out / "config.json").read_text(encoding="utf-8"))
        assert config["version"] == "9.9.9"
        assert config["wheel"] == "wheels/pypedigree-9.9.9-py3-none-any.whl"
        assert config["pyodide"] == build.PYODIDE_VERSION
        assert (out / config["wheel"]).exists()
        for name in build.PAGE_FILES:
            assert (out / name).exists(), name
        assert (out / "docs" / "index.html").exists()
        assert (out / ".nojekyll").exists()
        assert len(config["examples"]) == len(build.EXAMPLES)
        for example in config["examples"]:
            assert (out / example["file"]).exists()
            assert example["title"] and example["description"]
        styled = [e for e in config["examples"] if "style" in e]
        assert styled and all(isinstance(e["style"], dict) for e in styled)

    def test_build_replaces_an_old_site(self, tmp_path):
        wheel = tmp_path / "pypedigree-1.0.0-py3-none-any.whl"
        wheel.write_bytes(b"")
        site = tmp_path / "site"
        (site / "wheels").mkdir(parents=True)
        (site / "wheels" / "pypedigree-0.1.0-py3-none-any.whl").write_bytes(b"")
        build.build(site, wheel)
        assert [p.name for p in (site / "wheels").iterdir()] == [wheel.name]

    def test_examples_have_distinct_names(self):
        """A link names an example by its file name without the extension
        (``?example=quadrants``), so no two may share one."""
        names = [Path(str(e["file"])).stem for e in build.EXAMPLES]
        assert len(set(names)) == len(names)
        assert "quadrants" in names  # the name the help page links to

    def test_hidden_wins_over_display_rules(self):
        """The download buttons sit in a flex row; `hidden` must still hide them."""
        css = (WEB / "app.css").read_text(encoding="utf-8")
        assert re.search(r"\[hidden\]\s*\{\s*display:\s*none\s*!important;\s*\}", css)

    def test_script_only_uses_ids_of_the_page(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        ids = set(re.findall(r'id="([^"]+)"', html))
        used = set(re.findall(r'\$\("([^"]+)"\)', (WEB / "app.js").read_text(encoding="utf-8")))
        assert used and used <= ids, used - ids

    def test_form_fields_match_the_options_sent(self):
        html = (WEB / "index.html").read_text(encoding="utf-8")
        names = set(re.findall(r'name="([^"]+)"', html))
        assert set(web.STYLE_CHOICES) <= names
        assert {
            "family",
            "show_labels",
            "show_legend",
            "show_na_labels",
            "min_font_size",
            "style",
        } <= names
        for name, allowed in web.STYLE_CHOICES.items():
            select = re.search(rf'<select id="{name}".*?</select>', html, re.S)
            assert select is not None, name
            assert set(re.findall(r'value="([^"]+)"', select.group(0))) == set(allowed)
