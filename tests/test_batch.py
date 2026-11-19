"""One file per family: `render_families`, `family_file_names`, `--output-dir` (#49).

All pedigrees here are synthetic.
"""

from __future__ import annotations

import re
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import pytest
from click.testing import CliRunner

import pypedigree as pp
from pypedigree import batch
from pypedigree.cli import main
from pypedigree.errors import PedigreeBatchError, PedigreeDataWarning, PedigreeError
from pypedigree.phenotypes import resolve_phenotype_render

# Three families. IQ is a score over the file, but family G alone has only
# the values 1 and 2, which would read as a binary column; the sites differ
# from family to family.
COHORT = (
    "FID\tIID\tPAT\tMAT\tSEX\tPHENO\tb_autism\tIQ\tc_site\tPGS_autism\tPGS_x\n"
    "F\tfa\t0\t0\t1\t1\t1\t95\tParis\t0.2\t0.5\n"
    "F\tmo\t0\t0\t2\t1\t1\t110\tLyon\t-0.4\t-1.0\n"
    "F\tch\tfa\tmo\t1\t2\t2\t70\tParis\t1.8\t0.3\n"
    "G\tg1\t0\t0\t2\t1\t1\t1\tNice\t0.0\t0.1\n"
    "G\tg2\t0\t0\t1\t1\t1\t2\tNice\t0.1\t0.1\n"
    "H\th1\t0\t0\t1\t1\t1\t100\tParis\t1.0\t-0.2\n"
)


@pytest.fixture
def cohort() -> pp.Pedigree:
    return pp.parse_pedigree_from_string(COHORT)


class TestFileNames:
    def test_ids_kept_when_safe(self):
        assert batch.family_file_names(["FAM1", "fam-2", "F.3", "Été_4"]) == {
            "FAM1": "FAM1",
            "fam-2": "fam-2",
            "F.3": "F.3",
            "Été_4": "Été_4",
        }

    def test_other_characters_replaced(self):
        names = batch.family_file_names(["a/b", "c d", "e:f*g?", "..hidden.", "///"])
        assert names == {
            "a/b": "a_b",
            "c d": "c_d",
            "e:f*g?": "e_f_g_",
            "..hidden.": "hidden",
            "///": "___",
        }

    def test_empty_name(self):
        with pytest.warns(PedigreeDataWarning, match="'family_2'"):
            names = batch.family_file_names(["..", " ", "."])
        assert names == {
            "..": "family",
            " ": "_",
            ".": "family_2",
        }

    def test_windows_device_names(self):
        names = batch.family_file_names(["CON", "nul", "COM1", "con.x", "CONSOLE"])
        assert names == {
            "CON": "CON_",
            "nul": "nul_",
            "COM1": "COM1_",
            "con.x": "con.x_",
            "CONSOLE": "CONSOLE",
        }

    def test_long_names_are_cut(self):
        assert len(batch.family_file_names(["F" * 300])["F" * 300]) == 100

    def test_collisions_get_a_number_and_a_warning(self):
        with pytest.warns(PedigreeDataWarning) as record:
            names = batch.family_file_names(["F_1", "f/1", "F 1", "F_1_2"])
        assert [str(w.message) for w in record][0] == (
            "Family 'f/1' is written as 'f_1_2': 'f_1' is the file name of an earlier family."
        )
        assert len(record) == 3
        # Case is ignored: macOS and Windows file systems would merge them.
        assert names == {"F_1": "F_1", "f/1": "f_1_2", "F 1": "F_1_3", "F_1_2": "F_1_2_2"}

    def test_no_warning_without_collision(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            batch.family_file_names(["A", "B"])


class TestRenderFamilies:
    @pytest.mark.parametrize("fmt", ["png", "pdf", "svg"])
    def test_one_file_per_family(self, cohort, tmp_path, fmt):
        files = pp.render_families(cohort, tmp_path / "out", fmt=fmt)
        assert list(files) == ["F", "G", "H"]
        for fid, path in files.items():
            assert path == tmp_path / "out" / f"{fid}.{fmt}"
            assert path.stat().st_size > 0

    @pytest.mark.parametrize("script", [True, False])
    def test_interactive(self, cohort, tmp_path, script):
        files = pp.render_families(
            cohort, tmp_path, interactive=True, script=script, id_prefix="p-", responsive=True
        )
        svg = files["G"].read_text(encoding="utf-8")
        assert files["G"].suffix == ".svg"
        assert ("<script>" in svg) is script
        assert 'id="p-' in svg
        assert 'width="100%"' in svg
        assert "fa" not in re.findall(r'data-id="([^"]*)"', svg)

    def test_chosen_families_in_order(self, cohort, tmp_path):
        files = pp.render_families(cohort, tmp_path, family_ids=["H", "F"], fmt="svg")
        assert list(files) == ["H", "F"]
        assert sorted(p.name for p in tmp_path.iterdir()) == ["F.svg", "H.svg"]

    def test_unknown_family(self, cohort, tmp_path):
        with pytest.raises(KeyError, match="'Z'"):
            pp.render_families(cohort, tmp_path, family_ids=["F", "Z"])
        assert not list(tmp_path.iterdir())

    def test_unknown_format(self, cohort, tmp_path):
        with pytest.raises(ValueError, match="Unknown format 'jpg'"):
            pp.render_families(cohort, tmp_path, fmt="jpg")

    def test_invalid_style_draws_nothing(self, cohort, tmp_path):
        style = pp.PedigreeStyle(phenotype_colors={"autism": "not-a-color"})
        with pytest.raises(PedigreeError, match="not-a-color"):
            pp.render_families(cohort, tmp_path / "out", style=style)
        assert not (tmp_path / "out").exists()

    def test_circular(self, cohort, tmp_path):
        files = pp.render_families(cohort, tmp_path, circular=True)
        assert len(files) == 3

    def test_circular_interactive_warns(self, cohort, tmp_path):
        with pytest.warns(pp.PedigreeStyleWarning, match="circular"):
            pp.render_families(cohort, tmp_path, interactive=True, circular=True)

    def test_figures_are_closed(self, cohort, tmp_path):
        mine = plt.figure()
        before = plt.get_fignums()
        pp.render_families(cohort, tmp_path)
        assert plt.get_fignums() == before
        plt.close(mine)

    def test_a_failure_does_not_stop_the_others(self, cohort, tmp_path, monkeypatch):
        real = batch._save_figure

        def fail_on_g(pedigree, path, *args):
            if "G" in pedigree.families:
                raise RuntimeError("layout went wrong")
            real(pedigree, path, *args)

        monkeypatch.setattr(batch, "_save_figure", fail_on_g)
        with pytest.raises(PedigreeBatchError) as info:
            pp.render_families(cohort, tmp_path)
        error = info.value
        assert list(error.written) == ["F", "H"]
        assert error.failures == {"G": "layout went wrong"}
        assert "1 of 3 families could not be drawn (G: layout went wrong)" in str(error)
        assert sorted(p.name for p in tmp_path.iterdir()) == ["F.png", "H.png"]
        assert plt.get_fignums() == []


class TestSameAsTheWholeFile:
    """Each file draws its columns as the figure of the whole file does."""

    def test_modes_and_categories_from_the_whole_file(self, cohort):
        names = cohort.phenotype_names
        modes = cohort.phenotype_display_modes
        style = pp.resolve_style(None, names)
        whole = resolve_phenotype_render(cohort, names, modes, style)
        alone = resolve_phenotype_render(cohort.subset(["G"]).column_scope, names, modes, style)
        assert alone == whole
        assert whole.effective("IQ") == "gradient"
        assert set(whole.categorical_colors["site"]) == {"Paris", "Lyon", "Nice"}

    def test_every_file_has_the_same_legend(self, cohort, tmp_path):
        files = pp.render_families(cohort, tmp_path, interactive=True, script=False)
        legends = set()
        for path in files.values():
            svg = path.read_text(encoding="utf-8")
            legend = svg[svg.index('class="legend"') :]
            legends.add(re.sub(r'(x|y|cx|cy|x1|x2|y1|y2|points|transform|d)="[^"]*"', "", legend))
        assert len(legends) == 1
        assert "High" in legends.pop()  # IQ is a gradient in G's file too


class TestCli:
    @pytest.fixture
    def ped(self, tmp_path) -> Path:
        path = tmp_path / "cohort.ped"
        path.write_text(COHORT, encoding="utf-8")
        return path

    def test_output_dir(self, ped, tmp_path):
        out = tmp_path / "figures"
        result = CliRunner().invoke(main, [str(ped), "--output-dir", str(out)])
        assert result.exit_code == 0, result.output
        assert sorted(p.name for p in out.iterdir()) == ["F.png", "G.png", "H.png"]
        assert result.output.count("Saved: ") == 3

    def test_format_and_interactive(self, ped, tmp_path):
        out = tmp_path / "figures"
        for args, ext in ((["--format", "pdf"], "pdf"), (["--interactive"], "svg")):
            result = CliRunner().invoke(main, [str(ped), "--output-dir", str(out), *args])
            assert result.exit_code == 0, result.output
            assert (out / f"F.{ext}").exists()

    def test_with_family(self, ped, tmp_path):
        out = tmp_path / "figures"
        result = CliRunner().invoke(
            main, [str(ped), "--output-dir", str(out), "--family", "G", "--format", "svg"]
        )
        assert result.exit_code == 0, result.output
        assert [p.name for p in out.iterdir()] == ["G.svg"]

    @pytest.mark.parametrize("other", [["-o", "x.png"], ["--show"]])
    def test_conflicts(self, ped, tmp_path, other):
        result = CliRunner().invoke(main, [str(ped), "--output-dir", str(tmp_path), *other])
        assert result.exit_code == 2
        assert "--output-dir cannot be used with" in result.output

    def test_bad_id_prefix(self, ped, tmp_path):
        result = CliRunner().invoke(
            main,
            [str(ped), "--output-dir", str(tmp_path), "--interactive", "--id-prefix", "1x"],
        )
        assert result.exit_code == 2
        assert "--id-prefix" in result.output

    def test_failures_reported_and_exit_code(self, ped, tmp_path, monkeypatch):
        def fail_on_g(pedigree, path, *args):
            if "G" in pedigree.families:
                raise RuntimeError("layout went wrong")
            path.write_bytes(b"png")

        monkeypatch.setattr(batch, "_save_figure", fail_on_g)
        out = tmp_path / "figures"
        result = CliRunner().invoke(main, [str(ped), "--output-dir", str(out)])
        assert result.exit_code == 1
        assert result.output.count("Saved: ") == 2
        assert "Error: family 'G': layout went wrong" in result.output
        assert "1 of 3 families could not be drawn." in result.output
