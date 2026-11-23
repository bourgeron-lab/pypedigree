"""pTDT deviations computed from the parents' PGS (#52, second step).

(PGS_child - PGS_mid-parent) / SD(mid-parent), SD from the style (sqrt(1/2) by
default). All pedigrees here are synthetic.
"""

from __future__ import annotations

import math
import warnings

import pytest
from click.testing import CliRunner

from pypedigree import PedigreeStyle, parse_pedigree_from_string, render_interactive_svg
from pypedigree.cli import main
from pypedigree.errors import PedigreeDataWarning, PedigreeError
from pypedigree.transmission import ptdt_name, with_transmission

FAMILY = (
    "FID\tIID\tPAT\tMAT\tSEX\tPHENO\tb_autism\tPGS_autism\n"
    "F\tf\t0\t0\t1\t1\t1\t0.4\n"
    "F\tm\t0\t0\t2\t1\t1\t-0.2\n"
    "F\tc1\tf\tm\t1\t2\t2\t1.6\n"
    "F\tc2\tf\tm\t2\t1\t1\tNA\n"
    "F\tx\tf\t0\t1\t1\t1\t0.5\n"
)
ON = PedigreeStyle(pgs_transmission=True)


def _people(text: str = FAMILY, style: PedigreeStyle = ON):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PedigreeDataWarning)  # x's missing mother
        ped = parse_pedigree_from_string(text)
    out, computed = with_transmission(ped, style)
    return ped, out, computed, {i.sample_id: i for i in out.individuals()}


def test_name():
    assert ptdt_name("PGS_autism") == "pTDT_autism" and ptdt_name("PGS") == "pTDT"


def test_deviation_of_each_child():
    _, out, computed, people = _people()
    assert computed == {"pTDT_autism"}
    assert out.phenotype_names == ["autism", "PGS_autism", "pTDT_autism"]
    assert people["c1"].extra_phenotypes_raw["pTDT_autism"] == pytest.approx(1.5 / math.sqrt(0.5))
    assert people["c2"].extra_phenotypes_raw["pTDT_autism"] is None  # no PGS of its own
    assert people["x"].extra_phenotypes_raw["pTDT_autism"] is None  # one parent only
    assert people["f"].extra_phenotypes_raw["pTDT_autism"] is None  # founders: none


def test_sd_from_the_style():
    _, _, _, people = _people(style=PedigreeStyle(pgs_transmission=True, ptdt_midparent_sd=0.5))
    assert people["c1"].extra_phenotypes_raw["pTDT_autism"] == pytest.approx(3.0)


@pytest.mark.parametrize("sd", [0, -1, float("inf"), float("nan")])
def test_bad_sd(sd):
    with pytest.raises(PedigreeError, match="ptdt_midparent_sd"):
        _people(style=PedigreeStyle(pgs_transmission=True, ptdt_midparent_sd=sd))


def test_off_by_default_and_parsed_pedigree_untouched():
    ped, out, computed, _ = _people(style=PedigreeStyle())
    assert out is ped and not computed
    ped, out, _, _ = _people()
    assert "pTDT_autism" not in ped.phenotype_names
    assert "pTDT_autism" not in ped.families["F"].individuals["c1"].extra_phenotypes_raw


def test_a_parent_drawn_twice():
    """f has two partners: c3's father is f's duplicate, whose value is f's."""
    text = FAMILY + "F\tw\t0\t0\t2\t1\t1\t1.0\nF\tc3\tf\tw\t2\t1\t1\t0.7\n"
    _, _, _, people = _people(text)
    assert any(p.is_duplicate for p in people.values())
    assert people["c3"].extra_phenotypes_raw["pTDT_autism"] == pytest.approx(0.0)


def test_given_column_kept_with_a_warning():
    text = (
        FAMILY.replace("PGS_autism\n", "PGS_autism\tpTDT_autism\n")
        .replace("\n", "\tNA\n")
        .replace("PGS_autism\tpTDT_autism\tNA\n", "PGS_autism\tpTDT_autism\n")
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PedigreeDataWarning)
        ped = parse_pedigree_from_string(text)
    with pytest.warns(PedigreeDataWarning, match="drawn as given"):
        out, computed = with_transmission(ped, ON)
    assert out is ped and not computed


def test_drawn_without_the_mixed_curves_warning(recwarn):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PedigreeDataWarning)
        ped = parse_pedigree_from_string(FAMILY)
    svg = render_interactive_svg(ped, style=ON)
    assert "pTDT = +2.12" in svg
    assert not [w for w in recwarn if "both PGS and pTDT" in str(w.message)]


def test_cli_and_output_dir(tmp_path):
    ped = tmp_path / "f.ped"
    ped.write_text(FAMILY + FAMILY.split("\n", 1)[1].replace("F\t", "G\t"), encoding="utf-8")
    out = tmp_path / "f.svg"
    result = CliRunner().invoke(main, [str(ped), "--ptdt", "--interactive", "-o", str(out)])
    assert result.exit_code == 0, result.output
    assert out.read_text(encoding="utf-8").count("pTDT = +2.12") == 2
    outdir = tmp_path / "each"
    result = CliRunner().invoke(
        main, [str(ped), "--ptdt", "--interactive", "--output-dir", str(outdir)]
    )
    assert result.exit_code == 0, result.output
    files = sorted(outdir.glob("*.svg"))
    assert len(files) == 2 and all("pTDT = +2.12" in p.read_text(encoding="utf-8") for p in files)
