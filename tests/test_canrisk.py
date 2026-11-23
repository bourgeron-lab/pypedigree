"""CanRisk pedigree files (v1 to v4) read onto pypedigree's columns (#109).

The fixtures in ``tests/data/format/canrisk_*.txt`` are written by hand from
the CanRisk v3 and v4 specifications; every person is invented.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pytest
from click.testing import CliRunner

from pypedigree import (
    PedigreeStyle,
    parse_pedigree,
    parse_pedigree_from_string,
    render_interactive_svg,
    render_pedigree,
)
from pypedigree.canrisk import canrisk_to_text, is_canrisk
from pypedigree.cli import main
from pypedigree.errors import PedigreeDataWarning, PedigreeParseError
from pypedigree.models import Phenotype, Sex

DATA = Path(__file__).parent / "data" / "format"
V4 = DATA / "canrisk_v4.txt"

COLUMNS = [
    "FamID", "Name", "Target", "IndivID", "FathID", "MothID", "Sex", "MZtwin", "Dead", "Age",
    "Yob", "BC1", "BC2", "OC", "PRO", "PAN", "Ashkn", "BRCA1", "BRCA2", "PALB2", "ATM",
    "CHEK2", "BARD1", "RAD51D", "RAD51C", "BRIP1", "HOXB13", "ER:PR:HER2:CK14:CK56",
]  # fmt: skip


def _row(**values: object) -> str:
    cells = {c: "0" for c in COLUMNS}
    for gene in COLUMNS[17:27]:
        cells[gene] = "0:0"
    cells["ER:PR:HER2:CK14:CK56"] = "0:0:0:0:0"
    cells.update({k: str(v) for k, v in values.items()})
    return "\t".join(cells[c] for c in COLUMNS)


def _file(*rows: str, header: tuple[str, ...] = ()) -> str:
    return "\n".join(["##CanRisk 4.0", *header, "##" + "\t".join(COLUMNS), *rows])


def _people(path: Path = V4):
    ped = parse_pedigree(path)
    return ped, {i.sample_id: i for i in ped.individuals()}


class TestDetection:
    def test_first_line(self):
        assert is_canrisk(V4.read_text(encoding="utf-8"))
        assert is_canrisk("\n  ##canrisk 3\n")
        assert not is_canrisk("FID\tIID\n##CanRisk 4.0")

    @pytest.mark.parametrize("name", ["canrisk_v2.txt", "canrisk_v3.txt", "canrisk_v4.txt"])
    def test_versions(self, name):
        _, people = _people(DATA / name)
        assert people["t1"].is_proband
        assert people["t1"].notation.is_consultand

    def test_version_1_without_bard1(self):
        cols = [c for c in COLUMNS if c not in ("BARD1", "HOXB13")]
        full = dict(
            zip(COLUMNS, _row(FamID="FAM1", IndivID="a", Sex="F", BC1=44).split("\t"), strict=True)
        )
        row = "\t".join(full[c] for c in cols)
        text = "\n".join(["##CanRisk 1.0", "##" + "\t".join(cols), row])
        (person,) = parse_pedigree_from_string(text).individuals()
        assert person.extra_phenotypes_raw["breast_cancer"] == 2

    def test_unsupported_version(self):
        with pytest.raises(PedigreeParseError, match="unsupported CanRisk version"):
            parse_pedigree_from_string("##CanRisk 9.0\n##FamID\tIndivID\n")


class TestStructure:
    def test_ids_sex_parents_deaths_births(self):
        ped, people = _people()
        assert set(ped.families) == {"FAM1", "FAM2"}
        assert people["t1"].father_id == "f1" and people["t1"].mother_id == "m1"
        assert people["f1"].sex == Sex.MALE and people["m1"].sex == Sex.FEMALE
        assert people["m1"].is_deceased and not people["t1"].is_deceased
        assert people["t1"].date_of_birth == "1986"

    def test_names_are_not_used(self):
        text = canrisk_to_text(V4.read_text(encoding="utf-8"))
        assert "grandmother" not in text and "\tme\t" not in text

    def test_mz_twins(self):
        _, people = _people()
        assert [(t.twin_type.value, t.sibling_id) for t in people["s1"].twin_info] == [("MZ", "s2")]

    def test_target_is_consultand_and_affected_target_proband(self):
        _, people = _people()
        assert people["t1"].notation.is_consultand and people["t1"].is_proband
        assert people["t2"].notation.is_consultand and not people["t2"].is_proband
        assert not people["m1"].notation.is_consultand


class TestCancers:
    def test_only_cancers_someone_has(self):
        text = _file(
            _row(FamID="FAM1", IndivID="a", Sex="F", BC1=40, Target=1),
            _row(FamID="FAM1", IndivID="b", Sex="M"),
        )
        ped = parse_pedigree_from_string(text)
        assert ped.phenotype_names == ["breast_cancer"]

    def test_columns_and_affection(self):
        ped, people = _people()
        assert ped.phenotype_names[:4] == [
            "breast_cancer",
            "ovarian_cancer",
            "prostate_cancer",
            "pancreatic_cancer",
        ]
        assert people["t1"].extra_phenotypes_raw["breast_cancer"] == 2
        assert people["s1"].extra_phenotypes_raw["breast_cancer"] == 1
        assert people["gm1"].extra_phenotypes_raw["breast_cancer"] == 2  # AU
        assert people["t1"].phenotype == Phenotype.AFFECTED
        assert people["s1"].phenotype == Phenotype.UNAFFECTED

    def test_sex_limited_cancers_not_applicable(self):
        _, people = _people()
        assert "ovarian_cancer" in people["f1"].not_applicable
        assert "prostate_cancer" in people["m1"].not_applicable
        assert "prostate_cancer" not in people["f1"].not_applicable
        assert not {"breast_cancer", "pancreatic_cancer"} & people["f1"].not_applicable

    def test_an_age_given_is_kept_with_a_warning(self):
        text = _file(_row(FamID="FAM1", IndivID="a", Sex="M", OC=50, Target=1))
        with pytest.warns(PedigreeDataWarning, match="ovarian cancer; kept as affected"):
            (person,) = parse_pedigree_from_string(text).individuals()
        assert person.extra_phenotypes_raw["ovarian_cancer"] == 2
        assert not person.not_applicable


class TestAnnotations:
    def test_nsgc_order(self):
        _, people = _people()
        assert list(people["t1"].annotations.values()) == [
            "38y",
            "BC 35, BC2 37",
            "BRCA1 +",
            "ER− PR− HER2−",
            "Ashkenazi Jewish",
        ]
        assert list(people["m1"].annotations.values())[:3] == ["d. 61y", "OC 52", "BRCA1 +"]
        assert "BC ?" in people["gm1"].annotations.values()
        assert "PrC 61" in people["f1"].annotations.values()
        assert "HOXB13 + (het)" in people["f1"].annotations.values()
        assert "BRCA1 −" in people["b1"].annotations.values()

    def test_untested_genes_write_nothing(self):
        _, people = _people()
        assert not [a for a in people["s2"].annotations.values() if "BRCA" in a]

    def test_pathology_only_with_breast_cancer(self):
        text = _file(
            _row(FamID="FAM1", IndivID="a", Sex="F", OC=50, **{"ER:PR:HER2:CK14:CK56": "P:P:0:0:0"})
        )
        (person,) = parse_pedigree_from_string(text).individuals()
        assert not [a for a in person.annotations.values() if "ER" in a]


class TestPrs:
    def test_target_gets_the_z_scores(self):
        ped, people = _people()
        assert people["t1"].extra_phenotypes_raw["PGS_breast_cancer"] == 1.8
        assert people["t1"].extra_phenotypes_raw["PGS_ovarian_cancer"] == -0.6  # no space
        assert people["m1"].extra_phenotypes_raw["PGS_breast_cancer"] is None
        assert people["t2"].extra_phenotypes_raw["PGS_prostate_cancer"] == 0.4  # its own block
        assert people["t1"].extra_phenotypes_raw["PGS_prostate_cancer"] is None

    def test_drawn_on_the_z_scale_without_warning(self, recwarn):
        svg = render_interactive_svg(parse_pedigree(V4))
        assert "PGS breast cancer: z = +1.80" in svg
        assert not [w for w in recwarn if "z-scores" in str(w.message)]

    def test_prs_without_target_warns(self):
        text = _file(
            _row(FamID="FAM1", IndivID="a", Sex="F"), header=("##PRS_BC=alpha=0.4, zscore=1",)
        )
        with pytest.warns(PedigreeDataWarning, match="no target"):
            parse_pedigree_from_string(text)


class TestErrors:
    @pytest.mark.parametrize(
        ("change", "message"),
        [
            ({"Sex": "X"}, "Sex must be M or F"),
            ({"BC1": "old"}, "BC1 must be 0, an age or AU"),
            ({"BRCA1": "T:Q"}, "BRCA1 must be type:result"),
            ({"HOXB13": "T:P"}, "HOXB13 must be type:result"),
            ({"ER:PR:HER2:CK14:CK56": "P:P"}, "five of 0, N or P"),
            ({"Age": "x"}, "Age must be a whole number"),
        ],
    )
    def test_bad_values_name_the_line(self, change, message):
        text = _file(_row(**{"FamID": "FAM1", "IndivID": "a", "Sex": "F", **change}))
        with pytest.raises(PedigreeParseError, match=f"Line 3: .*{message}"):
            parse_pedigree_from_string(text)

    def test_row_length(self):
        with pytest.raises(PedigreeParseError, match="CanRisk fields"):
            parse_pedigree_from_string(_file("FAM1\ta\t0"))


class TestFilesCanRiskAccepts:
    """#119: files the CanRisk services read, as synthetic rows of the same shape."""

    def test_no_column_line_read_by_position(self):
        text = "##CanRisk 4.0\n" + _row(FamID="FAM1", IndivID="a", Sex="F", BC1=44)
        (person,) = parse_pedigree_from_string(text).individuals()
        assert person.sample_id == "a" and person.extra_phenotypes_raw["breast_cancer"] == 2

    @pytest.mark.parametrize(
        ("version", "drop"), [("2.0", ("HOXB13",)), ("1.0", ("HOXB13", "BARD1"))]
    )
    def test_no_column_line_older_formats(self, version, drop):
        cells = dict(
            zip(COLUMNS, _row(FamID="FAM1", IndivID="a", Sex="M", PRO=61).split("\t"), strict=True)
        )
        row = "\t".join(v for c, v in cells.items() if c not in drop)
        (person,) = parse_pedigree_from_string(f"##CanRisk {version}\n{row}").individuals()
        assert person.extra_phenotypes_raw["prostate_cancer"] == 2

    def test_tabs_and_spaces_mixed(self):
        row = _row(FamID="FAM1", IndivID="a", Sex="F", BC1=44)
        head, last = row.rsplit("\t", 1)
        (person,) = parse_pedigree_from_string(_file(head + "   " + last)).individuals()
        assert person.extra_phenotypes_raw["breast_cancer"] == 2

    def test_trailing_tab(self):
        (person,) = parse_pedigree_from_string(
            _file(_row(FamID="FAM1", IndivID="a", Sex="F") + "\t")
        ).individuals()
        assert person.sample_id == "a"

    def test_one_prs_line_before_two_families_warns(self):
        text = _file(
            _row(FamID="F1", IndivID="a", Sex="F", Target=1),
            _row(FamID="F2", IndivID="b", Sex="F", Target=1),
            header=("##PRS_BC=alpha=0.45,zscore=1.5",),
        )
        with pytest.warns(PedigreeDataWarning, match="targets of F2 get no PRS"):
            ped = parse_pedigree_from_string(text)
        people = {i.sample_id: i for i in ped.individuals()}
        assert people["a"].extra_phenotypes_raw["PGS_breast_cancer"] == 1.5
        assert people["b"].extra_phenotypes_raw["PGS_breast_cancer"] is None


class TestOutputs:
    def test_png(self):
        fig = render_pedigree(parse_pedigree(V4))
        try:
            texts = [t.get_text() for ax in fig.axes for t in ax.texts]
            assert "BRCA1 +" in texts and "Not applicable" not in texts  # N/A drawn white (#120)
        finally:
            plt.close(fig)

    def test_bands_too(self):
        svg = render_interactive_svg(
            parse_pedigree(V4), style=PedigreeStyle(phenotype_layout="bands")
        )
        assert "ovarian cancer: not applicable" in svg

    def test_cli(self, tmp_path):
        out = tmp_path / "canrisk.svg"
        result = CliRunner().invoke(main, [str(V4), "--interactive", "-o", str(out)])
        assert result.exit_code == 0, result.output
        assert "BRCA1 +" in out.read_text(encoding="utf-8")
