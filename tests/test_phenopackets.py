"""GA4GH Phenopackets v2 pedigrees: read a Family or a Pedigree, write a Pedigree (#108).

The fixtures ``tests/data/format/phenopacket_*.json`` are invented.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from pypedigree import (
    export_fam,
    export_phenopacket_pedigree,
    parse_pedigree,
    parse_pedigree_from_string,
    render_interactive_svg,
)
from pypedigree.cli import main
from pypedigree.errors import PedigreeDataWarning, PedigreeParseError
from pypedigree.models import Phenotype, Sex

DATA = Path(__file__).parent / "data" / "format"
FAMILY = DATA / "phenopacket_family.json"
PEDIGREE = DATA / "phenopacket_pedigree.json"
SNAKE = DATA / "phenopacket_snake_case.json"


def _people(path: Path):
    with pytest.warns(PedigreeDataWarning) if path == FAMILY else _nothing():
        ped = parse_pedigree(path)
    return ped, {i.sample_id: i for i in ped.individuals()}


class _nothing:
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


class TestRead:
    def test_bare_pedigree(self):
        _, people = _people(PEDIGREE)
        assert people["child"].father_id == "father"
        assert people["child"].phenotype == Phenotype.AFFECTED
        assert people["mother"].sex == Sex.FEMALE

    def test_family_proband_deaths_defaults(self):
        _, people = _people(FAMILY)
        assert people["II-1"].is_proband and not people["I-1"].is_proband
        assert people["I-1"].is_deceased and not people["I-2"].is_deceased
        assert people["II-2"].phenotype == Phenotype.UNKNOWN  # affectedStatus left out
        assert people["II-3"].sex == Sex.UNKNOWN  # sex left out
        assert people["I-1"].father_id is None

    def test_snake_case_and_enum_numbers(self):
        """Phenopackets numbers Sex FEMALE = 1, MALE = 2: the reverse of PED."""
        _, people = _people(SNAKE)
        assert people["a"].sex == Sex.MALE  # 2, PED's female
        assert people["b"].sex == Sex.FEMALE  # 1, PED's male
        assert people["a"].phenotype == Phenotype.UNAFFECTED
        assert people["b"].phenotype == Phenotype.AFFECTED
        assert people["c"].is_proband and people["c"].is_deceased
        assert people["c"].mother_id == "b"

    def test_declared_consanguinity_drawn_with_a_warning(self):
        with pytest.warns(PedigreeDataWarning, match="drawn with a double line anyway"):
            ped = parse_pedigree(FAMILY)
        union = ped.families["FAM1"].get_union_for_parents("I-1", "I-2")
        assert union is not None and union.is_consanguinous

    def test_consanguinity_found_in_the_structure_needs_no_warning(self, recwarn):
        persons = [
            ("g1", "0", "0", "MALE"),
            ("g2", "0", "0", "FEMALE"),
            ("a", "g1", "g2", "MALE"),
            ("b", "g1", "g2", "FEMALE"),
            ("x", "0", "0", "FEMALE"),
            ("y", "0", "0", "MALE"),
            ("c1", "a", "x", "MALE"),
            ("c2", "y", "b", "FEMALE"),
            ("p", "c1", "c2", "MALE"),
        ]
        data = {
            "proband": {"subject": {"id": "p"}},
            "consanguinousParents": True,
            "pedigree": {
                "persons": [
                    {
                        "familyId": "FAM1",
                        "individualId": i,
                        "paternalId": f,
                        "maternalId": m,
                        "sex": s,
                    }
                    for i, f, m, s in persons
                ]
            },
        }
        ped = parse_pedigree_from_string(json.dumps(data))
        assert ped.families["FAM1"].get_union_for_parents("c1", "c2").is_consanguinous
        assert not [w for w in recwarn if "consanguineous" in str(w.message)]

    @pytest.mark.parametrize(
        ("text", "message"),
        [
            ("{not json", "Invalid JSON"),
            ('{"cohort": []}', "Family .* or Pedigree"),
            ('{"persons": []}', "no persons"),
            ('{"persons": [{"sex": "MALE"}]}', "no individualId"),
            ('{"persons": [{"individualId": "a", "sex": "M"}]}', "unknown sex 'M'"),
            ('{"persons": [{"individualId": "a", "sex": 7}]}', "unknown sex number 7"),
            ('{"persons": [{"individualId": "a", "affectedStatus": "YES"}]}', "affectedStatus"),
        ],
    )
    def test_errors(self, text, message):
        with pytest.raises(PedigreeParseError, match=message):
            parse_pedigree_from_string(text)

    def test_svg_and_cli(self, tmp_path):
        svg = render_interactive_svg(parse_pedigree(PEDIGREE))
        assert svg.count('class="node"') == 3
        out = tmp_path / "family.png"
        result = CliRunner().invoke(main, [str(FAMILY), "-o", str(out)])
        assert result.exit_code == 0, result.output
        assert "double line anyway" in result.output
        assert out.exists()


class TestWrite:
    def test_round_trip(self, tmp_path):
        ped = parse_pedigree(DATA / "with_header.ped")
        out = export_phenopacket_pedigree(ped, tmp_path / "pedigree.json")
        data = json.loads(out.read_text(encoding="utf-8"))
        assert set(data) == {"persons"}
        assert {p["sex"] for p in data["persons"]} <= {"MALE", "FEMALE", "UNKNOWN_SEX"}
        back = parse_pedigree(out)
        fam_a, fam_b = tmp_path / "a.fam", tmp_path / "b.fam"
        export_fam(ped, fam_a)
        export_fam(back, fam_b)
        assert fam_a.read_text(encoding="utf-8") == fam_b.read_text(encoding="utf-8")

    def test_named_parent_without_a_row_is_written(self, tmp_path):
        """#108 comment: every paternalId / maternalId has its person."""
        text = "FID\tIID\tPAT\tMAT\tSEX\tPHENO\nF\tc\tdad\tmum\t1\t2\nF\tmum\t0\t0\t2\t1\n"
        with pytest.warns(PedigreeDataWarning):
            ped = parse_pedigree_from_string(text)
        out = export_phenopacket_pedigree(ped, tmp_path / "p.json")
        persons = {
            p["individualId"]: p for p in json.loads(out.read_text(encoding="utf-8"))["persons"]
        }
        assert set(persons) == {"c", "mum", "dad"}
        assert persons["dad"]["sex"] == "MALE" and persons["dad"]["affectedStatus"] == "MISSING"
        assert persons["dad"]["paternalId"] == persons["dad"]["maternalId"] == "0"
        named = {p["paternalId"] for p in persons.values()} | {
            p["maternalId"] for p in persons.values()
        }
        assert named - {"0"} <= set(persons)

    def test_defaults_written_out_and_no_phantoms(self, tmp_path):
        text = "FID\tIID\tPAT\tMAT\tSEX\tPHENO\nFAM1\tkid\tdad\t0\t0\t0\n"
        with pytest.warns(PedigreeDataWarning):
            ped = parse_pedigree_from_string(text)
        out = export_phenopacket_pedigree(ped, tmp_path / "p.json")
        persons = json.loads(out.read_text(encoding="utf-8"))["persons"]
        kid = next(p for p in persons if p["individualId"] == "kid")
        assert kid["sex"] == "UNKNOWN_SEX" and kid["affectedStatus"] == "MISSING"
        assert kid["maternalId"] == "0"
        assert not [p for p in persons if p["individualId"].startswith("_phantom_")]
