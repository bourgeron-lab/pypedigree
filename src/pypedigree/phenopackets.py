"""GA4GH Phenopackets v2 pedigrees (#108).

Rare-disease pipelines (Exomiser, LIRICAL, GA4GH registries) exchange
families as Phenopackets v2 (ISO 4454:2022) JSON. A ``Family`` message holds
the proband's phenopacket, the relatives' and a ``pedigree`` whose
``persons`` are PED rows written as JSON; a bare ``Pedigree`` is
``{"persons": [...]}``.

Reading: the persons are rewritten as a pedigree file with a header and
parsed as any other (`parser._parse_text`), so phantoms, validation and
consanguinity detection stay in one place. The proband comes from
``proband.subject.id``, deaths from ``subject.vitalStatus`` of every
phenopacket. proto3 JSON is accepted as protobuf writes it: camelCase or the
original snake_case names, fields at their default value left out, enums as
names or as numbers, read by the Phenopackets numbering (``Sex``:
``FEMALE = 1``, ``MALE = 2``, the reverse of PED), never as PED codes.

Writing: `export_phenopacket_pedigree` writes a ``Pedigree`` message from the
PED fields, without phantoms or duplicates, like `parser.export_fam`.

Standard library only (rule 1).
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import Any

from pypedigree.errors import PedigreeDataWarning, PedigreeParseError
from pypedigree.models import Pedigree, Phenotype, Sex

# Phenopackets v2 enums, by name and by number, as PED codes.
_SEX = {"UNKNOWN_SEX": "0", "FEMALE": "2", "MALE": "1", "OTHER_SEX": "0"}
_SEX_NUMBERS = {0: "UNKNOWN_SEX", 1: "FEMALE", 2: "MALE", 3: "OTHER_SEX"}
_AFFECTED = {"MISSING": "0", "UNAFFECTED": "1", "AFFECTED": "2"}
_AFFECTED_NUMBERS = {0: "MISSING", 1: "UNAFFECTED", 2: "AFFECTED"}
_VITAL_NUMBERS = {0: "UNKNOWN_STATUS", 1: "ALIVE", 2: "DECEASED"}


def load_phenopacket(text: str) -> dict[str, Any] | None:
    """The phenopacket JSON in ``text``, or None if it is not JSON at all.

    Raises:
        PedigreeParseError: If the text starts like JSON but is not valid
            JSON, or is not a ``Family`` or a ``Pedigree``.
    """
    if not text.lstrip().startswith("{"):
        return None
    try:
        data = json.loads(text)
    except json.JSONDecodeError as exc:
        raise PedigreeParseError(f"Invalid JSON: {exc}") from exc
    if not isinstance(data, dict) or not ("pedigree" in data or "persons" in data):
        raise PedigreeParseError(
            "JSON input must be a GA4GH Phenopackets v2 Family (with a 'pedigree') or "
            "Pedigree (with 'persons')"
        )
    return data


def _field(message: Any, camel: str, snake: str, default: Any = None) -> Any:
    """A field of a proto3 JSON message, by its camelCase or snake_case name."""
    if not isinstance(message, dict):
        return default
    if camel in message:
        return message[camel]
    return message.get(snake, default)


def _enum(value: Any, names: dict[str, str], numbers: dict[int, str], what: str, who: str) -> str:
    """An enum value, given by name or by its Phenopackets number."""
    if isinstance(value, bool):
        raise PedigreeParseError(f"Person {who!r}: {what} {value!r} is not a {what} value")
    if isinstance(value, int):
        if value not in numbers:
            raise PedigreeParseError(f"Person {who!r}: unknown {what} number {value}")
        value = numbers[value]
    if not isinstance(value, str) or value.upper() not in names:
        raise PedigreeParseError(
            f"Person {who!r}: unknown {what} {value!r} (expected {', '.join(names)})"
        )
    return names[value.upper()]


def _persons(data: dict[str, Any]) -> list[Any]:
    if "persons" in data:
        persons = data["persons"]
    else:
        persons = _field(data.get("pedigree"), "persons", "persons", [])
    if not isinstance(persons, list) or not persons:
        raise PedigreeParseError("The Phenopackets pedigree has no persons")
    return persons


def _phenopackets(data: dict[str, Any]) -> list[Any]:
    """The proband's phenopacket, then the relatives'."""
    packets = []
    if isinstance(data.get("proband"), dict):
        packets.append(data["proband"])
    relatives = data.get("relatives") or []
    if isinstance(relatives, list):
        packets.extend(p for p in relatives if isinstance(p, dict))
    return packets


def proband_id(data: dict[str, Any]) -> str | None:
    """``proband.subject.id`` of a ``Family``, if any."""
    subject = _field(data.get("proband"), "subject", "subject")
    value = _field(subject, "id", "id")
    return str(value) if value else None


def phenopacket_to_text(data: dict[str, Any]) -> str:
    """Rewrite a Phenopackets ``Family`` or ``Pedigree`` as a pedigree file with a header.

    Raises:
        PedigreeParseError: For a pedigree with no persons, a person without
            an ID, or an unknown sex or affected status.
    """
    deceased = set()
    for packet in _phenopackets(data):
        subject = _field(packet, "subject", "subject")
        vital = _field(subject, "vitalStatus", "vital_status")
        status = _field(vital, "status", "status", "UNKNOWN_STATUS")
        if isinstance(status, int) and not isinstance(status, bool):
            status = _VITAL_NUMBERS.get(status, "UNKNOWN_STATUS")
        if str(status).upper() == "DECEASED" and _field(subject, "id", "id"):
            deceased.add(str(_field(subject, "id", "id")))
    proband = proband_id(data)

    rows = [
        "\t".join(("FID", "IID", "Father", "Mother", "Sex", "Phenotype", "deceased", "proband"))
    ]
    for person in _persons(data):
        iid = _field(person, "individualId", "individual_id")
        if not iid:
            raise PedigreeParseError(f"A Phenopackets person has no individualId: {person!r}")
        iid = str(iid)
        fid = str(_field(person, "familyId", "family_id") or "FAM1")
        father = str(_field(person, "paternalId", "paternal_id") or "0")
        mother = str(_field(person, "maternalId", "maternal_id") or "0")
        sex = _enum(_field(person, "sex", "sex", 0), _SEX, _SEX_NUMBERS, "sex", iid)
        status = _enum(
            _field(person, "affectedStatus", "affected_status", 0),
            _AFFECTED,
            _AFFECTED_NUMBERS,
            "affectedStatus",
            iid,
        )
        cells = (
            fid,
            iid,
            father,
            mother,
            sex,
            status,
            "1" if iid in deceased else "0",
            "1" if iid == proband else "0",
        )
        if any("\t" in c or "\n" in c for c in cells):
            raise PedigreeParseError(f"Person {iid!r}: an ID contains a tab or a line break")
        rows.append("\t".join(cells))
    return "\n".join(rows)


def mark_declared_consanguinity(pedigree: Pedigree, data: dict[str, Any]) -> None:
    """Draw the proband's parents as consanguineous when the ``Family`` says so.

    pypedigree finds consanguinity from the structure; when the parents'
    common ancestors are not in the pedigree, their union is still drawn with
    a double line (NSGC allows it for known but undrawn consanguinity), with
    a warning saying why.
    """
    if not _field(data, "consanguinousParents", "consanguinous_parents", False):
        return
    pid = proband_id(data)
    for family in pedigree.families.values():
        proband = family.individuals.get(pid or "")
        if proband is None:
            continue
        union = None
        if proband.father_id and proband.mother_id:
            union = family.get_union_for_parents(proband.father_id, proband.mother_id)
        if union is None:
            warnings.warn(
                f"The Phenopacket declares the parents of {pid} consanguineous, but "
                f"{pid} has no two parents in the pedigree: nothing is drawn.",
                PedigreeDataWarning,
                stacklevel=3,
            )
        elif not union.is_consanguinous:
            union.is_consanguinous = True
            warnings.warn(
                f"The Phenopacket declares the parents of {pid} consanguineous, but the "
                "pedigree has none of their common ancestors: their union is drawn with a "
                "double line anyway.",
                PedigreeDataWarning,
                stacklevel=3,
            )
        return
    warnings.warn(
        "The Phenopacket declares the proband's parents consanguineous, but its proband "
        "is not in the pedigree: nothing is drawn.",
        PedigreeDataWarning,
        stacklevel=3,
    )


_SEX_NAMES = {Sex.MALE: "MALE", Sex.FEMALE: "FEMALE", Sex.UNKNOWN: "UNKNOWN_SEX"}
_AFFECTED_NAMES = {
    Phenotype.AFFECTED: "AFFECTED",
    Phenotype.UNAFFECTED: "UNAFFECTED",
    Phenotype.UNKNOWN: "MISSING",
}


def export_phenopacket_pedigree(pedigree: Pedigree, file_path: str | Path) -> Path:
    """Write a GA4GH Phenopackets v2 ``Pedigree`` message (``{"persons": [...]}``).

    From the PED fields, without phantom parents or duplicates, as
    `export_fam` does; unknown parents are ``"0"``, and ``UNKNOWN_SEX`` and
    ``MISSING`` are written out for readers that do not apply proto3
    defaults. A whole ``Family`` is not written: it needs phenopackets with
    clinical data and metadata pypedigree does not hold.

    Args:
        pedigree: The pedigree to export.
        file_path: Output file path (``.json``).

    Returns:
        The Path of the written file.
    """
    from pypedigree.parser import fam_rows  # the same rows as export_fam

    persons = [
        {
            "familyId": fid,
            "individualId": iid,
            "paternalId": father,
            "maternalId": mother,
            "sex": _SEX_NAMES.get(ind.sex, "UNKNOWN_SEX"),
            "affectedStatus": _AFFECTED_NAMES.get(ind.phenotype, "MISSING"),
        }
        for fid, iid, father, mother, ind in fam_rows(pedigree)
    ]
    path = Path(file_path)
    text = json.dumps({"persons": persons}, indent=2, ensure_ascii=False)
    path.write_text(text + "\n", encoding="utf-8", newline="\n")
    return path
