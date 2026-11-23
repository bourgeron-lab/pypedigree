"""Read CanRisk pedigree files (#109).

Cancer genetics services build families in the CanRisk tool (BOADICEA model)
and save them as CanRisk files: a first line ``##CanRisk <version>`` (1 to
4), ``##key=value`` header lines, a ``##FamID ...`` line naming the columns,
then one row per person. A file may hold several families, each with its own
header block. The columns are read by the names of that line, so versions 1
(no ``BARD1``), 2 and 3 (no ``HOXB13``) and 4 share this code.

The file is rewritten as a pedigree file with a header, then parsed as any
other (`parser._parse_text`), so phantoms, validation, twins and consanguinity
stay in one place:

- ``FamID``, ``IndivID``, ``FathID``, ``MothID``, ``Sex`` (M/F), ``Dead``,
  ``Yob`` (year of birth) and ``MZtwin`` (a code shared by identical twins)
  map onto the usual columns. ``Name`` is left out: in CanRisk it is often a
  real name or a role.
- ``Target`` (the person whose risk is computed) is the consultand; an
  affected target is also the proband.
- Each cancer someone has gets a binary column: breast (``BC1`` or ``BC2``),
  ovarian, prostate, pancreatic. Ovarian cancer does not apply to people
  assigned male at birth, prostate cancer to people assigned female at birth:
  those cells are "not applicable" (#107), unless the file gives an age,
  which is kept with a warning.
- Annotation lines, in NSGC order: age ("45y", "d. 61y"), cancers with
  their ages ("BC 40, BC2 52"), one line per gene tested ("BRCA1 +",
  "BRCA2 −", "HOXB13 + (het)"), breast pathology ("ER+ PR+ HER2−") and
  "Ashkenazi Jewish".
- ``##PRS_BC``, ``##PRS_OC`` and ``##PRS_PC`` (``alpha=…, zscore=…``) give the
  target's polygenic risk score as a z-score, in a ``PGS_`` column drawn on
  the z-score scale (the default, #106).

Other header keys (menarche, parity, BMI...) are risk factors, not pedigree
data, and are ignored. Only the standard library is used (rule 1).
"""

from __future__ import annotations

import re
import warnings
from collections.abc import Sequence
from dataclasses import dataclass, field

from pypedigree.errors import PedigreeDataWarning, PedigreeParseError

#: The first line of a CanRisk file starts with this (case-insensitive).
CANRISK_MAGIC = "##canrisk"
CANRISK_VERSIONS = ("1", "2", "3", "4")

MINUS = "−"

# Cancer columns: (pypedigree column, CanRisk columns, abbreviation of each in
# the annotation, the sex assigned at birth it applies to, or None for both).
_CANCERS: tuple[tuple[str, tuple[str, ...], tuple[str, ...], str | None], ...] = (
    ("breast_cancer", ("BC1", "BC2"), ("BC", "BC2"), None),
    ("ovarian_cancer", ("OC",), ("OC",), "F"),
    ("prostate_cancer", ("PRO",), ("PrC",), "M"),
    ("pancreatic_cancer", ("PAN",), ("PaC",), None),
)
_GENES = ("BRCA1", "BRCA2", "PALB2", "ATM", "CHEK2", "BARD1", "RAD51D", "RAD51C", "BRIP1")
_HOXB13 = "HOXB13"
_PATHOLOGY = "ER:PR:HER2:CK14:CK56"
_PATHOLOGY_MARKERS = ("ER", "PR", "HER2", "CK14", "CK56")
_PRS = {
    "prs_bc": "PGS_breast_cancer",
    "prs_oc": "PGS_ovarian_cancer",
    "prs_pc": "PGS_prostate_cancer",
}
_REQUIRED = ("FamID", "IndivID", "FathID", "MothID", "Sex")
# The columns of each format, in order, for a file without a ``##FamID`` line,
# which the CanRisk services read by position (#119): 28 fields in format 4,
# 27 in formats 2 and 3 (no HOXB13), 26 in format 1 (no BARD1 either).
_COLUMNS_V4 = (
    "FamID", "Name", "Target", "IndivID", "FathID", "MothID", "Sex", "MZtwin", "Dead", "Age",
    "Yob", "BC1", "BC2", "OC", "PRO", "PAN", "Ashkn", *_GENES, _HOXB13, _PATHOLOGY,
)  # fmt: skip
_DEFAULT_COLUMNS = {
    "4": _COLUMNS_V4,
    "3": tuple(c for c in _COLUMNS_V4 if c != _HOXB13),
    "2": tuple(c for c in _COLUMNS_V4 if c != _HOXB13),
    "1": tuple(c for c in _COLUMNS_V4 if c not in (_HOXB13, "BARD1")),
}
_ZSCORE = re.compile(r"zscore\s*=\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)", re.IGNORECASE)


def is_canrisk(text: str) -> bool:
    """Whether a text is a CanRisk file: its first line starts with ``##CanRisk``."""
    for line in text.splitlines():
        if line.strip():
            return line.strip().lower().startswith(CANRISK_MAGIC)
    return False


@dataclass
class _Person:
    line: int
    cells: dict[str, str]
    prs: dict[str, float] = field(default_factory=dict)

    def get(self, column: str) -> str:
        return self.cells.get(column, "0").strip() or "0"


def canrisk_to_text(text: str) -> str:
    """Rewrite a CanRisk file as a tab-separated pedigree file with a header.

    Args:
        text: The CanRisk file's content.

    Returns:
        The same families as pedigree file text, for `parser._parse_text`.

    Raises:
        PedigreeParseError: For an unsupported version, a missing column line
            or column, a row of the wrong length, or a value CanRisk does not
            allow (naming the line).
    """
    people = _read_blocks(text)
    if not people:
        raise PedigreeParseError("CanRisk file has no individuals")
    return _to_pedigree_text(people)


def _read_blocks(text: str) -> list[_Person]:
    """Every person of every family block, with their block's PRS."""
    people: list[_Person] = []
    columns: list[str] | None = None
    default_columns: list[str] | None = None
    prs: dict[str, float] = {}
    block: list[_Person] = []

    def close() -> None:
        _give_prs(block, prs)

    for number, raw in enumerate(text.splitlines(), start=1):
        line = raw.rstrip("\r\n")
        if not line.strip():
            continue
        stripped = line.strip()
        if stripped.lower().startswith(CANRISK_MAGIC):
            close()
            version = stripped[len(CANRISK_MAGIC) :].strip()
            if version.split(".")[0] not in CANRISK_VERSIONS:
                raise PedigreeParseError(
                    f"Line {number}: unsupported CanRisk version {version!r} "
                    f"(expected {', '.join(CANRISK_VERSIONS)})"
                )
            columns, prs, block = None, {}, []
            default_columns = list(_DEFAULT_COLUMNS.get(version.split(".")[0], ()))
            continue
        if stripped.startswith("##"):
            body = stripped[2:]
            words = body.split()
            if words and words[0].lower() == "famid":
                columns = words
                missing = [c for c in _REQUIRED if c not in columns]
                if missing:
                    raise PedigreeParseError(
                        f"Line {number}: CanRisk column line lacks {', '.join(missing)}"
                    )
            elif "=" in body:
                key, value = body.split("=", 1)
                if key.strip().lower() in _PRS:
                    prs[key.strip().lower()] = _prs_zscore(value, number)
            continue
        if columns is None:
            if not default_columns:
                raise PedigreeParseError(
                    f"Line {number}: a CanRisk row before the '##FamID ...' column line"
                )
            # No column line: the format's columns by position, as CanRisk reads them.
            columns = default_columns
        # Any run of whitespace separates fields, as the CanRisk services split
        # them (no field holds a space): tabs and spaces mixed in a row, or a
        # trailing tab from a spreadsheet, are read (#119).
        cells = line.split()
        if len(cells) != len(columns):
            raise PedigreeParseError(
                f"Line {number}: {len(cells)} CanRisk fields, the column line names {len(columns)}"
            )
        person = _Person(number, dict(zip(columns, (c.strip() for c in cells), strict=True)))
        people.append(person)
        block.append(person)
    close()
    return people


def _prs_zscore(value: str, line: int) -> float:
    """The z-score of a ``##PRS_..=alpha=…, zscore=…`` line."""
    match = _ZSCORE.search(value)
    if match is None:
        raise PedigreeParseError(f"Line {line}: CanRisk PRS line without a zscore: {value!r}")
    return float(match.group(1))


def _give_prs(block: Sequence[_Person], prs: dict[str, float]) -> None:
    """Attach a header block's PRS to the target of the family after it.

    As the CanRisk services read a file, the header lines belong to the family
    whose rows follow them; another family in the same block has none of its
    own, which a warning says (#119).
    """
    if not prs:
        return
    targets = [p for p in block if p.get("Target") == "1"]
    if not targets:
        warnings.warn(
            "A CanRisk family has a PRS line but no target (Target 1): the PRS is not drawn.",
            PedigreeDataWarning,
            stacklevel=4,
        )
        return
    first = block[0].get("FamID")
    owner = next((p for p in targets if p.get("FamID") == first), None)
    others = sorted({p.get("FamID") for p in targets if p.get("FamID") != first})
    if owner is not None:
        owner.prs.update(prs)
    if others or owner is None:
        warnings.warn(
            f"The CanRisk PRS lines before family {first} are that family's, as CanRisk reads "
            f"them: the targets of {', '.join(others)} get no PRS. Give each family its own "
            "'##CanRisk' header block with its PRS lines.",
            PedigreeDataWarning,
            stacklevel=4,
        )


def _sex(person: _Person) -> str:
    sex = person.get("Sex").upper()
    if sex not in ("M", "F"):
        raise PedigreeParseError(
            f"Line {person.line}: CanRisk Sex must be M or F, got {person.get('Sex')!r}"
        )
    return sex


def _cancer_age(person: _Person, column: str) -> str | None:
    """None when unaffected; the age at diagnosis, or "?" (AU, age unknown)."""
    value = person.get(column).upper()
    if value == "0":
        return None
    if value == "AU":
        return "?"
    if value.isdigit():
        return str(int(value))
    raise PedigreeParseError(
        f"Line {person.line}: CanRisk {column} must be 0, an age or AU, got {value!r}"
    )


def _gene_text(person: _Person, gene: str) -> str | None:
    """ "BRCA1 +", "BRCA2 −", "HOXB13 + (het)", or None when untested."""
    value = person.get(gene).upper()
    if value == "0":
        return None
    parts = value.split(":")
    allowed = ("0", "P", "N") if gene != _HOXB13 else ("0", "HET", "HOM", "N")
    if len(parts) != 2 or parts[0] not in ("0", "S", "T") or parts[1] not in allowed:
        raise PedigreeParseError(
            f"Line {person.line}: CanRisk {gene} must be type:result "
            f"(type 0, S or T; result {', '.join(allowed)}), got {value!r}"
        )
    kind, result = parts
    if kind == "0" or result == "0":
        return None
    if result == "N":
        return f"{gene} {MINUS}"
    if result in ("HET", "HOM"):
        return f"{gene} + ({result.lower()})"
    return f"{gene} +"


def _pathology_text(person: _Person) -> str | None:
    """ "ER+ PR+ HER2−" for the markers given, or None."""
    value = person.get(_PATHOLOGY).upper()
    if value == "0":
        return None
    parts = value.split(":")
    if len(parts) != len(_PATHOLOGY_MARKERS) or any(p not in ("0", "N", "P") for p in parts):
        raise PedigreeParseError(
            f"Line {person.line}: CanRisk {_PATHOLOGY} must be five of 0, N or P "
            f"separated by colons, got {value!r}"
        )
    shown = [
        f"{name}{'+' if p == 'P' else MINUS}"
        for name, p in zip(_PATHOLOGY_MARKERS, parts, strict=True)
        if p != "0"
    ]
    return " ".join(shown) or None


def _whole_number(person: _Person, column: str) -> int:
    value = person.get(column)
    if not value.isdigit():
        raise PedigreeParseError(
            f"Line {person.line}: CanRisk {column} must be a whole number, got {value!r}"
        )
    return int(value)


def _to_pedigree_text(people: list[_Person]) -> str:
    """The pedigree file text: structure, cancers, PGS, then annotations."""
    sexes = {id(p): _sex(p) for p in people}
    ages = {
        (id(p), column): _cancer_age(p, column)
        for p in people
        for _, columns, _, _ in _CANCERS
        for column in columns
        if column in p.cells
    }

    def affected(p: _Person, columns: Sequence[str]) -> bool:
        return any(ages.get((id(p), c)) is not None for c in columns)

    cancers = [c for c in _CANCERS if any(affected(p, c[1]) for p in people)]
    genes = [g for g in (*_GENES, _HOXB13) if any(_gene_text(p, g) for p in people)]
    pgs = [name for key, name in _PRS.items() if any(key in p.prs for p in people)]

    twins: dict[tuple[str, str], list[str]] = {}
    for p in people:
        code = p.get("MZtwin").upper()
        if code != "0":
            twins.setdefault((p.get("FamID"), code), []).append(p.get("IndivID"))

    header = [
        "FID",
        "IID",
        "Father",
        "Mother",
        "Sex",
        "Phenotype",
        "dob",
        "twins",
        "deceased",
        "proband",
        "consultand",
    ]
    header += [f"b_{name}" for name, _, _, _ in cancers]
    header += pgs
    header += ["a_age", "a_cancer"] + [f"a_{g}" for g in genes] + ["a_pathology", "a_ancestry"]

    rows = ["\t".join(header)]
    for p in people:
        sex = sexes[id(p)]
        fid, iid = p.get("FamID"), p.get("IndivID")
        dead = p.get("Dead") == "1"
        any_cancer = any(affected(p, c[1]) for c in _CANCERS)
        target = p.get("Target") == "1"
        yob = _whole_number(p, "Yob") if "Yob" in p.cells else 0
        code = p.get("MZtwin").upper()
        co_twins = [t for t in twins.get((fid, code), []) if t != iid] if code != "0" else []
        row = [
            fid,
            iid,
            p.get("FathID"),
            p.get("MothID"),
            "1" if sex == "M" else "2",
            "2" if any_cancer else "1",
            str(yob) if yob else "0",
            ",".join(f"MZ:{t}" for t in co_twins) or "0",
            "1" if dead else "0",
            "1" if target and any_cancer else "0",
            "1" if target else "0",
        ]
        for name, columns, _, applies_to in cancers:
            if affected(p, columns):
                if applies_to is not None and applies_to != sex:
                    warnings.warn(
                        f"Line {p.line}: {iid} (sex {sex} in CanRisk) has "
                        f"{name.replace('_', ' ')}; kept as affected.",
                        PedigreeDataWarning,
                        stacklevel=3,
                    )
                row.append("2")
            elif applies_to is not None and applies_to != sex:
                row.append("N/A")  # cannot apply (#107)
            else:
                row.append("1")
        for key, name in _PRS.items():
            if name in pgs:
                row.append(f"{p.prs[key]:g}" if key in p.prs else "NA")
        age = _whole_number(p, "Age") if "Age" in p.cells else 0
        row.append((f"d. {age}y" if dead else f"{age}y") if age else "")
        diagnoses = [
            f"{abbreviation} {ages[(id(p), column)]}"
            for _, columns, abbreviations, _ in _CANCERS
            for column, abbreviation in zip(columns, abbreviations, strict=True)
            if ages.get((id(p), column)) is not None
        ]
        row.append(", ".join(diagnoses))
        row += [_gene_text(p, g) or "" for g in genes]
        pathology = _pathology_text(p) if _PATHOLOGY in p.cells else None  # checked for all
        row.append((pathology or "") if affected(p, ("BC1", "BC2")) else "")
        row.append("Ashkenazi Jewish" if p.get("Ashkn") == "1" else "")
        rows.append("\t".join(row))
    return "\n".join(rows)
