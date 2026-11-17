#!/usr/bin/env python3
"""Generate the synthetic test fixtures under `complex/` and `phenotype/`.

Every pedigree this script writes is **fabricated**. Topologies are built from
the loop structures below and all phenotype values come from a seeded RNG, so
nothing here derives from a real family, cohort or participant. Gene symbols
used as annotation values are real, publicly documented symbols attached to
invented individuals; they say nothing about anyone.

The script is deterministic: the same seed always produces the same files, so
regenerating is a no-op unless the structure below is edited.

    uv run python tests/data/generate_fixtures.py

Structures are chosen to exercise specific renderer features, not to resemble
any particular study. See `tests/data/README.md`.
"""

from __future__ import annotations

import random
from pathlib import Path

DATA_DIR = Path(__file__).parent
SEED = 20260903

# Real gene symbols used as illustrative annotation labels on invented people.
GENE_LABELS = ["CHD8", "ARID1B", "ADNP", "SYNGAP1"]


# Column sets, one per fixture. Kept byte-identical to what the fixtures had
# before, so parsing behaviour and feature coverage are unchanged.
HEADER_EXTENDED = [
    "FID",
    "IID",
    "PAT",
    "MAT",
    "SEX",
    "PHENO",
    "DOB",
    "TWINS",
    "DECEASED",
    "PROBAND",
    "epilepsy",
    "IQ",
    "SRS",
    "PGS",
]
HEADER_TWO_BINARY_PGS = [
    "FID",
    "IID",
    "PAT",
    "MAT",
    "SEX",
    "PROBAND",
    "b_ASD",
    "b_ADHD",
    "IQ",
    "SRS",
    "PGS",
]
HEADER_ALL_MODES = [
    "FID",
    "IID",
    "PAT",
    "MAT",
    "SEX",
    "PROBAND",
    "b_Autism",
    "IQ",
    "g_SRS",
    "PGS_autism",
    "PGS_int",
    "PGS_ADHD",
    "a_Variant1",
    "a_Variant2",
]
HEADER_ANNOTATIONS = [
    "FID",
    "IID",
    "PAT",
    "MAT",
    "SEX",
    "PROBAND",
    "b_TSA",
    "IQ",
    "g_SRS",
    "PGS_autism",
    "PGS_IQ",
    "PGS_ADHD",
    "a_Variant1",
    "a_Variant2",
]


def _write(path: Path, header: list[str], rows: list[list[str]]) -> None:
    """Write a tab-separated pedigree file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(header)]
    lines += ["\t".join(str(f) for f in row) for row in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print(f"  wrote {path.relative_to(DATA_DIR.parent.parent)} ({len(rows)} individuals)")


def build_large_multi_union() -> tuple[list[str], list[list[str]]]:
    """A founder couple with thirteen children, most of whom marry out.

    Exercises: `b_` and `o_` prefixes, structural Deceased/Proband columns,
    `-9` no-data values, a crowded middle generation (over twenty
    individuals on one ring, which forces the circular layout to expand its
    base radius), and more than ten unions.

    The founders are IDs 1 and 2; everyone else is numbered in the order
    drawn. The shape is random: which of the thirteen children stay single
    (one or two), and how many children each couple has (1 to 4).
    """
    rng = random.Random(SEED + 202)
    header = [
        "FID",
        "ID",
        "Father",
        "Mother",
        "Sex",
        "b_ASD",
        "b_Dyslexia",
        "Deceased",
        "Proband",
        "o_Sleep disorder",
    ]
    fid = "FAM3"
    rows: list[list[str]] = []
    next_id = 1

    def add(pat: str, mat: str, sex: int) -> str:
        """Append one individual with fabricated phenotype values; return its ID."""
        nonlocal next_id
        sid = str(next_id)
        next_id += 1
        if rng.random() < 0.15:
            asd, syn = "-9", "-9"
        else:
            asd, syn = str(rng.choice([1, 1, 1, 2])), str(rng.choice([1, 1, 2]))
        rows.append([fid, sid, pat, mat, sex, asd, syn, 0, 0, rng.choice([1, 1, 2])])
        return sid

    father = add("0", "0", 1)
    mother = add("0", "0", 2)
    n_children = 13
    # A .ped file links partners only through their children, so a
    # childless couple cannot be written: an unmarried child stays single.
    single = set(rng.sample(range(n_children), k=rng.choice([1, 2])))
    grandchildren: list[str] = []
    for i in range(n_children):
        child_sex = rng.choice([1, 2])
        child = add(father, mother, child_sex)
        if i in single:
            continue
        spouse = add("0", "0", 2 if child_sex == 1 else 1)
        pat, mat = (child, spouse) if child_sex == 1 else (spouse, child)
        for _ in range(rng.choice([1, 1, 2, 2, 3, 4])):
            grandchildren.append(add(pat, mat, rng.choice([1, 2])))

    by_id = {r[1]: r for r in rows}
    # One proband among the grandchildren, affected; a few deaths in the
    # middle generation.
    proband = rng.choice(grandchildren)
    by_id[proband][5] = "2"
    by_id[proband][8] = 1
    middle = [r[1] for r in rows if r[2] == father]
    for sid in rng.sample(middle, k=2):
        by_id[sid][7] = 1
    return header, rows


def build_multigenerational() -> tuple[list[str], list[list[str]]]:
    """Five generations descending from one founder couple, with two loops.

    The gallery's large example: about seventy individuals, a first-cousin
    marriage and a marriage between two branches (second cousins), dizygotic
    and monozygotic twins, three phenotype bands (`b_`, `g_`, `cr_`),
    deceased and proband markers, a miscarriage, an adoption and annotations.
    Identifiers are cohort-style, generation then rank in the file (203 is
    the third individual of generation 2). They are not the clinical II-3
    numbering, which runs left to right across the drawing: the layout
    reorders sibships so the branches of a loop meet, so file order and
    drawing order differ.
    """
    rng = random.Random(SEED + 101)
    people: list[dict] = []

    def new(gen: int, sex: int, father: dict | None = None, mother: dict | None = None) -> dict:
        person = {"gen": gen, "sex": sex, "father": father, "mother": mother}
        person.update(deceased=0, proband=0, twins=[], pregnancy="0", adopted="0", note="0")
        people.append(person)
        return person

    def children(father: dict, mother: dict, sexes: list[int]) -> list[dict]:
        return [new(father["gen"] + 1, s, father, mother) for s in sexes]

    def marry_out(person: dict) -> tuple[dict, dict]:
        # The spouse comes right after the person in the file, so numbering
        # in file order runs roughly left to right across the drawing.
        spouse = new(person["gen"], 3 - person["sex"])
        people.pop()  # dicts compare by value: move it by position, not with remove()
        at = next(i for i, q in enumerate(people) if q is person)
        people.insert(at + 1, spouse)
        return (person, spouse) if person["sex"] == 1 else (spouse, person)

    def twins(a: dict, b: dict, kind: str) -> None:
        a["twins"].append((kind, b))
        b["twins"].append((kind, a))

    # I: the founders.
    f1, f2 = new(1, 1), new(1, 2)
    # II: four children, each married to someone from outside the family.
    gen2 = children(f1, f2, [1, 2, 1, 2])
    couples2 = [marry_out(p) for p in gen2]
    # III: 4, 3, 4 and 3 children; the third sibship ends with DZ twins.
    gen3 = [
        children(pat, mat, sexes)
        for (pat, mat), sexes in zip(
            couples2, [[1, 2, 2, 1], [1, 2, 1], [2, 1, 1, 2], [2, 1, 2]], strict=True
        )
    ]
    twins(gen3[2][2], gen3[2][3], "DZ")
    # The first two of each sibship marry out, the others more often than
    # not, and each couple has two to four children (IV). The twins stay
    # single: a spouse drawn beside one twin would sit inside their connector.
    couples3: list[tuple[dict, dict]] = []
    for sibship in gen3:
        for i, person in enumerate(sibship):
            if person["twins"]:
                continue
            if i < 2 or rng.random() < 0.6:
                couples3.append(marry_out(person))
    gen4 = {
        id(pat): children(pat, mat, [rng.choice([1, 2]) for _ in range(rng.randint(2, 4))])
        for pat, mat in couples3
    }
    sibships4 = list(gen4.values())
    twins(sibships4[2][0], sibships4[2][1], "MZ")
    sibships4[2][1]["sex"] = sibships4[2][0]["sex"]

    def descends_from(person: dict, ancestor: dict) -> bool:
        stack = [person]
        while stack:
            p = stack.pop()
            if p is ancestor:
                return True
            stack += [q for q in (p["father"], p["mother"]) if q is not None]
        return False

    def child_of(branch: dict, sex: int, skip: set[int]) -> dict:
        """A generation-IV descendant of `branch` (a generation-II person)."""
        for sibship in sibships4:
            for p in sibship:
                if id(p) not in skip and descends_from(p, branch) and not p["twins"]:
                    p["sex"] = sex
                    skip.add(id(p))
                    return p
        raise RuntimeError("no candidate")

    used: set[int] = set()
    # Loop 1, first cousins: two grandchildren of the first generation-II couple
    # through different children.
    a = child_of(gen3[0][0], 1, used)
    b = child_of(gen3[0][1], 2, used)
    # Loop 2, second cousins: descendants of the second and fourth branches.
    c = child_of(gen2[1], 1, used)
    d = child_of(gen2[3], 2, used)
    loop1, loop2 = (a, b), (c, d)
    # Two more couples in generation IV marry out.
    out4 = [marry_out(child_of(gen2[2], 1, used)), marry_out(child_of(gen2[0], 2, used))]
    # V: their children.
    gen5_loop1 = children(*loop1, [2, 1, 0, 1])
    gen5_loop2 = children(*loop2, [1, 2])
    for pat, mat in out4:
        children(pat, mat, [rng.choice([1, 2]) for _ in range(2)])

    # Markers and notation.
    for p in (f1, f2, gen2[0], couples2[2][1], gen3[3][2]):
        p["deceased"] = 1
        p["note"] = f"d. {rng.randint(55, 92)}y"
    gen5_loop1[0]["proband"] = 1
    gen5_loop1[2]["pregnancy"] = "SAB"
    gen5_loop1[2]["note"] = "10 wk"
    adoptee = sibships4[-1][-1]
    adoptee["adopted"] = "in"

    # Phenotypes: a recessive-looking condition, affected mostly below loops.
    loop_children = {id(p) for p in gen5_loop1 + gen5_loop2}
    for p in people:
        p["asd"] = (
            "2"
            if (id(p) in loop_children and rng.random() < 0.6)
            else ("2" if rng.random() < 0.08 else "1")
        )
        if p["gen"] <= 2 and rng.random() < 0.4:
            p["asd"] = "-9"
        affected = p["asd"] == "2"
        p["iq"] = (
            "-9"
            if p["gen"] <= 2
            else str(rng.randint(55, 85) if affected else rng.randint(85, 130))
        )
        p["carrier"] = "0"
    gen5_loop1[0]["asd"] = "2"
    gen5_loop1[0]["iq"] = "68"
    gen5_loop1[2]["asd"] = gen5_loop1[2]["iq"] = "-9"  # the miscarriage
    # Recessive story: the partners of each loop and one of the proband's
    # brothers are unaffected carriers.
    for p in (*loop1, *loop2, gen5_loop1[1]):
        if p["asd"] == "2":
            p["asd"], p["iq"] = "1", str(rng.randint(90, 125))
    for p in (*loop1, *loop2):
        p["carrier"] = "2"  # obligate: parents of an affected child
    for p in gen5_loop1 + gen5_loop2:
        if p["asd"] != "2" and p["pregnancy"] == "0":
            p["carrier"] = rng.choice(["0", "1"])
    adoptee["carrier"] = "0"
    gen5_loop1[1]["carrier"] = "1"

    # Number generation by generation, in file order.
    counters = [0] * 6
    for p in people:
        counters[p["gen"]] += 1
        p["id"] = f"{p['gen']}{counters[p['gen']]:02d}"

    header = [
        "FID",
        "IID",
        "Father",
        "Mother",
        "Sex",
        "Deceased",
        "Proband",
        "Twins",
        "b_ASD",
        "g_IQ",
        "cr_HEXA",
        "Pregnancy",
        "Adopted",
        "a_Note",
    ]
    rows = []
    for p in sorted(people, key=lambda q: q["gen"]):
        rows.append(
            [
                "FAM1",
                p["id"],
                p["father"]["id"] if p["father"] else "0",
                p["mother"]["id"] if p["mother"] else "0",
                p["sex"],
                p["deceased"],
                p["proband"],
                ",".join(f"{kind}:{q['id']}" for kind, q in p["twins"]) or "0",
                p["asd"],
                p["iq"],
                p["carrier"],
                p["pregnancy"],
                p["adopted"],
                p["note"],
            ]
        )
    return header, rows


def build_trio_with_siblings(
    fid: str,
    header: list[str],
    n_individuals: int,
    value_fn,
) -> tuple[list[str], list[list[str]]]:
    """Two parents and their children, `n_individuals` in total.

    `value_fn(rng, sid, is_child, proband=0)` returns the trailing phenotype
    fields, so each fixture decides its own column set.
    """
    rng = random.Random(SEED + len(fid))
    rows: list[list[str]] = []
    pat, mat = "002", "003"

    # Proband first, so ID 001 carries the arrow.
    rows.append([fid, "001", pat, mat, 1, *value_fn(rng, "001", True, proband=1)])
    rows.append([fid, pat, "0", "0", 1, *value_fn(rng, pat, False)])
    rows.append([fid, mat, "0", "0", 2, *value_fn(rng, mat, False)])
    for i in range(n_individuals - 3):
        sid = f"{i + 4:03d}"
        rows.append([fid, sid, pat, mat, 1 if i % 2 else 2, *value_fn(rng, sid, True)])
    return header, rows


def _maybe_missing(rng: random.Random, value: object, chance: float = 0.18) -> str:
    """Return the value, or the -9 no-data sentinel."""
    return "-9" if rng.random() < chance else str(value)


def _comma_decimal(value: float) -> str:
    """Format a score with a comma decimal separator (European locale)."""
    return f"{value:.10f}".replace(".", ",")


def build_quadrants() -> tuple[list[str], list[list[str]]]:
    """Three generations with the four conventional phenotypes.

    The gallery's example of the quadrant layout: autism, ADHD,
    intellectual disability and epilepsy as binary columns (each drawn in
    its conventional quadrant and color), plus an `o_` column for the
    individuals with sequencing data. The shape and every value are drawn
    at random from the seed.
    """
    rng = random.Random(SEED + 303)
    header = [
        "FID",
        "IID",
        "Father",
        "Mother",
        "Sex",
        "Deceased",
        "Proband",
        "b_Autism",
        "b_ADHD",
        "b_ID",
        "b_Epilepsy",
        "o_WGS",
    ]
    rows: list[list[str]] = []
    next_id = 1

    def add(pat: str, mat: str, sex: int, gen: int) -> str:
        nonlocal next_id
        sid = f"{gen}{next_id:02d}"
        next_id += 1
        if rng.random() < 0.08:
            values = ["-9"] * 4
        else:
            values = [
                str(rng.choice(choices))
                for choices in ([1, 1, 2], [1, 1, 2], [1, 1, 1, 2], [1, 1, 1, 1, 2])
            ]
        rows.append(["FAM1", sid, pat, mat, sex, 0, 0, *values, rng.choice([1, 2, 2])])
        return sid

    grandfather = add("0", "0", 1, 1)
    grandmother = add("0", "0", 2, 1)
    couples = []
    for _ in range(rng.choice([2, 3])):
        sex = rng.choice([1, 2])
        child = add(grandfather, grandmother, sex, 2)
        spouse = add("0", "0", 2 if sex == 1 else 1, 2)
        couples.append((child, spouse) if sex == 1 else (spouse, child))
    add(grandfather, grandmother, rng.choice([1, 2]), 2)  # stays single
    grandchildren = []
    for pat, mat in couples:
        for _ in range(rng.choice([2, 3, 4])):
            grandchildren.append(add(pat, mat, rng.choice([1, 2]), 3))

    by_id = {r[1]: r for r in rows}
    proband = rng.choice(grandchildren)
    by_id[proband][6] = 1
    by_id[proband][7] = "2"  # the proband has autism
    by_id[grandfather][5] = 1
    return header, rows


def main() -> None:
    print(f"Generating synthetic fixtures (seed {SEED})")

    header, rows = build_large_multi_union()
    _write(DATA_DIR / "complex/large_multi_union.ped", header, rows)

    # Extended 10-column format plus four phenotype columns, one a PGS.
    def extended_values(rng, sid, is_child, proband=0):
        return [
            rng.choice([1, 2]),  # PHENO
            0,  # DOB
            0,  # TWINS
            0,  # DECEASED
            proband,
            rng.choice([0, 1]),  # epilepsy
            rng.randint(75, 125),  # IQ
            rng.randint(40, 90),  # SRS
            f"{rng.uniform(-1.5, 1.5):.2f}",  # PGS
        ]

    header, rows = build_trio_with_siblings(
        "FAM1",
        HEADER_EXTENDED,
        7,
        extended_values,
    )
    _write(DATA_DIR / "complex/large_family.ped", header, rows)

    # Two binary columns and two continuous ones, with no-data values.
    def pgs_values(rng, sid, is_child, proband=0):
        return [
            proband,
            _maybe_missing(rng, rng.choice([1, 2])),  # b_ASD
            _maybe_missing(rng, rng.choice([1, 2])),  # b_ADHD
            _maybe_missing(rng, rng.randint(70, 120)),  # IQ
            _maybe_missing(rng, rng.randint(45, 155)),  # SRS
            f"{rng.uniform(-1.2, 1.2):.3f}",  # PGS
        ]

    header, rows = build_trio_with_siblings("FAM1", HEADER_TWO_BINARY_PGS, 7, pgs_values)
    _write(DATA_DIR / "complex/multi_phenotype_pgs.ped", header, rows)

    # Every display mode at once: binary, AUTO-gradient, explicit gradient,
    # three PGS columns and two annotation columns. Scores use comma decimals.
    def all_modes_values(gene_col_name):
        def _fn(rng, sid, is_child, proband=0):
            genes = rng.sample(GENE_LABELS, k=rng.choice([0, 1, 1, 2]))
            return [
                proband,
                rng.choice([1, 2]),  # binary
                _maybe_missing(rng, rng.randint(70, 120), 0.12),  # IQ (AUTO -> gradient)
                _maybe_missing(rng, rng.randint(40, 95), 0.12),  # g_SRS
                _comma_decimal(rng.uniform(-1, 1)),
                _comma_decimal(rng.uniform(-1, 1)),
                _comma_decimal(rng.uniform(-1, 1)),
                genes[0] if genes else "",
                genes[1] if len(genes) > 1 else "",
            ]

        del gene_col_name
        return _fn

    header, rows = build_trio_with_siblings(
        "FAM1",
        HEADER_ALL_MODES,
        6,
        all_modes_values("PGS_int"),
    )
    _write(DATA_DIR / "phenotype/test_all.ped", header, rows)

    header, rows = build_trio_with_siblings(
        "FAM1",
        HEADER_ANNOTATIONS,
        7,
        all_modes_values("PGS_IQ"),
    )
    _write(DATA_DIR / "phenotype/annotations.ped", header, rows)

    header, rows = build_multigenerational()
    _write(DATA_DIR / "complex/multigenerational.ped", header, rows)

    header, rows = build_quadrants()
    _write(DATA_DIR / "phenotype/quadrants.ped", header, rows)


if __name__ == "__main__":
    main()
