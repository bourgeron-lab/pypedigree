"""The row pass that lowers founders to their children (#163), kept fast.

`layout._lower_to_children` runs once per loop-breaking pass, so on a large
pedigree it ran 22 times and took about 45% of the layout (20,000 people,
1.0.0rc3). The optimised version (groups looked up once, an early return
when no line can get shorter, moves tried only from groups that can reach
such a line) must give exactly the rows of the plain algorithm below, which
is the 1.0.0rc3 one. All pedigrees here are synthetic.
"""

from __future__ import annotations

import sys
from collections import deque
from pathlib import Path

import pytest

from pypedigree import PedigreeStyle, parse_pedigree
from pypedigree.layout import (
    _assign_generations,
    _HierarchicalPlacer,
    _lower_to_children,
    _with_duplicate_partners,
)
from pypedigree.models import Family

sys.path.insert(0, str(Path(__file__).parent))
from conftest import (  # noqa: E402
    TEST_DATA_DIR,
    generate_clinical_pedigree,
    generate_realistic_pedigree,
)


def _reference(family: Family, generations: dict[str, int]) -> None:
    """The plain algorithm: every group tried as a start, every pass."""
    group = {sid: sid for sid in generations}

    def find(sid: str) -> str:
        while group[sid] != sid:
            group[sid] = group[group[sid]]
            sid = group[sid]
        return sid

    for union in family.unions:
        p1, p2 = union.partner_ids
        if p1 in group and p2 in group:
            group[find(p1)] = find(p2)
    members: dict[str, list[str]] = {}
    for sid in generations:
        members.setdefault(find(sid), []).append(sid)
    links = []
    for ind in family.individuals.values():
        sid = ind.sample_id
        if sid not in generations:
            continue
        for parent in (ind.father_id, ind.mother_id):
            if parent and parent in generations:
                links.append((parent, sid, 1))
        if ind.is_duplicate and ind.duplicate_of in generations:
            links.append((ind.duplicate_of, sid, 0))
    down: dict[str, list] = {}
    up: dict[str, list] = {}
    for link in links:
        down.setdefault(find(link[0]), []).append(link)
        up.setdefault(find(link[1]), []).append(link)

    def slack(link) -> int:
        return generations[link[1]] - generations[link[0]] - link[2]

    order = sorted(members, key=lambda g: generations[members[g][0]])
    for _ in range(len(links) + 1):
        moved = False
        for start in order:
            moving = {start}
            queue = deque([start])
            while queue:
                for link in down.get(queue.popleft(), []):
                    lower = find(link[1])
                    if lower not in moving and slack(link) == 0:
                        moving.add(lower)
                        queue.append(lower)
            shorter = [ln for g in moving for ln in down.get(g, []) if find(ln[1]) not in moving]
            longer = [ln for g in moving for ln in up.get(g, []) if find(ln[0]) not in moving]
            if len(shorter) <= len(longer):
                continue
            step = min(slack(ln) for ln in shorter)
            for g in moving:
                for sid in members[g]:
                    generations[sid] += step
            moved = True
        if not moved:
            break


def _families() -> list[Family]:
    families: list[Family] = []
    for path in sorted(TEST_DATA_DIR.rglob("*.ped")):
        try:
            families += list(parse_pedigree(path).families.values())
        except Exception:  # malformed-input fixtures
            continue
    families += [generate_realistic_pedigree(n, s) for s in range(4) for n in (40, 150)]
    families += [generate_clinical_pedigree(150, s, cousin_marriages=3) for s in range(4)]
    return families


FAMILIES = _families()


@pytest.mark.parametrize("family", FAMILIES, ids=lambda f: f"{f.family_id}-{len(f.individuals)}")
def test_same_rows_as_the_plain_algorithm(family):
    """At every loop-breaking pass, as `_place_hierarchically` sees the family."""
    style = PedigreeStyle()
    for _ in range(40):
        rows = _assign_generations(family)
        fast, plain = dict(rows), dict(rows)
        _lower_to_children(family, fast)  # idempotent on its own result
        assert fast == rows
        # The plain algorithm from the rows before lowering gives the same.
        _reference_family_rows(family, plain)
        assert plain == rows
        placer = _HierarchicalPlacer(family, rows, style, True)
        to_break = placer.loops_to_break(placer.run(), False)
        if not to_break:
            break
        family = _with_duplicate_partners(family, to_break)


def _reference_family_rows(family: Family, rows: dict[str, int]) -> None:
    """Recompute ``rows`` as `_assign_generations` does, lowering with the
    plain algorithm, in place."""
    import pypedigree.layout as layout

    real = layout._lower_to_children
    layout._lower_to_children = _reference
    try:
        rows.clear()
        rows.update(layout._assign_generations(family))
    finally:
        layout._lower_to_children = real


def test_the_reported_family_still_lowered():
    """The early return leaves alone only what cannot move (#163)."""
    family = parse_pedigree(TEST_DATA_DIR / "complex/married_in_founders.ped").get_family("FAM1")
    rows = _assign_generations(family)
    assert rows["21"] == rows["22"] == 1
