"""pTDT deviations computed from the parents' PGS (#52, second step).

The polygenic transmission disequilibrium test (Weiner et al. 2017, Nat Genet)
asks whether children inherited more of a polygenic score than expected from
their parents. Each child's deviation is

    (PGS_child - PGS_mid-parent) / SD(PGS_mid-parent)

N(0, 1) without over-transmission. With ``pgs_transmission`` (``--ptdt``),
pypedigree computes it for every child whose two parents have a value of a
``PGS_*`` column, into a ``pTDT_*`` column drawn on the z scale (#52, first
step). SD(mid-parent) is a cohort quantity one family cannot estimate: it is
``ptdt_midparent_sd``, √½ by default, the value for standardised scores and
unrelated parents. A child missing a parent, or a parent's value, gets no
deviation. Parents keep their PGS curve; a child with a deviation shows it
in a pTDT box beside its PGS box (`phenotypes.active_curves`, R6).

The deviations are computed over the whole file (`Pedigree.column_scope`) and
the pedigree is copied: the parsed one is not modified. No matplotlib (rule 1).
"""

from __future__ import annotations

import math
import warnings
from dataclasses import replace

from pypedigree.config import PedigreeStyle, is_pgs_column, is_ptdt_column
from pypedigree.errors import PedigreeDataWarning, PedigreeError
from pypedigree.models import Family, Individual, Pedigree


def ptdt_name(pgs_name: str) -> str:
    """The pTDT column computed from a PGS column: "PGS_autism" → "pTDT_autism"."""
    return "pTDT" + pgs_name[3:]


def check_transmission(style: PedigreeStyle) -> None:
    """Reject a ``ptdt_midparent_sd`` that is not a positive number.

    Raises:
        PedigreeError: If it is not positive and finite.
    """
    sd = style.ptdt_midparent_sd
    if not isinstance(sd, (int, float)) or isinstance(sd, bool) or not (0 < sd < math.inf):
        raise PedigreeError(f"ptdt_midparent_sd must be a positive number, got {sd!r}")


def _value(ind: Individual | None, name: str) -> float | None:
    if ind is None or ind.is_phantom:
        return None
    raw = ind.extra_phenotypes_raw.get(name)
    if raw is None or (isinstance(raw, float) and math.isnan(raw)):
        return None
    return float(raw)


def _deviations(family: Family, targets: dict[str, str], sd: float) -> Family:
    people = family.individuals

    def parent(pid: str | None) -> Individual | None:
        ind = people.get(pid) if pid else None
        if ind is not None and ind.is_duplicate and ind.duplicate_of:
            ind = people.get(ind.duplicate_of)
        return ind

    computed: dict[str, dict[str, float | None]] = {}
    for sid, ind in people.items():
        if ind.is_duplicate:
            continue
        father, mother = parent(ind.father_id), parent(ind.mother_id)
        values: dict[str, float | None] = {}
        for pgs, ptdt in targets.items():
            child, f, m = _value(ind, pgs), _value(father, pgs), _value(mother, pgs)
            values[ptdt] = None if None in (child, f, m) else (child - (f + m) / 2) / sd  # type: ignore[operator]
        computed[sid] = values
    new_people = {}
    for sid, ind in people.items():
        # A person drawn twice keeps the original's values.
        values = computed.get(ind.duplicate_of or "", {}) if ind.is_duplicate else computed[sid]
        new_people[sid] = replace(
            ind,
            extra_phenotypes_raw={**ind.extra_phenotypes_raw, **values},
            extra_phenotypes={**ind.extra_phenotypes, **values},
        )
    return replace(family, individuals=new_people)


def with_transmission(pedigree: Pedigree, style: PedigreeStyle) -> tuple[Pedigree, frozenset[str]]:
    """The pedigree with its children's pTDT deviations, and the columns computed.

    With ``style.pgs_transmission`` off, or no PGS column, the pedigree itself
    and no column. A ``pTDT_*`` column already in the file is kept, with a
    warning, and not computed.
    """
    if not style.pgs_transmission:
        return pedigree, frozenset()
    check_transmission(style)
    scope = pedigree.column_scope
    names = list(scope.phenotype_names)
    given = {n for n in names if is_ptdt_column(n)}
    targets = {}
    kept = []
    for name in names:
        if is_pgs_column(name):
            ptdt = ptdt_name(name)
            if ptdt in given:
                kept.append(ptdt)
            else:
                targets[name] = ptdt
    if kept:
        warnings.warn(
            f"pgs_transmission: {', '.join(kept)} is given in the file and is drawn as given, "
            "not computed from the parents.",
            PedigreeDataWarning,
            stacklevel=3,
        )
    if not targets:
        return pedigree, frozenset()
    order = []
    for name in names:
        order.append(name)
        if name in targets:
            order.append(targets[name])
    whole = replace(
        scope,
        families={
            fid: _deviations(f, targets, style.ptdt_midparent_sd)
            for fid, f in scope.families.items()
        },
        phenotype_names=order,
        source=None,
    )
    out = whole if pedigree.source is None else whole.subset(pedigree.families)
    return out, frozenset(targets.values())
