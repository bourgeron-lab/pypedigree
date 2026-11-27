"""Layout engine for pedigree diagrams.

Computes (x, y) coordinates for all individuals and connection lines.
No matplotlib dependency — produces pure numeric FamilyLayout objects.

The layout algorithm:
1. Assign generations (y-axis) via BFS from founders
2. Compute horizontal positions (x-axis): each union's descendants form a
   block, sibling blocks are packed against each other's contours, and
   every couple is centred over its children
3. Compute all connection geometry (relationship, descent, sibship, individual lines)
"""

from __future__ import annotations

import math
from collections import deque
from collections.abc import Mapping
from dataclasses import replace

from pypedigree import symbols
from pypedigree.config import PedigreeStyle, check_text_scale
from pypedigree.errors import PedigreeError
from pypedigree.models import (
    Adoption,
    DuplicateArc,
    EdgeLayout,
    EdgeType,
    Family,
    FamilyLayout,
    Individual,
    LineStyle,
    NodeLayout,
    Sex,
    TwinConnector,
    TwinType,
    Union,
)


def compute_layout(
    family: Family, style: PedigreeStyle | None = None, *, labels: bool = True
) -> FamilyLayout:
    """Compute the complete layout for a family.

    Args:
        family: The family to lay out.
        style: Visual style parameters. Uses defaults if None.
        labels: Whether labels will be drawn. Neighbours are then kept far
            enough apart that their labels (ID, notes and annotations) do
            not run into each other.

    Returns:
        FamilyLayout with all positions and connections computed.
    """
    if style is None:
        style = PedigreeStyle()
    check_labels(style)
    check_text_scale(style)

    if not family.individuals:
        return FamilyLayout(family_id=family.family_id)
    family = with_unknown_offspring(family)

    # Steps 1-2: Assign generations and horizontal positions. Loops whose
    # drawing would make sibship lines overlap are broken by drawing one
    # partner twice, so the family laid out may hold extra duplicate nodes.
    family, generations, x_positions = _place_hierarchically(family, style, labels)

    # Step 3: Increase generation gap for what sits below the symbols:
    # notes, annotations and childlessness markers (and, with labels
    # "both", the ID under the number).
    below = _room_below_labels(family, style, 1 if style.labels == "both" else 0)
    if below > 0:
        style = replace(style, generation_gap=style.generation_gap + below)

    # Step 4: Build node layouts
    nodes = _build_nodes(family, x_positions, generations, style)
    nodes = _numbered(nodes, family, generations, x_positions, style)

    # Step 5: Compute connection edges (twin groups excluded from individual lines)
    twin_groups = _find_twin_groups(family)
    edges = _compute_edges(family, x_positions, generations, style, twin_groups)

    # Step 6: Compute twin connectors
    twin_connectors = _compute_twin_connectors(family, x_positions, generations, style, twin_groups)

    edges += _childless_edges(nodes, style)
    edges += _couple_childless_edges(family, nodes, style)

    # Step 7: Compute duplicate arcs
    positions = {
        sid: (x, -generations.get(sid, 0) * style.generation_gap) for sid, x in x_positions.items()
    }
    duplicate_arcs = _route_duplicate_arcs(_compute_duplicate_arcs(family, positions, style), style)

    return FamilyLayout(
        family_id=family.family_id,
        nodes=nodes,
        edges=edges,
        twin_connectors=twin_connectors,
        duplicate_arcs=duplicate_arcs,
    )


def _assign_generations(family: Family) -> dict[str, int]:
    """Assign generation numbers to each individual.

    Founders (no parents) are generation 0. Children are placed one
    generation below their parents. Both partners in a union are
    aligned to the same generation.

    Args:
        family: The family to process.

    Returns:
        Map of sample_id to generation number (0 = top).
    """
    generations: dict[str, int] = {}

    # Build parent→children adjacency index for O(1) lookups
    parent_to_children: dict[str, list[Individual]] = {}
    for ind in family.individuals.values():
        if ind.father_id:
            parent_to_children.setdefault(ind.father_id, []).append(ind)
        if ind.mother_id:
            parent_to_children.setdefault(ind.mother_id, []).append(ind)

    # Start with founders at generation 0
    founders = family.get_founders()
    queue: deque[str] = deque()

    for founder in founders:
        generations[founder.sample_id] = 0
        queue.append(founder.sample_id)

    # BFS: assign children one generation below parents
    while queue:
        current_id = queue.popleft()
        current_gen = generations[current_id]

        for ind in parent_to_children.get(current_id, []):
            if ind.sample_id in generations:
                continue
            child_gen = current_gen + 1
            # If the other parent is already assigned, use max + 1
            other_parent = ind.mother_id if ind.father_id == current_id else ind.father_id
            if other_parent and other_parent in generations:
                parent_gen = max(current_gen, generations[other_parent])
                child_gen = parent_gen + 1

            generations[ind.sample_id] = child_gen
            queue.append(ind.sample_id)

    # Repeat the passes below until partners share a row, duplicates stand
    # no higher than the individual they repeat, and children sit below
    # their parents. The bound only matters for impossible constraints (a
    # partner who is also a descendant).
    for _ in range(len(family.individuals) + 1):
        before = dict(generations)

        # A duplicate node (drawn for a second partnership) has no parents
        # of its own: it stands at least on the row of the individual it
        # repeats, lower if their partner is from a younger generation.
        for ind in family.individuals.values():
            original = ind.duplicate_of
            if ind.is_duplicate and original in generations and ind.sample_id in generations:
                generations[ind.sample_id] = max(generations[ind.sample_id], generations[original])

        # Align partners: both should be at the same generation (use max
        # so children stay below both parents)
        for union in family.unions:
            p1, p2 = union.partner_ids
            if p1 in generations and p2 in generations:
                aligned_gen = max(generations[p1], generations[p2])
                generations[p1] = aligned_gen
                generations[p2] = aligned_gen

        # Re-propagate: ensure children are always below both parents
        changed = True
        while changed:
            changed = False
            for ind in family.individuals.values():
                if ind.sample_id not in generations:
                    continue
                parent_gen = -1
                if ind.father_id and ind.father_id in generations:
                    parent_gen = max(parent_gen, generations[ind.father_id])
                if ind.mother_id and ind.mother_id in generations:
                    parent_gen = max(parent_gen, generations[ind.mother_id])
                if parent_gen >= 0 and generations[ind.sample_id] <= parent_gen:
                    generations[ind.sample_id] = parent_gen + 1
                    changed = True

        if generations == before:
            break

    # Handle any unassigned individuals (shouldn't happen normally)
    max_gen = max(generations.values()) if generations else 0
    for sid in family.individuals:
        if sid not in generations:
            generations[sid] = max_gen + 1

    _lower_to_children(family, generations)
    return generations


def _lower_to_children(family: Family, generations: dict[str, int]) -> None:
    """Move people down when that shortens the lines joining them to the rest.

    Every founder starts on the top row, so the parents of someone who
    marries into a younger generation stay far above them, and the sibship
    line of their parents is drawn a row or more above them. This lowers
    such a set of people (founders, their partners, and the children and
    duplicates that must move with them to stay below their parents) as
    long as more lines get shorter than longer, in place.

    The constraints kept are the ones the passes above settle: partners
    share a row, children are below their parents, duplicates no higher
    than the individual they repeat. The total length of the lines only
    decreases, so this ends.

    Args:
        family: The family laid out.
        generations: Row of each individual, modified in place.
    """
    # Partners are drawn level: move each group of partners as one.
    group: dict[str, str] = {sid: sid for sid in generations}

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

    # Lines from an upper individual to a lower one, at least `gap` rows
    # apart: parent to child (1), original to duplicate (0).
    links: list[tuple[str, str, int]] = []
    for ind in family.individuals.values():
        sid = ind.sample_id
        if sid not in generations:
            continue
        for parent in (ind.father_id, ind.mother_id):
            if parent and parent in generations:
                links.append((parent, sid, 1))
        if ind.is_duplicate and ind.duplicate_of in generations:
            links.append((ind.duplicate_of, sid, 0))
    down: dict[str, list[tuple[str, str, int]]] = {}
    up: dict[str, list[tuple[str, str, int]]] = {}
    for link in links:
        down.setdefault(find(link[0]), []).append(link)
        up.setdefault(find(link[1]), []).append(link)

    def slack(link: tuple[str, str, int]) -> int:
        upper, lower, gap = link
        return generations[lower] - generations[upper] - gap

    order = sorted(members, key=lambda g: generations[members[g][0]])
    for _ in range(len(links) + 1):
        moved = False
        for start in order:
            # The groups that must move with this one: those hanging from
            # it by a line already as short as it can be.
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


class _Block:
    """A laid-out piece of the pedigree in its own coordinate frame.

    Keeps the position of every member and, for each generation, the
    leftmost and rightmost occupied x: the block's contour. Blocks are
    composed side by side by comparing contours generation by generation,
    so a subtree is packed as close to its neighbour as every one of its
    rows allows, and never interleaves with it.
    """

    __slots__ = ("hi", "lo", "xs")

    def __init__(self) -> None:
        self.xs: dict[str, float] = {}
        self.lo: dict[int, float] = {}
        self.hi: dict[int, float] = {}

    def add(self, sid: str, x: float, gen: int, left: float = 0.0, right: float = 0.0) -> None:
        """Put ``sid`` at ``x``, keeping extra room on their left and right."""
        self.xs[sid] = x
        self.lo[gen] = min(self.lo.get(gen, x - left), x - left)
        self.hi[gen] = max(self.hi.get(gen, x + right), x + right)

    def place(self, other: _Block, gap: float, target: float | None = None) -> float:
        """Merge ``other`` to the right of this block and return its offset.

        On every generation both blocks occupy, ``other`` lands at least
        ``gap`` to the right of this block. ``target`` is the preferred
        offset; it is honoured when the contours allow it.
        """
        shared = self.hi.keys() & other.lo.keys()
        if shared:
            offset = max(self.hi[g] + gap - other.lo[g] for g in shared)
            if target is not None:
                offset = max(offset, target)
        elif target is not None:
            offset = target
        elif self.xs and other.xs:
            offset = max(self.hi.values()) + gap - min(other.lo.values())
        else:
            offset = 0.0
        for sid, x in other.xs.items():
            self.xs[sid] = x + offset
        for g, x in other.lo.items():
            self.lo[g] = min(self.lo.get(g, x + offset), x + offset)
        for g, x in other.hi.items():
            self.hi[g] = max(self.hi.get(g, x + offset), x + offset)
        return offset


class _HierarchicalPlacer:
    """Horizontal placement of a family as nested, non-interleaving blocks.

    Each individual with parents in the pedigree is placed inside the
    block of their parents' union, together with their partner and their
    own descendants. Sibling blocks are packed left to right against each
    other's contours and every couple is centred over its children. As a
    result, the children of one union occupy a horizontal interval that no
    other sibship at the same generation enters, and sibship lines cannot
    overlap.

    Loops (both partners descend from people in the pedigree, as in
    consanguinity) cannot be drawn as a tree. The partner placed second is
    drawn next to the first one, on the side of their own family, and the
    sibship line of their parents stretches to reach them.
    """

    def __init__(
        self,
        family: Family,
        generations: dict[str, int],
        style: PedigreeStyle,
        labels: bool = True,
    ) -> None:
        self.family = family
        self.gen = generations
        self.style = style
        self.unions = [u for u in family.unions if all(p in generations for p in u.partner_ids)]
        self.unions_of: dict[str, list[int]] = {}
        self.parent_union: dict[str, int] = {}
        for ui, union in enumerate(self.unions):
            for pid in union.partner_ids:
                self.unions_of.setdefault(pid, []).append(ui)
            for cid in union.children_ids:
                if cid in generations:
                    self.parent_union.setdefault(cid, ui)
        self.duplicates_of: dict[str, list[str]] = {}
        for sid, ind in family.individuals.items():
            if ind.is_duplicate and ind.duplicate_of in generations:
                self.duplicates_of.setdefault(ind.duplicate_of, []).append(sid)
        self.attached = self._attached_to_lineage()
        self.outer_side = self._twin_outer_sides()
        self.claimed: set[str] = set()
        self.done_unions: set[int] = set()
        # Loop partner drawn next to their partner -> that partnership.
        self.attached_by: dict[str, int] = {}
        # Union whose lines a loop stretches -> (partnership, partner) to
        # duplicate if the stretched lines overlap another sibship.
        self.loops: dict[int, list[tuple[int, str]]] = {}
        self._lines: dict[str, frozenset[str]] = {}
        self.marker_room = self._marker_room()
        self.label_reach = self._label_reach() if labels else {}

    def _marker_room(self) -> dict[str, float]:
        """Distance each proband or consultand needs to their left neighbour.

        Their marker sits at the lower left of the symbol, level with the
        labels; with tight spacing or large text it would run into the
        label on the left (`symbols.marker_clearance`, for the longest
        label of the generation), or into a left neighbour's evaluation
        asterisk (the widest-reaching one of the generation, #143). Only
        those who need more than the sibling gap are listed.
        """
        style = self.style
        half = style.shape_size / 2
        mark_size = style.text_size(symbols.EVALUATION_MARK_FONT_SCALE)
        upp = symbols.units_per_point(style.text_scale)
        longest: dict[int, int] = {}
        # The asterisks of the generation, as far as they reach right of
        # their symbol, gap included (#143): any of them may be on the left.
        marks: dict[int, list[tuple[float, str]]] = {}
        for sid, g in self.gen.items():
            ind = self.family.individuals[sid]
            label = _get_node_label(sid, ind, style)
            longest[g] = max(longest.get(g, 0), len(label))
            if ind.notation.is_evaluated:
                brackets = style.line_width if ind.notation.adoption is not None else None
                sex = symbols.symbol_sex(ind.sex, ind.notation)
                reach = symbols.evaluation_mark_reach(
                    half, sex, mark_size, style.shape_line_width, brackets, upp
                )
                reach += symbols.EVALUATION_MARK_GAP_EM * mark_size * upp
                marks.setdefault(g, []).append((reach, sid))
        room: dict[str, float] = {}
        for sid, g in self.gen.items():
            ind = self.family.individuals[sid]
            if ind.is_proband or ind.notation.is_consultand:
                mark = max((r for r, other in marks.get(g, []) if other != sid), default=0.0)
                need = symbols.marker_clearance(
                    half, style.text_size(), longest[g], ind.is_proband, mark, upp
                )
                if need > style.sibling_gap:
                    room[sid] = need
        return room

    def _label_reach(self) -> dict[str, float]:
        """How far each label reaches from its symbol's centre, on each side.

        Half the widest line drawn under the symbol (ID, notes, annotations;
        `symbols.label_width`) plus half the space kept between two labels.
        Every labelled node is listed, short labels too: `_label_room` only
        counts what reaches beyond half the sibling gap, but a couple's gap
        adds both partners' reaches, so a long label next to a short one
        needs the short one's real reach, not nothing.
        """
        style = self.style
        size, small = style.text_size(), style.text_size(symbols.SMALL_TEXT_SCALE)
        upp = symbols.units_per_point(style.text_scale)
        pad = symbols.LABEL_GAP_EM * size * upp
        reach: dict[str, float] = {}
        for sid in self.gen:
            ind = self.family.individuals[sid]
            label = _get_node_label(sid, ind, style)
            below = symbols.notes(ind) + [v for v in ind.annotations.values() if v]
            if label:
                reach[sid] = (symbols.label_width(label, below, size, small, upp) + pad) / 2
        return reach

    def _label_room(self, sid: str) -> float:
        """Room beyond half the sibling gap that the label of ``sid`` needs."""
        return max(0.0, self.label_reach.get(sid, 0.0) - self.style.sibling_gap / 2)

    def _couple_gap(self, left: str, right: str, ui: int | None = None) -> float:
        """Distance between partners: the couple gap, or more for marker or labels.

        The labels need the sum of both partners' reaches, so they do not meet
        (#41). A couple with children also has its descent line in the middle
        of the gap, so each partner's text must stay on its side of it: the
        gap is then at least twice the larger reach (#130; a CanRisk panel's
        gene line under one partner crossed the line).
        """
        reach = self.label_reach
        left_reach, right_reach = reach.get(left, 0.0), reach.get(right, 0.0)
        labels = left_reach + right_reach
        if ui is not None and self.unions[ui].children_ids:
            labels = max(labels, 2 * max(left_reach, right_reach))
        return max(self.style.couple_gap, self.marker_room.get(right, 0.0), labels)

    def _add(self, block: _Block, sid: str, x: float, gen: int, pgs_right: bool = False) -> None:
        """Put ``sid`` in ``block``, with the room their marker and label need.

        Blocks keep the sibling gap between their contours, so the room on
        each side adds to it: two neighbours with long labels end up at
        least the sum of their labels' reaches apart. A right-hand partner,
        whose PGS box is drawn on their right (#78), keeps its room there:
        the sibling gap already holds the next symbol's box, on its left. A
        box moved further out than usual (past the evaluation asterisk, or a
        deceased line's end, #80) keeps the extra on its side.
        """
        label = self._label_room(sid)
        marker = max(0.0, self.marker_room.get(sid, 0.0) - self.style.sibling_gap)
        left, right = max(marker, label), label
        if self.style.pgs_box_room > 0:
            if pgs_right:
                right = max(right, self.style.pgs_box_room + self._pgs_extra(sid, "right"))
            else:
                left += self._pgs_extra(sid, "left")
        block.add(sid, x, gen, left, right)

    def _pgs_extra(self, sid: str, side: str) -> float:
        """How much further out than usual the PGS box of ``sid`` is drawn on
        ``side``: past a diamond's tip (#111), the evaluation asterisk (#78)
        or a deceased line's end (#80)."""
        style = self.style
        ind = self.family.individuals[sid]
        half = style.shape_size / 2
        size = style.text_size(symbols.EVALUATION_MARK_FONT_SCALE)
        sex = symbols.symbol_sex(ind.sex, ind.notation)
        near = symbols.pgs_box_near(
            half,
            sex,
            ind.notation.is_evaluated,
            size,
            side,
            struck=symbols.is_struck_through(ind),
            line_width=style.shape_line_width,
            bracket_width=style.line_width if ind.notation.adoption is not None else None,
            upp=symbols.units_per_point(style.text_scale),
        )
        # Measured from a square's box: a diamond reaches further (#111).
        return near - symbols.pgs_box_near(half, Sex.MALE, False, size, side)

    def _line(self, sid: str) -> frozenset[str]:
        """``sid`` and all their descendants."""
        line = self._lines.get(sid)
        if line is None:
            members = {sid}
            for ui in self.unions_of.get(sid, []):
                for cid in self.unions[ui].children_ids:
                    if cid in self.gen and cid not in members:
                        members |= self._line(cid)
            line = self._lines[sid] = frozenset(members)
        return line

    def _partner(self, ui: int, sid: str) -> str:
        p1, p2 = self.unions[ui].partner_ids
        return p2 if p1 == sid else p1

    def _twin_outer_sides(self) -> dict[str, str]:
        """Side on which each end twin of a twin group takes their partner.

        Twins are drawn adjacent (``_order_children``), and
        ``_children_group`` never moves them, so the first twin of a group
        is its leftmost and the last its rightmost. Their partners go on
        the outer side, so no spouse is drawn between twins and the twin
        connector does not span across one. Middle triplets have no outer
        side and keep the default.
        """
        sides: dict[str, str] = {}
        for union in self.unions:
            children = [c for c in union.children_ids if c in self.gen]
            if len(children) < 2:
                continue
            order = _order_children(self.family, children)
            for group, _tt in _find_twin_groups_in_list(self.family, order):
                members = [c for c in order if c in group]
                sides.setdefault(members[0], "left")
                sides.setdefault(members[-1], "right")
        return sides

    def _attached_to_lineage(self) -> set[str]:
        """Individuals reachable from someone with parents through partnerships.

        They are placed inside the block of that lineage (as partners, or
        as duplicates of a partner), so they must not start a block of
        their own.
        """
        links: dict[str, list[str]] = {}
        for union in self.unions:
            p1, p2 = union.partner_ids
            links.setdefault(p1, []).append(p2)
            links.setdefault(p2, []).append(p1)
        for original, dups in self.duplicates_of.items():
            for dup in dups:
                links.setdefault(original, []).append(dup)
                links.setdefault(dup, []).append(original)
        attached = set(self.parent_union)
        queue = deque(attached)
        while queue:
            for other in links.get(queue.popleft(), []):
                if other not in attached:
                    attached.add(other)
                    queue.append(other)
        return attached

    def run(self) -> dict[str, float]:
        layout = _Block()
        gap = self.style.sibling_gap
        roots: list[str] = []
        for union in self.unions:
            roots.extend(union.partner_ids[:1])
        roots.extend(self.family.individuals)
        for sid in roots:
            if sid in self.gen and sid not in self.claimed and sid not in self.attached:
                for block in self._place(sid):
                    layout.place(block, gap)
        # Whatever remains belongs to structures with no founder to start
        # from (e.g. a founder partnered with a descendant): start from the
        # highest generation.
        order = {sid: i for i, sid in enumerate(self.family.individuals)}
        for sid in sorted(self.gen, key=lambda s: (self.gen[s], order.get(s, 0))):
            if sid not in self.claimed:
                for block in self._place(sid):
                    layout.place(block, gap)
        if not layout.xs:
            return {}
        shift = min(layout.xs.values())
        return {sid: x - shift for sid, x in layout.xs.items()}

    def loops_to_break(
        self, x: dict[str, float], thorough: bool = False
    ) -> list[tuple[Union, str]]:
        """Loop partnerships whose drawing makes two sibship lines overlap.

        Returns (partnership, partner) pairs: drawing that partner a second
        time, next to the other partner, turns the loop into a tree. By
        default one loop per tangle (group of overlapping lines), so that
        no loop is broken that another break would have made unnecessary;
        ``thorough`` returns every loop stretching a tangled line.
        """
        if not self.loops:
            return []
        levels: dict[int, list[tuple[float, float, int]]] = {}
        for ui, union in enumerate(self.unions):
            kids = [x[c] for c in union.children_ids if c in x]
            if not kids:
                continue
            p1, p2 = union.partner_ids
            center = (x[p1] + x[p2]) / 2
            lo, hi = min(*kids, center), max(*kids, center)
            levels.setdefault(self.gen[p1], []).append((lo, hi, ui))
        width: dict[int, float] = {}
        root: dict[int, int] = {}

        def find(u: int) -> int:
            while root.setdefault(u, u) != u:
                root[u] = root[root[u]]
                u = root[u]
            return u

        for spans in levels.values():
            spans.sort()
            for lo, hi, ui in spans:
                width[ui] = hi - lo
            for i, (_lo, hi, ui) in enumerate(spans):
                for lo2, _hi2, uj in spans[i + 1 :]:
                    if lo2 >= hi - 1e-6:
                        break
                    root[find(ui)] = find(uj)
        # Lines that overlap, directly or through others, form a tangle.
        tangles: dict[int, list[int]] = {}
        for ui in root:
            tangles.setdefault(find(ui), []).append(ui)
        fixes: dict[tuple[int, str], None] = {}
        for members in tangles.values():
            stretched = [u for u in members if self.loops.get(u)]
            if not stretched:
                continue
            if thorough:
                for u in stretched:
                    fixes.update(dict.fromkeys(self.loops[u]))
            else:
                # Breaking the loop that stretches the widest line often
                # untangles the rest: break that one, and lay out again.
                widest = max(stretched, key=lambda u: (width[u], -u))
                fixes[self.loops[widest][0]] = None
        return [(self.unions[ui], pid) for ui, pid in fixes]

    def _place(self, sid: str) -> list[_Block]:
        """Lay out ``sid``, their partners and descendants.

        Returns the block, followed by the blocks of the individual's
        duplicates (drawn for second and later partnerships), which the
        caller packs immediately to the right so they stay close.
        """
        style = self.style
        self.claimed.add(sid)
        gen = self.gen[sid]

        left: tuple[str, int] | None = None
        right: tuple[str, int] | None = None
        remote: list[int] = []
        for ui in self.unions_of.get(sid, []):
            if ui in self.done_unions:
                continue
            self.done_unions.add(ui)
            other = self._partner(ui, sid)
            if other in self.claimed or self.gen[other] != gen:
                remote.append(ui)
                if self.gen[other] == gen:
                    self.loops.setdefault(ui, []).append((ui, other))
                continue
            if other in self.parent_union:
                # Loop: the partner's own family is not placed yet and will
                # be packed to the right, so put them on that side.
                side = "right"
            elif sid in self.outer_side:
                # An end twin: partner outside the twin group.
                side = self.outer_side[sid]
            else:
                side = "right" if self.unions[ui].partner_ids[0] == sid else "left"
            if side == "right" and right is None:
                right = (other, ui)
            elif side == "left" and left is None:
                left = (other, ui)
            elif right is None:
                right = (other, ui)
            elif left is None:
                left = (other, ui)
            else:
                remote.append(ui)
                continue
            self.claimed.add(other)
            if other in self.parent_union:
                self.attached_by[other] = ui

        block = _Block()
        gap, cg = style.sibling_gap, style.couple_gap
        # The partner on the right of each couple keeps room for their
        # proband marker from the one on the left, and both for their labels.
        cg_l = self._couple_gap(left[0], sid, left[1]) if left is not None else cg
        cg_r = self._couple_gap(sid, right[0], right[1]) if right is not None else cg
        c_left = c_right = None
        if left is not None:
            group, center = self._children_group(left[1])
            if center is not None:
                c_left = center + block.place(group, gap)
        if right is not None:
            group, center = self._children_group(right[1])
            if center is not None:
                target = None if c_left is None else c_left + max(cg_l, cg_r) - center
                c_right = center + block.place(group, gap, target)

        if c_left is not None and c_right is not None:
            x = (c_left + c_right) / 2
            x_left, x_right = 2 * c_left - x, 2 * c_right - x
        elif c_left is not None:
            x = c_left + cg_l / 2
            x_left, x_right = c_left - cg_l / 2, x + cg_r
        elif c_right is not None:
            x = c_right - cg_r / 2
            x_left, x_right = x - cg_l, c_right + cg_r / 2
        else:
            x = 0.0
            x_left, x_right = -cg_l, cg_r
        # The right-hand partner of a couple draws their PGS box on the
        # right (#78): ``sid`` with a partner on the left only, or the
        # partner on the right.
        self._add(block, sid, x, gen, pgs_right=left is not None and right is None)
        if left is not None:
            self._add(block, left[0], x_left, gen)
        if right is not None:
            self._add(block, right[0], x_right, gen, pgs_right=True)

        # Partnerships whose partner is drawn elsewhere: their children
        # still descend from this individual, so keep them underneath.
        for ui in remote:
            group, center = self._children_group(ui)
            if center is not None:
                block.place(group, gap, x - center)

        blocks = [block]
        for member in [sid] + [p[0] for p in (left, right) if p is not None]:
            for dup in self.duplicates_of.get(member, []):
                if dup in self.claimed:
                    continue
                partners = [self._partner(ui, dup) for ui in self.unions_of.get(dup, [])]
                if any(p in self.parent_union for p in partners):
                    continue  # placed next to that partner, inside their lineage
                blocks.extend(self._place(dup))
        return blocks

    def _children_group(self, ui: int) -> tuple[_Block, float | None]:
        """Lay out the children of a union side by side.

        Returns the group and the x its parents should be centred on: the
        midpoint of the children placed in it, or None when every child
        is drawn elsewhere.
        """
        union = self.unions[ui]
        children = [c for c in _order_children(self.family, union.children_ids) if c in self.gen]
        twin_groups = _find_twin_groups_in_list(self.family, children)
        twins = {sid for group, _tt in twin_groups for sid in group}

        # Loops: a child whose line leads to someone partnered with a
        # member of a family placed later goes last, and a child whose line
        # leads to someone already drawn next to their partner (in a family
        # placed earlier) goes first. Applied at every generation from the
        # common ancestors down, the two branches of a loop end up side by
        # side, and the line joining them crosses no other sibship.
        # A loop between two of these children's lines (cousins marrying)
        # closes inside this sibship: both ends go just before the children
        # leading out, so that no sibling marrying into a later family is
        # drawn between them.
        lines = {c: self._line(c) for c in children}
        inside = frozenset().union(*lines.values())

        def loop_partners(cid: str) -> list[str]:
            """Partners, not yet placed, of ``cid``'s line from other lineages."""
            line = lines[cid]
            return [
                q
                for d in line
                for ui2 in self.unions_of.get(d, [])
                if ui2 not in self.done_unions
                and (q := self._partner(ui2, d)) in self.parent_union
                and q not in self.claimed
                and q not in line
            ]

        movable = [c for c in children if c not in twins]
        first = [c for c in movable if any(d in self.attached_by for d in lines[c])]
        partners = {c: loop_partners(c) for c in movable if c not in first}
        last = [c for c, qs in partners.items() if any(q not in inside for q in qs)]
        within = [c for c, qs in partners.items() if qs and c not in last]
        middle = [c for c in children if c not in first and c not in within and c not in last]
        children = first + middle + within + last

        group = _Block()
        for cid in children:
            if cid in self.claimed:
                continue
            for block in self._place(cid):
                group.place(block, self.style.sibling_gap)
        for cid in union.children_ids:
            if cid not in group.xs and cid in self.attached_by:
                self.loops.setdefault(ui, []).append((self.attached_by[cid], cid))
        placed = [group.xs[c] for c in union.children_ids if c in group.xs]
        if not placed:
            return group, None
        return group, (min(placed) + max(placed)) / 2


_MAX_LOOP_BREAKING_PASSES = 40


def _place_hierarchically(
    family: Family, style: PedigreeStyle, labels: bool = True
) -> tuple[Family, dict[str, int], dict[str, float]]:
    """Assign generations and x positions, breaking loops that tangle lines.

    A loop partner is first drawn next to their partner, with their
    parents' sibship line stretched to reach them. Where that line would
    overlap another sibship, the partner is instead drawn twice, as the
    pedigree nomenclature allows (Bennett et al. 2008): once among their
    siblings and once as a duplicate next to their partner, the two joined
    by a dotted arc. The input family is not modified.

    Returns:
        The family as laid out (with any added duplicates), its generation
        assignments, and its x positions.
    """
    for attempt in range(_MAX_LOOP_BREAKING_PASSES):
        generations = _assign_generations(family)
        placer = _HierarchicalPlacer(family, generations, style, labels)
        x_positions = placer.run()
        thorough = attempt >= _MAX_LOOP_BREAKING_PASSES // 2
        to_duplicate = placer.loops_to_break(x_positions, thorough)
        if not to_duplicate:
            break
        family = _with_duplicate_partners(family, to_duplicate)
    return family, generations, x_positions


def _with_duplicate_partners(family: Family, items: list[tuple[Union, str]]) -> Family:
    """Copy of ``family`` where each (union, partner) uses a duplicate partner.

    The duplicate takes the partner's place in that union only, and
    becomes the parent of its children, as the parser does for second
    partnerships (see ``parser._create_duplicate_nodes``).
    """
    individuals = dict(family.individuals)
    unions = list(family.unions)
    for union, pid in items:
        index = next(i for i, u in enumerate(unions) if u is union)
        original = individuals[pid]
        other = union.partner_ids[1] if union.partner_ids[0] == pid else union.partner_ids[0]
        dup_id = f"_dup_{pid}_{other}"
        suffix = 2
        while dup_id in individuals:
            dup_id = f"_dup_{pid}_{other}_{suffix}"
            suffix += 1
        individuals[dup_id] = Individual(
            sample_id=dup_id,
            family_id=family.family_id,
            sex=original.sex,
            phenotype=original.phenotype,
            is_deceased=original.is_deceased,
            is_proband=False,
            is_duplicate=True,
            is_carrier=original.is_carrier,
            is_obligate_carrier=original.is_obligate_carrier,
            # Drawn identically to the original, phenotype bands included.
            extra_phenotypes=dict(original.extra_phenotypes),
            extra_phenotypes_raw=dict(original.extra_phenotypes_raw),
            extra_phenotypes_categorical=dict(original.extra_phenotypes_categorical),
            not_applicable=set(original.not_applicable),
            duplicate_of=original.duplicate_of or pid,
            notation=replace(original.notation, is_consultand=False),
            css_classes=list(original.css_classes),
            data_attributes=dict(original.data_attributes),
            is_inferred=original.is_inferred,
        )
        p1, p2 = union.partner_ids
        unions[index] = replace(
            union,
            partner_ids=(dup_id if p1 == pid else p1, dup_id if p2 == pid else p2),
            children_ids=list(union.children_ids),
        )
        for cid in union.children_ids:
            child = individuals.get(cid)
            if child is None:
                continue
            if child.father_id == pid:
                individuals[cid] = replace(child, father_id=dup_id)
            elif child.mother_id == pid:
                individuals[cid] = replace(child, mother_id=dup_id)
    return Family(family_id=family.family_id, individuals=individuals, unions=unions)


# Sample ID of the placeholder child; underscores are added if a real
# individual already uses it.
UNKNOWN_OFFSPRING_ID = "_unknown_offspring"


def with_unknown_offspring(family: Family) -> Family:
    """Add a "?" child under a couple that is alone in its family.

    A family holding only a man and a woman, both founders, with no child
    in the file (a trio whose child is missing, say), is drawn as their
    couple with a descent line to a "?", not as two unrelated symbols. The
    father stands on the left, as for every couple.

    Args:
        family: The family as parsed. Not modified.

    Returns:
        A copy with the placeholder child and the couple's union, or
        `family` itself when it does not have that shape.
    """
    people = list(family.individuals.values())
    if len(people) != 2 or family.unions:
        return family
    if any(not p.is_founder or p.is_phantom or p.is_duplicate for p in people):
        return family
    fathers = [p for p in people if p.sex == Sex.MALE]
    mothers = [p for p in people if p.sex == Sex.FEMALE]
    if len(fathers) != 1 or len(mothers) != 1:
        return family
    father, mother = fathers[0].sample_id, mothers[0].sample_id
    child_id = UNKNOWN_OFFSPRING_ID
    while child_id in family.individuals:
        child_id += "_"
    child = Individual(
        sample_id=child_id,
        family_id=family.family_id,
        father_id=father,
        mother_id=mother,
        is_unknown_offspring=True,
    )
    return replace(
        family,
        individuals={**family.individuals, child_id: child},
        unions=[Union(partner_ids=(father, mother), children_ids=[child_id])],
    )


def _get_node_label(sid: str, ind: Individual, style: PedigreeStyle | None = None) -> str:
    """Determine the display label for a node.

    - Phantom parents (generated IDs starting with "_phantom_"): no label
    - Duplicate nodes: show the original individual's name
    - Normal individuals: show their sample_id

    Args:
        sid: The node's sample_id.
        ind: The Individual object.
        style: Style parameters (for zero_pad_ids).

    Returns:
        The label string (may be empty for phantoms).
    """
    if sid.startswith("_phantom_") or ind.is_unknown_offspring:
        return ""
    label = ind.duplicate_of if (ind.is_duplicate and ind.duplicate_of) else sid
    if style and style.zero_pad_ids > 0:
        return label.zfill(style.zero_pad_ids)
    return label


def _build_nodes(
    family: Family,
    x_positions: dict[str, float],
    generations: dict[str, int],
    style: PedigreeStyle,
) -> list[NodeLayout]:
    """Build NodeLayout objects for all individuals.

    Args:
        family: The family data.
        x_positions: Horizontal positions.
        generations: Generation assignments.
        style: Style parameters.

    Returns:
        List of NodeLayout objects.
    """
    right_partners = _pgs_right_side(family, x_positions)

    nodes = []
    for sid, ind in family.individuals.items():
        x = x_positions.get(sid, 0.0)
        gen = generations.get(sid, 0)
        y = -gen * style.generation_gap  # Negative so gen 0 is at top

        nodes.append(
            NodeLayout(
                individual_id=sid,
                x=x,
                y=y,
                sex=symbols.symbol_sex(ind.sex, ind.notation),
                phenotype=ind.phenotype,
                is_deceased=symbols.is_struck_through(ind),
                is_proband=ind.is_proband,
                is_phantom=ind.is_phantom,
                is_inferred=ind.is_inferred,
                is_duplicate=ind.is_duplicate,
                is_carrier=ind.is_carrier,
                is_obligate_carrier=ind.is_obligate_carrier,
                label=_get_node_label(sid, ind, style),
                extra_phenotypes=dict(ind.extra_phenotypes),
                extra_phenotypes_raw=dict(ind.extra_phenotypes_raw),
                extra_phenotypes_categorical=dict(ind.extra_phenotypes_categorical),
                not_applicable=set(ind.not_applicable),
                pgs_side="right" if sid in right_partners else "left",
                annotations=[v for v in ind.annotations.values() if v],
                notation=ind.notation,
                notes=symbols.notes(ind),
                marker_depth=symbols.childless_depth(ind.notation),
                is_unknown_offspring=ind.is_unknown_offspring,
                css_classes=list(ind.css_classes),
                data_attributes=dict(ind.data_attributes),
            )
        )

    return nodes


# Room below each line of notes or annotations, beyond its step.
_LINE_BELOW_MARGIN = 0.05


def _pgs_right_side(family: Family, positions: Mapping[str, float]) -> set[str]:
    """Who draws their PGS box on the right (#78): the right-hand partner of
    a couple, unless they also have a partner on their right, so that no box
    sits on a line between partners."""
    right: set[str] = set()
    left: set[str] = set()
    for union in family.unions:
        p1, p2 = union.partner_ids
        if p1 in positions and p2 in positions:
            lo, hi = (p1, p2) if positions[p1] <= positions[p2] else (p2, p1)
            left.add(lo)
            right.add(hi)
    return right - left


def check_labels(style: PedigreeStyle) -> None:
    """Reject an unknown ``labels`` mode (#125) before anything is laid out.

    Raises:
        PedigreeError: If ``style.labels`` is not one of `symbols.LABEL_MODES`.
    """
    if style.labels not in symbols.LABEL_MODES:
        raise PedigreeError(
            f"Unknown labels {style.labels!r}: expected one of {', '.join(symbols.LABEL_MODES)}"
        )


def _numbered(
    nodes: list[NodeLayout],
    family: Family,
    generations: Mapping[str, int],
    x_positions: Mapping[str, float],
    style: PedigreeStyle,
) -> list[NodeLayout]:
    """Label the nodes with their NSGC numbers (#125), with ``labels``
    "generation" or "both".

    Every labelled node is numbered (a named parent with no row too; not the
    unnamed phantoms nor the unknown-offspring placeholder), by its place in
    the drawing; a duplicate takes its original's number. With "both" the ID
    is written right under the number, as the first note ("SB", "AMAB"), in
    their smaller text.
    """
    if style.labels == "id":
        return nodes
    numbered = [
        sid
        for sid, ind in family.individuals.items()
        if not ind.is_duplicate and _get_node_label(sid, ind, style)
    ]
    numbers = symbols.generation_numbers(generations, x_positions, numbered)
    out = []
    for node in nodes:
        ind = family.individuals.get(node.individual_id)
        original = ind.duplicate_of if ind is not None and ind.is_duplicate else None
        number = numbers.get(original or node.individual_id)
        if number is None:
            out.append(node)
            continue
        notes = [node.label, *node.notes] if style.labels == "both" else node.notes
        out.append(replace(node, label=number, number=number, notes=notes))
    return out


def _room_below_labels(family: Family, style: PedigreeStyle, extra_lines: int = 0) -> float:
    """Extra generation gap for what sits below the symbols.

    Each line of notes or annotations takes its step below the line above
    (`symbols.below_label_offsets`, at the style's text scale: 0.25, more for text
    enlarged by ``min_font_size``, #66) plus `_LINE_BELOW_MARGIN`, and a
    childlessness marker its own depth, for the individual who needs the
    most room. `extra_lines` more for everyone (the ID under an NSGC number,
    #125).
    """
    size, small = style.text_size(), style.text_size(symbols.SMALL_TEXT_SCALE)

    def text_room(lines: int) -> float:
        if not lines:
            return 0.0
        upp = symbols.units_per_point(style.text_scale)
        offsets = symbols.below_label_offsets(size, [small] * lines, upp)
        return offsets[-1] + lines * _LINE_BELOW_MARGIN

    return max(
        (
            text_room(symbols.text_lines_below(ind) + extra_lines)
            + symbols.childless_depth(ind.notation)
            for ind in family.individuals.values()
        ),
        default=0.0,
    )


def _individual_line_style(family: Family, child_id: str) -> LineStyle:
    """Dashed for a child adopted in: the parents drawn are adoptive."""
    child = family.individuals.get(child_id)
    if child is not None and child.notation.adoption == Adoption.IN:
        return LineStyle.DASHED
    return LineStyle.SOLID


def _couple_childless_edges(
    family: Family, nodes: list[NodeLayout], style: PedigreeStyle, *, outward: bool = False
) -> list[EdgeLayout]:
    """A couple's childlessness marker (#104): from the middle of the line
    between the partners, down (outward in the circular layout), one or two
    bars."""
    half = style.shape_size / 2
    at = {n.individual_id: (n.x, n.y) for n in nodes}
    edges = []
    for union in family.unions:
        if union.childlessness is None or union.children_ids:
            continue
        a, b = union.partner_ids
        if a not in at or b not in at:
            continue
        x = (at[a][0] + at[b][0]) / 2
        y = (at[a][1] + at[b][1]) / 2
        direction = -math.pi / 2
        if outward and (x or y):
            direction = math.atan2(y, x)
        for start, end in symbols.couple_childless_segments(
            x, y, half, union.childlessness, direction
        ):
            edges.append(EdgeLayout(start=start, end=end, edge_type=EdgeType.CHILDLESS))
    return edges


def _childless_edges(nodes: list[NodeLayout], style: PedigreeStyle) -> list[EdgeLayout]:
    """Childlessness markers: a stem below the symbol and one or two bars.

    The marker leaves the symbol on the side of its label: below, or
    outward in the circular layout.
    """
    half = style.shape_size / 2
    edges = []
    for node in nodes:
        kind = node.notation.childlessness
        if kind is None or node.notation.childless_partner is not None:
            continue  # a couple's marker hangs from their line (#104)
        direction = node.label_angle if node.label_angle is not None else -math.pi / 2
        for start, end in symbols.childless_segments(
            node.x, node.y, node.sex, half, kind, direction
        ):
            edges.append(EdgeLayout(start=start, end=end, edge_type=EdgeType.CHILDLESS))
    return edges


def _compute_edges(
    family: Family,
    x_positions: dict[str, float],
    generations: dict[str, int],
    style: PedigreeStyle,
    twin_groups: list[tuple[set[str], TwinType]] | None = None,
) -> list[EdgeLayout]:
    """Compute all connection edges for the family.

    Edge types:
    - RELATIONSHIP: horizontal line between partners
    - DESCENT: vertical line from relationship midpoint downward
    - SIBSHIP: horizontal line at the bottom of descent connecting siblings
    - INDIVIDUAL: vertical lines from sibship line to each child

    Args:
        family: The family data.
        x_positions: Horizontal positions.
        generations: Generation assignments.
        style: Style parameters.

    Returns:
        List of EdgeLayout objects.
    """
    edges: list[EdgeLayout] = []
    half_shape = style.shape_size / 2

    def top_of(child_id: str, child_gen: int) -> float:
        # A child pulled to a lower generation (partner of someone from a
        # younger generation) is reached on their own row.
        row = max(generations.get(child_id, child_gen), child_gen)
        return -row * style.generation_gap + half_shape

    for union in family.unions:
        p1, p2 = union.partner_ids
        if p1 not in x_positions or p2 not in x_positions:
            continue

        p1_x = x_positions[p1]
        p2_x = x_positions[p2]
        gen = generations.get(p1, 0)
        y = -gen * style.generation_gap

        # Ensure p1 is left, p2 is right
        left_x = min(p1_x, p2_x)
        right_x = max(p1_x, p2_x)

        # 1. Relationship line (horizontal between partners)
        line_style = LineStyle.SOLID
        if union.is_consanguinous:
            line_style = LineStyle.DOUBLE
        elif union.is_broken:
            line_style = LineStyle.BROKEN

        edges.append(
            EdgeLayout(
                start=(left_x + half_shape, y),
                end=(right_x - half_shape, y),
                edge_type=EdgeType.RELATIONSHIP,
                line_style=line_style,
                has_descent=bool(union.children_ids),
            )
        )

        if not union.children_ids:
            continue

        union_children_set = set(union.children_ids)

        # 2. Line of descent (vertical from relationship midpoint down)
        couple_center_x = (left_x + right_x) / 2
        child_gen = gen + 1
        child_y = -child_gen * style.generation_gap

        # Descent goes from relationship line down to sibship level
        descent_bottom_y = child_y + half_shape + style.generation_gap * 0.3

        # For consanguinous unions, start below the bottom double-line
        descent_top_y = y
        if union.is_consanguinous:
            descent_top_y = y - style.double_line_gap

        edges.append(
            EdgeLayout(
                start=(couple_center_x, descent_top_y),
                end=(couple_center_x, descent_bottom_y),
                edge_type=EdgeType.DESCENT,
            )
        )

        # Identify twin members in this union's children
        twin_member_ids: set[str] = set()
        # Map child_id -> set of twin group members (for finding convergence x)
        child_twin_group: dict[str, set[str]] = {}
        if twin_groups:
            for group_ids, _tt in twin_groups:
                for cid in group_ids:
                    if cid in union_children_set:
                        twin_member_ids.add(cid)
                        child_twin_group[cid] = group_ids

        # Non-twin children get regular individual lines
        non_twin_children = [cid for cid in union.children_ids if cid not in twin_member_ids]

        # Build "connection points" on the sibship line:
        # - Non-twin children: their x position
        # - Twin groups: their convergence x (center of twin group)
        connection_xs: list[float] = []
        for cid in non_twin_children:
            connection_xs.append(x_positions.get(cid, 0.0))

        seen_twin_groups: set[frozenset[str]] = set()
        if twin_groups:
            for group_ids, _tt in twin_groups:
                group_in_union = group_ids & union_children_set
                if not group_in_union:
                    continue
                key = frozenset(group_in_union)
                if key in seen_twin_groups:
                    continue
                seen_twin_groups.add(key)
                twin_xs = [x_positions.get(s, 0.0) for s in group_in_union]
                convergence_x = (min(twin_xs) + max(twin_xs)) / 2
                connection_xs.append(convergence_x)

        if len(connection_xs) <= 1:
            # Single connection point: one child, or one twin group.
            # Pedigree nomenclature uses orthogonal lines only. When the
            # child could not be centred under its parents, bridge the
            # offset with a short horizontal segment at sibship level, as a
            # multi-child sibship line does, instead of a diagonal. The same
            # segment connects a twin group whose convergence point is off
            # the descent line, which would otherwise leave a gap.
            if connection_xs and not math.isclose(connection_xs[0], couple_center_x):
                edges.append(
                    EdgeLayout(
                        start=(min(connection_xs[0], couple_center_x), descent_bottom_y),
                        end=(max(connection_xs[0], couple_center_x), descent_bottom_y),
                        edge_type=EdgeType.SIBSHIP,
                    )
                )
            if non_twin_children:
                child_id = non_twin_children[0]
                child_x = x_positions.get(child_id, couple_center_x)
                edges.append(
                    EdgeLayout(
                        start=(child_x, descent_bottom_y),
                        end=(child_x, top_of(child_id, child_gen)),
                        edge_type=EdgeType.INDIVIDUAL,
                        line_style=_individual_line_style(family, child_id),
                    )
                )
            # If only twins (no non-twin children), the twin connector draws
            # the V-lines down from its convergence point.
        else:
            # Multiple connection points: sibship line + individual lines
            # Include descent point (couple center) so the sibship bar
            # always connects to the vertical descent line, even when
            # children were shifted right to avoid overlap.
            leftmost = min(min(connection_xs), couple_center_x)
            rightmost = max(max(connection_xs), couple_center_x)

            # 3. Sibship line (horizontal at descent_bottom_y)
            edges.append(
                EdgeLayout(
                    start=(leftmost, descent_bottom_y),
                    end=(rightmost, descent_bottom_y),
                    edge_type=EdgeType.SIBSHIP,
                )
            )

            # 4. Individual lines for non-twin children only
            for child_id in non_twin_children:
                child_x = x_positions.get(child_id, 0.0)
                edges.append(
                    EdgeLayout(
                        start=(child_x, descent_bottom_y),
                        end=(child_x, top_of(child_id, child_gen)),
                        edge_type=EdgeType.INDIVIDUAL,
                        line_style=_individual_line_style(family, child_id),
                    )
                )

    return edges


def _order_children(family: Family, children_ids: list[str]) -> list[str]:
    """Order children by DOB (oldest first) with twins forced adjacent.

    Sorting priority:
    1. DOB-based ordering: oldest (smallest DOB string) on the left,
       youngest on the right. Individuals without DOB sort after those
       with DOB, preserving their original order.
    2. Twin adjacency: after DOB sorting, twins are grouped together.
       The group is placed at the position of the first twin in the
       DOB-sorted order.

    Args:
        family: The family data (for accessing Individual.date_of_birth
            and twin_info).
        children_ids: The children to order.

    Returns:
        Reordered list of children IDs.
    """
    if len(children_ids) <= 1:
        return list(children_ids)

    # Build twin groups for these children
    twin_groups = _find_twin_groups_in_list(family, children_ids)
    twin_member_to_group: dict[str, int] = {}
    for gi, (group_ids, _tt) in enumerate(twin_groups):
        for sid in group_ids:
            twin_member_to_group[sid] = gi

    # Sort by DOB: individuals with DOB come first, sorted ascending.
    # Individuals without DOB maintain original order at the end.
    def sort_key(sid: str) -> tuple[int, str, int]:
        ind = family.individuals.get(sid)
        if ind and ind.date_of_birth:
            return (0, ind.date_of_birth, 0)
        # No DOB: use original position to maintain stability
        return (1, "", children_ids.index(sid))

    sorted_ids = sorted(children_ids, key=sort_key)

    # Now group twins together: for each twin group, place all members
    # at the position of the first member in sorted order
    if not twin_groups:
        return sorted_ids

    result: list[str] = []
    placed: set[str] = set()

    for sid in sorted_ids:
        if sid in placed:
            continue

        if sid in twin_member_to_group:
            gi = twin_member_to_group[sid]
            group_ids, _tt = twin_groups[gi]
            # Place all twins in this group together, preserving their
            # relative DOB order within the group
            group_sorted = [s for s in sorted_ids if s in group_ids and s not in placed]
            result.extend(group_sorted)
            placed.update(group_sorted)
        else:
            result.append(sid)
            placed.add(sid)

    return result


def _find_twin_groups_in_list(
    family: Family, children_ids: list[str]
) -> list[tuple[set[str], TwinType]]:
    """Find twin groups among a specific list of children.

    Uses union-find to merge overlapping twin declarations into groups.

    Args:
        family: The family data.
        children_ids: Children to search among.

    Returns:
        List of (group_member_ids, twin_type) tuples.
    """
    children_set = set(children_ids)
    # Build adjacency from twin_info
    edges: list[tuple[str, str, TwinType]] = []
    for cid in children_ids:
        ind = family.individuals.get(cid)
        if not ind:
            continue
        for ti in ind.twin_info:
            if ti.sibling_id in children_set:
                edges.append((cid, ti.sibling_id, ti.twin_type))

    if not edges:
        return []

    # Union-find
    parent: dict[str, str] = {}

    def find(x: str) -> str:
        while parent.get(x, x) != x:
            parent[x] = parent.get(parent[x], parent[x])
            x = parent[x]
        return x

    def unite(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    group_types: dict[str, TwinType] = {}
    for a, b, tt in edges:
        unite(a, b)
        # Use the first twin type encountered for the group
        root = find(a)
        if root not in group_types:
            group_types[root] = tt

    # Collect groups
    groups: dict[str, set[str]] = {}
    for a, b, _tt in edges:
        root = find(a)
        groups.setdefault(root, set()).add(a)
        groups.setdefault(root, set()).add(b)

    return [(members, group_types[root]) for root, members in groups.items()]


def _find_twin_groups(family: Family) -> list[tuple[set[str], TwinType]]:
    """Find all twin groups in a family.

    Scans all individuals' twin_info to discover connected twin groups.

    Args:
        family: The family data.

    Returns:
        List of (group_member_ids, twin_type) tuples.
    """
    all_ids = list(family.individuals.keys())
    return _find_twin_groups_in_list(family, all_ids)


def _compute_twin_connectors(
    family: Family,
    x_positions: dict[str, float],
    generations: dict[str, int],
    style: PedigreeStyle,
    twin_groups: list[tuple[set[str], TwinType]],
) -> list[TwinConnector]:
    """Compute twin connector geometry for all twin groups.

    For each twin group, computes a convergence point on the sibship line
    and V-lines from there down to each twin's shape.

    The convergence point is at the horizontal center of the twin group,
    at the sibship line y-level (descent_bottom_y).

    Args:
        family: The family data.
        x_positions: Horizontal positions.
        generations: Generation assignments.
        style: Style parameters.
        twin_groups: Twin groups from _find_twin_groups().

    Returns:
        List of TwinConnector objects.
    """
    if not twin_groups:
        return []

    half_shape = style.shape_size / 2
    connectors: list[TwinConnector] = []

    for group_ids, twin_type in twin_groups:
        # All twins in a group should share the same generation
        gens = {generations.get(sid, 0) for sid in group_ids}
        if not gens:
            continue
        child_gen = min(gens)
        child_y = -child_gen * style.generation_gap

        # Convergence at sibship level (same formula as _compute_edges)
        descent_bottom_y = child_y + half_shape + style.generation_gap * 0.3

        # Convergence point: centered above twin group at sibship level
        twin_xs = [x_positions.get(sid, 0.0) for sid in group_ids]
        convergence_x = (min(twin_xs) + max(twin_xs)) / 2
        convergence_point = (convergence_x, descent_bottom_y)

        # Child positions: top-center of each twin's shape
        child_positions = [
            (x_positions.get(sid, 0.0), child_y + half_shape)
            for sid in sorted(group_ids, key=lambda s: x_positions.get(s, 0.0))
        ]

        connectors.append(
            TwinConnector(
                convergence_point=convergence_point,
                child_positions=child_positions,
                twin_type=twin_type,
            )
        )

    return connectors


# Largest bulge of a duplicate arc, as a fraction of the generation gap.
_MAX_ARC_BULGE = 0.4


def _compute_duplicate_arcs(
    family: Family,
    positions: dict[str, tuple[float, float]],
    style: PedigreeStyle,
) -> list[DuplicateArc]:
    """Compute dotted arc positions connecting duplicate nodes to their originals.

    For each duplicate individual, creates a DuplicateArc from the original
    node's position to the duplicate's position.

    Args:
        family: The family data.
        positions: (x, y) of every individual, from either layout.
        style: Style parameters.

    Returns:
        List of DuplicateArc objects.
    """
    arcs: list[DuplicateArc] = []

    # Centre of the tree (non-duplicate nodes only), so each arc can bow
    # away from the bulk of the pedigree.
    placed = [
        positions[s] for s, i in family.individuals.items() if not i.is_duplicate and s in positions
    ]
    if placed:
        center_x = sum(p[0] for p in placed) / len(placed)
        center_y = sum(p[1] for p in placed) / len(placed)
    else:
        center_x = 0.0
        center_y = 0.0

    for sid, ind in family.individuals.items():
        if not ind.is_duplicate or not ind.duplicate_of:
            continue

        original_id = ind.duplicate_of
        if original_id not in positions or sid not in positions:
            continue

        orig_x, orig_y = positions[original_id]
        dup_x, dup_y = positions[sid]

        # Choose curvature sign so the arc bows AWAY from the tree centre.
        #
        # Matplotlib arc3 convention: the control point offset from the
        # midpoint is (rad*dy, -rad*dx) — i.e. rad > 0 pushes the arc
        # in the direction of the RIGHT perpendicular (dy, -dx) of A→B.
        #
        # We pick the sign so that offset points away from the tree
        # centre: compute the right-perpendicular of A→B and check if
        # it aligns with the outward vector (midpoint → away from centre).
        dx = dup_x - orig_x
        dy = dup_y - orig_y
        mid_x = (orig_x + dup_x) / 2.0
        mid_y = (orig_y + dup_y) / 2.0
        # Right-perpendicular of (dx, dy) — matches matplotlib's arc3
        rperp_x, rperp_y = dy, -dx
        # Vector from tree centre to arc midpoint
        out_x = mid_x - center_x
        out_y = mid_y - center_y
        dot = rperp_x * out_x + rperp_y * out_y
        # arc3 bows by |rad| * length / 2 at its midpoint. Cap that bulge so
        # a long arc stays within the gap between two generations instead
        # of sweeping over several rows.
        length = math.hypot(dx, dy)
        magnitude = 0.3
        if length > 0:
            magnitude = min(magnitude, 2 * _MAX_ARC_BULGE * style.generation_gap / length)
        curvature = magnitude if dot >= 0 else -magnitude

        arcs.append(
            DuplicateArc(
                from_pos=(orig_x, orig_y),
                to_pos=(dup_x, dup_y),
                curvature=curvature,
            )
        )

    return arcs


# Routing of the dotted lines between a person's two copies (#149): where the
# line leaves the top of a copy, in half symbol sizes from its centre (off its
# own vertical line, and short of a deceased line's end at the corner); the
# share of the band above the row kept clear at its bottom and top; the
# radius of the line's corners; and the room kept between two lines' spans.
_ARC_EXIT = 0.5
_ARC_BAND_MARGIN = 0.2
_ARC_CORNER = 0.1
_ARC_LANE_PAD = 0.1


def _route_duplicate_arcs(arcs: list[DuplicateArc], style: PedigreeStyle) -> list[DuplicateArc]:
    """Route each dotted line through the band above its row (#149).

    The two copies of a person are always in the same row, often far apart;
    a curve from centre to centre ran through every symbol in between. The
    line now leaves the top of each copy `_ARC_EXIT` half-sizes towards the
    other, rises into the band between the symbols' tops and their sibship
    lines (`generation_gap` / 2 above the row), which holds no symbol, runs
    across and comes down, with rounded corners. Lines of a row whose spans
    overlap get lanes at different heights, the shorter ones lower, so a
    line nested in another stays inside it.
    """
    half = style.shape_size / 2
    rows: dict[float, list[int]] = {}
    for i, arc in enumerate(arcs):
        rows.setdefault(round(arc.from_pos[1], 6), []).append(i)
    routed = list(arcs)
    for row_y, members in rows.items():
        spans = {}
        for i in members:
            (x1, _), (x2, _) = arcs[i].from_pos, arcs[i].to_pos
            side = 1.0 if x2 >= x1 else -1.0
            spans[i] = (x1 + side * half * _ARC_EXIT, x2 - side * half * _ARC_EXIT)
        lanes: list[list[tuple[float, float]]] = []
        lane_of: dict[int, int] = {}
        for i in sorted(members, key=lambda j: abs(spans[j][1] - spans[j][0])):
            lo, hi = sorted(spans[i])
            for k, taken in enumerate(lanes):
                if all(hi + _ARC_LANE_PAD < a or b + _ARC_LANE_PAD < lo for a, b in taken):
                    taken.append((lo, hi))
                    lane_of[i] = k
                    break
            else:
                lanes.append([(lo, hi)])
                lane_of[i] = len(lanes) - 1
        bottom = row_y + half
        band = style.generation_gap / 2 - half
        low, high = bottom + band * _ARC_BAND_MARGIN, bottom + band * (1 - _ARC_BAND_MARGIN)
        for i in members:
            k = lane_of[i]
            lane = low + (high - low) * (k + 1) / (len(lanes) + 1)
            routed[i] = replace(arcs[i], waypoints=_arc_route(*spans[i], bottom, lane))
    return routed


def _arc_route(xa: float, xb: float, bottom: float, lane: float) -> list[tuple[float, float]]:
    """Up from (xa, bottom), across at ``lane``, down to (xb, bottom), with
    quarter-circle corners of radius `_ARC_CORNER` (less when short)."""
    side = 1.0 if xb >= xa else -1.0
    r = min(_ARC_CORNER, (lane - bottom) / 2, abs(xb - xa) / 2)
    steps = 6
    points = [(xa, bottom)]
    for k in range(steps + 1):  # first corner: from going up to going across
        a = math.pi / 2 * k / steps
        points.append((xa + side * r * (1 - math.cos(a)), lane - r + r * math.sin(a)))
    for k in range(steps + 1):  # second corner: from going across to going down
        a = math.pi / 2 * k / steps
        points.append((xb - side * r + side * r * math.sin(a), lane - r * (1 - math.cos(a))))
    points.append((xb, bottom))
    return points


# =====================================================================
# Circular / radial layout
# =====================================================================


def compute_layout_circular(
    family: Family,
    style: PedigreeStyle | None = None,
) -> FamilyLayout:
    """Compute a circular (radial) layout for a family.

    Founders are placed at the centre and descendants radiate outward
    in concentric rings, one ring per generation.

    Experimental: its structure (rings, order, duplicates) is stable, but how
    it is drawn (ring spacing, labels, lines between rings) may change in a
    1.x release. The hierarchical layout (`compute_layout`) is the reference.

    The angular order comes from the hierarchical layout: each individual's
    horizontal position becomes an angle, so the circular layout inherits
    its guarantees (siblings stay together, sibships of different couples
    never interleave, partners sit side by side, tangled loops are drawn
    with duplicates) whatever the number of founder couples.

    Args:
        family: The family to lay out.
        style: Visual style parameters. Uses defaults if None.

    Returns:
        FamilyLayout with all positions and connections computed.
    """
    if style is None:
        style = PedigreeStyle()

    if not family.individuals:
        return FamilyLayout(family_id=family.family_id)
    family = with_unknown_offspring(family)

    # Steps 1-2: Generations and left-to-right order, from the hierarchical
    # placement (the family may gain duplicate nodes there).
    family, generations, x_positions = _place_hierarchically(family, style, labels=False)
    check_labels(style)
    check_text_scale(style)

    # Step 3: Angles and ring radii
    angles, radii = _circular_positions(family, generations, x_positions, style)

    # Step 4: Polar → Cartesian
    xy: dict[str, tuple[float, float]] = {
        sid: (radii[sid] * math.cos(angles[sid]), radii[sid] * math.sin(angles[sid]))
        for sid in family.individuals
    }

    # Step 5: Build node layouts; NSGC numbers (#125) follow the hierarchical
    # order the angles come from, so a family gets the same in both layouts.
    nodes = _build_circular_nodes(family, xy, angles, radii, generations, style)
    nodes = _numbered(nodes, family, generations, x_positions, style)

    # Step 6: Edges and twin connectors
    twin_groups = _find_twin_groups(family)
    edges, twin_connectors = _compute_circular_edges(
        family, xy, angles, radii, generations, style, twin_groups
    )
    edges += _childless_edges(nodes, style)
    edges += _couple_childless_edges(family, nodes, style, outward=True)

    return FamilyLayout(
        family_id=family.family_id,
        nodes=nodes,
        edges=edges,
        twin_connectors=twin_connectors,
        duplicate_arcs=_compute_duplicate_arcs(family, xy, style),
    )


def _circular_positions(
    family: Family,
    generations: dict[str, int],
    x_positions: dict[str, float],
    style: PedigreeStyle,
) -> tuple[dict[str, float], dict[str, float]]:
    """Turn hierarchical positions into angles and ring radii.

    The hierarchical x axis is wrapped around the circle, counter-clockwise
    and centred at the bottom, so the left of the pedigree is on the left.
    A large pedigree uses the whole circle, with one sibling gap left open
    at the top so its first and last individuals do not touch; a small one
    only the arc it needs. Angles are kept continuous (no wrap at ±π), so
    the children of a union always span one unbroken arc.

    A lone founder couple sits at the centre, ``couple_gap`` apart (#151).
    Every other generation gets a
    ring; rings are ``style.radial_gap`` apart, and the first descendant
    ring has ``style.base_radius``, or more if one ring would otherwise be
    too crowded for its nodes to be ``1.8 × shape_size`` apart.

    Returns:
        Tuple of (angles, radii) mapping sample_id to radians and radius.
    """
    min_arc = style.shape_size * 1.8
    xs = list(x_positions.values())
    x_mid = (min(xs) + max(xs)) / 2
    turn = (max(xs) - min(xs)) + style.sibling_gap
    # Radians per unit of x: the whole circle for a large pedigree, but no
    # more than needed for neighbours on the first ring to be min_arc
    # apart, so a small family stays grouped at the bottom instead of being
    # spread around the circle (two siblings would face each other).
    closest = min(style.sibling_gap, style.couple_gap)
    scale = min(2 * math.pi / turn, min_arc / (style.base_radius * closest))
    angles = {sid: 3 * math.pi / 2 + scale * (x - x_mid) for sid, x in x_positions.items()}

    min_gen = min(generations.values())
    founders = [sid for sid, g in generations.items() if g == min_gen]
    centre_couple = len(founders) == 2 and any(
        set(u.partner_ids) == set(founders) for u in family.unions
    )
    radii: dict[str, float] = {}
    if centre_couple:
        left, right = sorted(founders, key=lambda s: x_positions[s])
        angles[left], angles[right] = math.pi, 0.0
        # As far apart as any couple, whose gap holds their PGS boxes and
        # labels (#151); 0.7 symbol sizes each made their boxes overlap.
        radii[left] = radii[right] = max(style.shape_size * 0.7, style.couple_gap / 2)
    first_ring = min_gen + 1 if centre_couple else min_gen

    # Angular room on each ring decides how far out the rings must start.
    by_ring: dict[int, list[float]] = {}
    for sid, g in generations.items():
        if sid not in radii:
            by_ring.setdefault(g, []).append(angles[sid])
    start = style.base_radius
    for g, ring_angles in by_ring.items():
        ring_angles.sort()
        gaps = [b - a for a, b in zip(ring_angles, ring_angles[1:], strict=False)]
        gaps.append(ring_angles[0] + 2 * math.pi - ring_angles[-1])
        smallest = min(gaps)
        if smallest > 0:
            start = max(start, min_arc / smallest - (g - first_ring) * style.radial_gap)
    for sid, g in generations.items():
        if sid not in radii:
            radii[sid] = start + (g - first_ring) * style.radial_gap
    return angles, radii


def _build_circular_nodes(
    family: Family,
    xy: dict[str, tuple[float, float]],
    angles: dict[str, float],
    radii: dict[str, float],
    generations: dict[str, int],
    style: PedigreeStyle,
) -> list[NodeLayout]:
    """Build NodeLayout objects for the circular layout.

    Args:
        family: The family data.
        xy: Cartesian (x, y) positions.
        angles: Angular positions (radians).
        radii: Radial positions.
        generations: Generation assignments.
        style: Style parameters.

    Returns:
        List of NodeLayout objects.
    """
    # By x on the page, not by angle: angles grow right to left at the top of
    # the circle, and the centre couple's left partner is at angle π (#151).
    right_partners = _pgs_right_side(family, {sid: p[0] for sid, p in xy.items()})

    nodes = []
    for sid, ind in family.individuals.items():
        x, y = xy.get(sid, (0.0, 0.0))
        # Shapes stay upright: a square turned to face outward reads as the
        # unknown-sex diamond, and bands would change direction around the
        # circle. The label goes radially outward instead of below.
        r = radii.get(sid, 0.0)
        theta = angles.get(sid, 0.0)
        label_angle = theta if r > style.shape_size else None
        nodes.append(
            NodeLayout(
                individual_id=sid,
                x=x,
                y=y,
                sex=symbols.symbol_sex(ind.sex, ind.notation),
                phenotype=ind.phenotype,
                is_deceased=symbols.is_struck_through(ind),
                is_proband=ind.is_proband,
                is_phantom=ind.is_phantom,
                is_inferred=ind.is_inferred,
                is_duplicate=ind.is_duplicate,
                is_carrier=ind.is_carrier,
                is_obligate_carrier=ind.is_obligate_carrier,
                label=_get_node_label(sid, ind, style),
                extra_phenotypes=dict(ind.extra_phenotypes),
                extra_phenotypes_raw=dict(ind.extra_phenotypes_raw),
                extra_phenotypes_categorical=dict(ind.extra_phenotypes_categorical),
                not_applicable=set(ind.not_applicable),
                label_angle=label_angle,
                pgs_side="right" if sid in right_partners else "left",
                annotations=[v for v in ind.annotations.values() if v],
                notation=ind.notation,
                notes=symbols.notes(ind),
                marker_depth=symbols.childless_depth(ind.notation),
                is_unknown_offspring=ind.is_unknown_offspring,
                css_classes=list(ind.css_classes),
                data_attributes=dict(ind.data_attributes),
            )
        )
    return nodes


def _circular_sex(family: Family, sid: str) -> Sex:
    """The shape drawn for ``sid`` (a duplicate's is its original's)."""
    ind = family.individuals[sid]
    return symbols.symbol_sex(ind.sex, ind.notation)


def _compute_circular_edges(
    family: Family,
    xy: dict[str, tuple[float, float]],
    angles: dict[str, float],
    radii: dict[str, float],
    generations: dict[str, int],
    style: PedigreeStyle,
    twin_groups: list[tuple[set[str], TwinType]] | None = None,
) -> tuple[list[EdgeLayout], list[TwinConnector]]:
    """Compute connection edges and twin connectors for the circular layout.

    Uses proper pedigree structure:
    - RELATIONSHIP: straight line between partners
    - DESCENT: line from the couple midpoint to the middle of the sibship arc
    - SIBSHIP: arc connecting siblings on the same ring
    - INDIVIDUAL: short radial line from sibship arc to each child
    - Twins: a V from one point of the sibship arc to the inner edge of
      each twin, as in the hierarchical layout

    Angles are continuous (see ``_circular_positions``), so the children of
    a union span the arc [min angle, max angle] without wrapping.

    Args:
        family: The family data.
        xy: Cartesian positions.
        angles: Angular positions.
        radii: Radial positions.
        generations: Generation assignments.
        style: Style parameters.
        twin_groups: Twin groups, from ``_find_twin_groups``.

    Returns:
        Tuple of (edges, twin connectors).
    """
    edges: list[EdgeLayout] = []
    twin_connectors: list[TwinConnector] = []
    half = style.shape_size / 2
    # The sibship arc sits this far inside the children's ring: clear of
    # their shapes by a third of the gap between rings, so twin lines form a
    # proper V instead of running nearly along the arc.
    sib_offset = half + 0.3 * style.radial_gap

    for union in family.unions:
        p1, p2 = union.partner_ids
        if p1 not in xy or p2 not in xy:
            continue

        x1, y1 = xy[p1]
        x2, y2 = xy[p2]

        # 1. Relationship line between partners
        line_style = LineStyle.SOLID
        if union.is_consanguinous:
            line_style = LineStyle.DOUBLE
        elif union.is_broken:
            line_style = LineStyle.BROKEN

        children = [cid for cid in union.children_ids if cid in xy]
        start, end = (x1, y1), (x2, y2)
        if line_style == LineStyle.BROKEN:
            # From outline to outline, so the separation marks (#141), placed
            # along the line, fall between the partners and not inside them.
            direction = math.atan2(y2 - y1, x2 - x1)
            r1 = symbols.shape_reach(_circular_sex(family, p1), half, direction)
            r2 = symbols.shape_reach(_circular_sex(family, p2), half, direction + math.pi)
            ux, uy = math.cos(direction), math.sin(direction)
            start, end = (x1 + r1 * ux, y1 + r1 * uy), (x2 - r2 * ux, y2 - r2 * uy)
        edges.append(
            EdgeLayout(
                start=start,
                end=end,
                edge_type=EdgeType.RELATIONSHIP,
                line_style=line_style,
                has_descent=bool(children),
            )
        )

        if not children:
            continue

        mid_x = (x1 + x2) / 2
        mid_y = (y1 + y2) / 2

        # Single child: one line straight from the couple midpoint to the
        # child, with no sibship arc to hang it from. Typed DESCENT to match
        # what the hierarchical layout emits for the same structure, so
        # callers filtering on EdgeType see parent-child links in both.
        if len(children) == 1:
            edges.append(
                EdgeLayout(start=(mid_x, mid_y), end=xy[children[0]], edge_type=EdgeType.DESCENT)
            )
            continue

        sib_r = max(radii[cid] for cid in children) - sib_offset

        # Twin groups hang from one point of the arc (their mean angle);
        # the other children each from their own angle.
        in_twin_group: set[str] = set()
        groups: list[tuple[list[str], TwinType]] = []
        for group, twin_type in twin_groups or []:
            members = [cid for cid in children if cid in group]
            if len(members) >= 2:
                in_twin_group.update(members)
                groups.append((members, twin_type))
        group_angles = [sum(angles[m] for m in members) / len(members) for members, _ in groups]
        connections = sorted(
            [angles[cid] for cid in children if cid not in in_twin_group] + group_angles
        )
        min_a, max_a = connections[0], connections[-1]
        mid_a = (min_a + max_a) / 2

        # 2. DESCENT: couple midpoint → middle of the sibship arc
        edges.append(
            EdgeLayout(
                start=(mid_x, mid_y),
                end=(sib_r * math.cos(mid_a), sib_r * math.sin(mid_a)),
                edge_type=EdgeType.DESCENT,
            )
        )

        # 3. SIBSHIP: arc joining the connection points. A union whose only
        # children are one twin group has a single point, and no arc.
        if len(connections) > 1:
            arc_wps = _arc_waypoints(min_a, max_a, sib_r, n_points=max(8, len(children) * 4))
            if len(arc_wps) >= 2:
                edges.append(
                    EdgeLayout(
                        start=arc_wps[0],
                        end=arc_wps[-1],
                        edge_type=EdgeType.SIBSHIP,
                        waypoints=arc_wps[1:-1],
                    )
                )

        # 4. Twins: one V per twin group, hanging from the arc
        for (members, twin_type), a in zip(groups, group_angles, strict=True):
            twin_connectors.append(
                TwinConnector(
                    convergence_point=(sib_r * math.cos(a), sib_r * math.sin(a)),
                    child_positions=[
                        (
                            (radii[m] - half) * math.cos(angles[m]),
                            (radii[m] - half) * math.sin(angles[m]),
                        )
                        for m in members
                    ],
                    twin_type=twin_type,
                )
            )

        # 5. INDIVIDUAL: short radial line from sibship arc to each child
        for cid in children:
            if cid in in_twin_group:
                continue
            ca = angles[cid]
            edges.append(
                EdgeLayout(
                    start=(sib_r * math.cos(ca), sib_r * math.sin(ca)),
                    end=xy[cid],
                    edge_type=EdgeType.INDIVIDUAL,
                    line_style=_individual_line_style(family, cid),
                )
            )

    return edges, twin_connectors


def _arc_waypoints(
    angle_start: float,
    angle_end: float,
    radius: float,
    n_points: int = 8,
    inset: float = 0.0,
) -> list[tuple[float, float]]:
    """Generate waypoints along a circular arc.

    Args:
        angle_start: Start angle in radians.
        angle_end: End angle in radians.
        radius: Arc radius.
        n_points: Number of intermediate points.
        inset: Shrink arc at each end by this distance (shape half-size).

    Returns:
        List of (x, y) waypoints along the arc.
    """
    if radius <= 0:
        return []

    # Inset in angular terms
    inset_angle = inset / radius if radius > 0 else 0.0

    a0 = min(angle_start, angle_end) + inset_angle
    a1 = max(angle_start, angle_end) - inset_angle

    if a1 <= a0:
        # Arc too small after inset
        mid = (angle_start + angle_end) / 2
        return [(radius * math.cos(mid), radius * math.sin(mid))]

    waypoints = []
    for i in range(n_points):
        t = i / max(n_points - 1, 1)
        a = a0 + t * (a1 - a0)
        waypoints.append((radius * math.cos(a), radius * math.sin(a)))

    return waypoints
