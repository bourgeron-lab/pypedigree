"""Interactive SVG renderer for pedigree diagrams.

Generates standalone SVG files with embedded JavaScript for:
- Tooltips on hover (sample_id, sex, phenotype, DOB, extra phenotypes)
- Click to highlight connected relatives (parents, children, siblings)
- Click on background to reset
- Mouse wheel zoom and drag to pan

Uses the same FamilyLayout pipeline as the matplotlib renderer, and the
same phenotype and carrier decisions (`pypedigree.phenotypes`,
`pypedigree.carriers`), but writes SVG directly, without matplotlib: colors
are read by `pypedigree.colors`.
"""

from __future__ import annotations

import html
import json
import math
import re
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

from pypedigree import symbols
from pypedigree.carriers import (
    CARRIER,
    LEGACY_CARRIER_KEY,
    OBLIGATE_CARRIER,
    CarrierScheme,
    carrier_columns,
    carrier_state,
    check_notation,
    resolve_carrier_scheme,
)
from pypedigree.config import (
    DEFAULT_PGS_COLOR,
    PedigreeStyle,
    apply_outline_spacing,
    apply_pgs_spacing,
    check_text_scale,
    is_curve_column,
    resolve_style,
)
from pypedigree.layout import compute_layout
from pypedigree.models import (
    DuplicateArc,
    EdgeLayout,
    Family,
    FamilyLayout,
    Individual,
    LineStyle,
    NodeLayout,
    Pedigree,
    Phenotype,
    PhenotypeDisplayMode,
    Pregnancy,
    Sex,
    TwinConnector,
    TwinType,
)
from pypedigree.phenotypes import (
    FALLBACK_COLOR,
    KEY_ROW_KINDS,
    LEG_CAPTION_COLOR,
    LEG_CAPTION_GAP,
    LEG_CATEGORY_COLOR,
    LEG_RING_INSET_PT,
    LEG_RING_WIDTH_PT,
    LEG_SEX_SYMBOL_GAP,
    LEG_SWATCH_EDGE,
    LEG_UNIT_X_EM,
    LEG_UNIT_Y_EM,
    LEG_X0,
    NA_FILL,
    NA_TEXT_COLOR,
    OUTLINE_WIDTH_SCALE,
    PGS_DOT_PT,
    PGS_MARKER_LINE_PT,
    PGS_VALUE_GAP_EM,
    PGS_VALUE_HALO_PT,
    PGS_VALUE_LINE_EM,
    QUADRANT_FONT_SCALE,
    SEPARATOR_COLOR,
    UNKNOWN_FILL,
    LegendRow,
    PhenotypeRender,
    active_curves,
    all_bands_na,
    band_cell,
    band_room,
    bin_text,
    binned_row_geometry,
    box_tail_positions,
    cell_texts,
    check_band_text,
    check_curve_options,
    check_gradient_bins,
    check_legend,
    check_min_font_size,
    check_pgs_scale,
    check_phenotype_layout,
    curve_spacing,
    curve_value_text,
    display_name,
    fit_text,
    kind_columns,
    legend_binary_captions,
    legend_gradient_width,
    legend_layout,
    legend_min_right,
    legend_row_centers,
    legend_rows,
    legend_sex_items,
    legend_slot_count,
    letter_box,
    outline_colors,
    outline_columns,
    outline_reach,
    outline_ring_offsets,
    pgs_box_left,
    pgs_caption_above,
    pgs_caption_offset,
    pgs_curve,
    pgs_curve_base,
    pgs_curve_x,
    pgs_dot_size,
    pgs_icon_geometry,
    pgs_marker_values,
    pgs_row_layout,
    pgs_values_width,
    pgs_values_x,
    pregnancy_loss_fill,
    quadrant_box,
    quadrant_cell,
    quadrant_center,
    quadrant_key_geometry,
    quadrant_room,
    resolve_phenotype_render,
    spread_dots,
    spread_labels,
    text_color_on,
    to_hex,
    validate_style_colors,
    warn_asymptomatic_sections,
    warn_pgs_scale,
    with_outline_rings,
    with_quadrants,
)
from pypedigree.phenotypes import box_width as pgs_box_w
from pypedigree.transmission import with_transmission


def render_interactive_svg(
    pedigree: Pedigree,
    style: PedigreeStyle | None = None,
    output: str | Path | None = None,
    show_labels: bool = True,
    show_legend: bool = True,
    *,
    id_prefix: str = "",
    responsive: bool = False,
) -> str:
    """Render pedigree as interactive SVG with tooltips and click highlighting.

    The drawing matches `render_pedigree`: phenotype bands, outline colors,
    PGS curves, annotations, carrier notation and the legend all come from
    the same per-figure decisions (`pypedigree.phenotypes`,
    `pypedigree.carriers`). Several families are drawn side by side.

    The SVG can be opened on its own or inlined in an HTML page next to
    other SVGs: its root has the class ``pypedigree``, its CSS only applies
    under that class, and its script only acts on its own root. Opened on
    its own, the mouse wheel zooms; inlined in a page, only Ctrl + wheel
    (or a trackpad pinch) zooms, so the page still scrolls.

    Args:
        pedigree: The parsed pedigree data.
        style: Visual style parameters. Uses defaults if None.
        output: If provided, write the SVG to this path.
        show_labels: Whether to draw sample_id labels under shapes.
        show_legend: Whether to draw the phenotype and carrier legend.
        id_prefix: Prepended to every ``id`` in the SVG (and to the
            ``url(#...)`` and ``href="#..."`` references to them), so that
            several figures can share one page. Letters, digits, ``_``,
            ``-`` and ``.``, not starting with a digit, ``-`` or ``.``.
            Empty by default: ids are then ``main-group``, ``tooltip``...
        responsive: Write ``width="100%"`` and no height instead of a size
            in pixels: the figure then fills the width of its container,
            with the proportions of its ``viewBox``.

    Returns:
        The SVG content as a string.
    """
    return _render(pedigree, style, output, show_labels, show_legend, id_prefix, responsive, True)


def render_svg(
    pedigree: Pedigree,
    style: PedigreeStyle | None = None,
    output: str | Path | None = None,
    show_labels: bool = True,
    show_legend: bool = True,
    *,
    id_prefix: str = "",
    responsive: bool = False,
) -> str:
    """Render pedigree as a static SVG: no script, no tooltip.

    The same drawing as `render_interactive_svg`, for print and for hosts
    that add their own behaviour. Each individual keeps its ``<title>``
    (shown as a native tooltip by browsers) and its ``data-*`` attributes;
    the ``dimmed`` class still fades a node or an edge.

    Args:
        pedigree: The parsed pedigree data.
        style: Visual style parameters. Uses defaults if None.
        output: If provided, write the SVG to this path.
        show_labels: Whether to draw sample_id labels under shapes.
        show_legend: Whether to draw the phenotype and carrier legend.
        id_prefix: Prepended to every ``id`` in the SVG, as in
            `render_interactive_svg`.
        responsive: Write ``width="100%"`` and no height, as in
            `render_interactive_svg`.

    Returns:
        The SVG content as a string.
    """
    return _render(pedigree, style, output, show_labels, show_legend, id_prefix, responsive, False)


def _render(
    pedigree: Pedigree,
    style: PedigreeStyle | None,
    output: str | Path | None,
    show_labels: bool,
    show_legend: bool,
    id_prefix: str,
    responsive: bool,
    interactive: bool,
) -> str:
    check_id_prefix(id_prefix)
    wrap = _SvgWrap(id_prefix=id_prefix, responsive=responsive, interactive=interactive)
    if style is None:
        style = PedigreeStyle()
    validate_style_colors(style)
    # The children's pTDT deviations, from their parents' PGS (#52).
    pedigree, _computed = with_transmission(pedigree, style)
    names = list(pedigree.phenotype_names or [])
    modes = pedigree.phenotype_display_modes
    # Same style resolution as render_pedigree, so both outputs share colors
    # and spacing. Each renderer applies the PGS spacing exactly once.
    style = resolve_style(style, names, modes)
    check_curve_options(style)
    pgs_cols = [n for n in names if is_curve_column(n, style)]
    warn_pgs_scale(pedigree.column_scope.individuals(), pgs_cols, style)
    # The only place the PGS spacing is applied (see `apply_pgs_spacing`):
    # room for the most kinds of curve one individual has (R6).
    style = apply_pgs_spacing(
        style, names, **curve_spacing(pedigree.individuals(), pgs_cols, style)
    )
    check_notation(style.carrier_notation)
    families = list(pedigree.families.values())
    # A subset is drawn as in the whole file (`Pedigree.column_scope`).
    scope = pedigree.column_scope
    carriers = resolve_carrier_scheme(
        scope.individuals(), carrier_columns(names, modes), style.carrier_notation
    )
    render = resolve_phenotype_render(scope, names, modes, style)
    warn_asymptomatic_sections(pedigree.individuals(), carriers)
    check_phenotype_layout(style)
    check_pgs_scale(style)
    check_min_font_size(style)
    check_text_scale(style)
    check_band_text(style)
    check_legend(style)
    check_gradient_bins(style)
    render = with_quadrants(render, render.band_names(names), carriers, style)
    # Several o_ borders (#46): room around every symbol, applied once here.
    render = with_outline_rings(render, pedigree.individuals(), names)
    style = apply_outline_spacing(style, outline_reach(render.outline_rings, style))
    # Bands need colors; PGS curves need none (see renderer._render_layout).
    banded = bool(
        names and render.effective_modes and (style.phenotype_colors or render.curve_names(names))
    )
    ctx = _SvgContext(
        style=style,
        carriers=carriers,
        render=render,
        display_modes=modes,
        band_names=render.band_names(names) if banded else [],
        outline_names=render.outline_names(names) if banded else [],
        pgs_names=render.curve_names(names) if banded else [],
        banded=banded,
        show_labels=show_labels,
    )

    if not families:
        return _write(_wrap_svg("", (0.0, 0.0, 100.0, 100.0), style, wrap), output)

    parts: list[str] = []
    if carriers.keys:
        parts.append(_svg_carrier_patterns(carriers, style))

    extents: list[_Extent] = []
    cursor = 0.0
    titled = len(families) > 1
    for index, family in enumerate(families):
        layout = _flip_vertical(compute_layout(family, style, labels=show_labels))
        x0, y0, x1, y1 = _family_extent(layout, ctx)
        if index > 0:
            # Side by side, left to right, like render_pedigree's subplots.
            dx = cursor - x0

            def shift(x: float, y: float, dx: float = dx) -> tuple[float, float]:
                return (x + dx, y)

            layout = _map_layout(layout, shift)
            x0, x1 = x0 + dx, x1 + dx
        if titled:
            title_fs = (style.text_size() + 2) * symbols.units_per_point(style.text_scale)
            parts.append(
                _svg_text(
                    (x0 + x1) / 2,
                    y0 - title_fs * 0.5,
                    f"Family: {family.family_id}",
                    title_fs,
                    style,
                    anchor="middle",
                    baseline="auto",
                    cls="family-title",
                )
            )
            y0 -= title_fs * 2
        parts.append(_render_family_svg(layout, family, _build_adjacency(family), ctx, f"f{index}"))
        extents.append((x0, y0, x1, y1))
        cursor = x1 + style.shape_size * 3

    bounds = _union(extents)
    show_pheno_legend = (
        show_legend and names and (style.phenotype_colors or render.curve_names(names))
    )
    show_carrier_legend = show_legend and (carriers.has_legacy_carriers or bool(carriers.keys))
    if show_pheno_legend or show_carrier_legend:
        rows = legend_rows(names, style, render, modes, carriers)
        if rows:
            legend, legend_extent = _svg_legend(
                rows, style, carriers, bounds[0], bounds[3] + style.shape_size / 2
            )
            parts.append(legend)
            bounds = _union([bounds, legend_extent])

    return _write(_wrap_svg("\n".join(parts), bounds, style, wrap), output)


# Font sizes and line widths, given in points in PedigreeStyle, are scaled
# to layout units by `symbols.units_per_point(style.text_scale)` (#145):
# `symbols.UNITS_PER_POINT` at the default `text_scale`.

_Extent = tuple[float, float, float, float]  # x0, y0, x1, y1 (y down)


@dataclass(frozen=True)
class _SvgContext:
    """Per-figure decisions every node of the SVG is drawn from."""

    style: PedigreeStyle
    carriers: CarrierScheme
    render: PhenotypeRender
    display_modes: Mapping[str, PhenotypeDisplayMode]
    band_names: list[str]
    outline_names: list[str]
    pgs_names: list[str]
    banded: bool
    show_labels: bool


def _write(svg: str, output: str | Path | None) -> str:
    if output is not None:
        Path(output).write_text(svg, encoding="utf-8")
    return svg


def _union(extents: Sequence[_Extent]) -> _Extent:
    return (
        min(e[0] for e in extents),
        min(e[1] for e in extents),
        max(e[2] for e in extents),
        max(e[3] for e in extents),
    )


def _map_layout(
    layout: FamilyLayout, f: Callable[[float, float], tuple[float, float]]
) -> FamilyLayout:
    """Apply a point transformation to every coordinate of a layout."""

    def p(pt: tuple[float, float]) -> tuple[float, float]:
        return f(pt[0], pt[1])

    nodes = []
    for n in layout.nodes:
        x, y = f(n.x, n.y)
        nodes.append(replace(n, x=x, y=y))
    return replace(
        layout,
        nodes=nodes,
        edges=[
            replace(e, start=p(e.start), end=p(e.end), waypoints=[p(w) for w in e.waypoints])
            for e in layout.edges
        ],
        twin_connectors=[
            replace(
                tc,
                convergence_point=p(tc.convergence_point),
                child_positions=[p(c) for c in tc.child_positions],
            )
            for tc in layout.twin_connectors
        ],
        duplicate_arcs=[
            replace(
                a,
                from_pos=p(a.from_pos),
                to_pos=p(a.to_pos),
                waypoints=[p(w) for w in a.waypoints],
            )
            for a in layout.duplicate_arcs
        ],
    )


def _flip_vertical(layout: FamilyLayout) -> FamilyLayout:
    """Mirror a layout top to bottom for SVG coordinates.

    The layout puts generation 0 at y = 0 and descendants at negative y
    (matplotlib's axis points up); SVG's y axis points down, so without
    this the founders were drawn at the bottom. `_svg_duplicate_arc` bows
    to the opposite side of matplotlib's arc3 for the same curvature, which
    the mirroring compensates exactly.
    """
    return _map_layout(layout, lambda x, y: (x, -y))


def _family_extent(layout: FamilyLayout, ctx: _SvgContext) -> _Extent:
    """Bounding box of everything drawn for one family."""
    pt = symbols.units_per_point(ctx.style.text_scale)  # layout units per point (#145)
    style = ctx.style
    half = style.shape_size / 2
    ring = outline_reach(ctx.render.outline_rings, style)  # several o_ borders (#46)
    reach = half * 1.2 + ring  # the diamond is the widest shape
    xs = [n.x for n in layout.nodes]
    ys = [n.y for n in layout.nodes]
    for e in layout.edges:
        for x, y in (e.start, e.end, *e.waypoints):
            xs.append(x)
            ys.append(y)
    if not xs:
        return (0.0, 0.0, 1.0, 1.0)
    x0, x1 = min(xs) - reach, max(xs) + reach
    y0, y1 = min(ys) - reach, max(ys) + reach
    if ctx.pgs_names:
        # Each PGS box, left of its symbol or right of it (#78), and its
        # values beyond it.
        values = pgs_values_width(layout.nodes, ctx.pgs_names, style)
        box_w = pgs_box_w(style)
        for n in layout.nodes:
            for index in range(len(active_curves(n, ctx.pgs_names, style))):
                box_left = pgs_box_left(n, style, _node_ring(n, ctx), index)
                x0 = min(x0, box_left - values)
                x1 = max(x1, box_left + box_w + values)
    if ctx.show_labels and layout.nodes:
        label_bottom = max(
            n.y
            + half
            + style.label_offset
            + n.marker_depth
            + style.text_size() * pt * 1.2
            + (_label_line_offsets(n, style) or [0.0])[-1]
            for n in layout.nodes
        )
        y1 = max(y1, label_bottom)
        # A label wider than its symbol reaches beyond it on both sides.
        for n in layout.nodes:
            if n.label and n.label_angle is None:
                half_label = (
                    symbols.label_width(
                        n.label,
                        n.notes + n.annotations,
                        style.text_size(),
                        style.text_size(symbols.SMALL_TEXT_SCALE),
                    )
                    / 2
                )
                x0 = min(x0, n.x - half_label)
                x1 = max(x1, n.x + half_label)
    return (x0, y0, x1, y1)


def _build_adjacency(family: Family) -> dict[str, dict[str, list[str]]]:
    """Build adjacency data for each individual (parents, children, siblings).

    Returns:
        Dict mapping sample_id -> {"parents": [...], "children": [...], "siblings": [...]}.
    """
    adj: dict[str, dict[str, list[str]]] = {}

    for sid in family.individuals:
        adj[sid] = {"parents": [], "children": [], "siblings": []}

    for ind in family.individuals.values():
        if ind.father_id and ind.father_id in family.individuals:
            adj[ind.sample_id]["parents"].append(ind.father_id)
        if ind.mother_id and ind.mother_id in family.individuals:
            adj[ind.sample_id]["parents"].append(ind.mother_id)

    for union in family.unions:
        p1, p2 = union.partner_ids
        for cid in union.children_ids:
            if cid in adj:
                if p1 in family.individuals:
                    adj[p1]["children"].append(cid)
                if p2 in family.individuals:
                    adj[p2]["children"].append(cid)
        # Siblings: all children in same union
        for i, c1 in enumerate(union.children_ids):
            for c2 in union.children_ids[i + 1 :]:
                if c1 in adj and c2 in adj:
                    adj[c1]["siblings"].append(c2)
                    adj[c2]["siblings"].append(c1)

    return adj


def _render_family_svg(
    layout: FamilyLayout,
    family: Family,
    adjacency: dict[str, dict[str, list[str]]],
    ctx: _SvgContext,
    prefix: str,
) -> str:
    """Render one family's layout as SVG elements.

    Args:
        prefix: Makes the element ids of this family unique in the document
            (sample ids may repeat across families, and need not be valid
            XML ids).
    """
    pt = symbols.units_per_point(ctx.style.text_scale)  # layout units per point (#145)
    parts: list[str] = []
    style = ctx.style

    for edge in layout.edges:
        parts.append(_svg_edge(edge, style))
    for tc in layout.twin_connectors:
        parts.append(_svg_twin_connector(tc, style))
    for arc in layout.duplicate_arcs:
        parts.append(_svg_duplicate_arc(arc, style))

    empty: dict[str, list[str]] = {"parents": [], "children": [], "siblings": []}
    for i, node in enumerate(layout.nodes):
        if node.is_unknown_offspring:
            # A placeholder child: a "?", with no tooltip or highlighting.
            qx, qy = symbols.unknown_offspring_anchor(node.x, node.y, style.shape_size / 2, up=-1.0)
            parts.append(
                _svg_text(
                    qx,
                    qy,
                    symbols.UNKNOWN_OFFSPRING_TEXT,
                    style.text_size(symbols.UNKNOWN_OFFSPRING_FONT_SCALE) * pt,
                    style,
                    color=style.outline_color,
                    baseline="hanging",
                    cls="unknown-offspring",
                )
            )
            continue
        ind = family.individuals.get(node.individual_id)
        adj = adjacency.get(node.individual_id, empty)
        parts.append(_svg_node(node, ind, adj, ctx, f"{prefix}n{i}"))

    return "\n".join(parts)


# SVG drawing of each carrier pattern: (rotation in degrees, lines, dots).
# Lines are "v" (vertical) and/or "h" (horizontal) before rotation.
_SVG_PATTERNS: dict[str, tuple[int, str, str]] = {
    "diagonal": (45, "v", ""),
    "antidiagonal": (-45, "v", ""),
    "vertical": (0, "v", ""),
    "horizontal": (0, "h", ""),
    "crossed": (45, "vh", ""),
    "grid": (0, "vh", ""),
    "dotted": (0, "", "dot"),
    "circles": (0, "", "ring"),
}


def _svg_pattern_id(pattern: str, obligate: bool) -> str:
    return f"carrier-{pattern}-{'sparse' if obligate else 'dense'}"


def _svg_carrier_patterns(carriers: CarrierScheme, style: PedigreeStyle) -> str:
    """``<pattern>`` definitions for every carrier source of the figure."""
    out = ["<defs>"]
    stroke = html.escape(style.fill_color)
    for key in carriers.keys:
        pattern = carriers.pattern(key)
        for obligate in (False, True):
            if obligate and key not in carriers.obligate_keys:
                continue
            size = style.shape_size * (0.25 if obligate else 0.125)
            angle, lines, dots = _SVG_PATTERNS[pattern]
            body = ""
            paint = f'stroke="{stroke}" stroke-width="0.02"'
            centre = f'cx="{size / 2:.3f}" cy="{size / 2:.3f}"'
            if "v" in lines:
                body += f'<line x1="0" y1="0" x2="0" y2="{size:.3f}" {paint}/>'
            if "h" in lines:
                body += f'<line x1="0" y1="0" x2="{size:.3f}" y2="0" {paint}/>'
            if dots == "dot":
                body += f'<circle {centre} r="{size / 8:.3f}" fill="{stroke}"/>'
            elif dots == "ring":
                body += (
                    f'<circle cx="{size / 2:.3f}" cy="{size / 2:.3f}" r="{size / 4:.3f}" '
                    f'fill="none" stroke="{stroke}" stroke-width="0.015"/>'
                )
            out.append(
                f'<pattern id="{_svg_pattern_id(pattern, obligate)}" patternUnits="userSpaceOnUse" '
                f'width="{size:.3f}" height="{size:.3f}" patternTransform="rotate({angle})">'
                f"{body}</pattern>"
            )
    out.append("</defs>")
    return "".join(out)


# Root class of every SVG: the CSS below only applies under it, and hosts
# can select pypedigree figures with it.
ROOT_CLASS = "pypedigree"

_ID_PREFIX_RE = re.compile(r"[A-Za-z_][A-Za-z0-9_.-]*")


def check_id_prefix(id_prefix: str) -> None:
    """Raise ValueError unless ``id_prefix`` can start an XML id.

    Empty, or a letter or ``_`` followed by letters, digits, ``_``, ``-``
    or ``.``: such ids need no escaping in CSS attribute selectors or in
    ``url(#...)`` references.
    """
    if id_prefix and not _ID_PREFIX_RE.fullmatch(id_prefix):
        raise ValueError(
            f"id_prefix {id_prefix!r} is not a valid id start: use letters, digits, "
            "'_', '-' or '.', starting with a letter or '_'."
        )


# Every id the renderer writes, and every reference to one, is in an
# attribute. Text and titles are escaped (html.escape turns '"' into
# '&quot;'), so these patterns never match inside them.
_ID_ATTR_RE = re.compile(r'(\sid=")')
_URL_REF_RE = re.compile(r'(="url\(#)')
_HREF_REF_RE = re.compile(r'(\shref="#)')


def _prefix_ids(markup: str, id_prefix: str) -> str:
    """Prepend ``id_prefix`` to every id in ``markup`` and every reference to one."""
    if not id_prefix:
        return markup
    for pattern in (_ID_ATTR_RE, _URL_REF_RE, _HREF_REF_RE):
        markup = pattern.sub(lambda m: m.group(1) + id_prefix, markup)
    return markup


@dataclass(frozen=True)
class _SvgWrap:
    """How the drawing is wrapped: ids, size, and whether it is interactive."""

    id_prefix: str = ""
    responsive: bool = False
    interactive: bool = True


_PX_PER_UNIT = 60.0
_MAX_PX_W = 1600
_MAX_PX_H = 1200
_MIN_PX_W = 600

# Scoped under the root class, so an inlined figure styles nothing else on
# the page.
_STATIC_CSS = f"""
.{ROOT_CLASS} .node, .{ROOT_CLASS} .edge {{ transition: opacity 0.3s; }}
.{ROOT_CLASS} .node.dimmed {{ opacity: 0.15; }}
.{ROOT_CLASS} .edge.dimmed {{ opacity: 0.1; }}"""

_INTERACTIVE_CSS = f"""
.{ROOT_CLASS} .node {{ cursor: pointer; }}
.{ROOT_CLASS} .node:hover .shape {{ filter: brightness(0.9); }}
.{ROOT_CLASS} .tooltip-box {{ pointer-events: none; }}
.{ROOT_CLASS} text {{ user-select: none; }}"""


def _wrap_svg(
    content: str, bounds: _Extent, style: PedigreeStyle, wrap: _SvgWrap | None = None
) -> str:
    """Wrap SVG content in a complete SVG document, interactive or not."""
    if wrap is None:
        wrap = _SvgWrap()
    min_x, min_y, max_x, max_y = bounds
    pad = style.margin
    vb_x = min_x - pad
    vb_y = min_y - pad
    vb_w = max_x - min_x + 2 * pad
    vb_h = max_y - min_y + 2 * pad

    # Default display size: 60 px per layout unit, scaled down to fit in
    # 1600 x 1200 and up to be at least 600 px wide (unless that would be
    # taller than 1200). One scale for both axes, so the size has the
    # drawing's proportions; clamping each axis on its own letterboxed wide
    # pedigrees with bands of white above and below.
    if wrap.responsive:
        size = 'width="100%"'
    else:
        fit = min(_MAX_PX_W / vb_w, _MAX_PX_H / vb_h)
        scale = min(max(_PX_PER_UNIT, _MIN_PX_W / vb_w), fit)
        size = f'width="{round(vb_w * scale)}" height="{round(vb_h * scale)}"'

    css = _STATIC_CSS + (_INTERACTIVE_CSS if wrap.interactive else "")
    tooltip = (
        """
<!-- Tooltip -->
<g id="tooltip" class="tooltip-box" visibility="hidden">
  <rect id="tooltip-bg" rx="0.06" ry="0.06" fill="rgba(0,0,0,0.85)" />
  <text id="tooltip-text" fill="white" font-size="0.28" />
</g>"""
        if wrap.interactive
        else ""
    )
    body = f"""<defs>
<marker id="arrowhead" viewBox="0 0 10 10" refX="9" refY="5"
        markerWidth="5" markerHeight="5" orient="auto-start-reverse">
<path d="M0,0 L10,5 L0,10 z" fill="{style.line_color}"/></marker>
<style>{css}
</style>
</defs>
<g id="main-group">
{content}
</g>{tooltip}"""
    body = _prefix_ids(body, wrap.id_prefix)
    script = _interactive_js(wrap.id_prefix) if wrap.interactive else ""

    return f"""<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" class="{ROOT_CLASS}"
     viewBox="{vb_x:.2f} {vb_y:.2f} {vb_w:.2f} {vb_h:.2f}"
     {size}
     style="background: white; font-family: {style.font_family};">
{body}
{script}
</svg>"""


# Shift from the requested anchor to the alphabetic baseline, in font sizes
# (cap height of a typical sans-serif font is about 0.7 em).
_BASELINE_SHIFT = {"hanging": 0.75, "central": 0.35, "auto": 0.0}


def _svg_text(
    x: float,
    y: float,
    text: str,
    size: float,
    style: PedigreeStyle,
    *,
    anchor: str = "middle",
    baseline: str = "central",
    color: str = "black",
    weight: str = "",
    italic: bool = False,
    cls: str = "",
    halo: float = 0.0,
) -> str:
    """One `<text>` element; `size` is in layout units.

    `baseline` says where `y` falls on the text: "hanging" (top),
    "central" (middle) or "auto" (baseline). It is turned into an explicit
    offset rather than `dominant-baseline`, which several SVG renderers
    (librsvg, older Inkscape) do not honour. ``halo`` (layout units) draws a
    white outline under the glyphs.
    """
    y += size * _BASELINE_SHIFT[baseline]
    extra = ""
    if weight:
        extra += f' font-weight="{weight}"'
    if italic:
        extra += ' font-style="italic"'
    if halo:
        extra += (
            f' stroke="white" stroke-width="{halo:.3f}" stroke-linejoin="round"'
            ' paint-order="stroke"'
        )
    if cls:
        extra += f' class="{cls}"'
    return (
        f'<text x="{x:.3f}" y="{y:.3f}" font-size="{size:.3f}" '
        f'text-anchor="{anchor}" fill="{html.escape(color)}"{extra}>{html.escape(text)}</text>'
    )


_GENDER_LABELS = {"man": "Man", "woman": "Woman", "nonbinary": "Non-binary"}
_PREGNANCY_LABELS = {
    "P": "Pregnancy",
    "SAB": "Miscarriage (SAB)",
    "TOP": "Termination of pregnancy (TOP)",
    "ECT": "Ectopic pregnancy (ECT)",
    "SB": "Stillbirth (SB)",
}
_CHILDLESS_LABELS = {"choice": "No children by choice", "infertility": "Infertility"}
_ROLE_LABELS = {"D": "Gamete donor", "S": "Surrogate"}


def _tooltip(
    node: NodeLayout,
    ind: Individual | None,
    carriers: CarrierScheme,
    style: PedigreeStyle,
    outline_names: Sequence[str] = (),
) -> str:
    """Text of a node's tooltip, lines separated by " | "."""
    lines = [f"ID: {node.individual_id}"]
    if node.number:
        lines.append(f"Number: {node.number}")  # NSGC numbering (#125)
    sex = ind.sex if ind is not None else node.sex
    lines.append(f"Sex: {({Sex.MALE: 'Male', Sex.FEMALE: 'Female'}).get(sex, 'Unknown')}")
    notation = node.notation
    if notation.gender is not None:
        lines.append(f"Gender: {_GENDER_LABELS[notation.gender.value]}")
    pheno_label = {
        Phenotype.AFFECTED: "Affected",
        Phenotype.UNAFFECTED: "Unaffected",
    }.get(node.phenotype, "Unknown")
    lines.append(f"Phenotype: {pheno_label}")

    if node.is_carrier:
        lines.append("Status: Carrier")
    if node.is_obligate_carrier:
        lines.append("Status: Obligate carrier")
    for key in carriers.keys:
        state = carrier_state(node, key)
        if key != LEGACY_CARRIER_KEY and state in (CARRIER, OBLIGATE_CARRIER):
            shown = display_name(key, style)
            lines.append(carriers.label(key, state == OBLIGATE_CARRIER, shown))
        elif state is None and key in node.not_applicable:
            lines.append(f"{display_name(key, style)}: not applicable")  # #121

    if notation.pregnancy is not None:
        lines.append(_PREGNANCY_LABELS[notation.pregnancy.value])
    if notation.adoption is not None:
        lines.append(f"Adopted {notation.adoption.value}")
    if notation.childlessness is not None:
        label = _CHILDLESS_LABELS[notation.childlessness.value]
        if notation.childless_partner:
            label += f" (with {notation.childless_partner})"  # a couple (#104)
        lines.append(label)
    if notation.reproductive_role is not None:
        lines.append(_ROLE_LABELS[notation.reproductive_role.value])
    if notation.count:
        lines.append(f"{notation.count} individuals")  # #104
    if notation.is_asymptomatic:
        lines.append("Asymptomatic or presymptomatic carrier")
    if notation.is_consultand:
        lines.append("Consultand")
    if node.is_inferred:
        lines.append("Not sampled (inferred)")

    if ind is not None and ind.date_of_birth:
        lines.append(f"DOB: {ind.date_of_birth}")

    borders = outline_columns(node, outline_names)
    if len(borders) > 1:
        # Which border is which, inner first (#46).
        lines.append("Borders: " + ", ".join(display_name(n, style) for n in borders))

    raw = node.extra_phenotypes_raw
    categorical = node.extra_phenotypes_categorical
    for name in dict.fromkeys([*raw, *categorical]):
        if name in carriers.keys:
            continue
        category = categorical.get(name)
        value = raw.get(name)
        shown = display_name(name, style)
        if name in node.not_applicable:
            lines.append(f"{shown}: not applicable")  # #107
        elif category is not None:
            lines.append(f"{shown}: {category}")
        elif value is not None and math.isnan(value):
            lines.append(f"{shown}: NA")
        elif value is not None and is_curve_column(name, style):
            lines.append(f"{shown}: {curve_value_text(name, value, style)}")
        elif value is not None:
            # With clinical thresholds, the step it falls in (#57).
            step = bin_text(name, value, style)
            lines.append(f"{shown}: {value:g}" + (f" ({step})" if step else ""))
    return " | ".join(lines)


def _svg_node(
    node: NodeLayout,
    ind: Individual | None,
    adj: dict[str, list[str]],
    ctx: _SvgContext,
    uid: str,
) -> str:
    """Generate SVG for one node with data attributes for interactivity.

    Mirrors `renderer._render_layout`: a node with phenotype data is drawn
    with bands (and its PGS curve) when the figure has any, otherwise with
    its affection status and carrier sections.
    """
    pt = symbols.units_per_point(ctx.style.text_scale)  # layout units per point (#145)
    style = ctx.style
    half = style.shape_size / 2
    # data-x/data-y: the symbol's centre and data-size its size, in the
    # coordinates of #main-group (the viewBox's user units before zoom and
    # pan), so a host can draw an overlay in the same system. The host's own
    # classes and data-* attributes (css_class / data_ columns) follow.
    classes = " ".join(["node", *node.css_classes])
    extra = "".join(
        f' data-{name}="{html.escape(value)}"' for name, value in node.data_attributes.items()
    )
    number = f' data-number="{html.escape(node.number)}"' if node.number else ""
    parts = [
        f'<g class="{html.escape(classes)}" data-id="{html.escape(node.individual_id)}" '
        f'data-x="{node.x + 0.0:.3f}" data-y="{node.y + 0.0:.3f}" '
        f'data-size="{style.shape_size:.3f}"{extra}{number} '
        f'data-parents="{html.escape(",".join(adj.get("parents", [])))}" '
        f'data-children="{html.escape(",".join(adj.get("children", [])))}" '
        f'data-siblings="{html.escape(",".join(adj.get("siblings", [])))}">'
    ]

    banded = ctx.banded and bool(node.extra_phenotypes) and not node.is_phantom
    triangle = symbols.is_pregnancy_loss(node.notation)
    if triangle:
        names = ctx.band_names if banded else []
        fill = pregnancy_loss_fill(node, names, style, ctx.render, ctx.display_modes)
        parts.append(_svg_pregnancy_loss(node, style, fill))
    elif banded and ctx.band_names:
        parts.append(_svg_banded_body(node, ctx, uid))
    else:
        parts.append(_svg_plain_body(node, ctx, uid))

    if not triangle:
        parts.append(_svg_legacy_carrier(node, ctx, uid))
    parts.append(_svg_notation(node, ctx, banded and bool(ctx.band_names)))

    sw = style.shape_line_width * pt
    ring = _node_ring(node, ctx) if banded else 0.0  # around several o_ borders (#46)
    if node.is_deceased and not triangle:
        o = (half + ring) * symbols.DECEASED_OVERSHOOT + ring
        parts.append(
            f'<line x1="{node.x - half - o:.3f}" y1="{node.y + half + o:.3f}" '
            f'x2="{node.x + half + o:.3f}" y2="{node.y - half - o:.3f}" '
            f'stroke="{style.outline_color}" stroke-width="{sw:.3f}"/>'
        )

    if node.is_proband or node.notation.is_consultand:
        if triangle:
            tip_x, tip_y = symbols.triangle_points(node.x, node.y, half, up=-1.0)[1]
        else:
            tip_x, tip_y = node.x - half - ring, node.y + half + ring
        size = style.text_size() * pt
        ax, ay = symbols.marker_points(tip_x, tip_y, size, up=-1.0)
        if node.is_proband:
            parts.append(
                _svg_text(
                    ax - symbols.MARKER_LETTER_GAP_EM * size,
                    ay,
                    "P",
                    size,
                    style,
                    anchor="end",
                    color=style.line_color,
                    weight="bold",
                )
            )
        parts.append(
            f'<line x1="{ax:.3f}" y1="{ay:.3f}" '
            f'x2="{tip_x:.3f}" y2="{tip_y:.3f}" '
            f'stroke="{style.line_color}" stroke-width="{sw:.3f}" '
            f'marker-end="url(#arrowhead)"/>'
        )

    if banded and ctx.pgs_names:
        parts.append(_svg_pgs(node, ctx))

    if ctx.show_labels and node.label:
        parts.append(_svg_labels(node, style))

    tip = _tooltip(node, ind, ctx.carriers, style, ctx.outline_names)
    parts.append(f"<title>{html.escape(tip)}</title>")
    parts.append("</g>")
    return "\n".join(p for p in parts if p)


def _svg_pregnancy_loss(node: NodeLayout, style: PedigreeStyle, fill: str) -> str:
    """Triangle of a miscarriage, ectopic pregnancy or termination.

    Mirrors `renderer._draw_pregnancy_loss`: filled as decided by
    ``phenotypes.pregnancy_loss_fill``, struck through for a termination.
    """
    pt = symbols.units_per_point(style.text_scale)  # layout units per point (#145)
    half = style.shape_size / 2
    points = symbols.triangle_points(node.x, node.y, half, up=-1.0)
    sw = style.shape_line_width * pt
    out = [
        f'<polygon class="shape" points="{" ".join(f"{x:.3f},{y:.3f}" for x, y in points)}" '
        f'fill="{html.escape(fill)}" stroke="{html.escape(style.outline_color)}" '
        f'stroke-width="{sw:.3f}"/>'
    ]
    if node.notation.pregnancy == Pregnancy.TERMINATION:
        (x1, y1), (x2, y2) = symbols.termination_stroke(node.x, node.y, half, up=-1.0)
        out.append(
            f'<line class="termination" x1="{x1:.3f}" y1="{y1:.3f}" x2="{x2:.3f}" '
            f'y2="{y2:.3f}" stroke="{html.escape(style.outline_color)}" '
            f'stroke-width="{style.line_width * pt:.3f}"/>'
        )
    return "\n".join(out)


def _svg_notation(node: NodeLayout, ctx: _SvgContext, banded: bool) -> str:
    """Adoption brackets, the documented-evaluation asterisk, and the letter (P, D or S).

    Mirrors `renderer._draw_notation`.
    """
    pt = symbols.units_per_point(ctx.style.text_scale)  # layout units per point (#145)
    style, carriers = ctx.style, ctx.carriers
    half = style.shape_size / 2
    ring = _node_ring(node, ctx) if banded else 0.0
    if ring:
        # Round the outer border of several (#46).
        half = symbols.grown_half(node.sex, half, ring)
    triangle = symbols.is_pregnancy_loss(node.notation)
    out = []
    if node.notation.adoption is not None:
        for line in symbols.bracket_lines(node.x, node.y, half, node.sex, triangle):
            pts = " ".join(f"{x:.3f},{y:.3f}" for x, y in line)
            out.append(
                f'<polyline class="adoption-bracket" points="{pts}" fill="none" '
                f'stroke="{html.escape(style.line_color)}" '
                f'stroke-width="{style.line_width * pt:.3f}"/>'
            )
    if node.notation.is_evaluated:
        mark_pt = style.text_size(symbols.EVALUATION_MARK_FONT_SCALE)
        width = style.shape_line_width
        if banded and outline_colors(node, ctx.outline_names, style):
            width *= OUTLINE_WIDTH_SCALE
        brackets = style.line_width if node.notation.adoption is not None else None
        ex, ey = symbols.evaluation_mark_position(
            node.x,
            node.y,
            half,
            node.sex,
            mark_pt,
            width,
            pt,
            triangle,
            up=-1.0,
            bracket_width=brackets,
        )
        size = mark_pt * pt
        out.append(
            _svg_text(
                ex,
                # "*" sits high in its line: put its ink's middle on ey (#110).
                ey + symbols.EVALUATION_MARK_INK_MIDDLE_EM * size,
                symbols.EVALUATION_MARK,
                size,
                style,
                baseline="auto",
                color=style.line_color,
                weight="bold",
                cls="documented-evaluation",
            )
        )
    if node.notation.is_asymptomatic and not triangle:
        # A vertical line through the symbol, past its top and bottom (#104).
        (x1, y1), (x2, y2) = symbols.asymptomatic_line(
            node.x, node.y, style.shape_size / 2, node.sex
        )
        out.append(
            f'<line class="asymptomatic" x1="{x1:.3f}" y1="{y1:.3f}" x2="{x2:.3f}" '
            f'y2="{y2:.3f}" stroke="{html.escape(style.outline_color)}" '
            f'stroke-width="{style.shape_line_width * pt:.3f}"/>'
        )
    letter = symbols.inner_letter(node.notation)
    if letter is not None and not triangle:
        # At the centre, on a white backing over bands, quadrants or
        # carrier sections (#105).
        divided = (
            banded
            or bool(carriers.marks(node))
            or (carriers.legacy and (node.is_carrier or node.is_obligate_carrier))
        )
        color = style.line_color
        size = style.text_size() * pt
        if divided:
            # Wider for a number of two digits or more (#104).
            half_w = size * (0.55 + symbols.LETTER_CLEAR_CHAR_EM * (len(letter) - 1))
            out.append(
                f'<rect class="symbol-letter-backing" x="{node.x - half_w:.3f}" '
                f'y="{node.y - size * 0.6:.3f}" width="{2 * half_w:.3f}" '
                f'height="{size * 1.2:.3f}" rx="{size * 0.25:.3f}" fill="white" '
                f'fill-opacity="0.9"/>'
            )
        elif node.phenotype == Phenotype.AFFECTED:
            color = text_color_on(style.fill_color)
        out.append(_svg_text(node.x, node.y, letter, size, style, color=color, cls="symbol-letter"))
    return "\n".join(out)


def _svg_plain_body(node: NodeLayout, ctx: _SvgContext, uid: str) -> str:
    """Shape filled by affection status, divided into carrier sections.

    With the NSGC 2022 notation, someone with carrier results gets a symbol
    divided into vertical sections (see pypedigree.carriers).
    """
    pt = symbols.units_per_point(ctx.style.text_scale)  # layout units per point (#145)
    style, carriers = ctx.style, ctx.carriers
    half = style.shape_size / 2
    legacy_carrier = carriers.legacy and (node.is_carrier or node.is_obligate_carrier)
    affected = node.phenotype == Phenotype.AFFECTED and not legacy_carrier
    marks = carriers.marks(node)
    sections: list[str] = []
    if marks:
        if affected:
            sections.append(html.escape(style.fill_color))
        sections.extend(f"url(#{_svg_pattern_id(mk.pattern, mk.obligate)})" for mk in marks)
    fill = style.fill_color if affected and not sections else "white"

    stroke = style.outline_color
    sw = style.shape_line_width * pt
    dash = 'stroke-dasharray="0.1 0.05"' if node.is_phantom or node.is_inferred else ""
    out = [
        f'{_svg_shape_outline(node, half, cls="shape")} fill="{html.escape(fill)}" '
        f'stroke="{html.escape(stroke)}" stroke-width="{sw:.3f}" {dash}/>'
    ]
    if sections:
        out.append(_svg_sections(node, uid, half, sections, stroke, sw, dash))
    return "\n".join(out)


def _svg_banded_body(node: NodeLayout, ctx: _SvgContext, uid: str) -> str:
    """Horizontal phenotype bands, then one band per carrier source.

    Mirrors `renderer._draw_banded_node`, from the same `band_cell`
    decisions.
    """
    pt = symbols.units_per_point(ctx.style.text_scale)  # layout units per point (#145)
    style, carriers = ctx.style, ctx.carriers
    half = style.shape_size / 2
    cx, cy = node.x, node.y
    keys = list(carriers.keys)
    n_bands = len(ctx.band_names) + len(keys)
    band_h = 2 * half / n_bands
    top = cy - half
    left, width = cx - half * 1.5, half * 3
    clip = f"bands-{uid}"
    clip_attr = f'clip-path="url(#{clip})"'
    sep = f'stroke="{SEPARATOR_COLOR}" stroke-width="{0.5 * pt:.3f}" {clip_attr}'
    band_fs = style.text_size(symbols.BAND_TEXT_SCALE)

    # The letter at the centre, which band text keeps clear of (#105).
    letter = letter_box(node.notation, style.text_size())

    def fitted(texts: tuple[str, ...], yc: float, size: float, color: str) -> str | None:
        """The text of a band, fitted to its width at the SVG's scale (#55)."""
        room = band_room(node.sex, half, yc - cy)
        fit = fit_text(
            texts,
            size,
            style,
            pt,
            room,
            band_h if n_bands else 2 * half,
            letter=letter,
            offset=(0.0, yc - cy),
        )
        if fit is None:
            return None
        x = cx + fit.dx
        return _svg_text(x, yc, fit.text, fit.size * pt, style, color=color, weight="bold")

    def separator(y: float) -> str:
        return f'<line x1="{left:.3f}" y1="{y:.3f}" x2="{left + width:.3f}" y2="{y:.3f}" {sep}/>'

    out = [
        f'<defs><clipPath id="{clip}">{_svg_shape_outline(node, half)}/></clipPath></defs>',
        f'{_svg_shape_outline(node, half)} fill="white" stroke="none"/>',
    ]
    band_names = ctx.band_names
    if ctx.render.quadrants:
        quadrants, all_na = _svg_quadrants(node, ctx, clip_attr)
        out.extend(quadrants)
        band_names, keys = [], []
    else:
        all_na = not keys and all_bands_na(node, band_names)

    for i, name in enumerate(band_names):
        cell = band_cell(node, name, style, ctx.render, ctx.display_modes)
        y = top + i * band_h
        out.append(
            f'<rect class="band" x="{left:.3f}" y="{y:.3f}" width="{width:.3f}" '
            f'height="{band_h:.3f}" fill="{to_hex(cell.color)}" {clip_attr}/>'
        )
        if i > 0:
            out.append(
                f'<line x1="{left:.3f}" y1="{y:.3f}" x2="{left + width:.3f}" y2="{y:.3f}" {sep}/>'
            )
        if cell.text is not None and not (cell.is_na and all_na):
            texts = cell_texts(cell, style.band_text)
            if (svg := fitted(texts, y + band_h / 2, band_fs, cell.text_color)) is not None:
                out.append(svg)

    for j, key in enumerate(keys):
        y = top + (len(band_names) + j) * band_h
        state = carrier_state(node, key)
        not_applicable = state is None and key in node.not_applicable  # #121
        geometry = f'x="{left:.3f}" y="{y:.3f}" width="{width:.3f}" height="{band_h:.3f}"'
        fill = NA_FILL if state is None and not not_applicable else "white"
        out.append(f'<rect class="band" {geometry} fill="{fill}" {clip_attr}/>')
        if state in (CARRIER, OBLIGATE_CARRIER):
            pid = _svg_pattern_id(carriers.pattern(key), state == OBLIGATE_CARRIER)
            out.append(f'<rect {geometry} fill="url(#{pid})" {clip_attr}/>')
        out.append(
            f'<line x1="{left:.3f}" y1="{y:.3f}" x2="{left + width:.3f}" y2="{y:.3f}" {sep}/>'
        )
        na_band = state is None and not not_applicable and style.show_na_labels
        if na_band and (svg := fitted(("NA",), y + band_h / 2, band_fs, NA_TEXT_COLOR)):
            out.append(svg)

    if all_na and style.show_na_labels:
        na = fit_text(
            ("NA",),
            style.text_size(symbols.NA_TEXT_SCALE),
            style,
            pt,
            band_room(node.sex, half, 0.0),
            2 * half,
            letter=letter,
        )
        if na is not None:
            size = na.size * pt
            na_x = cx + na.dx
            out.append(
                _svg_text(na_x, cy, na.text, size, style, color=NA_TEXT_COLOR, weight="bold")
            )

    borders = outline_colors(node, ctx.outline_names, style)
    border = borders[0] if borders else None
    sw = style.shape_line_width * pt * (OUTLINE_WIDTH_SCALE if border else 1.0)
    dash = 'stroke-dasharray="0.1 0.05"' if node.is_phantom or node.is_inferred else ""
    out.append(
        f'{_svg_shape_outline(node, half, cls="shape")} fill="none" '
        f'stroke="{html.escape(border or style.outline_color)}" stroke-width="{sw:.3f}" {dash}/>'
    )
    # The other borders around it, outwards, in column order (#46).
    offsets = outline_ring_offsets(len(borders), style)
    for color, offset in zip(borders[1:], offsets, strict=True):
        grown = symbols.grown_half(node.sex, half, offset * pt)
        out.append(
            f'{_svg_shape_outline(node, grown, cls="outline-ring")} fill="none" '
            f'stroke="{html.escape(color)}" stroke-width="{sw:.3f}" {dash}/>'
        )
    return "\n".join(out)


def _node_ring(node: NodeLayout, ctx: _SvgContext) -> float:
    """Room the borders after the first take around this symbol (#46)."""
    return outline_reach(len(outline_colors(node, ctx.outline_names, ctx.style)), ctx.style)


def _svg_quadrants(node: NodeLayout, ctx: _SvgContext, clip_attr: str) -> tuple[list[str], bool]:
    """The symbol cut in four, one quadrant per phenotype or carrier source.

    Mirrors `renderer._draw_quadrants`, from the same `quadrant_cell` and
    `quadrant_box` decisions (with the y axis pointing down).

    Returns:
        The elements, and whether every quadrant has no data.
    """
    pt = symbols.units_per_point(ctx.style.text_scale)  # layout units per point (#145)
    style, carriers = ctx.style, ctx.carriers
    half = style.shape_size / 2
    cx, cy = node.x, node.y
    columns = [n for n in ctx.render.quadrants if n not in carriers.keys]
    has_carriers = any(n in carriers.keys for n in ctx.render.quadrants)
    all_na = not has_carriers and all_bands_na(node, columns)
    fs = style.text_size(QUADRANT_FONT_SCALE)
    # The letter at the centre, which quadrant text keeps clear of (#105).
    letter = letter_box(node.notation, style.text_size())
    out: list[str] = []
    for name, position in ctx.render.quadrants.items():
        x, y, w, h = quadrant_box(position, cx, cy, half, up=-1.0)
        geometry = f'x="{x:.3f}" y="{y:.3f}" width="{w:.3f}" height="{h:.3f}"'
        texts: tuple[str, ...] = ()
        color = NA_TEXT_COLOR
        if name in carriers.keys:
            state = carrier_state(node, name)
            no_data = state is None and name not in node.not_applicable  # N/A: #121
            fill = NA_FILL if no_data else "white"
            out.append(f'<rect class="band" {geometry} fill="{fill}" {clip_attr}/>')
            if state in (CARRIER, OBLIGATE_CARRIER):
                pid = _svg_pattern_id(carriers.pattern(name), state == OBLIGATE_CARRIER)
                out.append(f'<rect {geometry} fill="url(#{pid})" {clip_attr}/>')
            if no_data and style.show_na_labels:
                texts = ("NA",)
        else:
            cell = quadrant_cell(node, name, style, ctx.render, ctx.display_modes)
            out.append(f'<rect class="band" {geometry} fill="{to_hex(cell.color)}" {clip_attr}/>')
            if cell.text is not None and not (cell.is_na and all_na):
                texts, color = cell_texts(cell, style.band_text), cell.text_color
        room = quadrant_room(node.sex, half, position)
        tx, ty = quadrant_center(position, cx, cy, half, up=-1.0)
        fit = fit_text(texts, fs, style, pt, room, half, letter=letter, offset=(tx - cx, ty - cy))
        if fit is not None:
            size = fit.size * pt
            out.append(_svg_text(tx, ty, fit.text, size, style, color=color, weight="bold"))
    return out, all_na


def _svg_legacy_carrier(node: NodeLayout, ctx: _SvgContext, uid: str) -> str:
    """Legacy notation: left half filled for a carrier, central dot for an obligate one."""
    style = ctx.style
    if not ctx.carriers.legacy:
        return ""
    half = style.shape_size / 2
    out = []
    if node.is_carrier:
        clip = f"half-{uid}"
        out.append(
            f'<defs><clipPath id="{clip}">{_svg_shape_outline(node, half)}/></clipPath></defs>'
        )
        out.append(
            f'<rect class="carrier-half-fill" x="{node.x - half * 1.5:.3f}" '
            f'y="{node.y - half * 1.5:.3f}" width="{half * 1.5:.3f}" height="{half * 3:.3f}" '
            f'fill="{html.escape(style.fill_color)}" clip-path="url(#{clip})"/>'
        )
    if node.is_obligate_carrier:
        out.append(
            f'<circle class="obligate-carrier-dot" cx="{node.x:.3f}" cy="{node.y:.3f}" '
            f'r="{style.carrier_dot_radius:.3f}" fill="{html.escape(style.fill_color)}"/>'
        )
    return "\n".join(out)


def _svg_pgs(node: NodeLayout, ctx: _SvgContext) -> str:
    """An individual's curve boxes beside its symbol, one per kind of curve (R6).

    Mirrors `renderer._draw_pgs_ridge`: each curve runs bottom to top, with
    each score's position on [0, 1] marked by a thin solid line and a dot (#158).
    """
    style = ctx.style
    columns = kind_columns(ctx.pgs_names, style)
    out = []
    for index, (kind, markers) in enumerate(active_curves(node, ctx.pgs_names, style)):
        tails = box_tail_positions(columns[kind], style)
        out.append(_svg_curve_box(node, ctx, index, kind, markers, tails))
    return "\n".join(out)


def _svg_curve_box(
    node: NodeLayout,
    ctx: _SvgContext,
    index: int,
    caption: str,
    active: list[tuple[str, float]],
    tails: tuple[float, ...],
) -> str:
    """One box of `_svg_pgs`, ``index`` boxes out from the symbol."""
    pt = symbols.units_per_point(ctx.style.text_scale)  # layout units per point (#145)
    style = ctx.style
    half = style.shape_size / 2
    box_w, box_h = pgs_box_w(style), half * 2
    # Clear of several o_ borders (#46); right of a right-hand partner (#78).
    left = pgs_box_left(node, style, _node_ring(node, ctx), index)
    bottom = node.y + half
    thin = 0.5 * pt

    # Peak to the left of a baseline near the right edge (#126).
    curve = [(pgs_curve_x(left, box_w, d), bottom - t * box_h) for t, d in pgs_curve(50)]
    base_x = pgs_curve_base(left, box_w)
    pts = " ".join(f"{x:.3f},{y:.3f}" for x, y in curve)
    area = f"{base_x:.3f},{bottom:.3f} {pts} {base_x:.3f},{bottom - box_h:.3f}"
    out = [
        '<g class="pgs">',
        f'<rect x="{left:.3f}" y="{bottom - box_h:.3f}" width="{box_w:.3f}" height="{box_h:.3f}" '
        f'rx="0.02" fill="white" fill-opacity="0.5" stroke="{SEPARATOR_COLOR}" '
        f'stroke-width="{thin:.3f}"/>',
        f'<polygon points="{area}" fill="{FALLBACK_COLOR}" fill-opacity="0.08"/>',
        f'<polyline points="{pts}" fill="none" stroke="#aaaaaa" stroke-width="{0.6 * pt:.3f}"/>',
    ]
    for tail in tails:
        ty = bottom - tail * box_h
        out.append(
            f'<line class="pgs-tail" x1="{left:.3f}" y1="{ty:.3f}" x2="{left + box_w:.3f}" '
            f'y2="{ty:.3f}" stroke="#aaaaaa" stroke-width="{thin:.3f}" '
            'stroke-dasharray="0.02 0.03"/>'
        )
    ys = [bottom - value * box_h for _name, value in active]
    # Dots closer than their diameter side by side, so none hides another (#54).
    # Within the box, so no dot reaches the symbol (#89), and smaller when the
    # box cannot hold them all side by side (#147).
    dot = pgs_dot_size(len(active), box_w, pt)
    offsets = spread_dots(ys, dot, box_w)
    for (name, _value), y, dx in zip(active, ys, offsets, strict=True):
        color = html.escape(style.pgs_colors.get(name, DEFAULT_PGS_COLOR))
        out.append(
            f'<line class="pgs-marker" x1="{left:.3f}" y1="{y:.3f}" x2="{left + box_w:.3f}" '
            f'y2="{y:.3f}" stroke="{color}" stroke-width="{PGS_MARKER_LINE_PT * pt:.3f}"/>'
        )
        cx = left + box_w / 2 + dx
        out.append(f'<circle cx="{cx:.3f}" cy="{y:.3f}" r="{dot / 2:.3f}" fill="{color}"/>')
    if style.show_pgs_values:
        out.extend(_svg_pgs_values(node, active, ys, left, bottom - box_h, bottom, style))
    above = pgs_caption_above(node)
    caption_size = style.text_size(symbols.PGS_CAPTION_SCALE)
    gap = pgs_caption_offset(half, caption_size, pt)
    out.append(
        _svg_text(
            left + box_w / 2,
            bottom - box_h - gap if above else bottom + gap,
            caption,
            caption_size * pt,
            style,
            baseline="auto" if above else "hanging",
            color="#666666",
        )
    )
    out.append("</g>")
    return "\n".join(out)


def _svg_pgs_values(
    node: NodeLayout,
    active: list[tuple[str, float]],
    ys: list[float],
    left: float,
    top: float,
    bottom: float,
    style: PedigreeStyle,
) -> list[str]:
    """Each PGS value left of the curve's box, at its marker's height (#54)."""
    pt = symbols.units_per_point(style.text_scale)  # layout units per point (#145)
    texts = pgs_marker_values(node, [name for name, _ in active], style)
    size = style.text_size(symbols.PGS_CAPTION_SCALE)
    step = PGS_VALUE_LINE_EM * size * pt
    placed = spread_labels(ys, step, top + step / 2, bottom - step / 2)
    box_w = pgs_box_w(style)
    # The far side of the box from the symbol (#78).
    x, side = pgs_values_x(node, left, box_w, PGS_VALUE_GAP_EM * size * pt)
    anchor = "start" if side == "left" else "end"
    out = []
    for (name, _value), y in zip(active, placed, strict=True):
        if name in texts:
            color = style.pgs_colors.get(name, DEFAULT_PGS_COLOR)
            out.append(
                _svg_text(
                    x,
                    y,
                    texts[name],
                    size * pt,
                    style,
                    anchor=anchor,
                    color=color,
                    cls="pgs-value",
                    halo=PGS_VALUE_HALO_PT * pt,
                )
            )
    return out


def _label_line_offsets(node: NodeLayout, style: PedigreeStyle) -> list[float]:
    """Where each note and annotation starts below the label's top (#66)."""
    pt = symbols.units_per_point(style.text_scale)  # layout units per point (#145)
    small = style.text_size(symbols.SMALL_TEXT_SCALE)
    lines = len(node.notes) + len(node.annotations)
    return symbols.below_label_offsets(style.text_size(), [small] * lines, pt)


def _svg_labels(node: NodeLayout, style: PedigreeStyle) -> str:
    """The ID label under a symbol, its notes ("SB", "AMAB"...), then its annotations."""
    pt = symbols.units_per_point(style.text_scale)  # layout units per point (#145)
    y = node.y + style.shape_size / 2 + style.label_offset + node.marker_depth
    out = [
        _svg_text(
            node.x,
            y,
            node.label,
            style.text_size() * pt,
            style,
            baseline="hanging",
            color=style.line_color,
        )
    ]
    below = [(text, style.line_color, False, "note") for text in node.notes]
    below += [(text, "#666666", True, "annotation") for text in node.annotations]
    for (text, color, italic, cls), offset in zip(
        below, _label_line_offsets(node, style), strict=True
    ):
        out.append(
            _svg_text(
                node.x,
                y + offset,
                text,
                style.text_size(symbols.SMALL_TEXT_SCALE) * pt,
                style,
                baseline="hanging",
                color=color,
                italic=italic,
                cls=cls,
            )
        )
    return "\n".join(out)


def _svg_shape_outline(node: NodeLayout, half: float, cls: str = "") -> str:
    """Geometry of an individual's shape, without paint attributes.

    The element is left open, so callers append attributes and ``/>``.
    """
    attr = f' class="{cls}"' if cls else ""
    if node.sex == Sex.MALE:
        return (
            f'<rect{attr} x="{node.x - half:.3f}" y="{node.y - half:.3f}" '
            f'width="{half * 2:.3f}" height="{half * 2:.3f}"'
        )
    if node.sex == Sex.FEMALE:
        return f'<circle{attr} cx="{node.x:.3f}" cy="{node.y:.3f}" r="{half:.3f}"'
    r = half * 1.2
    points = (
        f"{node.x:.3f},{node.y - r:.3f} {node.x + r:.3f},{node.y:.3f} "
        f"{node.x:.3f},{node.y + r:.3f} {node.x - r:.3f},{node.y:.3f}"
    )
    return f'<polygon{attr} points="{points}"'


def _svg_sections(
    node: NodeLayout,
    uid: str,
    half: float,
    fills: list[str],
    stroke: str,
    sw: float,
    dash: str,
) -> str:
    """Equal vertical sections clipped to the shape, then the outline on top."""
    clip_id = f"sections-{uid}"
    reach = half * 1.2 if node.sex not in (Sex.MALE, Sex.FEMALE) else half
    width = 2 * reach / len(fills)
    out = [f'<defs><clipPath id="{clip_id}">{_svg_shape_outline(node, half)}/></clipPath></defs>']
    for i, fill in enumerate(fills):
        left = node.x - reach + i * width
        out.append(
            f'<rect x="{left:.3f}" y="{node.y - reach:.3f}" width="{width:.3f}" '
            f'height="{reach * 2:.3f}" fill="{fill}" clip-path="url(#{clip_id})"/>'
        )
        if i > 0:
            out.append(
                f'<line x1="{left:.3f}" y1="{node.y - reach:.3f}" x2="{left:.3f}" '
                f'y2="{node.y + reach:.3f}" stroke="{stroke}" stroke-width="{sw / 2:.3f}" '
                f'clip-path="url(#{clip_id})"/>'
            )
    out.append(
        f'{_svg_shape_outline(node, half)} fill="none" stroke="{stroke}" '
        f'stroke-width="{sw:.3f}" {dash}/>'
    )
    return "".join(out)


# Average advance width of a character, in font sizes. Python cannot measure
# text, so the legend border is sized from this estimate.
_CHAR_WIDTH_EM = 0.55


def _svg_legend(
    rows: Sequence[LegendRow],
    style: PedigreeStyle,
    carriers: CarrierScheme,
    x0: float,
    y0: float,
) -> tuple[str, _Extent]:
    """The legend panel, its top-left corner at (x0, y0).

    Same rows (`legend_rows`) and same geometry (`legend_metrics`, in units
    of the font size) as the matplotlib legend.

    Returns:
        The SVG group and its bounding box.
    """
    pt = symbols.units_per_point(style.text_scale)  # layout units per point (#145)
    binary_headers = any(r.kind in ("binary", "binary_shared") for r in rows)
    rows, m = legend_layout(rows, style)
    ux = LEG_UNIT_X_EM * m.font * pt
    uy = LEG_UNIT_Y_EM * m.font * pt
    fs = m.fs * pt
    body: list[str] = []
    right = legend_min_right(style)  # rightmost content, in legend units (#79)

    def X(u: float) -> float:
        return x0 + u * ux

    def Y(v: float) -> float:
        # Legend units have y pointing up from the bottom of the panel.
        return y0 + (m.total_h - v) * uy

    def text(u: float, v: float, s: str, size: float, anchor: str = "start", **kw: Any) -> None:
        nonlocal right
        width = len(s) * _CHAR_WIDTH_EM * size / ux
        end = {"start": u + width, "middle": u + width / 2}.get(anchor, u)
        right = max(right, end)
        body.append(_svg_text(X(u), Y(v), s, size, style, anchor=anchor, **kw))

    def swatch(
        u: float,
        vc: float,
        w: float,
        h: float,
        fill: str,
        stroke: str = LEG_SWATCH_EDGE,
        width_pt: float = 0.5,
    ) -> None:
        stroke_attr = (
            f'stroke="{html.escape(stroke)}" stroke-width="{width_pt * pt:.3f}"'
            if stroke != "none"
            else 'stroke="none"'
        )
        body.append(
            f'<rect x="{X(u):.3f}" y="{Y(vc + h / 2):.3f}" width="{w * ux:.3f}" '
            f'height="{h * uy:.3f}" fill="{html.escape(fill)}" {stroke_attr}/>'
        )

    def caption(u: float, v: float, s: str, anchor: str = "start", **kw: Any) -> None:
        kw.setdefault("color", LEG_CAPTION_COLOR)
        size = m.caption_fs(kw.pop("size_factor", 1.0)) * pt
        text(u, v, s, size, anchor, baseline="hanging", **kw)

    # Header: title and sex symbols, each if the style has it (#67).
    if m.title:
        text(
            LEG_X0,
            m.title_y,
            m.title,
            m.title_fs * pt,
            baseline="hanging",
            weight="bold",
        )
    if m.sex_y is not None:
        v = m.sex_y
        side = m.marker_size * pt
        for u, shape, label in legend_sex_items(style, m):
            if shape == "square":
                body.append(
                    f'<rect x="{X(u) - side / 2:.3f}" y="{Y(v) - side / 2:.3f}" '
                    f'width="{side:.3f}" height="{side:.3f}" fill="white" stroke="black" '
                    f'stroke-width="{pt:.3f}"/>'
                )
            else:
                body.append(
                    f'<circle cx="{X(u):.3f}" cy="{Y(v):.3f}" r="{side / 2:.3f}" fill="white" '
                    f'stroke="black" stroke-width="{pt:.3f}"/>'
                )
            text(u + LEG_SEX_SYMBOL_GAP, v, label, fs)
    y_top = m.content_top

    # The quadrant key first (and its note, #67), so the captions sit above
    # the swatch rows.
    keys = [r for r in rows if r.kind in KEY_ROW_KINDS]
    for row, vc in zip(keys, legend_row_centers(keys, m, y_top), strict=True):
        if row.kind == "key_note":
            size = m.caption_fs() * pt
            text(LEG_X0, vc, row.label, size, baseline="central", color=LEG_CATEGORY_COLOR)
            continue
        g = quadrant_key_geometry(row, m, vc)
        sx, sy, sw, sh = g.square
        square = (
            f'<rect x="{X(sx):.3f}" y="{Y(sy + sh / 2):.3f}" width="{sw * ux:.3f}" '
            f'height="{sh * uy:.3f}"'
        )
        shapes = [("legend-key-square", square, (sx, sy + sh / 2))]
        if g.circle is not None:
            ccx, ccy, rx, _ry = g.circle
            circle = f'<circle cx="{X(ccx):.3f}" cy="{Y(ccy):.3f}" r="{rx * ux:.3f}"'
            shapes.append(("legend-key-circle", circle, (ccx - rx, ccy + sh / 2)))
        for shape_id, outline, (u0, v_top) in shapes:
            body.append(f'<defs><clipPath id="{shape_id}">{outline}/></clipPath></defs>')
            for cell in row.quadrants:
                qx = X(u0 + (sw / 2 if cell.position.endswith("right") else 0.0))
                qy = Y(v_top - (0.0 if cell.position.startswith("top") else sh / 2))
                geometry = (
                    f'x="{qx:.3f}" y="{qy:.3f}" width="{sw / 2 * ux:.3f}" '
                    f'height="{sh / 2 * uy:.3f}" clip-path="url(#{shape_id})"'
                )
                body.append(f'<rect {geometry} fill="{html.escape(cell.color)}"/>')
                if cell.pattern:
                    pid = _svg_pattern_id(cell.pattern, False)
                    body.append(f'<rect {geometry} fill="url(#{pid})"/>')
            body.append(f'{outline} fill="none" stroke="black" stroke-width="{pt:.3f}"/>')
        for u, v, label, anchor in g.labels:
            text(u, v, label, fs, anchor)
    y_top -= legend_slot_count(keys) * m.row_h
    rows = [r for r in rows if r.kind not in KEY_ROW_KINDS]

    if binary_headers:
        for i, s in enumerate(legend_binary_captions(style)):
            caption(m.swatch_x(i) + m.sw_w / 2, y_top, s, "middle")
        y_top -= m.binary_header_h

    for row, vc in zip(rows, legend_row_centers(rows, m, y_top), strict=True):
        if row.kind == "gradient":
            bar_w = legend_gradient_width(m)  # as wide as the other bars (#95)
            step = bar_w / len(row.colors)
            for k, rgb in enumerate(row.colors):
                swatch(LEG_X0 + k * step, vc, step * 1.02, m.sw_h, to_hex(rgb), stroke="none")
            swatch(LEG_X0, vc, bar_w, m.sw_h, "none")
            cap_v = vc - m.sw_h / 2 - LEG_CAPTION_GAP
            caption(LEG_X0, cap_v, row.low)
            caption(LEG_X0 + bar_w, cap_v, row.high, "end")
        elif row.kind == "binned":
            # A gradient's steps (`gradient_bins`, #57), each with its limits.
            # The steps' names under their limits (#81).
            steps = binned_row_geometry(row, m)
            names = row.bin_names or ("",) * len(row.bins)
            for k, (rgb, limits, step_name) in enumerate(
                zip(row.colors, row.bins, names, strict=True)
            ):
                u = LEG_X0 + k * steps.step_w
                swatch(u, vc, steps.step_w * 1.01, steps.bar_h, to_hex(rgb), stroke="none")
                v = vc - steps.bar_h / 2 - LEG_CAPTION_GAP
                caption(u + steps.step_w / 2, v, limits, "middle")
                if step_name and row.names_inside:
                    # In the step, in a color that reads on it (#90).
                    text(
                        u + steps.step_w / 2,
                        vc,
                        step_name,
                        steps.name_size * pt,
                        "middle",
                        baseline="central",
                        color=text_color_on(to_hex(rgb)),
                    )
                elif step_name:
                    caption(
                        u + steps.step_w / 2,
                        v - steps.name_dy,
                        step_name,
                        "middle",
                        color=LEG_CATEGORY_COLOR,
                    )
            swatch(LEG_X0, vc, steps.step_w * len(row.colors), steps.bar_h, "none")
        elif row.kind == "categorical":
            u = LEG_X0
            for value, color in row.categories.items():
                swatch(u, vc, m.sw_w * 0.8, m.sw_h, color)
                caption(
                    u + m.sw_w * 0.4,
                    vc - m.sw_h / 2 - LEG_CAPTION_GAP,
                    value,
                    "middle",
                    color=LEG_CATEGORY_COLOR,
                    size_factor=0.8,
                )
                u += m.sw_w + m.sw_gap * 0.5
        elif row.kind == "binary":
            for k, fill in enumerate((row.color, "white", UNKNOWN_FILL)):
                swatch(m.swatch_x(k), vc, m.sw_w, m.sw_h, fill)
        elif row.kind == "binary_shared":
            # The "affected" swatch cut into one stripe per column of the key.
            stripe = m.sw_w / len(row.colors)
            for k, rgb in enumerate(row.colors):
                swatch(m.swatch_x(0) + k * stripe, vc, stripe, m.sw_h, to_hex(rgb), stroke="none")
            swatch(m.swatch_x(0), vc, m.sw_w, m.sw_h, "none")
            for k, fill in ((1, "white"), (2, UNKNOWN_FILL)):
                swatch(m.swatch_x(k), vc, m.sw_w, m.sw_h, fill)
        elif row.kind == "outline":
            swatch(LEG_X0, vc, m.sw_w, m.sw_h, "white", stroke=row.color, width_pt=2.5)
        elif row.kind == "asymptomatic":
            # A white swatch with the carrier's vertical line (#104).
            swatch(LEG_X0, vc, m.sw_w, m.sw_h, "white", stroke=style.outline_color, width_pt=0.5)
            reach = m.sw_h / 2 * (1 + symbols.ASYMPTOMATIC_OVERSHOOT)
            u = LEG_X0 + m.sw_w / 2
            body.append(
                f'<line class="asymptomatic" x1="{X(u):.3f}" y1="{Y(vc - reach):.3f}" '
                f'x2="{X(u):.3f}" y2="{Y(vc + reach):.3f}" '
                f'stroke="{html.escape(style.outline_color)}" stroke-width="{1.0 * pt:.3f}"/>'
            )
        elif row.kind == "outline_rings":
            # Two nested borders, inner first (#46).
            inner, outer = (to_hex(c) for c in row.colors)
            du, dv = LEG_RING_INSET_PT * pt / ux, LEG_RING_INSET_PT * pt / uy
            swatch(LEG_X0, vc, m.sw_w, m.sw_h, "white", stroke=outer, width_pt=LEG_RING_WIDTH_PT)
            swatch(
                LEG_X0 + du,
                vc,
                m.sw_w - 2 * du,
                m.sw_h - 2 * dv,
                "none",
                stroke=inner,
                width_pt=LEG_RING_WIDTH_PT,
            )
        elif row.kind == "pgs":
            # Upright, as beside a symbol in the figure (#159).
            icon = pgs_icon_geometry(row, m, vc)
            left, bottom, w, h = icon.box
            body.append(
                f'<rect class="pgs-legend-box" x="{X(left):.3f}" y="{Y(bottom + h):.3f}" '
                f'width="{w * ux:.3f}" height="{h * uy:.3f}" fill="white" '
                f'stroke="{SEPARATOR_COLOR}" stroke-width="{0.5 * pt:.3f}"/>'
            )
            pts = " ".join(f"{X(u):.3f},{Y(v):.3f}" for u, v in icon.curve)
            body.append(
                f'<polygon points="{X(icon.base_x):.3f},{Y(bottom):.3f} {pts} '
                f'{X(icon.base_x):.3f},{Y(bottom + h):.3f}" fill="{FALLBACK_COLOR}" '
                f'fill-opacity="0.08"/>'
            )
            body.append(
                f'<polyline class="pgs-legend-curve" points="{pts}" fill="none" '
                f'stroke="#aaaaaa" stroke-width="{0.6 * pt:.3f}"/>'
            )
            for v in icon.tails:
                body.append(
                    f'<line class="pgs-legend-tail" x1="{X(left):.3f}" y1="{Y(v):.3f}" '
                    f'x2="{X(left + w):.3f}" y2="{Y(v):.3f}" stroke="#aaaaaa" '
                    f'stroke-width="{0.5 * pt:.3f}" stroke-dasharray="0.02 0.03"/>'
                )
            for u, v, label, anchor in icon.labels:
                text(
                    u,
                    v,
                    label,
                    icon.size * pt,
                    anchor=anchor,
                    baseline="central",
                    color=LEG_CATEGORY_COLOR,
                )
        elif row.kind == "legacy_carrier":
            big = 9 * m.scale * pt
            body.append(
                f'<rect x="{X(0.9) - big / 2:.3f}" y="{Y(vc) - big / 2:.3f}" width="{big:.3f}" '
                f'height="{big:.3f}" fill="white" stroke="black" stroke-width="{pt:.3f}"/>'
            )
            if row.obligate:
                body.append(
                    f'<circle cx="{X(0.9):.3f}" cy="{Y(vc):.3f}" r="{1.5 * m.scale * pt:.3f}" '
                    f'fill="black"/>'
                )
            else:
                body.append(
                    f'<rect x="{X(0.9) - big / 2:.3f}" y="{Y(vc) - big / 2:.3f}" '
                    f'width="{big / 2:.3f}" height="{big:.3f}" fill="black"/>'
                )
        else:  # carrier: patterned swatch
            pid = _svg_pattern_id(row.pattern, row.obligate)
            swatch(LEG_X0, vc, m.sw_w, m.sw_h, "white", stroke=style.fill_color)
            swatch(LEG_X0, vc, m.sw_w, m.sw_h, f"url(#{pid})", stroke="none")

        swatch_rows = ("gradient", "categorical", "binary", "binary_shared")
        # One column for the names of swatch and bar rows (#88).
        label_u = m.name_x if row.kind in (*swatch_rows, "binned") else m.icon_text_x
        if row.kind != "pgs":
            text(label_u, vc, row.label, fs, baseline="central")
            continue
        # The caption, then the PGS key: each score's name after a line of its
        # color with a dot, as it is marked in the figure (#56, #158), on the
        # lines below when too long for one (#68).
        layout = pgs_row_layout(row, label_u, m)
        for k, line in enumerate(layout.captions):
            text(layout.text_x, vc - k * m.line_h, line, fs, baseline="central")
        for item in layout.items:
            v = vc - item.line * m.line_h
            color = html.escape(item.color)
            body.append(
                f'<line class="pgs-key-line" x1="{X(item.dash_x0):.3f}" y1="{Y(v):.3f}" '
                f'x2="{X(item.dash_x1):.3f}" y2="{Y(v):.3f}" stroke="{color}" '
                f'stroke-width="{PGS_MARKER_LINE_PT * pt:.3f}"/>'
            )
            body.append(
                f'<circle class="pgs-key-dot" cx="{X((item.dash_x0 + item.dash_x1) / 2):.3f}" '
                f'cy="{Y(v):.3f}" r="{PGS_DOT_PT * m.scale * pt / 2:.3f}" fill="{color}"/>'
            )
            text(item.text_x, v, item.name, fs, baseline="central")

    width = (right + LEG_X0 - 0.1) * ux
    border = (
        f'<rect x="{X(0.1):.3f}" y="{Y(m.total_h - 0.1):.3f}" width="{width:.3f}" '
        f'height="{(m.total_h - 0.2) * uy:.3f}" rx="{0.1 * uy:.3f}" fill="white" '
        f'stroke="{LEG_SWATCH_EDGE}" stroke-width="{0.8 * pt:.3f}"/>'
    )
    svg = '<g class="legend">\n' + "\n".join([border, *body]) + "\n</g>"
    return svg, (x0, y0, X(0.1) + width, y0 + m.total_h * uy)


def _svg_edge(edge: EdgeLayout, style: PedigreeStyle) -> str:
    """Generate SVG line for an edge."""
    pt = symbols.units_per_point(style.text_scale)  # layout units per point (#145)
    x1, y1 = edge.start
    x2, y2 = edge.end
    sw = style.line_width * pt

    if edge.line_style == LineStyle.DOUBLE:
        gap = style.double_line_gap
        if y1 == y2:  # Horizontal
            return (
                f'<g class="edge">'
                f'<line x1="{x1:.3f}" y1="{y1 + gap:.3f}" '
                f'x2="{x2:.3f}" y2="{y2 + gap:.3f}" '
                f'stroke="{style.line_color}" stroke-width="{sw:.3f}"/>'
                f'<line x1="{x1:.3f}" y1="{y1 - gap:.3f}" '
                f'x2="{x2:.3f}" y2="{y2 - gap:.3f}" '
                f'stroke="{style.line_color}" stroke-width="{sw:.3f}"/>'
                f"</g>"
            )
        else:
            return (
                f'<g class="edge">'
                f'<line x1="{x1 + gap:.3f}" y1="{y1:.3f}" '
                f'x2="{x2 + gap:.3f}" y2="{y2:.3f}" '
                f'stroke="{style.line_color}" stroke-width="{sw:.3f}"/>'
                f'<line x1="{x1 - gap:.3f}" y1="{y1:.3f}" '
                f'x2="{x2 - gap:.3f}" y2="{y2:.3f}" '
                f'stroke="{style.line_color}" stroke-width="{sw:.3f}"/>'
                f"</g>"
            )

    if edge.line_style == LineStyle.BROKEN:
        return _svg_separation(edge, style)

    dash = ""
    if edge.line_style == LineStyle.DASHED:
        dash = 'stroke-dasharray="0.1 0.05"'

    return (
        f'<line class="edge" x1="{x1:.3f}" y1="{y1:.3f}" '
        f'x2="{x2:.3f}" y2="{y2:.3f}" '
        f'stroke="{style.line_color}" stroke-width="{sw:.3f}" {dash}/>'
    )


def _svg_separation(edge: EdgeLayout, style: PedigreeStyle) -> str:
    """A separated couple's relationship line, stopped round its two slashes (#141).

    Mirrors `renderer._draw_separation`, with y pointing down.
    """
    pt = symbols.units_per_point(style.text_scale)  # layout units per point (#145)
    line, slashes = symbols.separation_marks(
        edge.start, edge.end, style.shape_size / 2, edge.has_descent, up=-1.0
    )
    stroke = f'stroke="{html.escape(style.line_color)}" stroke-width="{style.line_width * pt:.3f}"'
    out = ['<g class="edge">']
    for cls, segments in (("relationship", line), ("separation-mark", slashes)):
        for (x0, y0), (x1, y1) in segments:
            out.append(
                f'<line class="{cls}" x1="{x0:.3f}" y1="{y0:.3f}" '
                f'x2="{x1:.3f}" y2="{y1:.3f}" {stroke}/>'
            )
    out.append("</g>")
    return "\n".join(out)


def _svg_twin_connector(tc: TwinConnector, style: PedigreeStyle) -> str:
    """Generate SVG for twin connector lines."""
    pt = symbols.units_per_point(style.text_scale)  # layout units per point (#145)
    cx, cy = tc.convergence_point
    sw = style.line_width * pt
    parts = ['<g class="edge">']

    for child_x, child_y in tc.child_positions:
        parts.append(
            f'<line x1="{cx:.3f}" y1="{cy:.3f}" '
            f'x2="{child_x:.3f}" y2="{child_y:.3f}" '
            f'stroke="{style.line_color}" stroke-width="{sw:.3f}"/>'
        )

    if tc.twin_type == TwinType.MONOZYGOTIC and len(tc.child_positions) >= 2:
        frac = 0.5
        bar_points = []
        for child_x, child_y in tc.child_positions:
            bx = cx + frac * (child_x - cx)
            by = cy + frac * (child_y - cy)
            bar_points.append((bx, by))
        parts.append(
            f'<line x1="{bar_points[0][0]:.3f}" y1="{bar_points[0][1]:.3f}" '
            f'x2="{bar_points[-1][0]:.3f}" y2="{bar_points[-1][1]:.3f}" '
            f'stroke="{style.line_color}" stroke-width="{sw:.3f}"/>'
        )

    if tc.twin_type == TwinType.UNKNOWN and len(tc.child_positions) >= 2:
        # Mirrors `renderer._draw_twin_connector` (#139): y points down here.
        mark_pt = style.text_size(symbols.TWIN_MARK_FONT_SCALE)
        dx, dy = symbols.twin_mark_offset(
            tc.convergence_point, tc.child_positions, mark_pt, style.line_width
        )
        size = mark_pt * pt
        qy = cy + dy * pt + symbols.TWIN_MARK_INK_MIDDLE_EM * size
        parts.append(
            _svg_text(
                cx + dx * pt,
                qy,
                symbols.TWIN_MARK,
                size,
                style,
                baseline="auto",
                color=style.line_color,
                weight="bold",
            )
        )

    parts.append("</g>")
    return "\n".join(parts)


def _svg_duplicate_arc(arc: DuplicateArc, style: PedigreeStyle) -> str:
    """Generate SVG for duplicate arc (dotted curve)."""
    pt = symbols.units_per_point(style.text_scale)  # layout units per point (#145)
    x1, y1 = arc.from_pos
    x2, y2 = arc.to_pos
    sw = style.line_width * pt
    if arc.waypoints:
        # The hierarchical layout's route, above the row (#149).
        points = " ".join(f"{x:.3f},{y:.3f}" for x, y in arc.waypoints)
        return (
            f'<polyline class="duplicate-link" points="{points}" fill="none" '
            f'stroke="{style.line_color}" stroke-width="{sw:.3f}" stroke-dasharray="0.05 0.05"/>'
        )

    # Compute control point for quadratic bezier
    mid_x = (x1 + x2) / 2
    mid_y = (y1 + y2) / 2
    dx = x2 - x1
    dy = y2 - y1
    # Perpendicular offset for curvature
    ctrl_x = mid_x - dy * arc.curvature
    ctrl_y = mid_y + dx * arc.curvature

    return (
        f'<path class="edge" d="M {x1:.3f} {y1:.3f} Q {ctrl_x:.3f} {ctrl_y:.3f} '
        f'{x2:.3f} {y2:.3f}" '
        f'fill="none" stroke="{style.line_color}" stroke-width="{sw:.3f}" '
        f'stroke-dasharray="0.05 0.05"/>'
    )


_INTERACTIVE_JS = """<script>
<![CDATA[
(function() {
  const P = __ID_PREFIX__;
  // Act on our own root only: the <svg> holding this script. Fallback: the
  // last pypedigree SVG parsed so far, which is ours while the page loads.
  const me = document.currentScript;
  let svg = me && me.closest ? me.closest('svg') : null;
  if (!svg) {
    const all = document.querySelectorAll('svg.pypedigree');
    svg = all[all.length - 1];
  }
  if (!svg) return;
  const byId = function(id) { return svg.querySelector('[id="' + P + id + '"]'); };
  const mainGroup = byId('main-group');
  const tooltip = byId('tooltip');
  const tooltipBg = byId('tooltip-bg');
  const tooltipText = byId('tooltip-text');
  const nodes = svg.querySelectorAll('.node');
  const edges = svg.querySelectorAll('.edge');
  let activeNode = null;

  // Click on node: highlight relatives
  nodes.forEach(node => {
    node.addEventListener('click', function(e) {
      e.stopPropagation();
      const id = this.dataset.id;
      const parents = (this.dataset.parents || '').split(',').filter(Boolean);
      const children = (this.dataset.children || '').split(',').filter(Boolean);
      const siblings = (this.dataset.siblings || '').split(',').filter(Boolean);
      const connected = new Set([id, ...parents, ...children, ...siblings]);

      if (activeNode === id) {
        // Toggle off
        nodes.forEach(n => n.classList.remove('dimmed'));
        edges.forEach(e => e.classList.remove('dimmed'));
        activeNode = null;
      } else {
        nodes.forEach(n => {
          if (connected.has(n.dataset.id)) {
            n.classList.remove('dimmed');
          } else {
            n.classList.add('dimmed');
          }
        });
        edges.forEach(e => e.classList.add('dimmed'));
        activeNode = id;
      }
    });

    // Hover: show tooltip
    node.addEventListener('mouseenter', function(e) {
      const title = this.querySelector('title');
      if (!title) return;
      const text = title.textContent;
      const lines = text.split(' | ');

      // Anchor above the node, in the SVG's own coordinates, so the tooltip
      // follows the node through zoom and pan.
      const bbox = this.getBBox();
      const toSvg = svg.getScreenCTM().inverse().multiply(this.getScreenCTM());
      const anchor = new DOMPoint(bbox.x + bbox.width / 2, bbox.y).matrixTransform(toSvg);

      tooltip.removeAttribute('transform');
      tooltipText.textContent = '';
      lines.forEach((line, i) => {
        const tspan = document.createElementNS('http://www.w3.org/2000/svg', 'tspan');
        tspan.textContent = line;
        tspan.setAttribute('x', 0);
        tspan.setAttribute('dy', i === 0 ? '0' : '0.32');
        tooltipText.appendChild(tspan);
      });

      const tbox = tooltipText.getBBox();
      const pad = 0.1;
      tooltipBg.setAttribute('x', tbox.x - pad);
      tooltipBg.setAttribute('y', tbox.y - pad);
      tooltipBg.setAttribute('width', tbox.width + 2 * pad);
      tooltipBg.setAttribute('height', tbox.height + 2 * pad);

      // Centre over the node, kept inside the visible area.
      const vb = svg.viewBox.baseVal;
      let left = anchor.x - tbox.width / 2;
      left = Math.min(left, vb.x + vb.width - tbox.width - 2 * pad);
      left = Math.max(left, vb.x + pad);
      let top = anchor.y - 0.15 - tbox.height - pad;
      if (top < vb.y + pad) top = vb.y + pad;
      tooltip.setAttribute('transform',
        'translate(' + (left - tbox.x) + ',' + (top - tbox.y) + ')');
      tooltip.setAttribute('visibility', 'visible');
    });

    node.addEventListener('mouseleave', function() {
      tooltip.setAttribute('visibility', 'hidden');
    });
  });

  // Click on background: reset
  svg.addEventListener('click', function() {
    nodes.forEach(n => n.classList.remove('dimmed'));
    edges.forEach(e => e.classList.remove('dimmed'));
    activeNode = null;
  });

  // Zoom with the mouse wheel. Opened on its own (the SVG is the document),
  // the wheel zooms. Inlined in a page, only Ctrl/Cmd + wheel zooms (a
  // trackpad pinch arrives as Ctrl + wheel), so the page keeps scrolling.
  let scale = 1;
  let panX = 0, panY = 0;
  let isPanning = false;
  let startX, startY;
  const standalone = svg === document.documentElement;

  function updateTransform() {
    mainGroup.setAttribute('transform',
      'translate(' + panX + ',' + panY + ') scale(' + scale + ')');
  }

  svg.addEventListener('wheel', function(e) {
    if (!standalone && !e.ctrlKey && !e.metaKey) return;
    e.preventDefault();
    const factor = e.deltaY > 0 ? 0.9 : 1.1;
    scale *= factor;
    scale = Math.max(0.1, Math.min(10, scale));
    updateTransform();
  }, {passive: false});

  // Pan with the middle button or Ctrl + drag. Pointer capture keeps the
  // drag on this SVG, with no listener on the window.
  svg.addEventListener('pointerdown', function(e) {
    if (e.button === 1 || (e.button === 0 && e.ctrlKey)) {
      isPanning = true;
      startX = e.clientX;
      startY = e.clientY;
      svg.setPointerCapture(e.pointerId);
      e.preventDefault();
    }
  });

  svg.addEventListener('pointermove', function(e) {
    if (!isPanning) return;
    const vb = svg.viewBox.baseVal;
    const rect = svg.getBoundingClientRect();
    const sx = vb.width / rect.width;
    const sy = vb.height / rect.height;
    panX += (e.clientX - startX) * sx;
    panY += (e.clientY - startY) * sy;
    startX = e.clientX;
    startY = e.clientY;
    updateTransform();
  });

  function stopPanning(e) {
    if (!isPanning) return;
    isPanning = false;
    if (svg.hasPointerCapture(e.pointerId)) svg.releasePointerCapture(e.pointerId);
  }
  svg.addEventListener('pointerup', stopPanning);
  svg.addEventListener('pointercancel', stopPanning);
})();
]]>
</script>"""


def _interactive_js(id_prefix: str) -> str:
    """The interactive script, looking up its elements under ``id_prefix``."""
    return _INTERACTIVE_JS.replace("__ID_PREFIX__", json.dumps(id_prefix))
