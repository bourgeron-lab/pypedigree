"""Phenotype display decisions shared by the matplotlib and SVG renderers.

This module decides *what* each phenotype column looks like: its effective
display mode, the color and text of every band, and the rows of the legend.
The renderers only turn these decisions into matplotlib artists or SVG
elements, so a PNG and an interactive SVG of the same pedigree always agree.

Colors are parsed by `pypedigree.colors`, which reads them as
``matplotlib.colors`` does without importing matplotlib, so a style accepts
the same color names in both outputs and this module needs no matplotlib.
"""

from __future__ import annotations

import math
import warnings
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field, replace
from statistics import NormalDist
from typing import Any, Protocol

from pypedigree import symbols
from pypedigree.carriers import CarrierScheme
from pypedigree.colors import parse_rgb, rgb_to_hex
from pypedigree.config import DEFAULT_PGS_COLOR, PedigreeStyle, is_pgs_column, known_phenotype
from pypedigree.errors import PedigreeError, PedigreeStyleWarning
from pypedigree.models import (
    Individual,
    NodeLayout,
    Pedigree,
    Phenotype,
    PhenotypeDisplayMode,
    Sex,
)

RGB = tuple[float, ...]
Color = str | RGB

# Distinct colors for categorical phenotype values.
CATEGORICAL_COLORS = [
    "#e74c3c",  # red
    "#3498db",  # blue
    "#2ecc71",  # green
    "#f39c12",  # orange
    "#9b59b6",  # purple
    "#1abc9c",  # teal
    "#e67e22",  # dark orange
    "#34495e",  # dark grey
    "#e91e63",  # pink
    "#00bcd4",  # cyan
]

UNKNOWN_FILL = "#d9d9d9"  # binary 0 / missing value
# "No data" (-9 / NA) is grey like "unknown"; its "NA" text tells them apart.
NA_FILL = UNKNOWN_FILL
FALLBACK_COLOR = "#999999"  # column without a phenotype_colors entry
SEPARATOR_COLOR = "#cccccc"  # thin line between bands
NA_TEXT_COLOR = "#555555"  # "NA" written on NA_FILL (contrast 5:1)
DARK_TEXT = "#333333"  # band text on a light background
GRADIENT_STEPS = 20  # swatches making up a legend gradient bar

# Shape of the PGS distribution curve: a gaussian of this width on [0, 1].
PGS_CURVE_SIGMA = 0.18
# The curve is a normal density over [0, 1] with this sigma, centred on 0.5,
# so it spans 0.5 / PGS_CURVE_SIGMA (about 2.8) standard deviations each way:
# a z-score is placed at 0.5 + z * PGS_CURVE_SIGMA, clipped to the box.
PGS_SCALES = ("minmax", "zscore", "percentile")
PGS_TAIL_SD = 2.0  # tails marked on a fixed scale
_PERCENTILE_LIMITS = (0.01, 99.99)  # a percentile of 0 or 100 has no z


# --- Colors -----------------------------------------------------------------


def to_rgb(color: str) -> tuple[float, float, float]:
    """Convert a color string to (r, g, b) floats in [0, 1].

    Accepts the colors matplotlib accepts (`pypedigree.colors`): hex
    (``#e74c3c``), a named color (``grey``), a grey level (``"0.5"``); and,
    when matplotlib is installed, its ``xkcd:`` names and ``C0`` shorthands.

    Raises:
        PedigreeError: If the color is not recognised.
    """
    try:
        return parse_rgb(color)
    except (ValueError, TypeError) as exc:
        raise PedigreeError(
            f"Invalid color {color!r}: expected a hex value like '#e74c3c' "
            f"or a named color like 'grey'"
        ) from exc


def to_hex(color: Color) -> str:
    """Hex form (``#rrggbb``) of a color string or an RGB tuple."""
    if isinstance(color, str):
        return rgb_to_hex(parse_rgb(color))
    return rgb_to_hex(color)


def validate_style_colors(style: PedigreeStyle) -> None:
    """Check every color in a style before any drawing starts.

    Reports all offending values at once, so a style file with several
    typos does not have to be fixed one render at a time.

    Raises:
        PedigreeError: If any color is not recognised, listing each one
            with the style option it came from.
    """
    candidates: list[tuple[str, object]] = [
        ("fill_color", style.fill_color),
        ("outline_color", style.outline_color),
        ("line_color", style.line_color),
    ]
    candidates += [(f"phenotype_colors[{k!r}]", v) for k, v in style.phenotype_colors.items()]
    candidates += [(f"pgs_colors[{k!r}]", v) for k, v in style.pgs_colors.items()]
    for name, conf in style.gradient_diverging.items():
        if not isinstance(conf, dict):
            continue
        for key in ("color_low", "color_high"):
            if key in conf:
                candidates.append((f"gradient_diverging[{name!r}][{key!r}]", conf[key]))

    bad: list[str] = []
    for option, value in candidates:
        try:
            # Any value: parse_rgb raising is how a bad one is detected.
            parse_rgb(value)
        except (ValueError, TypeError):
            bad.append(f"{option} = {value!r}")

    if bad:
        raise PedigreeError(
            "Invalid color(s) in style: "
            + "; ".join(bad)
            + ". Expected hex values like '#e74c3c' or named colors like 'grey'."
        )


def blend_with_white(rgb: Sequence[float], alpha: float) -> RGB:
    """Blend a color toward white by a factor (0=white, 1=full color)."""
    return tuple(1.0 - alpha * (1.0 - c) for c in rgb)


def normalize_to_unit(value: float, low: float, high: float) -> float:
    """Map a value onto [0, 1] against an explicit range, clamping outliers.

    Returns:
        The position in [0, 1]; 0.5 for a degenerate range.
    """
    if high == low:
        return 0.5
    return max(0.0, min(1.0, (value - low) / (high - low)))


def diverging_color_at(conf: Mapping[str, Any], t: float) -> RGB:
    """Color of a diverging (tricolor) gradient at a normalized position.

    Shared by the bands and the legend bar so the two can never drift
    apart. The gradient runs ``color_low`` → white at the midpoint →
    ``color_high``.

    Args:
        conf: A ``gradient_diverging`` entry (``min``, ``mid``, ``max``,
            ``color_low``, ``color_high``; each optional).
        t: Position in [0, 1] between ``min`` and ``max``.
    """
    d_min = conf.get("min", 0)
    d_mid = conf.get("mid", 50)
    d_max = conf.get("max", 100)
    c_low = to_rgb(conf.get("color_low", "#e74c3c"))
    c_high = to_rgb(conf.get("color_high", "#2ecc71"))

    mid_t = normalize_to_unit(d_mid, d_min, d_max)
    if t <= mid_t:
        # Below the midpoint: full color_low at t=0, white at the midpoint.
        frac = t / mid_t if mid_t > 0 else 0.0
        return blend_with_white(c_low, 1.0 - frac)
    # Above the midpoint: white at the midpoint, full color_high at t=1.
    frac = (t - mid_t) / (1.0 - mid_t) if mid_t < 1 else 0.0
    return blend_with_white(c_high, frac)


def diverging_color_for_value(conf: Mapping[str, Any], value: float) -> RGB:
    """Color of a diverging gradient for a raw phenotype value."""
    t = normalize_to_unit(value, conf.get("min", 0), conf.get("max", 100))
    return diverging_color_at(conf, t)


def text_color_on(background: Color) -> str:
    """Dark or white text, whichever reads on a band of this color."""
    if background == "white":
        return DARK_TEXT
    r, g, b = to_rgb(background) if isinstance(background, str) else background[:3]
    brightness = 0.299 * r + 0.587 * g + 0.114 * b
    return "white" if brightness < 0.55 else DARK_TEXT


def pgs_curve(n_points: int) -> list[tuple[float, float]]:
    """The PGS distribution curve as (position, density) pairs on [0, 1].

    The highest sample has density 1. Both renderers draw the curve beside
    each node and in the legend from these points.
    """
    ts = [i / (n_points - 1) for i in range(n_points)]
    dens = [math.exp(-0.5 * ((t - 0.5) / PGS_CURVE_SIGMA) ** 2) for t in ts]
    peak = max(dens)
    return [(t, d / peak) for t, d in zip(ts, dens, strict=True)]


# --- Per-figure resolution --------------------------------------------------


@dataclass(frozen=True)
class PhenotypeRender:
    """How each phenotype column is drawn, resolved once per figure.

    Bands and legend read the same resolution, so they can never disagree,
    and nothing is carried over between successive renders in one process.

    Attributes:
        effective_modes: Column name → "binary", "gradient", "categorical",
            "outline", "carrier" or "pgs". AUTO columns are resolved from
            their values.
        categorical_colors: Column name → {category value → hex color}, for
            categorical columns only.
        quadrants: With the "quadrants" layout, each phenotype column or
            carrier source → its quadrant ("top-left", ...), in drawing
            order. Empty with the "bands" layout, or when the figure has
            more columns than quadrants.
        outline_rings: The most colored borders one individual of the
            drawn families has, one per ``o_`` column they are affected in
            (`with_outline_rings`, #46). Above one, every symbol is spaced
            and labelled for that many.
    """

    effective_modes: dict[str, str] = field(default_factory=dict)
    categorical_colors: dict[str, dict[str, str]] = field(default_factory=dict)
    quadrants: dict[str, str] = field(default_factory=dict)
    outline_rings: int = 0

    def effective(self, name: str) -> str:
        """Effective display mode of a column ("binary" if unknown)."""
        return self.effective_modes.get(name, "binary")

    def categorical_color(self, name: str, value: str) -> str | None:
        """Color assigned to a category value, or None if unassigned."""
        return self.categorical_colors.get(name, {}).get(value)

    def band_names(self, names: Sequence[str]) -> list[str]:
        """Columns drawn as horizontal bands inside the symbol, in order."""
        return [n for n in names if self.effective(n) not in ("pgs", "outline", "carrier")]

    def outline_names(self, names: Sequence[str]) -> list[str]:
        """Columns drawn as a colored symbol border (``o_`` prefix)."""
        return [n for n in names if self.effective(n) == "outline"]

    def pgs_names(self, names: Sequence[str]) -> list[str]:
        """Columns drawn on the PGS distribution curve."""
        return [n for n in names if self.effective(n) == "pgs"]


def is_na_value(value: float | None) -> bool:
    """Whether a phenotype value carries no data (missing or NaN)."""
    return value is None or (isinstance(value, float) and math.isnan(value))


def detect_auto_mode(name: str, individuals: Sequence[Individual]) -> str:
    """Resolve an AUTO phenotype column from all of its values.

    The decision is taken per column rather than per individual, so a
    column never mixes binary and gradient bands across the diagram.

    Returns:
        "categorical" if any value is a string, "gradient" if any numeric
        value falls outside the 0/1/2 coding, "binary" otherwise.
    """
    for ind in individuals:
        if ind.extra_phenotypes_categorical.get(name) is not None:
            return "categorical"
    for ind in individuals:
        raw = ind.extra_phenotypes_raw.get(name)
        if is_na_value(raw) or raw in (0.0, 1.0, 2.0):
            continue
        return "gradient"
    return "binary"


def assign_categorical_colors(
    name: str,
    individuals: Sequence[Individual],
    style: PedigreeStyle,
) -> dict[str, str]:
    """Assign a distinct color to each category value of one column.

    Values are taken in file order, so the mapping is reproducible for a
    given pedigree. The column's own ``phenotype_colors`` entry is used
    for the first category and then excluded from the palette, otherwise
    two categories could be handed the same color (the palettes overlap).

    Returns:
        Mapping of category value → color.
    """
    values: list[str] = []
    for ind in individuals:
        cat_val = ind.extra_phenotypes_categorical.get(name)
        if cat_val is not None and cat_val not in values:
            values.append(cat_val)

    base = style.phenotype_colors.get(name)
    # Skip colors used by other phenotype columns, so a category swatch is
    # not confused with the band of another column on the same individual.
    taken = {c for n, c in style.phenotype_colors.items() if n != name}
    palette = [c for c in CATEGORICAL_COLORS if c != base and c not in taken]
    if len(palette) < len(values):
        # Too many categories to stay clear of the other columns.
        palette = [c for c in CATEGORICAL_COLORS if c != base]
    if base:
        palette.insert(0, base)

    return {val: palette[i % len(palette)] for i, val in enumerate(values)}


def resolve_phenotype_render(
    pedigree: Pedigree,
    phenotype_names: Sequence[str],
    display_modes: Mapping[str, PhenotypeDisplayMode],
    style: PedigreeStyle,
) -> PhenotypeRender:
    """Resolve how every phenotype column of a pedigree is drawn.

    Args:
        pedigree: The pedigree whose values decide: the one rendered, or for a
            subset the whole pedigree (`Pedigree.column_scope`).
        phenotype_names: Ordered phenotype column names.
        display_modes: Modes declared by column prefixes.
        style: Resolved style.
    """
    individuals = list(pedigree.individuals())

    effective_modes: dict[str, str] = {}
    categorical_colors: dict[str, dict[str, str]] = {}

    for name in phenotype_names:
        mode = display_modes.get(name, PhenotypeDisplayMode.AUTO)
        if is_pgs_column(name):
            effective_modes[name] = "pgs"
            continue
        if mode == PhenotypeDisplayMode.OUTLINE:
            effective_modes[name] = "outline"
            continue
        if mode == PhenotypeDisplayMode.CARRIER:
            effective_modes[name] = "carrier"
            continue

        if mode == PhenotypeDisplayMode.CATEGORICAL:
            effective = "categorical"
        elif mode in (PhenotypeDisplayMode.GRADIENT, PhenotypeDisplayMode.GRADIENT_INVERTED):
            effective = "gradient"
        elif mode == PhenotypeDisplayMode.BINARY:
            effective = "binary"
        elif name in style.gradient_diverging or name in style.gradient_range:
            # An explicit gradient configuration implies a continuous column.
            effective = "gradient"
        else:
            effective = detect_auto_mode(name, individuals)
        effective_modes[name] = effective

        if effective == "categorical":
            categorical_colors[name] = assign_categorical_colors(name, individuals, style)

    return PhenotypeRender(effective_modes, categorical_colors)


# --- Names shown (#56) --------------------------------------------------------

LEGEND_MODES = ("auto", "compact", "full")


def display_name(name: str, style: PedigreeStyle) -> str:
    """The name a column is shown with: its `phenotype_labels` entry, or its
    name with underscores as spaces."""
    label = style.phenotype_labels.get(name)
    return label if label is not None else name.replace("_", " ")


def pgs_trait_name(name: str, style: PedigreeStyle) -> str:
    """A PGS column's name in the legend's key: its label, or its name
    without "PGS_" ("PGS_autism" -> "autism")."""
    label = style.phenotype_labels.get(name)
    if label is not None:
        return label
    trait = name[4:] if name.upper().startswith("PGS_") else name
    return trait.replace("_", " ")


def check_legend(style: PedigreeStyle) -> None:
    """Reject an unknown ``legend`` or an invalid ``phenotype_labels``.

    Raises:
        PedigreeError: If ``style.legend`` is not one of `LEGEND_MODES`, or
            a label is not a string.
    """
    if style.legend not in LEGEND_MODES:
        raise PedigreeError(
            f"Unknown legend {style.legend!r}: expected one of {', '.join(LEGEND_MODES)}"
        )
    bad = [k for k, v in style.phenotype_labels.items() if not isinstance(v, str)]
    if bad:
        raise PedigreeError(f"phenotype_labels must map names to text; not for {bad!r}")


# --- Bands -------------------------------------------------------------------


class _HasPhenotypes(Protocol):
    """A NodeLayout, as seen by the band decisions."""

    extra_phenotypes: dict[str, float | None]
    extra_phenotypes_raw: dict[str, float | None]
    extra_phenotypes_categorical: dict[str, str | None]


@dataclass(frozen=True)
class BandCell:
    """One phenotype band of one individual.

    Attributes:
        color: Fill color: a color string or an RGB tuple.
        text: Text written in the band, or None.
        is_na: Whether the value is "no data" (-9 / NA); the band is grey
            (``NA_FILL``, as for "unknown") and reads "NA" unless every band
            of the symbol is NA, in which case a single centred "NA"
            replaces them. With ``style.show_na_labels`` off there is no
            "NA" at all and ``text`` is None.
        short: The value alone, for a gradient band whose ``text`` is
            "name: value"; None when ``text`` is already as short (#55).
    """

    color: Color
    text: str | None
    is_na: bool
    short: str | None = None

    @property
    def text_color(self) -> str:
        """Color of the band's text."""
        if self.is_na:
            return NA_TEXT_COLOR
        return text_color_on(self.color)


def band_cell(
    node: _HasPhenotypes,
    name: str,
    style: PedigreeStyle,
    render: PhenotypeRender,
    display_modes: Mapping[str, PhenotypeDisplayMode],
) -> BandCell:
    """Color and text of one phenotype band of one individual."""
    color_hex = style.phenotype_colors.get(name, FALLBACK_COLOR)
    norm_val = node.extra_phenotypes.get(name)
    raw_val = node.extra_phenotypes_raw.get(name)
    cat_val = node.extra_phenotypes_categorical.get(name)
    mode = display_modes.get(name, PhenotypeDisplayMode.AUTO)
    inverted = mode == PhenotypeDisplayMode.GRADIENT_INVERTED
    # Effective display type is resolved per column, not per individual,
    # so every band of a column is drawn the same way as its legend row.
    effective = render.effective(name)

    if isinstance(norm_val, float) and math.isnan(norm_val):
        return BandCell(NA_FILL, "NA" if style.show_na_labels else None, True)

    color: Color
    text: str | None = None
    if effective == "categorical":
        cat_color = render.categorical_color(name, cat_val) if cat_val is not None else None
        color = cat_color if cat_color is not None else UNKNOWN_FILL
        text = cat_val
    elif effective == "binary":
        if norm_val is None:
            color = UNKNOWN_FILL
        elif norm_val == 1.0:
            color = color_hex
        else:
            color = "white"
        if style.show_phenotype_labels and norm_val is not None:
            text = "YES" if norm_val == 1.0 else "NO"
    else:
        rgb = to_rgb(color_hex)
        rv = raw_val if raw_val is not None else 0
        if norm_val is None and raw_val is None:
            color = UNKNOWN_FILL
        elif name in style.gradient_diverging:
            color = diverging_color_for_value(style.gradient_diverging[name], rv)
        elif name in style.gradient_range:
            g_min, g_max = style.gradient_range[name]
            t = normalize_to_unit(rv, g_min, g_max)
            color = blend_with_white(rgb, 1.0 - t if inverted else t)
        else:
            nv = norm_val if norm_val is not None else 0
            color = blend_with_white(rgb, 1.0 - nv if inverted else nv)
        if raw_val is not None:
            # Gradient values are always written, shortened or smaller when
            # they do not fit (`fit_text`).
            shown = display_name(name, style)
            return BandCell(color, f"{shown}: {raw_val:g}", False, f"{raw_val:g}")
    return BandCell(color, text, False)


def pregnancy_loss_fill(
    node: NodeLayout,
    band_names: Sequence[str],
    style: PedigreeStyle,
    render: PhenotypeRender,
    display_modes: Mapping[str, PhenotypeDisplayMode],
) -> str:
    """Fill of a pregnancy-loss triangle, which is too small for bands.

    In a figure with phenotype bands, the color the individual's first band
    would have (grey for no data, without the "NA" the triangle has no room
    for); otherwise the affection status, as for any symbol.
    """
    if band_names:
        cell = band_cell(node, band_names[0], style, render, display_modes)
        return to_hex(cell.color)
    return style.fill_color if node.phenotype == Phenotype.AFFECTED else "white"


def all_bands_na(node: _HasPhenotypes, band_names: Sequence[str]) -> bool:
    """Whether every phenotype band of a symbol is "no data"."""
    return all(is_na_value(node.extra_phenotypes.get(n)) for n in band_names)


# --- Quadrants ----------------------------------------------------------------

PHENOTYPE_LAYOUTS = ("auto", "quadrants", "bands")
QUADRANT_POSITIONS = ("top-left", "top-right", "bottom-left", "bottom-right")
# Conventional quadrant of the known phenotypes (see config.known_phenotype).
KNOWN_PHENOTYPE_QUADRANTS = {
    "autism": "top-left",
    "adhd": "top-right",
    "intellectual disability": "bottom-left",
    "epilepsy": "bottom-right",
}
# A quadrant box reaches this far from the centre, relative to half the
# shape size, so it covers the diamond too; it is clipped to the shape.
QUADRANT_REACH = 1.5
# Quadrant text, relative to the style's font size (bands use 0.75).
QUADRANT_FONT_SCALE = 0.65


def check_phenotype_layout(style: PedigreeStyle) -> None:
    """Reject an unknown layout or quadrant name before anything is drawn.

    Raises:
        PedigreeError: If `phenotype_layout` is not "auto", "quadrants" or "bands",
            or `phenotype_quadrants` names an unknown quadrant or puts two
            columns in the same one.
    """
    if style.phenotype_layout not in PHENOTYPE_LAYOUTS:
        raise PedigreeError(
            f"Unknown phenotype_layout {style.phenotype_layout!r}: "
            f"expected one of {', '.join(PHENOTYPE_LAYOUTS)}"
        )
    for name, position in style.phenotype_quadrants.items():
        if position not in QUADRANT_POSITIONS:
            raise PedigreeError(
                f"phenotype_quadrants: {name!r} is given {position!r}; "
                f"expected one of {', '.join(QUADRANT_POSITIONS)}"
            )
    positions = list(style.phenotype_quadrants.values())
    shared = sorted({p for p in positions if positions.count(p) > 1})
    if shared:
        raise PedigreeError(f"phenotype_quadrants puts several columns in {', '.join(shared)}")


def assign_quadrants(columns: Sequence[str], style: PedigreeStyle) -> dict[str, str]:
    """The quadrant of each column, for the "quadrants" and "auto" layouts.

    A column listed in `style.phenotype_quadrants` takes that quadrant; a
    known phenotype its conventional one (autism top left, ADHD top right,
    intellectual disability bottom left, epilepsy bottom right) when free;
    the others the free quadrants in reading order.

    Args:
        columns: Phenotype band columns, then carrier sources, in order.
        style: Resolved style.

    Returns:
        Column → quadrant, in the order of `columns`; empty with the
        "bands" layout, without columns, or with more than four: the figure
        then uses bands, silently with "auto" (the default), with a
        `PedigreeStyleWarning` when "quadrants" was asked for.
    """
    if style.phenotype_layout == "bands" or not columns:
        return {}
    if len(columns) > len(QUADRANT_POSITIONS):
        if style.phenotype_layout == "auto":
            return {}
        warnings.warn(
            f"phenotype_layout 'quadrants' holds at most {len(QUADRANT_POSITIONS)} "
            f"columns and this figure has {len(columns)} ({', '.join(columns)}): "
            "drawn as bands",
            PedigreeStyleWarning,
            stacklevel=3,
        )
        return {}
    assigned = {n: style.phenotype_quadrants[n] for n in columns if n in style.phenotype_quadrants}
    free = [p for p in QUADRANT_POSITIONS if p not in assigned.values()]
    for name in columns:
        key = known_phenotype(name)
        position = KNOWN_PHENOTYPE_QUADRANTS.get(key) if key else None
        if name not in assigned and position in free:
            assigned[name] = position
            free.remove(position)
    for name in columns:
        if name not in assigned:
            assigned[name] = free.pop(0)
    return {name: assigned[name] for name in columns}


def with_quadrants(
    render: PhenotypeRender,
    band_names: Sequence[str],
    carriers: CarrierScheme,
    style: PedigreeStyle,
) -> PhenotypeRender:
    """`render` with its quadrant assignment, once the carrier sources are known.

    With no phenotype column, nothing is divided into quadrants: both
    renderers then draw the affection status and the carrier sections on
    the whole symbol, as with bands, so the legend has no quadrant key
    either.
    """
    if not band_names:
        return replace(render, quadrants={})
    return replace(render, quadrants=assign_quadrants([*band_names, *carriers.keys], style))


def quadrant_box(
    position: str, cx: float, cy: float, half: float, up: float = 1.0
) -> tuple[float, float, float, float]:
    """A quadrant's box, to be clipped to the symbol: (x, y, width, height).

    (x, y) is the corner with the smallest coordinates in the drawing's own
    axes, so the SVG renderer, whose y axis points down, passes ``up=-1``.
    """
    reach = half * QUADRANT_REACH
    x = cx - reach if position.endswith("left") else cx
    upper_in_drawing = position.startswith("top") == (up > 0)
    y = cy if upper_in_drawing else cy - reach
    return x, y, reach, reach


def quadrant_center(
    position: str, cx: float, cy: float, half: float, up: float = 1.0
) -> tuple[float, float]:
    """Where a quadrant's text goes: halfway between the centre and the corner."""
    dx = -half / 2 if position.endswith("left") else half / 2
    dy = half / 2 if position.startswith("top") else -half / 2
    return cx + dx, cy + up * dy


def quadrant_cell(
    node: _HasPhenotypes,
    name: str,
    style: PedigreeStyle,
    render: PhenotypeRender,
    display_modes: Mapping[str, PhenotypeDisplayMode],
) -> BandCell:
    """Color and text of one phenotype quadrant: a band cell, shorter text.

    A quadrant is half as wide as a band, so a gradient value is written
    alone ("103"), without its column name, which the legend key gives.
    """
    cell = band_cell(node, name, style, render, display_modes)
    if cell.is_na or render.effective(name) != "gradient":
        return cell
    raw = node.extra_phenotypes_raw.get(name)
    return BandCell(cell.color, None if raw is None else f"{raw:g}", False)


# --- Text that fits its band or quadrant (#55) --------------------------------

BAND_TEXT_MODES = ("auto", "name_value", "value", "none")
# With band_text "auto", a text too wide for its place is drawn smaller, down
# to this fraction of `font_size` (or to `min_font_size` when it is set),
# before it is left out.
FIT_FLOOR_SCALE = 0.5
_FIT_PAD_EM = 0.15  # space kept on each side of a fitted text, in ems
_FIT_HEIGHT_EM = 0.8  # height of a line of digits and capitals, in ems


@dataclass(frozen=True)
class FittedText:
    """The text a band or quadrant shows, and its size in points."""

    text: str
    size: float


# Width available to a centred line of text of a given height (both in
# layout units): see `band_room` and `quadrant_room`.
Room = Callable[[float], float]


def check_band_text(style: PedigreeStyle) -> None:
    """Reject an unknown ``band_text`` before anything is drawn.

    Raises:
        PedigreeError: If ``style.band_text`` is not one of `BAND_TEXT_MODES`.
    """
    if style.band_text not in BAND_TEXT_MODES:
        raise PedigreeError(
            f"Unknown band_text {style.band_text!r}: expected one of {', '.join(BAND_TEXT_MODES)}"
        )


def cell_texts(cell: BandCell, mode: str) -> tuple[str, ...]:
    """The texts a band may show, the preferred first, under ``band_text``.

    "NA" is always offered (``show_na_labels`` decides on it); "none" drops
    the values, "value" writes them without the column name, "name_value"
    keeps the name, and "auto" offers both, the name first.
    """
    if cell.text is None:
        return ()
    if cell.is_na:
        return (cell.text,)
    short = cell.short or cell.text
    return {
        "none": (),
        "value": (short,),
        "name_value": (cell.text,),
    }.get(mode, tuple(dict.fromkeys((cell.text, short))))


def band_room(sex: Sex, half: float, dy: float) -> Room:
    """Width inside a symbol for a text centred ``dy`` from its centre."""
    return lambda h: 2 * symbols.shape_half_width_at(sex, half, abs(dy) + h / 2)


def quadrant_room(sex: Sex, half: float, position: str) -> Room:
    """Width for a text at `quadrant_center`, between the symbol's edge and its axis."""
    dx, dy = half / 2, half / 2

    def room(h: float) -> float:
        edge = symbols.shape_half_width_at(sex, half, dy + h / 2)
        return 2 * max(0.0, min(dx, edge - dx))

    return room


def fit_text(
    texts: Sequence[str],
    size: float,
    style: PedigreeStyle,
    pt: float,
    room: Room,
    height: float,
) -> FittedText | None:
    """The text and size a band or quadrant shows, so that it fits (#55).

    With ``band_text`` "auto", the first text that fits at ``size`` is
    taken; else the last (shortest) is drawn smaller, down to
    `min_font_size`, or `FIT_FLOOR_SCALE` times `font_size` without it;
    else nothing. Other modes draw their text at ``size`` as it comes. The
    renderers call it with their own scale: the SVG with its fixed one, the
    PNG when the figure is drawn, so at the same scale both choose the same.

    Args:
        texts: Candidates, preferred first (`cell_texts`).
        size: The size wanted, in points.
        style: The style (``band_text``, ``font_size``, ``min_font_size``).
        pt: Layout units per point at the drawing's scale.
        room: Width available for a text of a given height (`band_room`,
            `quadrant_room`), in layout units.
        height: Height of the band or quadrant, in layout units.

    Returns:
        The text and its size, or None when nothing is written.
    """
    if not texts:
        return None
    if style.band_text != "auto":
        return FittedText(texts[0], size)
    if style.min_font_size is not None:
        floor = min(size, style.min_font_size)
    else:
        floor = min(size, style.font_size * FIT_FLOOR_SCALE)

    def fits(text: str, s: float) -> bool:
        h = s * pt * _FIT_HEIGHT_EM
        width = symbols.text_width(text, s, bold=True) / symbols.UNITS_PER_POINT
        return h <= height and (width + 2 * _FIT_PAD_EM * s) * pt <= room(h)

    for text in texts:
        if fits(text, size):
            return FittedText(text, size)
    text = texts[-1]
    if not fits(text, floor):
        return None
    low, high = floor, size
    for _ in range(20):
        mid = (low + high) / 2
        low, high = (mid, high) if fits(text, mid) else (low, mid)
    return FittedText(text, low)


#: Width of a colored border (``o_`` column), in multiples of the shape's line width.
OUTLINE_WIDTH_SCALE = 2.5
#: White gap between two colored borders, in points (#46).
OUTLINE_RING_GAP_PT = 1.0
#: More borders than this make a symbol much larger: a warning says so.
MAX_QUIET_OUTLINE_RINGS = 3


def outline_columns(node: _HasPhenotypes, outline_names: Sequence[str]) -> list[str]:
    """The ``o_`` columns this individual is affected in, in column order."""
    return [n for n in outline_names if node.extra_phenotypes.get(n) == 1.0]


def outline_colors(
    node: _HasPhenotypes, outline_names: Sequence[str], style: PedigreeStyle
) -> list[str]:
    """One border color per ``o_`` column this individual is affected in.

    The first is drawn on the symbol's outline, the others around it,
    outwards, in column order (#46).
    """
    return [
        style.phenotype_colors.get(n, FALLBACK_COLOR) for n in outline_columns(node, outline_names)
    ]


def outline_color(
    node: _HasPhenotypes, outline_names: Sequence[str], style: PedigreeStyle
) -> str | None:
    """Border color from the first ``o_`` column this individual is affected in."""
    colors = outline_colors(node, outline_names, style)
    return colors[0] if colors else None


def outline_ring_offsets(count: int, style: PedigreeStyle) -> list[float]:
    """How far each border after the first lies outside the symbol, in points.

    Measured from the outline to the border's centre line: one border's
    width and a white gap (`OUTLINE_RING_GAP_PT`) per border.
    """
    step = style.shape_line_width * OUTLINE_WIDTH_SCALE + OUTLINE_RING_GAP_PT
    return [k * step for k in range(1, count)]


def outline_reach(count: int, style: PedigreeStyle) -> float:
    """Room the borders after the first take around a symbol, in layout units.

    Converted with `symbols.UNITS_PER_POINT`, the scale of the SVG; the PNG
    draws at that scale or larger (see `renderer._OutlineRings`).
    """
    offsets = outline_ring_offsets(count, style)
    return offsets[-1] * symbols.UNITS_PER_POINT if offsets else 0.0


def with_outline_rings(
    render: PhenotypeRender, nodes: Iterable[_HasPhenotypes], names: Sequence[str]
) -> PhenotypeRender:
    """`render` with the most borders one individual of `nodes` has (#46).

    Warns when that is more than `MAX_QUIET_OUTLINE_RINGS`: every symbol then
    takes much more room.
    """
    outline_names = render.outline_names(names)
    if not outline_names:
        return render
    most, who = 0, ""
    for node in nodes:
        count = len(outline_columns(node, outline_names))
        if count > most:
            most = count
            who = str(getattr(node, "sample_id", "") or getattr(node, "individual_id", ""))
    if most > MAX_QUIET_OUTLINE_RINGS:
        warnings.warn(
            f"{who} has {most} colored borders (one per o_ column they are affected in): "
            f"above {MAX_QUIET_OUTLINE_RINGS}, symbols take much more room",
            PedigreeStyleWarning,
            stacklevel=3,
        )
    return replace(render, outline_rings=most)


def check_pgs_scale(style: PedigreeStyle) -> None:
    """Reject an unknown ``pgs_scale`` before anything is drawn.

    Raises:
        PedigreeError: If ``style.pgs_scale`` is not one of `PGS_SCALES`.
    """
    if style.pgs_scale not in PGS_SCALES:
        raise PedigreeError(
            f"Unknown pgs_scale {style.pgs_scale!r}: expected one of {', '.join(PGS_SCALES)}"
        )


def check_min_font_size(style: PedigreeStyle) -> None:
    """Reject a ``min_font_size`` that is not a positive number, before drawing.

    Raises:
        PedigreeError: If ``style.min_font_size`` is set and not a number
            above zero.
    """
    value = style.min_font_size
    if value is None:
        return
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not value > 0:
        raise PedigreeError(
            f"min_font_size must be a number of points above zero, or None; got {value!r}"
        )


def pgs_zscore(raw: float, scale: str) -> float:
    """The z-score a PGS value stands for, on the "zscore" or "percentile" scale."""
    if scale == "percentile":
        low, high = _PERCENTILE_LIMITS
        return NormalDist().inv_cdf(min(max(raw, low), high) / 100)
    return raw


def pgs_position(z: float) -> float:
    """Position of a z-score on the curve's [0, 1] axis, clipped to the box."""
    return min(max(0.5 + z * PGS_CURVE_SIGMA, 0.0), 1.0)


def pgs_tail_positions(scale: str) -> tuple[float, ...]:
    """Where the +/-2 SD tails are marked on the curve: none on "minmax"."""
    if scale == "minmax":
        return ()
    return (pgs_position(-PGS_TAIL_SD), pgs_position(PGS_TAIL_SD))


def pgs_caption_above(node: NodeLayout) -> bool:
    """Whether the "PGS" caption of a curve goes above it rather than below.

    A proband's or a consultand's marker sits at the symbol's lower left,
    just under the curve, where the caption would otherwise be written.
    """
    return node.is_proband or node.notation.is_consultand


def active_pgs(
    node: _HasPhenotypes, pgs_names: Sequence[str], scale: str = "minmax"
) -> list[tuple[str, float]]:
    """PGS columns with a value for this individual, with their [0, 1] position.

    On "minmax" the position is the value normalised over the pedigree; on
    "zscore" and "percentile" it is the value's z-score on the drawn
    normal curve (`pgs_position`), the same for every pedigree.
    """
    out = []
    for name in pgs_names:
        if scale == "minmax":
            value = node.extra_phenotypes.get(name)
            if not is_na_value(value):
                out.append((name, value or 0.0))
            continue
        raw = node.extra_phenotypes_raw.get(name)
        if not is_na_value(raw) and raw is not None:
            out.append((name, pgs_position(pgs_zscore(raw, scale))))
    return out


def _ordinal(n: int) -> str:
    suffix = "th" if 10 <= n % 100 <= 20 else {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th")
    return f"{n}{suffix}"


def pgs_value_text(raw: float, scale: str) -> str:
    """A PGS value as the tooltip writes it.

    "z = +1.10, 86th percentile" on "zscore", "86th percentile, z = +1.08"
    on "percentile", the value alone on "minmax".
    """
    if scale == "minmax":
        return f"{raw:g}"
    z = pgs_zscore(raw, scale)
    percentile = round(NormalDist().cdf(z) * 100)
    percentile = min(max(percentile, 1), 99)
    if scale == "percentile":
        return f"{_ordinal(round(min(max(raw, 0), 100)))} percentile, z = {z:+.2f}"
    return f"z = {z:+.2f}, {_ordinal(percentile)} percentile"


# --- PGS values and markers (#54) --------------------------------------------

MINUS = "\u2212"  # the typographic minus of the values written beside the curve
PGS_DOT_PT = 3.0  # diameter of a PGS marker's dot, in points
PGS_DOT_SPREAD = 1.5  # distance between the centres of dots set side by side, in diameters
PGS_BOX_PAD = 0.02  # how far the curve's rounded box reaches past its edge, in layout units
PGS_VALUE_GAP_EM = 0.4  # between a value and the curve's box, in ems of the value
PGS_VALUE_LINE_EM = 1.15  # smallest distance between two values' centres, in ems
PGS_VALUE_HALO_PT = 2.0  # white outline around a value, so a line behind it does not cross it


def pgs_marker_text(raw: float, scale: str) -> str:
    """A PGS value as written beside its marker (``show_pgs_values``).

    Signed with one decimal on "zscore" ("+2.1", "\u22120.4"), the rounded
    percentile on "percentile" ("98"), three significant digits on
    "minmax" ("0.712").
    """
    if scale == "zscore":
        text = f"{raw:+.1f}"
    elif scale == "percentile":
        text = f"{round(raw):d}"
    else:
        text = f"{raw:.3g}"
    return text.replace("-", MINUS)


def pgs_marker_values(node: _HasPhenotypes, pgs_names: Sequence[str], scale: str) -> dict[str, str]:
    """The text of each PGS value an individual has, by column."""
    out = {}
    for name in pgs_names:
        raw = node.extra_phenotypes_raw.get(name)
        if raw is not None and not is_na_value(raw):
            out[name] = pgs_marker_text(raw, scale)
    return out


def pgs_values_width(
    nodes: Iterable[_HasPhenotypes], pgs_names: Sequence[str], style: PedigreeStyle
) -> float:
    """Room the values need on the left of the curves, in layout units.

    The widest value of these individuals plus the gap to the box, measured
    at the SVG's scale as labels are (#41); 0 without ``show_pgs_values``.
    """
    if not style.show_pgs_values or not pgs_names:
        return 0.0
    size = style.text_size(symbols.PGS_CAPTION_SCALE)
    widths = [
        symbols.text_width(text, size)
        for ind in nodes
        for text in pgs_marker_values(ind, pgs_names, style.pgs_scale).values()
    ]
    if not widths:
        return 0.0
    return max(widths) + PGS_VALUE_GAP_EM * size * symbols.UNITS_PER_POINT + PGS_BOX_PAD


def spread_dots(ys: Sequence[float], dot: float) -> list[float]:
    """Horizontal offsets that keep the markers' dots apart (#54).

    Dots closer than ``dot`` (their diameter) vertically form a group, set
    side by side, centred on the curve's axis, `PGS_DOT_SPREAD` diameters
    apart. A dot alone keeps offset 0.

    Args:
        ys: Height of each dot, in any unit.
        dot: The dots' diameter, in the same unit.

    Returns:
        The offset of each dot, in the order of ``ys``, in that unit.
    """
    order = sorted(range(len(ys)), key=lambda i: ys[i])
    offsets = [0.0] * len(ys)
    group: list[int] = []
    for i in [*order, None]:
        if i is not None and (not group or ys[i] - ys[group[-1]] < dot):
            group.append(i)
            continue
        for k, j in enumerate(group):
            offsets[j] = (k - (len(group) - 1) / 2) * dot * PGS_DOT_SPREAD
        group = [] if i is None else [i]
    return offsets


def spread_labels(ys: Sequence[float], step: float, low: float, high: float) -> list[float]:
    """Heights for labels at ``ys``, at least ``step`` apart, kept in [low, high].

    Labels too close are set ``step`` apart around the mean of their
    heights, in order; a group reaching past ``low`` or ``high`` is moved
    back inside (when it fits).

    Returns:
        The height of each label, in the order of ``ys``.
    """
    order = sorted(range(len(ys)), key=lambda i: ys[i])
    groups: list[list[int]] = []

    def placed(group: list[int]) -> list[float]:
        centre = sum(ys[i] for i in group) / len(group)
        first = centre - step * (len(group) - 1) / 2
        span = step * (len(group) - 1)
        first = min(max(first, low), max(low, high - span))
        return [first + k * step for k in range(len(group))]

    for i in order:
        groups.append([i])
        while len(groups) > 1 and placed(groups[-2])[-1] + step > placed(groups[-1])[0]:
            last = groups.pop()
            groups[-1].extend(last)
    out = [0.0] * len(ys)
    for group in groups:
        for i, y in zip(group, placed(group), strict=True):
            out[i] = y
    return out


# --- Legend ------------------------------------------------------------------

# Legend panel geometry, in the legend's own coordinates (y pointing up). One unit is a
# fixed multiple of the base font size (`style.font_size`), horizontally and
# vertically, so the panel has the same physical size and proportions whatever
# the size of the pedigree above it, and its swatches grow with its text.
# These are internal layout constants, not user-tunable style.
LEG_UNIT_X_EM = 4.2  # one horizontal unit, in font sizes
LEG_UNIT_Y_EM = 2.0  # one vertical unit, in font sizes
LEG_X0 = 0.4  # left margin for swatches and icons
LEG_GRAD_W = 2.6  # width of a gradient bar
LEG_HEADER_H = 1.8  # title + male/female symbol row
LEG_BINARY_HEADER_H = 0.5  # Affected / Unaffected / Unknown captions
LEG_ROW_H = 1.2  # one phenotype or carrier row
LEG_PAD = 0.8  # room below the last row for its captions
LEG_REFERENCE_FONT = 8.0  # font size the marker sizes below are tuned for
LEG_TEXT_SCALE = 0.85  # legend row text, as a fraction of the legend's font
LEG_SWATCH_EDGE = "#cccccc"
LEG_CAPTION_COLOR = "#999999"
LEG_CATEGORY_COLOR = "#666666"
LEG_KEY_SLOTS = 2  # rows taken by the quadrant key
LEG_KEY_SIDE = 1.5  # height of the key's square and circle, in vertical units
LEG_KEY_GAP = 0.12  # between the key's symbols and their labels
LEG_CHAR_WIDTH_EM = 0.55  # average glyph width, to make room for labels
LEG_SHARED_BINARY_LABEL = "Phenotypes in the key (their colors)"
LEG_OUTLINE_RINGS_LABEL = "Several borders: one per column, inner first"
LEG_RING_WIDTH_PT = 1.6  # each border of the "outline_rings" icon
LEG_RING_INSET_PT = 2.2  # its inner border, inside the outer one
LEG_KEY_DASH_EM = 1.4  # length of a PGS score's dash in the legend's key, in ems
LEG_KEY_GAP_EM = 0.3  # between a dash and its score's name
LEG_KEY_SPACE_EM = 0.9  # between one score's name and the next dash


@dataclass(frozen=True)
class LegendMetrics:
    """Derived geometry of the legend panel.

    Lengths are in legend units (see `LEG_UNIT_X_EM`); font and marker sizes
    are in points. All of them scale with the legend's font (`font`), never
    with the figure.

    Attributes:
        n_rows: Number of phenotype/carrier rows (at least 1).
        row_h: Height of one row.
        total_h: Height of the whole panel.
        scale: Marker size relative to the reference font size.
        sw_h: Swatch height.
        sw_w: Swatch width.
        sw_gap: Horizontal gap between swatches.
        text_x: Label x for rows built from three swatches.
        icon_text_x: Label x for rows built from a single icon.
        fs: Row font size.
        label_fs: Caption font size (column headers, scale ends).
        marker_size: Size of plot-marker icons.
        unit_x: One horizontal legend unit, in inches.
        unit_y: One vertical legend unit, in inches.
        font: The font size the legend is built on: `style.font_size`, or
            more so that its row text reaches `style.min_font_size`.
        title_fs: Size of the "Legend" title.
        min_fs: `style.min_font_size`, or 0.
    """

    n_rows: int
    row_h: float
    total_h: float
    scale: float
    sw_h: float
    sw_w: float
    sw_gap: float
    text_x: float
    icon_text_x: float
    fs: float
    label_fs: float
    marker_size: float
    unit_x: float
    unit_y: float
    font: float = 8.0
    title_fs: float = 7.8
    min_fs: float = 0.0

    def caption_fs(self, size_factor: float = 1.0) -> float:
        """Size of a caption at ``size_factor`` times `label_fs`, never below `min_fs`."""
        return max(self.label_fs * size_factor, self.min_fs)

    def row_center(self, index: int, y_top: float) -> float:
        """Vertical centre of row `index`, counting down from `y_top`."""
        return y_top - index * self.row_h - self.row_h / 2

    def swatch_x(self, index: int) -> float:
        """Left edge of the `index`-th swatch in a multi-swatch row."""
        return LEG_X0 + index * (self.sw_w + self.sw_gap)


def legend_metrics(
    n_pheno_rows: int, style: PedigreeStyle, *, binary_headers: bool = False
) -> LegendMetrics:
    """Compute the legend panel geometry for a given number of rows.

    Rows keep a constant height: the panel grows with its row count, and
    the figure grows with the panel (see `_add_legend_axes`).

    Args:
        n_pheno_rows: Number of phenotype and carrier rows to fit.
        style: Style supplying the base font size and family.
        binary_headers: Whether the Affected / Unaffected / Unknown captions
            are drawn above the rows.

    Returns:
        The derived metrics.
    """
    n_rows = max(n_pheno_rows, 1)
    row_h = LEG_ROW_H
    sw_h = row_h * 0.45
    sw_w = 0.8
    sw_gap = 0.35

    # Three-swatch rows put their label just past the last swatch.
    text_x = LEG_X0 + 3 * sw_w + 2 * sw_gap + 0.3
    # Single-icon rows align on the widest icon: PGS curve (ends at 1.3),
    # outline swatch, or carrier marker (~1.2).
    icon_right = max(1.3, LEG_X0 + sw_w, 1.2)
    icon_text_x = icon_right + 0.3

    # With a minimum text size the whole legend grows until its row text
    # reaches it (#51); its smaller captions are then raised to it too.
    min_fs = style.min_font_size or 0.0
    font = max(style.font_size, min_fs / LEG_TEXT_SCALE)
    fs = font * LEG_TEXT_SCALE
    em = font / 72  # inches
    header_h = LEG_HEADER_H + (LEG_BINARY_HEADER_H if binary_headers else 0.0)

    return LegendMetrics(
        n_rows=n_rows,
        row_h=row_h,
        total_h=header_h + n_rows * row_h + LEG_PAD,
        scale=font / LEG_REFERENCE_FONT,
        sw_h=sw_h,
        sw_w=sw_w,
        sw_gap=sw_gap,
        text_x=text_x,
        icon_text_x=icon_text_x,
        fs=fs,
        label_fs=max(fs * 0.75, min_fs),
        marker_size=font,
        unit_x=LEG_UNIT_X_EM * em,
        unit_y=LEG_UNIT_Y_EM * em,
        font=font,
        title_fs=max(fs + 1, font * 0.7),
        min_fs=min_fs,
    )


@dataclass(frozen=True)
class LegendRow:
    """One row of the legend panel.

    Attributes:
        kind: "binary", "binary_shared" (the binary columns of the quadrant
            key, compact legend), "gradient", "categorical", "outline",
            "outline_rings" (how to read several borders, #46), "pgs",
            "legacy_carrier" (half-fill / dot) or "carrier" (patterned).
        label: Row text.
        color: Swatch color ("binary", "outline").
        colors: Gradient bar, left to right ("gradient"); the columns'
            colors, side by side in the "affected" swatch ("binary_shared");
            the two borders of the icon, inner first ("outline_rings").
        low: Caption under the left end of a gradient bar.
        high: Caption under the right end of a gradient bar.
        categories: Category value → color ("categorical").
        pgs: PGS column → marker color ("pgs").
        pattern: Carrier pattern name ("carrier").
        obligate: Obligate-carrier row ("carrier", "legacy_carrier").
        quadrants: Cells of the quadrant key ("quadrants"), one per column.
        key: The PGS scores' names and colors, written after the label
            ("pgs", #56).
    """

    kind: str
    label: str
    color: str = ""
    colors: tuple[RGB, ...] = ()
    low: str = ""
    high: str = ""
    categories: Mapping[str, str] = field(default_factory=dict)
    pgs: Mapping[str, str] = field(default_factory=dict)
    pattern: str = ""
    obligate: bool = False
    quadrants: tuple[QuadrantKeyCell, ...] = ()
    key: tuple[tuple[str, str], ...] = ()

    @property
    def slots(self) -> int:
        """Number of legend rows this row takes (the quadrant key takes two)."""
        return LEG_KEY_SLOTS if self.kind == "quadrants" else 1


@dataclass(frozen=True)
class QuadrantKeyCell:
    """One quadrant of the legend key: which column it shows, and how.

    Attributes:
        position: "top-left", "top-right", "bottom-left" or "bottom-right".
        label: Column name, or the carrier source's label.
        color: Fill (the column's color; white for a carrier source).
        pattern: Carrier pattern name, for a carrier source.
    """

    position: str
    label: str
    color: str
    pattern: str = ""


@dataclass(frozen=True)
class QuadrantKeyGeometry:
    """Where the quadrant key's parts go, in legend units (y up).

    The key is a square and a circle, each cut in four; left-quadrant labels
    sit left of the square, right-quadrant labels right of the circle.

    Attributes:
        square: (x, y_center, width, height) of the square.
        circle: (x_center, y_center, x_radius, y_radius) of the circle; the
            radii differ in legend units so it is round on the page.
        labels: (x, y, text, anchor) per label; anchor is "end" or "start".
    """

    square: tuple[float, float, float, float]
    circle: tuple[float, float, float, float]
    labels: tuple[tuple[float, float, str, str], ...]


def legend_row_centers(rows: Sequence[LegendRow], m: LegendMetrics, y_top: float) -> list[float]:
    """Vertical centre of each row, counting down from `y_top`."""
    centers = []
    slot = 0
    for row in rows:
        centers.append(y_top - (slot + row.slots / 2) * m.row_h)
        slot += row.slots
    return centers


def legend_slot_count(rows: Sequence[LegendRow]) -> int:
    """Rows of the panel the legend rows take (pass to `legend_metrics`)."""
    return sum(r.slots for r in rows)


def quadrant_key_geometry(row: LegendRow, m: LegendMetrics, y_center: float) -> QuadrantKeyGeometry:
    """Lay out the quadrant key row centred on `y_center`."""
    side_y = LEG_KEY_SIDE
    side_x = side_y * m.unit_y / m.unit_x  # square on the page
    char_w = LEG_CHAR_WIDTH_EM * (m.fs / 72) / m.unit_x
    left = [c for c in row.quadrants if c.position.endswith("left")]
    left_w = max((len(c.label) for c in left), default=0) * char_w
    square_x = LEG_X0 + (left_w + LEG_KEY_GAP if left else 0.0)
    circle_cx = square_x + side_x * 1.6 + side_x / 2
    labels = []
    for cell in row.quadrants:
        v = y_center + (side_y / 4 if cell.position.startswith("top") else -side_y / 4)
        if cell.position.endswith("left"):
            labels.append((square_x - LEG_KEY_GAP, v, cell.label, "end"))
        else:
            labels.append((circle_cx + side_x / 2 + LEG_KEY_GAP, v, cell.label, "start"))
    return QuadrantKeyGeometry(
        square=(square_x, y_center, side_x, side_y),
        circle=(circle_cx, y_center, side_x / 2, side_y / 2),
        labels=tuple(labels),
    )


@dataclass(frozen=True)
class LegendKeyItem:
    """One score of the PGS row's key: a dash of its color, then its name.

    Positions are in legend units (see `LEG_UNIT_X_EM`).
    """

    dash_x0: float
    dash_x1: float
    text_x: float
    name: str
    color: str


def legend_key_items(row: LegendRow, label_x: float, m: LegendMetrics) -> list[LegendKeyItem]:
    """Where the PGS row's key goes, after its label written from ``label_x`` (#56).

    Text is measured in DejaVu Sans at the row's size (`symbols.text_width`),
    so both renderers set the dashes and names at the same places.
    """

    def units(points: float) -> float:
        return points / 72 / m.unit_x

    em = units(m.fs)
    x = label_x + units(symbols.text_width(row.label, m.fs) / symbols.UNITS_PER_POINT)
    items = []
    for name, color in row.key:
        x += LEG_KEY_SPACE_EM * em
        dash_x1 = x + LEG_KEY_DASH_EM * em
        text_x = dash_x1 + LEG_KEY_GAP_EM * em
        items.append(LegendKeyItem(x, dash_x1, text_x, name, color))
        x = text_x + units(symbols.text_width(name, m.fs) / symbols.UNITS_PER_POINT)
    return items


def pgs_legend_description(pgs_names: Sequence[str], scale: str = "minmax") -> str:
    """Caption of the PGS row, from its column names and the scale.

    A single column named "PGS" is described in words; otherwise the caption
    ends with a colon, and the row's key (`LegendRow.key`) follows it: each
    score's name after a dash of its color (#56).
    """
    if len(pgs_names) == 1 and pgs_names[0].upper() == "PGS":
        caption, keyed = "PGS distribution (vertical) with individual score", False
    else:
        caption, keyed = "PGS distribution", True
    if scale == "zscore":
        caption += ", z-score, \u00b12 SD marked"
    elif scale == "percentile":
        caption += ", population percentile, \u00b12 SD marked"
    return f"{caption}:" if keyed else caption


def legend_rows(
    phenotype_names: Sequence[str],
    style: PedigreeStyle,
    render: PhenotypeRender,
    display_modes: Mapping[str, PhenotypeDisplayMode],
    carriers: CarrierScheme,
) -> list[LegendRow]:
    """The rows of the legend, in drawing order.

    With the quadrant layout, a key showing which quadrant holds which
    column comes first. Then one row per band column, one per outline
    column, one shared row for all PGS columns, and the carrier rows: the
    half-fill and dot rows with the legacy notation, then one per carrier
    source and per obligate-carrier source.
    """
    rows: list[LegendRow] = []
    names = list(phenotype_names)
    steps = [g / (GRADIENT_STEPS - 1) for g in range(GRADIENT_STEPS)]
    binary_in_key: list[tuple[int, str]] = []  # binary rows of key columns

    for name in render.band_names(names):
        color = style.phenotype_colors.get(name, FALLBACK_COLOR)
        effective = render.effective(name)
        if name in style.gradient_diverging:
            conf = style.gradient_diverging[name]
            rows.append(
                LegendRow(
                    "gradient",
                    display_name(name, style),
                    colors=tuple(diverging_color_at(conf, t) for t in steps),
                    low=str(conf.get("min", 0)),
                    high=str(conf.get("max", 100)),
                )
            )
        elif effective == "gradient":
            inverted = display_modes.get(name) == PhenotypeDisplayMode.GRADIENT_INVERTED
            rgb = to_rgb(color)
            rows.append(
                LegendRow(
                    "gradient",
                    display_name(name, style),
                    colors=tuple(blend_with_white(rgb, 1.0 - t if inverted else t) for t in steps),
                    low="Low",
                    high="High",
                )
            )
        elif effective == "categorical":
            rows.append(
                LegendRow(
                    "categorical",
                    display_name(name, style),
                    categories=render.categorical_colors.get(name, {}),
                )
            )
        else:
            if name in render.quadrants:
                binary_in_key.append((len(rows), color))
            rows.append(LegendRow("binary", display_name(name, style), color=color))

    # Compact legend (#56): with quadrants, the key already says which color
    # is which column, so two or more binary columns share one row.
    compact = style.legend == "compact" or (style.legend == "auto" and bool(render.quadrants))
    if compact and len(binary_in_key) >= 2:
        first = binary_in_key[0][0]
        shared = LegendRow(
            "binary_shared",
            LEG_SHARED_BINARY_LABEL,
            colors=tuple(to_rgb(c) for _i, c in binary_in_key),
        )
        merged = {i for i, _c in binary_in_key}
        rows = [r for i, r in enumerate(rows) if i not in merged]
        rows.insert(first, shared)

    for name in render.outline_names(names):
        rows.append(
            LegendRow(
                "outline",
                f"Outline: {display_name(name, style)} (colored border when affected)",
                color=style.phenotype_colors.get(name, FALLBACK_COLOR),
            )
        )
    outline_names = render.outline_names(names)
    if render.outline_rings > 1 and len(outline_names) > 1:
        # Some individual has several borders: say how to read them (#46).
        rows.append(
            LegendRow(
                "outline_rings",
                LEG_OUTLINE_RINGS_LABEL,
                colors=tuple(
                    to_rgb(style.phenotype_colors.get(n, FALLBACK_COLOR)) for n in outline_names[:2]
                ),
            )
        )

    pgs_names = render.pgs_names(names)
    if pgs_names:
        colors = {n: style.pgs_colors.get(n, DEFAULT_PGS_COLOR) for n in pgs_names}
        plain = len(pgs_names) == 1 and pgs_names[0].upper() == "PGS"
        rows.append(
            LegendRow(
                "pgs",
                pgs_legend_description(pgs_names, style.pgs_scale),
                pgs=colors,
                key=()
                if plain
                else tuple((pgs_trait_name(n, style), colors[n]) for n in pgs_names),
            )
        )

    if carriers.legacy:
        if carriers.has_legacy_plain:
            rows.append(LegendRow("legacy_carrier", "Carrier (half-filled = heterozygous carrier)"))
        if carriers.has_legacy_obligate:
            rows.append(
                LegendRow(
                    "legacy_carrier",
                    "Obligate carrier (center dot = inferred from pedigree)",
                    obligate=True,
                )
            )

    for key in carriers.keys:
        for obligate in (False, True):
            if obligate and key not in carriers.obligate_keys:
                continue
            rows.append(
                LegendRow(
                    "carrier",
                    carriers.label(key, obligate, display_name(key, style)),
                    pattern=carriers.pattern(key),
                    obligate=obligate,
                )
            )

    if render.quadrants:
        # The key comes first: which quadrant holds which column.
        cells = []
        for name, position in render.quadrants.items():
            if name in carriers.keys:
                cells.append(
                    QuadrantKeyCell(
                        position,
                        carriers.label(name, False, display_name(name, style)),
                        "white",
                        carriers.pattern(name),
                    )
                )
            elif render.effective(name) == "categorical":
                cells.append(
                    QuadrantKeyCell(position, display_name(name, style), LEG_CATEGORY_COLOR)
                )
            else:
                cells.append(
                    QuadrantKeyCell(
                        position,
                        display_name(name, style),
                        style.phenotype_colors.get(name, FALLBACK_COLOR),
                    )
                )
        rows.insert(0, LegendRow("quadrants", "", quadrants=tuple(cells)))
    return rows
