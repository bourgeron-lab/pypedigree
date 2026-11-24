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

import bisect
import math
import re
import warnings
from collections.abc import Callable, Collection, Iterable, Mapping, Sequence
from dataclasses import dataclass, field, replace
from difflib import get_close_matches
from statistics import NormalDist
from typing import Any, Protocol

from pypedigree import symbols
from pypedigree.carriers import CarrierScheme
from pypedigree.colors import parse_rgb, rgb_to_hex
from pypedigree.config import (
    DEFAULT_PGS_COLOR,
    PedigreeStyle,
    curve_trait,
    is_curve_column,
    is_pgs_column,
    is_ptdt_column,
    known_phenotype,
)
from pypedigree.errors import PedigreeDataWarning, PedigreeError, PedigreeStyleWarning
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
# "Not applicable" (#107): the band or quadrant is white, as for unaffected,
# with no text; the SVG tooltip says "not applicable". The grey diagonal drawn
# in 0.10.0 read as the deceased line on every sex-limited cancer (#120).
NOT_APPLICABLE_FILL = "white"
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


# --- Clinical thresholds (#57) -----------------------------------------------


def bin_index(value: float, cutoffs: Sequence[float]) -> int:
    """Which step of ``gradient_bins`` a value falls in: 0 below the first
    cut-off; a value at a cut-off goes to the step above it."""
    return bisect.bisect_right(list(cutoffs), value)


def bin_color(rgb: Sequence[float], index: int, cutoffs: Sequence[float], inverted: bool) -> RGB:
    """Color of step ``index``: white to the column's color in even steps,
    the other way for a ``gi_`` column (low = atypical, as its gradient)."""
    t = index / len(cutoffs)  # len(cutoffs) + 1 steps
    return blend_with_white(rgb, 1.0 - t if inverted else t)


def bin_labels(cutoffs: Sequence[float]) -> tuple[str, ...]:
    """The steps' limits, as the legend writes them under each step.

    "< 60", "60–65", "66–75", "≥ 76" for the cut-offs 60, 66, 76: with
    whole-number cut-offs a step ends one below the next; otherwise it is
    written up to the next cut-off ("1.5–2.5"), which it excludes.
    """
    whole = all(float(c).is_integer() for c in cutoffs)
    labels = [f"< {cutoffs[0]:g}"]
    for low, high in zip(cutoffs, cutoffs[1:], strict=False):
        if whole and high - low == 1:
            labels.append(f"{low:g}")
        else:
            labels.append(f"{low:g}\u2013{high - 1 if whole else high:g}")
    labels.append(f"\u2265 {cutoffs[-1]:g}")
    return tuple(labels)


def bin_names(name: str, style: PedigreeStyle) -> tuple[str, ...]:
    """The steps' names (#81), "" for an unnamed step; () without any."""
    cutoffs = style.gradient_bins.get(name)
    names = style.gradient_bin_labels.get(name)
    if not cutoffs or not names or len(names) != len(cutoffs) + 1:
        return ()
    shown = tuple(n or "" for n in names)
    return shown if any(shown) else ()


BIN_LABEL_POSITIONS = ("below", "inside")


def bin_label_position(name: str, style: PedigreeStyle) -> str:
    """Where the legend writes a column's step names (#90): "below" or "inside"."""
    position = style.gradient_bin_label_position
    if isinstance(position, dict):
        return position.get(name, "below")
    return position


def bin_text(name: str, value: float, style: PedigreeStyle) -> str | None:
    """The step a value falls in, as the legend writes it, with its name if
    it has one ("\u2265 76, severe", #81), or None without bins."""
    cutoffs = style.gradient_bins.get(name)
    if not cutoffs:
        return None
    step = bin_index(value, cutoffs)
    names = bin_names(name, style)
    text = bin_labels(cutoffs)[step]
    return f"{text}, {names[step]}" if names and names[step] else text


def check_gradient_bins(style: PedigreeStyle) -> None:
    """Reject ``gradient_bins`` that are not lists of ascending numbers, and
    ``gradient_bin_labels`` that do not name each step of such a column (#81).

    Raises:
        PedigreeError: Naming each column whose cut-offs are not a non-empty
            list of strictly increasing numbers, or a column whose step names
            are not one text (or None) per step of its cut-offs.
    """
    bad = []
    for name, cutoffs in style.gradient_bins.items():
        numbers = isinstance(cutoffs, list | tuple) and all(
            isinstance(c, int | float) and not isinstance(c, bool) and math.isfinite(c)
            for c in cutoffs
        )
        if (
            not numbers
            or not cutoffs
            or any(b <= a for a, b in zip(cutoffs, cutoffs[1:], strict=False))
        ):
            bad.append(f"{name!r}: {cutoffs!r}")
    if bad:
        raise PedigreeError(
            "gradient_bins must map each column to ascending cut-offs, e.g. "
            f'{{"SRS": [60, 66, 76]}}; not for {", ".join(bad)}'
        )
    # Where the names go (#90): "below" or "inside", for all or by column.
    position = style.gradient_bin_label_position
    positions = list(position.values()) if isinstance(position, dict) else [position]
    if not all(p in BIN_LABEL_POSITIONS for p in positions):
        raise PedigreeError(
            'gradient_bin_label_position must be "below" or "inside", or a column -> '
            f"position mapping; got {position!r}"
        )
    # Step names (#81): one per step, for a column with cut-offs.
    for name, names in style.gradient_bin_labels.items():
        limits = style.gradient_bins.get(name)
        if limits is None:
            raise PedigreeError(
                f"gradient_bin_labels names the steps of {name!r}, which has no gradient_bins"
            )
        if not isinstance(names, list | tuple) or len(names) != len(limits) + 1:
            raise PedigreeError(
                f"gradient_bin_labels[{name!r}] must give one name per step: "
                f"{len(limits) + 1} for {len(limits)} cut-off(s), from the lowest step up "
                f'(None or "" for an unnamed step); got {names!r}'
            )
        if not all(n is None or isinstance(n, str) for n in names):
            raise PedigreeError(
                f"gradient_bin_labels[{name!r}] must hold text or None; got {names!r}"
            )


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


# The curve in its box: values grow upward, the bell's baseline this far
# (a fraction of the box's width) inside the box's right edge, its peak
# `PGS_CURVE_REACH` of the width to the left (#126). Turned a quarter turn
# clockwise it reads as the usual bell, peak up, positive values on the right,
# as the legend's curve; before 0.11.0 it was mirrored (baseline on the left).
PGS_CURVE_MARGIN = 0.1
PGS_CURVE_REACH = 0.8


def pgs_curve_base(box_left: float, box_width: float) -> float:
    """x of the PGS curve's baseline in its box: near the right edge (#126)."""
    return box_left + box_width * (1 - PGS_CURVE_MARGIN)


def pgs_curve_x(box_left: float, box_width: float, density: float) -> float:
    """x of the PGS curve at a density (0 to 1): from the baseline, leftward (#126)."""
    return pgs_curve_base(box_left, box_width) - density * box_width * PGS_CURVE_REACH


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
        asymptomatic: Whether some individual is an asymptomatic or
            presymptomatic carrier (#104): a vertical line through the symbol,
            which the legend explains and which quadrants would hide.
    """

    effective_modes: dict[str, str] = field(default_factory=dict)
    categorical_colors: dict[str, dict[str, str]] = field(default_factory=dict)
    quadrants: dict[str, str] = field(default_factory=dict)
    outline_rings: int = 0
    asymptomatic: bool = False

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
        if is_curve_column(name):  # a PGS or a pTDT column (#52)
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
        elif (
            name in style.gradient_diverging
            or name in style.gradient_range
            or name in style.gradient_bins
        ):
            # An explicit gradient configuration implies a continuous column.
            effective = "gradient"
        else:
            effective = detect_auto_mode(name, individuals)
        effective_modes[name] = effective

        if effective == "categorical":
            categorical_colors[name] = assign_categorical_colors(name, individuals, style)

    asymptomatic = any(ind.notation.is_asymptomatic for ind in individuals)
    return PhenotypeRender(effective_modes, categorical_colors, asymptomatic=asymptomatic)


# --- Names shown (#56) --------------------------------------------------------

LEGEND_MODES = ("auto", "compact", "full", "minimal")


def display_name(name: str, style: PedigreeStyle) -> str:
    """The name a column is shown with: its `phenotype_labels` entry, or its
    name with underscores as spaces."""
    label = style.phenotype_labels.get(name)
    return label if label is not None else name.replace("_", " ")


def pgs_trait_name(name: str, style: PedigreeStyle) -> str:
    """A PGS or pTDT column's name in the legend's key: its label, or its
    name without "PGS_" / "pTDT_" ("PGS_autism" -> "autism")."""
    label = style.phenotype_labels.get(name)
    if label is not None:
        return label
    trait = curve_trait(name)
    # Only a prefix followed by "_" goes: "PGSautism" keeps its name.
    trait = trait[1:] if trait.startswith("_") and len(trait) > 1 else name
    return trait.replace("_", " ")


def check_legend(style: PedigreeStyle) -> None:
    """Reject an unknown ``legend``, an invalid legend option or ``phenotype_labels``.

    Raises:
        PedigreeError: If ``style.legend`` is not one of `LEGEND_MODES`,
            ``legend_title`` is not text or None, ``legend_sex_symbols`` not a
            bool, ``legend_spacing`` outside `LEG_SPACING_RANGE`, or a label
            is not a string.
    """
    if style.legend not in LEGEND_MODES:
        raise PedigreeError(
            f"Unknown legend {style.legend!r}: expected one of {', '.join(LEGEND_MODES)}"
        )
    if style.legend_title is not None and not isinstance(style.legend_title, str):
        raise PedigreeError(f"legend_title must be text or null, not {style.legend_title!r}")
    if not isinstance(style.legend_sex_symbols, bool):
        raise PedigreeError(
            f"legend_sex_symbols must be true or false, not {style.legend_sex_symbols!r}"
        )
    low, high = LEG_SPACING_RANGE
    spacing = style.legend_spacing
    if (
        isinstance(spacing, bool)
        or not isinstance(spacing, int | float)
        or not low <= spacing <= high
    ):
        raise PedigreeError(
            f"legend_spacing must be a number from {low} to {high}, not {spacing!r}"
        )
    if style.legend_bar_width not in BAR_WIDTHS:
        raise PedigreeError(
            f'legend_bar_width must be "equal" or "own", not {style.legend_bar_width!r}'
        )
    width = style.legend_max_width
    if width is not None:
        valid = (
            isinstance(width, str)
            and re.fullmatch(r"\s*([0-9]*\.?[0-9]+)\s*(cm|mm|in)\s*", width) is not None
            and legend_max_width_inches(style) != 0
        ) or (not isinstance(width, bool | str) and isinstance(width, int | float) and width > 0)
        if not valid:
            raise PedigreeError(
                "legend_max_width must be a width in centimetres, or text with its unit "
                f'("85mm", "8.5cm", "3.3in"), not {width!r}'
            )
    unknown = [k for k in style.legend_text if k not in LEGEND_TEXT_DEFAULTS]
    if unknown:
        hints = []
        for key in unknown:
            close = get_close_matches(str(key), list(LEGEND_TEXT_DEFAULTS), n=1)
            hints.append(f"{key!r}" + (f" (did you mean {close[0]!r}?)" if close else ""))
        raise PedigreeError(
            f"Unknown legend_text key(s): {', '.join(hints)}; expected "
            + ", ".join(LEGEND_TEXT_DEFAULTS)
        )
    not_text = [k for k, v in style.legend_text.items() if not isinstance(v, str)]
    if not_text:
        raise PedigreeError(f"legend_text must map roles to text; not for {not_text!r}")
    bad = [k for k, v in style.phenotype_labels.items() if not isinstance(v, str)]
    if bad:
        raise PedigreeError(f"phenotype_labels must map names to text; not for {bad!r}")


# --- Bands -------------------------------------------------------------------


class _HasPhenotypes(Protocol):
    """A NodeLayout, as seen by the band decisions."""

    extra_phenotypes: dict[str, float | None]
    extra_phenotypes_raw: dict[str, float | None]
    extra_phenotypes_categorical: dict[str, str | None]
    not_applicable: set[str]


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
        not_applicable: Whether the column does not apply to the individual
            (N/A, #107): white, as unaffected, with no text.
    """

    color: Color
    text: str | None
    is_na: bool
    short: str | None = None
    not_applicable: bool = False

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
    if name in node.not_applicable:
        return BandCell(NOT_APPLICABLE_FILL, None, False, not_applicable=True)
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
        elif style.gradient_bins.get(name):
            cutoffs = style.gradient_bins[name]
            color = bin_color(rgb, bin_index(rv, cutoffs), cutoffs, inverted)
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


def band_has_no_data(node: _HasPhenotypes, name: str) -> bool:
    """Whether a symbol has no data for a column: -9 / NA, or no value at all.

    A 0 ("unknown") is a value, drawn grey without "NA". The parser stores a
    -9 as NaN, a 0 of a binary column as None, and every categorical value
    as None with its category in ``extra_phenotypes_categorical`` (None when
    missing); an individual a separate phenotype file does not list has no
    entry.
    """
    if name in node.not_applicable:
        return False  # "not applicable" is said, not missing (#107)
    if name not in node.extra_phenotypes:
        return True
    value = node.extra_phenotypes[name]
    if isinstance(value, float) and math.isnan(value):
        return True
    return name in node.extra_phenotypes_categorical and (
        node.extra_phenotypes_categorical[name] is None
    )


def all_bands_na(node: _HasPhenotypes, band_names: Sequence[str]) -> bool:
    """Whether a symbol has no data in any of its bands (`band_has_no_data`).

    Then a single centred "NA" is written instead of one per band. A symbol
    with a 0 ("unknown") or a category in one of them does not qualify. A
    band that does not apply (N/A, #107) is no data either, but keeps its
    stroke; a symbol whose bands are all N/A gets no "NA".
    """
    missing = [band_has_no_data(node, n) for n in band_names]
    applicable = [n not in node.not_applicable for n in band_names]
    return any(missing) and all(m or not a for m, a in zip(missing, applicable, strict=True))


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
    either. Under "auto", a single source (one band column, no carrier
    column) fills the whole symbol too, as NSGC draws one condition: in one
    quadrant, an affected individual looked three quarters unaffected (#108).
    """
    if not band_names:
        return replace(render, quadrants={})
    if style.phenotype_layout == "auto" and len(band_names) + len(carriers.keys) == 1:
        return replace(render, quadrants={})
    if render.asymptomatic and style.phenotype_layout != "bands":
        # The asymptomatic carrier's line runs where the quadrants meet (#104):
        # "auto" draws bands instead; "quadrants" keeps them, with a warning.
        if style.phenotype_layout == "auto":
            return replace(render, quadrants={})
        warnings.warn(
            "Some individuals are asymptomatic carriers: their vertical line runs along the "
            "line between quadrants. Use phenotype_layout 'bands' (or 'auto') to tell them apart.",
            PedigreeStyleWarning,
            stacklevel=2,
        )
    return replace(render, quadrants=assign_quadrants([*band_names, *carriers.keys], style))


def warn_asymptomatic_sections(individuals: Iterable[Any], carriers: CarrierScheme) -> None:
    """Warn when an asymptomatic carrier's line crosses carrier sections (#104).

    The NSGC 2022 carrier sections also split a symbol vertically; the line
    overshoots the symbol so it can still be told apart, but the reader
    should know.
    """
    both = [
        ind.sample_id
        for ind in individuals
        if ind.notation.is_asymptomatic and not carriers.legacy and carriers.marks(ind)
    ]
    if both:
        shown = ", ".join(both[:5]) + (", ..." if len(both) > 5 else "")
        warnings.warn(
            f"{len(both)} asymptomatic carrier(s) ({shown}) also have carrier sections: the "
            "vertical line is drawn across them.",
            PedigreeDataWarning,
            stacklevel=2,
        )


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
    if cell.is_na or cell.not_applicable or render.effective(name) != "gradient":
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
    """The text a band or quadrant shows, its size in points, and how far
    right of the place asked for it goes (layout units), to keep clear of
    the symbol's letter (#105)."""

    text: str
    size: float
    dx: float = 0.0


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
    *,
    letter: float | LetterBox | None = None,
    offset: tuple[float, float] = (0.0, 0.0),
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
        letter: The font size of the letter written at the symbol's centre
            ("P", "D", "S"), in points, if it has one (#105). A text that
            would run into it goes right of it (a band's, centred) or gets
            less room (a quadrant's).
        offset: Where the text is centred, from the symbol's centre, in
            layout units (only needed with ``letter``).

    Returns:
        The text and its size, or None when nothing is written.
    """
    if not texts:
        return None

    def place(h: float) -> tuple[float, float]:
        """Room for a text ``h`` high, and how far right it moves."""
        return clear_of_letter(room(h), h, offset, letter, pt)

    if style.band_text != "auto":
        return FittedText(texts[0], size, place(size * pt * _FIT_HEIGHT_EM)[1])
    if style.min_font_size is not None:
        floor = min(size, style.min_font_size)
    else:
        floor = min(size, style.font_size * FIT_FLOOR_SCALE)

    def fits(text: str, s: float) -> bool:
        h = s * pt * _FIT_HEIGHT_EM
        width = symbols.text_width(text, s, bold=True) / symbols.UNITS_PER_POINT
        return h <= height and (width + 2 * _FIT_PAD_EM * s) * pt <= place(h)[0]

    for text in texts:
        if fits(text, size):
            return FittedText(text, size, place(size * pt * _FIT_HEIGHT_EM)[1])
    text = texts[-1]
    if not fits(text, floor):
        return None
    low, high = floor, size
    for _ in range(20):
        mid = (low + high) / 2
        low, high = (mid, high) if fits(text, mid) else (low, mid)
    return FittedText(text, low, place(low * pt * _FIT_HEIGHT_EM)[1])


@dataclass(frozen=True)
class LetterBox:
    """The text at a symbol's centre that band text keeps clear of (#105):
    its font size in points, and how many characters it has ("P", "12")."""

    size: float
    chars: int = 1


def letter_box(notation: Any, size: float) -> LetterBox | None:
    """The box of a symbol's inner text (`symbols.inner_letter`), or None."""
    text = symbols.inner_letter(notation)
    return LetterBox(size, len(text)) if text else None


def clear_of_letter(
    width: float,
    h: float,
    offset: tuple[float, float],
    letter: float | LetterBox | None,
    pt: float,
) -> tuple[float, float]:
    """Room for a band's or quadrant's text, kept clear of the inner letter (#105).

    The letter ("P", "D", "S") is at the symbol's centre with its backing,
    `symbols.LETTER_CLEAR_HALF_W_EM` by `LETTER_CLEAR_HALF_H_EM` ems of its
    font. A text whose line reaches that height goes, if centred on the
    symbol (a band), right of the letter, centred in what is left of the
    band; if centred off the axis (a quadrant), it keeps its place with
    less room.

    Args:
        width: Room the text has without the letter, centred on ``offset``.
        h: The text's height, in layout units.
        offset: The text's centre from the symbol's centre.
        letter: The letter's font size in points, or None without a letter.
        pt: Layout units per point at the drawing's scale.

    Returns:
        The room, and how far right of ``offset`` the text goes.
    """
    if letter is None:
        return width, 0.0
    box = letter if isinstance(letter, LetterBox) else LetterBox(letter)
    ox, oy = offset
    em = box.size * pt
    half_w = (symbols.LETTER_CLEAR_HALF_W_EM + symbols.LETTER_CLEAR_CHAR_EM * (box.chars - 1)) * em
    half_h = symbols.LETTER_CLEAR_HALF_H_EM * em
    gap = symbols.LETTER_CLEAR_GAP_EM * em
    if abs(oy) >= h / 2 + half_h + gap:
        return width, 0.0
    if abs(ox) < 1e-9:
        side = max(0.0, width / 2 - half_w - gap)
        return side, half_w + gap + side / 2
    return max(0.0, min(width, 2 * (abs(ox) - half_w - gap))), 0.0


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

    Converted with `symbols.units_per_point` of the style, the scale of the
    SVG; the PNG draws at that scale or larger (see `renderer._OutlineRings`).
    """
    offsets = outline_ring_offsets(count, style)
    return offsets[-1] * symbols.units_per_point(style.text_scale) if offsets else 0.0


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
    just under a curve on that side, where the caption would otherwise be
    written; a curve on the right (#78) keeps its caption below.
    """
    return (node.is_proband or node.notation.is_consultand) and node.pgs_side != "right"


def pgs_box_left(node: NodeLayout, style: PedigreeStyle, ring: float = 0.0) -> float:
    """Left edge of a node's PGS box, in layout units.

    Left of the symbol, or right of it for a right-hand partner
    (`NodeLayout.pgs_side`, #78), `symbols.pgs_box_near` from its centre,
    plus ``ring``, the room of its outer borders (#46). Both renderers
    place the box with this.
    """
    half = style.shape_size / 2
    mark = style.text_size(symbols.EVALUATION_MARK_FONT_SCALE)
    near = symbols.pgs_box_near(
        half,
        node.sex,
        node.notation.is_evaluated,
        mark,
        node.pgs_side,
        struck=node.is_deceased,
        ring=ring,
        line_width=style.shape_line_width * (OUTLINE_WIDTH_SCALE if ring else 1.0),
        bracket_width=style.line_width if node.notation.adoption is not None else None,
        upp=symbols.units_per_point(style.text_scale),
    )
    if node.pgs_side == "right":
        return node.x + near
    return node.x - near - box_width(style)


def pgs_values_x(node: NodeLayout, box_left: float, box_w: float, gap: float) -> tuple[float, str]:
    """Where a curve's values are written (#54): x and "left" / "right", the
    side of the text at x. On the far side of the box from the symbol: left
    of a box on the left, right of a box on the right (#78). ``gap`` is the
    space from the box, in layout units."""
    if node.pgs_side == "right":
        return box_left + box_w + PGS_BOX_PAD + gap, "left"
    return box_left - PGS_BOX_PAD - gap, "right"


# Between a curve's box and its "PGS" caption: a fraction of the symbol's
# half-size, or more when the caption needs it (#66): room for a marker dot
# on the box's edge, then this many ems of the caption.
PGS_CAPTION_OFFSET = 0.15
PGS_CAPTION_GAP_EM = 0.3


def pgs_caption_offset(half: float, size: float, units_per_point: float) -> float:
    """How far the "PGS" caption is from its curve's box, in layout units.

    The larger of `PGS_CAPTION_OFFSET` times the symbol's half-size and,
    converted with ``units_per_point``, half a marker dot (a score at either
    end of the curve puts its dot on the box's edge) plus `PGS_CAPTION_GAP_EM`
    of the caption's size. At the default sizes the first is larger, so
    nothing moves; a caption enlarged by ``min_font_size`` keeps its gap.

    Args:
        half: Half the symbol's size, in layout units.
        size: The caption's font size, in points.
        units_per_point: Layout units per point at the scale drawn
            (`symbols.UNITS_PER_POINT` in the SVG).
    """
    needed = (PGS_DOT_PT / 2 + PGS_CAPTION_GAP_EM * size) * units_per_point
    return max(half * PGS_CAPTION_OFFSET, needed)


def column_scale(name: str, scale: str) -> str:
    """The scale a curve column is placed on: ``pgs_scale``, or "zscore" for
    a pTDT column, whose deviations are N(0, 1) by construction (#52)."""
    return "zscore" if is_ptdt_column(name) else scale


PTDT_CAPTION = "pTDT"


def curve_caption(active: Sequence[tuple[str, float]]) -> str:
    """The caption of an individual's curve: "pTDT" when it shows pTDT
    deviations (`active_pgs`), "PGS" otherwise."""
    return PTDT_CAPTION if active and is_ptdt_column(active[0][0]) else "PGS"


def active_pgs(
    node: _HasPhenotypes, pgs_names: Sequence[str], scale: str = "minmax"
) -> list[tuple[str, float]]:
    """Curve columns with a value for this individual, with their [0, 1] position.

    On "minmax" the position is the value normalised over the pedigree; on
    "zscore" and "percentile" it is the value's z-score on the drawn
    normal curve (`pgs_position`), the same for every pedigree. A pTDT
    column is always on the z scale (#52).

    An individual's curve shows one kind of score: its pTDT deviations if
    it has any, otherwise its PGS (`warn_mixed_curves` reports individuals
    with both), so the parents' PGS and the children's pTDT can share a
    figure, each curve captioned by `curve_caption`.
    """
    ptdt = [n for n in pgs_names if is_ptdt_column(n)]
    if ptdt:
        shown = _active_curves(node, ptdt, scale)
        if shown:
            return shown
    return _active_curves(node, [n for n in pgs_names if not is_ptdt_column(n)], scale)


def warn_mixed_curves(
    individuals: Iterable[_HasPhenotypes],
    names: Sequence[str],
    computed: Collection[str] = (),
) -> None:
    """Warn when some individual has both PGS and pTDT values (#52): its
    curve shows the pTDT deviations only. Columns pypedigree computed from
    the parents (`computed`, ``pgs_transmission``) are left out: every child
    has both, by design."""
    ptdt = [n for n in names if is_ptdt_column(n) and n not in computed]
    pgs = [n for n in names if is_curve_column(n) and not is_ptdt_column(n)]
    if not ptdt or not pgs:
        return

    def has(ind: _HasPhenotypes, cols: list[str]) -> bool:
        raw = ind.extra_phenotypes_raw
        return any(n in raw and not is_na_value(raw[n]) for n in cols)

    both = [getattr(i, "sample_id", "?") for i in individuals if has(i, ptdt) and has(i, pgs)]
    if both:
        shown = ", ".join(both[:5]) + (", ..." if len(both) > 5 else "")
        warnings.warn(
            f"{len(both)} individual(s) have both PGS and pTDT values ({shown}): their curve "
            "shows the pTDT deviations only. Give the PGS to the parents and the pTDT "
            "deviation to the children.",
            PedigreeDataWarning,
            stacklevel=2,
        )


# On "zscore" (the default since 0.10.0, #106), a PGS column whose values do
# not look like z-scores draws a warning: any value beyond this, ...
ZSCORE_MAX_ABS = 5.0
# ... every value within this of 0, with at least this many values not all 0
# (raw sums such as PLINK's, around 1e-4; #117: checked in a trio too), ...
ZSCORE_TINY = 0.05
ZSCORE_TINY_MIN_VALUES = 2
# ... values spread over less than this standard deviation (equal values say
# nothing), with at least this many values, ...
ZSCORE_MIN_SD = 0.05
ZSCORE_SD_MIN_VALUES = 5
# ... or a mean this far from 0 over at least this many values (a quartet of
# raw scores around 3 drew from +2.4 to +3.6 SD without a word, #117).
ZSCORE_MAX_MEAN = 2.0
ZSCORE_MEAN_MIN_VALUES = 3


def warn_pgs_scale(individuals: Iterable[_HasPhenotypes], names: Sequence[str], scale: str) -> None:
    """Warn when a PGS column drawn on "zscore" does not look like z-scores (#106, #117).

    On "zscore" a value is placed on the standard normal and clipped to the
    curve's ends, so raw scores would all land at one end (or, if tiny, at
    the centre) without a word. pTDT columns are always z-scores and are not
    checked. Callers pass the whole file's individuals (`Pedigree.column_scope`),
    so every family drawn from it gets the same message; duplicates (a person
    drawn twice) are counted once.
    """
    if scale != "zscore":
        return
    people = [i for i in individuals if not getattr(i, "is_duplicate", False)]
    odd = []
    for name in names:
        if not is_pgs_column(name) or is_ptdt_column(name):
            continue
        values = [
            float(raw)
            for ind in people
            if (raw := ind.extra_phenotypes_raw.get(name)) is not None and not is_na_value(raw)
        ]
        if not values:
            continue
        mean = sum(values) / len(values)
        sd = math.sqrt(sum((v - mean) ** 2 for v in values) / len(values))
        largest = max(values, key=abs)
        if abs(largest) > ZSCORE_MAX_ABS:
            odd.append(f"{name} (a value of {largest:g})")
        elif len(values) >= ZSCORE_TINY_MIN_VALUES and 0 < abs(largest) < ZSCORE_TINY:
            odd.append(f"{name} (every value within ±{ZSCORE_TINY:g})")
        elif len(values) >= ZSCORE_SD_MIN_VALUES and 0 < sd < ZSCORE_MIN_SD:
            odd.append(f"{name} (spread {sd:.2g})")
        elif len(values) >= ZSCORE_MEAN_MIN_VALUES and abs(mean) > ZSCORE_MAX_MEAN:
            odd.append(f"{name} (mean {mean:.2f})")
    if odd:
        warnings.warn(
            f"PGS column(s) {', '.join(odd)} do not look like z-scores, but are drawn on the "
            "z-score scale (pgs_scale 'zscore', the default). For raw scores, use pgs_scale "
            "'minmax' (--pgs-scale minmax).",
            PedigreeDataWarning,
            stacklevel=2,
        )


def _active_curves(
    node: _HasPhenotypes, pgs_names: Sequence[str], scale: str
) -> list[tuple[str, float]]:
    out = []
    for name in pgs_names:
        if column_scale(name, scale) == "minmax":
            value = node.extra_phenotypes.get(name)
            if not is_na_value(value):
                out.append((name, value or 0.0))
            continue
        raw = node.extra_phenotypes_raw.get(name)
        if not is_na_value(raw) and raw is not None:
            out.append((name, pgs_position(pgs_zscore(raw, column_scale(name, scale)))))
    return out


def _ordinal(n: int) -> str:
    suffix = "th" if 10 <= n % 100 <= 20 else {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th")
    return f"{n}{suffix}"


def curve_value_text(name: str, raw: float, scale: str) -> str:
    """A curve column's value as the tooltip writes it: `pgs_value_text`, or
    "pTDT = +2.12" for a pTDT deviation, which has no percentile (#52)."""
    if is_ptdt_column(name):
        return f"{PTDT_CAPTION} = {raw:+.2f}".replace("-", MINUS)
    return pgs_value_text(raw, scale)


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
# A score's line across its box (#158): solid, as thin as the grey curve, so it
# reads as one value (dashes looked like an interval); the dot gives its color.
PGS_MARKER_LINE_PT = 0.6
PGS_DOT_SPREAD = 1.5  # distance between the centres of dots set side by side, in diameters
PGS_DOT_MIN_PT = 1.5  # the smallest a dot is drawn to fit its box (#147), in points
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
            out[name] = pgs_marker_text(raw, column_scale(name, scale))
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
    upp = symbols.units_per_point(style.text_scale)
    widths = [
        symbols.text_width(text, size, upp=upp)
        for ind in nodes
        for text in pgs_marker_values(
            ind, [n for n, _ in active_pgs(ind, pgs_names, style.pgs_scale)], style.pgs_scale
        ).values()
    ]
    if not widths:
        return 0.0
    return max(widths) + PGS_VALUE_GAP_EM * size * upp + PGS_BOX_PAD


def pgs_dot_count(pgs_names: Sequence[str]) -> int:
    """The most dots one individual's box can hold: the columns of the more
    numerous kind, PGS or pTDT (`active_pgs` gives one kind per individual)."""
    ptdt = sum(1 for n in pgs_names if is_ptdt_column(n))
    return max(ptdt, len(pgs_names) - ptdt)


def pgs_box_width(pgs_names: Sequence[str], style: PedigreeStyle) -> tuple[float, float]:
    """Width of the PGS boxes, and how much wider than usual, in layout units (#147).

    `symbols.PGS_BOX_WIDTH` half-sizes, or wider to hold all of an
    individual's dots side by side at full size, `PGS_DOT_SPREAD` diameters
    apart (`spread_dots`), at the style's text scale
    (`symbols.units_per_point`): with several scores the usual box made them
    cover each other when drawn small. Each renderer passes it to
    `config.apply_pgs_spacing`, which keeps it in the style (`box_width`).
    """
    usual = style.shape_size / 2 * symbols.PGS_BOX_WIDTH
    count = pgs_dot_count(pgs_names)
    if count < 2:
        return usual, 0.0
    dot = PGS_DOT_PT * symbols.units_per_point(style.text_scale)
    width = max(usual, dot * (1 + (count - 1) * PGS_DOT_SPREAD))
    return width, width - usual


def box_width(style: PedigreeStyle) -> float:
    """Width of the PGS boxes, in layout units: set by `apply_pgs_spacing`
    when wider than usual (#147), else `symbols.PGS_BOX_WIDTH` half-sizes."""
    return style.pgs_box_width or style.shape_size / 2 * symbols.PGS_BOX_WIDTH


def pgs_dot_size(count: int, width: float, pt: float) -> float:
    """Diameter of an individual's PGS dots, in layout units (#147).

    `PGS_DOT_PT`, or smaller when the box (inner ``width``) cannot hold all
    ``count`` of the individual's dots side by side, `PGS_DOT_SPREAD`
    diameters apart, as `spread_dots` sets them when their scores coincide:
    drawn small, full-size dots covered each other. Never below
    `PGS_DOT_MIN_PT`. Both renderers size the dots with this at the scale
    they draw at (``pt`` layout units per point).

    Args:
        count: How many dots the individual has (one per score).
        width: The box's inner width, in layout units.
        pt: Layout units per point at the drawing's scale.
    """
    full = PGS_DOT_PT * pt
    if count < 2:
        return full
    fit = width / (1 + (count - 1) * PGS_DOT_SPREAD)
    return min(full, max(PGS_DOT_MIN_PT * pt, fit))


def spread_dots(ys: Sequence[float], dot: float, width: float | None = None) -> list[float]:
    """Horizontal offsets that keep the markers' dots apart (#54).

    Dots closer than ``dot`` (their diameter) vertically form a group, set
    side by side, centred on the curve's axis, `PGS_DOT_SPREAD` diameters
    apart, or closer when the group would pass the box (#89): within
    ``width``, its dots may then overlap a little, each keeping its color
    and its line. A dot alone keeps offset 0.

    Args:
        ys: Height of each dot, in any unit.
        dot: The dots' diameter, in the same unit.
        width: The box's inner width, in the same unit, that the dots stay
            within; None for no limit.

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
        step = dot * PGS_DOT_SPREAD
        if width is not None and len(group) > 1:
            step = min(step, max(0.0, (width - dot) / (len(group) - 1)))
        for k, j in enumerate(group):
            offsets[j] = (k - (len(group) - 1) / 2) * step
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
LEG_SEX_SYMBOL_GAP = 0.25  # between a sex symbol and its "= Male" text
LEG_SEX_ITEM_GAP = 0.6  # between "= Male" and the female symbol, when the text is long
LEG_OUTLINE_RINGS_LABEL = "Several borders: one per column, inner first"
LEG_RING_WIDTH_PT = 1.6  # each border of the "outline_rings" icon
LEG_RING_INSET_PT = 2.2  # its inner border, inside the outer one
LEG_KEY_DASH_EM = 1.4  # a PGS score's line in the legend's key (dot in its middle), in ems
LEG_KEY_GAP_EM = 0.3  # between a dash and its score's name
LEG_KEY_SPACE_EM = 0.9  # between one score's name and the next dash
# The PGS row's icon (#159): the box drawn beside a symbol, upright as in the
# figure, in ems of the legend's font, with "high" above, "low" below and the
# ±2 SD lines named on its right, `LEG_PGS_TEXT_GAP_EM` (ems of that text) away.
LEG_PGS_BOX_W_EM = 1.45
LEG_PGS_BOX_H_EM = 3.4
LEG_PGS_TEXT_GAP_EM = 0.3
LEG_PGS_ICON_GAP = 0.3  # between the icon's labels and the row's text, in legend units
# The header (#67), before `legend_spacing`: the margin above it (smaller with
# neither title nor sex symbols), the title's line, and the sex symbols' row
# (above its centre when it comes first, below it).
LEG_TOP_BARE = 0.5
LEG_TITLE_H = 1.0
LEG_SEX_ABOVE = 0.4
LEG_SEX_BELOW = 0.8
LEG_DEFAULT_TITLE = "Legend"
LEG_MINIMAL_SPACING = 0.7  # `legend_spacing` with legend "minimal", unless set
LEG_SPACING_RANGE = (0.6, 3.0)  # below 0.6 a gradient's end captions reach the next row
# "minimal" (#67): the binary columns of the quadrant key are read from this
# caption, in a row under the key, instead of a row of swatches.
LEG_KEY_NOTE_LABEL = "Colored: affected \u00b7 white: unaffected \u00b7 grey: unknown"
# A row is never shorter than its swatch, the gap to the captions under it
# (a gradient's ends, category names) and those captions, plus this margin,
# so a small `legend_spacing` cannot make them reach the next row.
LEG_CAPTION_GAP = 0.05
LEG_ROW_MARGIN = 0.08
LEG_LINE_EM = 1.17  # a line of text's height, in ems (DejaVu Sans)
# Below the last row's captions (#77): the frame's inset from the panel's edge,
# then a gap, so they end inside the frame.
LEG_FRAME_INSET = 0.1
LEG_BOTTOM_PAD = 0.1
CAPTIONED_ROW_KINDS = ("gradient", "binned", "categorical")  # captions under the swatch
# Rows that write their name after a swatch or a bar, in one column (#88).
NAMED_ROW_KINDS = ("gradient", "binned", "categorical", "binary", "binary_shared")
LEG_WRAP_LINE_EM = 1.25  # between the lines of a wrapped PGS row (#68), in ems of its text
LEG_WIDTH_UNITS = {"cm": 1 / 2.54, "mm": 1 / 25.4, "in": 1.0}  # `legend_max_width` -> inches
KEY_ROW_KINDS = ("quadrants", "key_note")  # drawn first, above the binary captions

# Words a wrapped legend row keeps together (#68): joined by a no-break space.
NBSP = "\u00a0"

# The legend's fixed wording, by role (#79): ``legend_text`` overrides any of
# it. "pgs" and "ptdt" have no fixed default: their caption is built from the
# scale (`pgs_legend_description`, `ptdt_legend_description`).
LEGEND_TEXT_DEFAULTS: dict[str, str] = {
    "key_note": LEG_KEY_NOTE_LABEL,
    "shared_binary": LEG_SHARED_BINARY_LABEL,
    "affected": "Affected",
    "unaffected": "Unaffected",
    "unknown": "Unknown",
    "pgs": "",
    "ptdt": "",
    "male": "Male",
    "female": "Female",
    "low": "Low",
    "high": "High",
    # Beside the PGS / pTDT icon (#159).
    "pgs_high": "high",
    "pgs_low": "low",
    "pgs_plus_2sd": f"+2{NBSP}SD",
    "pgs_minus_2sd": f"{MINUS}2{NBSP}SD",
    # Still accepted, unused since 0.10.1: N/A is drawn white with no legend row (#120).
    "not_applicable": "Not applicable",
    "asymptomatic": "Asymptomatic or presymptomatic carrier",
}


def legend_text(style: PedigreeStyle, role: str) -> str:
    """The legend's wording for ``role``: ``legend_text``, else its default."""
    return style.legend_text.get(role) or LEGEND_TEXT_DEFAULTS[role]


def legend_binary_captions(style: PedigreeStyle) -> tuple[str, str, str]:
    """The Affected / Unaffected / Unknown captions above the swatch rows."""
    return (
        legend_text(style, "affected"),
        legend_text(style, "unaffected"),
        legend_text(style, "unknown"),
    )


def legend_sex_items(style: PedigreeStyle, m: LegendMetrics) -> list[tuple[float, str, str]]:
    """The sex row: (symbol x, "square" / "circle", text) for male, then female.

    The female symbol keeps its place unless a longer "male" text (#79)
    needs the room; the text follows its symbol at `LEG_SEX_SYMBOL_GAP`.
    """
    male = f"= {legend_text(style, 'male')}"
    female = f"= {legend_text(style, 'female')}"
    male_end = 0.6 + LEG_SEX_SYMBOL_GAP + _legend_text_width(male, m.fs, m)
    return [(0.6, "square", male), (max(2.5, male_end + LEG_SEX_ITEM_GAP), "circle", female)]


def legend_min_right(style: PedigreeStyle) -> float:
    """The least width of the legend's content, in legend units: the frame
    hugs a "minimal" legend's content (#79), others keep a minimum."""
    return LEG_X0 if style.legend == "minimal" else LEG_X0 + 4.0


def legend_title(style: PedigreeStyle) -> str:
    """The legend's title ("" for none): ``legend_title``, none for "minimal" by default."""
    title = style.legend_title
    if style.legend == "minimal" and title == LEG_DEFAULT_TITLE:
        return ""
    return title or ""


def legend_sex_symbols(style: PedigreeStyle) -> bool:
    """Whether the legend has its "= Male / = Female" row (not by default with "minimal")."""
    return style.legend_sex_symbols and style.legend != "minimal"


def legend_spacing(style: PedigreeStyle) -> float:
    """The legend's spacing factor: ``legend_spacing``, 0.7 for "minimal" by default."""
    if style.legend == "minimal" and style.legend_spacing == 1.0:
        return LEG_MINIMAL_SPACING
    return float(style.legend_spacing)


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
        spacing: Factor on the row height and the gaps (``legend_spacing``).
        title: The title ("" for none).
        title_y: Top of the title.
        sex_y: Centre of the "= Male / = Female" row, or None without it.
        content_top: Top of the first row (the quadrant key's, or the binary
            captions above the swatch rows), below the header.
        binary_header_h: Height of the Affected / Unaffected / Unknown captions.
        line_h: Distance between two lines of a wrapped row (#68).
        bar_w: The width every gradient and binned bar takes, the widest
            one's (``legend_bar_width`` "equal", #95), or 0: each its own.
            Set by `legend_layout`.
        name_x: Where the rows with a swatch or a bar (binary, gradient,
            binned, categorical) write their name: `text_x`, or past the
            widest binned bar, one column for the whole legend (#88). Set by
            `legend_layout`.
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
    spacing: float = 1.0
    title: str = LEG_DEFAULT_TITLE
    title_y: float = 0.0
    sex_y: float | None = None
    content_top: float = 0.0
    binary_header_h: float = LEG_BINARY_HEADER_H
    line_h: float = 0.5
    name_x: float = 0.0
    bar_w: float = 0.0

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
    n_pheno_rows: int,
    style: PedigreeStyle,
    *,
    binary_headers: bool = False,
    extra_lines: int = 0,
    bottom_captions: int = 0,
    bottom_lines: int = 1,
    extra_h: float = 0.0,
) -> LegendMetrics:
    """Compute the legend panel geometry for a given number of rows.

    Rows keep a constant height: the panel grows with its row count, and
    the figure grows with the panel (see `_add_legend_axes`). The header
    holds the title and the sex symbols, each optional (#67); the row
    height and every gap are multiplied by the spacing factor.

    Args:
        n_pheno_rows: Number of phenotype and carrier rows to fit
            (`legend_slot_count`).
        style: Style supplying the base font size and family.
        binary_headers: Whether the Affected / Unaffected / Unknown captions
            are drawn above the rows.
        extra_lines: Lines added by wrapped rows, beyond their first (#68).
        bottom_captions: Lines of captions under the last row's swatch
            (`legend_caption_lines`); the panel grows below it by as much as
            they would pass the frame (#77).
        bottom_lines: Lines the last row takes (its extra lines already hold
            part of its captions, #81).
        extra_h: Room rows take beyond their slots (`LegendRow.extra_h`,
            #90, #94); the part below holds the same part of their captions.

    Returns:
        The derived metrics.
    """
    n_rows = max(n_pheno_rows, 1)
    spacing = legend_spacing(style)
    sw_h = LEG_ROW_H * spacing * 0.45
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
    label_fs = max(fs * 0.75, min_fs)
    # Never shorter than a swatch and its captions (only reached below 1.0).
    caption_h = LEG_LINE_EM * label_fs / (LEG_UNIT_Y_EM * font)
    row_h = max(LEG_ROW_H * spacing, sw_h + LEG_CAPTION_GAP + caption_h + LEG_ROW_MARGIN)

    # The header, from the top (y up from the panel's bottom, so its height
    # is counted first and placed once the total is known).
    title = legend_title(style)
    sex = legend_sex_symbols(style)
    top = (LEG_PAD if title or sex else LEG_TOP_BARE) * spacing
    title_h = LEG_TITLE_H * spacing if title else 0.0
    sex_above = LEG_SEX_ABOVE * spacing if sex and not title else 0.0
    sex_below = LEG_SEX_BELOW * spacing if sex else 0.0
    header_h = top + title_h + sex_above + sex_below
    binary_header_h = LEG_BINARY_HEADER_H * spacing
    line_h = LEG_WRAP_LINE_EM * fs / (LEG_UNIT_Y_EM * font)
    # The last row's captions hang below its slot: room for them inside the
    # frame, only when they need it, so other legends keep their height.
    hang = (
        sw_h / 2
        + LEG_CAPTION_GAP
        + (caption_h + (bottom_captions - 1) * line_h if bottom_captions else 0.0)
        - row_h / 2
        - (bottom_lines - 1) * line_h
    )
    bottom = max(0.0, hang + LEG_FRAME_INSET + LEG_BOTTOM_PAD) if bottom_captions else 0.0
    total_h = (
        header_h
        + (binary_header_h if binary_headers else 0.0)
        + n_rows * row_h
        + extra_lines * line_h
        + extra_h
        + bottom
    )
    title_y = total_h - top
    sex_y = title_y - title_h - sex_above if sex else None

    return LegendMetrics(
        n_rows=n_rows,
        row_h=row_h,
        total_h=total_h,
        scale=font / LEG_REFERENCE_FONT,
        sw_h=sw_h,
        sw_w=sw_w,
        sw_gap=sw_gap,
        text_x=text_x,
        icon_text_x=icon_text_x,
        fs=fs,
        label_fs=label_fs,
        marker_size=font,
        unit_x=LEG_UNIT_X_EM * em,
        unit_y=LEG_UNIT_Y_EM * em,
        font=font,
        title_fs=max(fs + 1, font * 0.7),
        min_fs=min_fs,
        spacing=spacing,
        title=title,
        title_y=title_y,
        sex_y=sex_y,
        content_top=total_h - header_h,
        binary_header_h=binary_header_h,
        line_h=line_h,
        name_x=text_x,
    )


@dataclass(frozen=True)
class LegendRow:
    """One row of the legend panel.

    Attributes:
        kind: "binary", "binary_shared" (the binary columns of the quadrant
            key, compact legend), "gradient", "binned" (a gradient cut into
            steps by ``gradient_bins``, #57), "categorical", "outline",
            "outline_rings" (how to read several borders, #46), "pgs",
            "legacy_carrier" (half-fill / dot) or "carrier" (patterned).
        label: Row text.
        color: Swatch color ("binary", "outline").
        colors: Gradient bar, left to right ("gradient"); the columns'
            colors, side by side in the "affected" swatch ("binary_shared");
            the two borders of the icon, inner first ("outline_rings").
        low: Caption under the left end of a gradient bar; below the icon
            of a "pgs" row (#159).
        high: Caption under the right end of a gradient bar; above the icon
            of a "pgs" row.
        categories: Category value → color ("categorical").
        pgs: PGS column → marker color ("pgs").
        pattern: Carrier pattern name ("carrier").
        obligate: Obligate-carrier row ("carrier", "legacy_carrier").
        quadrants: Cells of the quadrant key ("quadrants"), one per column.
        key: The PGS scores' names and colors, written after the label
            ("pgs", #56).
        single_symbol: The quadrant key shows one symbol, not a square and
            a circle ("quadrants", legend "minimal", #67).
        bins: The steps' limits ("binned", clinical thresholds, #57), under
            the swatches of `colors`.
        bin_names: The steps' names, under their limits ("" unnamed, () when
            none is named, #81).
        names_inside: The names are written in the steps, not under the
            limits (``gradient_bin_label_position`` "inside", #90).
        extra_h: Room the row takes beyond its slot, in legend units: what a
            bar taller than a swatch adds (#90), half above it and half below
            (#94). Set by `legend_layout`.
        wrap_right: Where a "pgs" row's text must end, in legend units, or
            None: one line whatever its length (#68). Set by `legend_layout`.
        lines: Lines of text the row takes (more than one when wrapped).
        tails: The ±2 SD lines of a "pgs" row's icon: (position on the
            curve's [0, 1] axis, name), none on the "minmax" scale (#159).

    A "key_note" row is a caption under the quadrant key, saying how to read
    its binary columns (legend "minimal").
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
    single_symbol: bool = False
    bins: tuple[str, ...] = ()
    bin_names: tuple[str, ...] = ()
    names_inside: bool = False
    extra_h: float = 0.0
    wrap_right: float | None = None
    lines: int = 1
    tails: tuple[tuple[float, str], ...] = ()

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

    The key is a square and a circle, each cut in four, or the square alone
    (`LegendRow.single_symbol`); left-quadrant labels sit left of the square,
    right-quadrant labels right of the last symbol.

    Attributes:
        square: (x, y_center, width, height) of the square.
        circle: (x_center, y_center, x_radius, y_radius) of the circle; the
            radii differ in legend units so it is round on the page. None
            when the key has one symbol.
        labels: (x, y, text, anchor) per label; anchor is "end" or "start".
    """

    square: tuple[float, float, float, float]
    circle: tuple[float, float, float, float] | None
    labels: tuple[tuple[float, float, str, str], ...]


def legend_row_centers(rows: Sequence[LegendRow], m: LegendMetrics, y_top: float) -> list[float]:
    """Vertical centre of each row, counting down from `y_top`.

    A row of several lines (a wrapped PGS row, #68) has the centre of its
    first line, where its icon goes; its other lines follow `m.line_h` apart.
    A row with `LegendRow.extra_h` (a bar taller than a swatch, #90) takes it
    half above and half below: its centre moves down by half (#94), so the
    bar keeps the same gap from the row above as a swatch does.
    """
    centers = []
    y = y_top
    for row in rows:
        first = m.row_h / 2 if row.lines > 1 else row.slots * m.row_h / 2
        centers.append(y - row.extra_h / 2 - first)
        y -= row.slots * m.row_h + (row.lines - 1) * m.line_h + row.extra_h
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
    right_edge = square_x + side_x if row.single_symbol else circle_cx + side_x / 2
    labels = []
    for cell in row.quadrants:
        v = y_center + (side_y / 4 if cell.position.startswith("top") else -side_y / 4)
        if cell.position.endswith("left"):
            labels.append((square_x - LEG_KEY_GAP, v, cell.label, "end"))
        else:
            labels.append((right_edge + LEG_KEY_GAP, v, cell.label, "start"))
    return QuadrantKeyGeometry(
        square=(square_x, y_center, side_x, side_y),
        circle=None if row.single_symbol else (circle_cx, y_center, side_x / 2, side_y / 2),
        labels=tuple(labels),
    )


@dataclass(frozen=True)
class LegendKeyItem:
    """One score of the PGS row's key: a line of its color with a dot in its
    middle, as the score is marked in the figure (#158), then its name.

    Positions are in legend units (see `LEG_UNIT_X_EM`).
    """

    dash_x0: float
    dash_x1: float
    text_x: float
    name: str
    color: str
    line: int = 0  # the row's line it is on, 0 for the first (#68)


@dataclass(frozen=True)
class PgsRowLayout:
    """The text of a PGS row: its caption's lines, then its key (#68).

    Attributes:
        text_x: Where its text starts: past the icon's labels (#159).
        captions: The caption, one entry per line (one when it fits).
        items: The key's scores, each with the line it is on.
        lines: Lines the row takes.
    """

    text_x: float
    captions: tuple[str, ...]
    items: tuple[LegendKeyItem, ...]
    lines: int


def _legend_text_width(text: str, size: float, m: LegendMetrics) -> float:
    """Width of `text` at `size` points, in horizontal legend units (DejaVu Sans)."""
    return symbols.text_width(text, size) / symbols.UNITS_PER_POINT / 72 / m.unit_x


def _wrap_words(text: str, x: float, right: float, m: LegendMetrics) -> list[str]:
    """`text` from `x` in lines ending before `right`, broken between words."""
    lines: list[str] = []
    for word in text.split(" "):
        if lines and _legend_text_width(f"{lines[-1]} {word}", m.fs, m) <= right - x:
            lines[-1] = f"{lines[-1]} {word}"
        else:
            lines.append(word)
    return lines


def pgs_row_layout(row: LegendRow, label_x: float, m: LegendMetrics) -> PgsRowLayout:
    """Where a PGS row's caption and key go, its text starting at ``label_x``,
    or further right when the icon's labels reach past it (#159).

    On one line when it ends before `LegendRow.wrap_right` (or that is
    None): the caption, then each score's name after a line of its color
    (#56). Otherwise (#68) the caption on the first line, wrapped between
    words if it is itself too long, and the key on the lines below,
    starting under the caption's text and wrapped between scores, never
    inside a name. Text is measured in DejaVu Sans at the row's size
    (`symbols.text_width`), so both renderers break at the same places.
    """
    em = m.fs / 72 / m.unit_x  # one em of the row's text, in legend units
    label_x = max(label_x, pgs_icon_geometry(row, m, 0.0).right + LEG_PGS_ICON_GAP)

    def place(x: float, name: str, color: str, line: int) -> LegendKeyItem:
        dash_x1 = x + LEG_KEY_DASH_EM * em
        text_x = dash_x1 + LEG_KEY_GAP_EM * em
        return LegendKeyItem(x, dash_x1, text_x, name, color, line)

    def end(item: LegendKeyItem) -> float:
        return item.text_x + _legend_text_width(item.name, m.fs, m)

    # One line, as before #68.
    items: list[LegendKeyItem] = []
    x = label_x + _legend_text_width(row.label, m.fs, m)
    for name, color in row.key:
        items.append(place(x + LEG_KEY_SPACE_EM * em, name, color, 0))
        x = end(items[-1])
    right = row.wrap_right
    if right is None or x <= right:
        return PgsRowLayout(label_x, (row.label,), tuple(items), 1)

    captions = _wrap_words(row.label, label_x, right, m)
    line = len(captions)
    items = []
    x = label_x
    for name, color in row.key:
        item = place(x, name, color, line)
        if x > label_x and end(item) > right:
            line += 1
            item = place(label_x, name, color, line)
        items.append(item)
        x = end(item) + LEG_KEY_SPACE_EM * em
    return PgsRowLayout(label_x, tuple(captions), tuple(items), line + 1 if items else line)


def legend_key_items(row: LegendRow, label_x: float, m: LegendMetrics) -> list[LegendKeyItem]:
    """Where the PGS row's key goes, after its label written from ``label_x`` (#56).

    See `pgs_row_layout`, which also gives the caption's lines.
    """
    return list(pgs_row_layout(row, label_x, m).items)


@dataclass(frozen=True)
class PgsIconGeometry:
    """The PGS row's icon (#159): the box as drawn beside a symbol, upright.

    Positions are in legend units (y up). The curve's peak is on the left,
    high values at the top, as in the figure (#126); "high" is written above
    the box, "low" below it, and each ±2 SD line's name on its right.

    Attributes:
        box: (left, bottom, width, height) of the box.
        curve: The curve's points, bottom to top.
        base_x: x of the curve's baseline (the shaded area runs to it).
        tails: y of each ±2 SD line, across the box.
        labels: (x, y, text, anchor) per label, ``y`` its centre; anchor is
            "middle" (centred on x) or "start".
        size: Size of the labels, in points.
        right: Where the icon and its labels end.
        height: Height of the icon with "high" and "low".
    """

    box: tuple[float, float, float, float]
    curve: tuple[tuple[float, float], ...]
    base_x: float
    tails: tuple[float, ...]
    labels: tuple[tuple[float, float, str, str], ...]
    size: float
    right: float
    height: float


def pgs_icon_geometry(row: LegendRow, m: LegendMetrics, y_center: float) -> PgsIconGeometry:
    """Where the parts of a "pgs" row's icon go (#159).

    ``y_center`` is the row's centre from `legend_row_centers` (its first
    line's, when the row is wrapped): the icon is centred on the row's text.
    Sizes are in ems of the legend's font, so the icon scales with it as the
    swatches do; its labels are captions (`LegendMetrics.caption_fs`).
    """
    w = LEG_PGS_BOX_W_EM / LEG_UNIT_X_EM
    h = LEG_PGS_BOX_H_EM / LEG_UNIT_Y_EM
    size = m.caption_fs()
    line = LEG_LINE_EM * size / (LEG_UNIT_Y_EM * m.font)
    gap_y = LEG_PGS_TEXT_GAP_EM * size / (LEG_UNIT_Y_EM * m.font)
    gap_x = LEG_PGS_TEXT_GAP_EM * size / 72 / m.unit_x
    left = LEG_X0
    bottom = y_center - (row.lines - 1) * m.line_h / 2 - h / 2
    curve = tuple((pgs_curve_x(left, w, d), bottom + t * h) for t, d in pgs_curve(30))

    def width(text: str) -> float:
        return _legend_text_width(text, size, m)

    labels: list[tuple[float, float, str, str]] = []
    right = left + w
    # "high" and "low" centred on the box, never past the frame's left edge.
    for text, y in (
        (row.high, bottom + h + gap_y + line / 2),
        (row.low, bottom - gap_y - line / 2),
    ):
        if text:
            x = max(left + w / 2, LEG_FRAME_INSET + LEG_BOTTOM_PAD + width(text) / 2)
            labels.append((x, y, text, "middle"))
            right = max(right, x + width(text) / 2)
    tails = []
    for position, text in row.tails:
        y = bottom + position * h
        tails.append(y)
        labels.append((left + w + gap_x, y, text, "start"))
        right = max(right, left + w + gap_x + width(text))
    return PgsIconGeometry(
        box=(left, bottom, w, h),
        curve=curve,
        base_x=pgs_curve_base(left, w),
        tails=tuple(tails),
        labels=tuple(labels),
        size=size,
        right=right,
        height=h + 2 * (gap_y + line),
    )


def caption_overhang(row: LegendRow, m: LegendMetrics) -> float:
    """How far the captions under a row's swatch (a gradient's ends, category
    names, a binned row's limits) reach into the row below it."""
    n = legend_caption_lines(row)
    if not n:
        return 0.0
    caption_h = LEG_LINE_EM * m.label_fs / (LEG_UNIT_Y_EM * m.font)
    first = m.row_h / 2 if row.lines > 1 else row.slots * m.row_h / 2
    reach = m.sw_h / 2 + LEG_CAPTION_GAP + caption_h + (n - 1) * m.line_h
    return max(0.0, reach - first - (row.lines - 1) * m.line_h - row.extra_h / 2)


def pgs_row_extra_h(row: LegendRow, m: LegendMetrics, above: LegendRow | None = None) -> float:
    """Room a "pgs" row takes beyond its text (#159): its icon with "high" and
    "low", the margin a swatch has above and below it in its row, and, twice
    (the icon is centred), what the captions of the row ``above`` reach into
    it, so "high" keeps from them the distance a swatch keeps from its row's top."""
    text_h = m.row_h + (row.lines - 1) * m.line_h
    margin = m.row_h - m.sw_h
    hang = caption_overhang(above, m) if above is not None else 0.0
    return max(0.0, pgs_icon_geometry(row, m, 0.0).height + margin + 2 * hang - text_h)


LEG_BIN_GAP_EM = 0.6  # between two steps' limits, in ems of the caption
LEG_LABEL_GAP = 0.3  # between a bar and its row's label, in legend units
LEG_BIN_INSIDE_PAD_EM = 0.3  # above and below a step's name written inside it (#90), in ems
LEG_BIN_INSIDE_SIDE_EM = 0.5  # on each side of it (#94), in ems


@dataclass(frozen=True)
class BinnedGeometry:
    """Where a "binned" row's steps and label go, in legend units (#57).

    The steps share a bar as wide as a gradient's, or wider when their
    limits need it. ``label_x`` is the least x its name needs, past the
    bar; the names of all rows start at `LegendMetrics.name_x`, the
    largest of these (#88).
    """

    step_w: float
    label_x: float
    name_dy: float = 0.0  # from a step's limits down to its name: a row line (#81, #94)
    bar_h: float = 0.0  # the bar's height: a swatch's, or the names' inside it (#90)
    name_size: float = 0.0  # the names' font size, in points


def binned_row_geometry(row: LegendRow, m: LegendMetrics) -> BinnedGeometry:
    """Lay out a "binned" row: each step as wide as its widest caption,
    limits or name (#81), the names one caption line below the limits; or,
    with the names inside the steps (#90), each step as wide as its name
    plus a margin on each side, and the bar as tall as a line of them plus
    a margin, when that is taller than a swatch."""
    size = m.caption_fs()
    em = size / 72 / m.unit_x  # one em of the captions, horizontally
    gap = LEG_BIN_GAP_EM * em
    line = LEG_LINE_EM * size / (LEG_UNIT_Y_EM * m.font)
    inside = row.names_inside and bool(row.bin_names)
    below = () if inside else row.bin_names
    widest = max((_legend_text_width(b, size, m) for b in [*row.bins, *below]), default=0.0)
    step_w = max(LEG_GRAD_W / max(len(row.colors), 1), widest + gap)
    if m.bar_w:
        # The legend's common bar width (#95), never narrower than its steps.
        step_w = max(step_w, m.bar_w / max(len(row.colors), 1))
    bar_h = m.sw_h
    if inside:
        names = max((_legend_text_width(n, size, m) for n in row.bin_names), default=0.0)
        step_w = max(step_w, names + 2 * LEG_BIN_INSIDE_SIDE_EM * em)
        pad_v = LEG_BIN_INSIDE_PAD_EM * size / (LEG_UNIT_Y_EM * m.font)
        bar_h = max(bar_h, line + 2 * pad_v)
    bar_end = LEG_X0 + step_w * len(row.colors)
    # Names below go one row line (`m.line_h`) under the limits: the line the
    # row takes for them, so the next row keeps the usual gap (#94).
    return BinnedGeometry(step_w, max(m.text_x, bar_end + LEG_LABEL_GAP), m.line_h, bar_h, size)


def legend_row_right(row: LegendRow, m: LegendMetrics) -> float:
    """Where a row's content ends, in legend units (DejaVu Sans widths).

    Used to wrap the PGS row at the widest other row (#68).
    """
    if row.kind == "quadrants":
        g = quadrant_key_geometry(row, m, 0.0)
        sx, _sy, sw, _sh = g.square
        right = sx + sw if g.circle is None else g.circle[0] + g.circle[2]
        for x, _y, text, anchor in g.labels:
            if anchor == "start":
                right = max(right, x + _legend_text_width(text, m.fs, m))
        return right
    if row.kind == "key_note":
        return LEG_X0 + _legend_text_width(row.label, m.caption_fs(), m)
    if row.kind in NAMED_ROW_KINDS:
        right = m.name_x + _legend_text_width(row.label, m.fs, m)
        if row.kind == "categorical":
            right = max(right, LEG_X0 + len(row.categories) * (m.sw_w + m.sw_gap * 0.5))
        return right
    if row.kind == "pgs":
        layout = pgs_row_layout(row, m.icon_text_x, m)
        ends = [layout.text_x + _legend_text_width(c, m.fs, m) for c in layout.captions]
        ends += [i.text_x + _legend_text_width(i.name, m.fs, m) for i in layout.items]
        return max(ends)
    return m.icon_text_x + _legend_text_width(row.label, m.fs, m)


def legend_max_width_inches(style: PedigreeStyle) -> float | None:
    """``legend_max_width`` in inches: a number is in centimetres, or text
    with its unit ("85mm", "8.5cm", "3.3in"); None when unset."""
    value = style.legend_max_width
    if value is None:
        return None
    if isinstance(value, str):
        match = re.fullmatch(r"\s*([0-9]*\.?[0-9]+)\s*(cm|mm|in)\s*", value)
        assert match is not None  # checked by `check_legend`
        return float(match.group(1)) * LEG_WIDTH_UNITS[match.group(2)]
    return float(value) * LEG_WIDTH_UNITS["cm"]


def legend_layout(
    rows: Sequence[LegendRow], style: PedigreeStyle
) -> tuple[list[LegendRow], LegendMetrics]:
    """The legend's rows, PGS row wrapped (#68), and the panel's metrics.

    The PGS row wraps where the widest other row ends, or, with
    ``legend_max_width``, where the panel reaches that width (never before
    the widest other row: they do not wrap). Its height then counts each
    extra line. Both renderers call this, so they wrap at the same places.
    """
    # A binned row with named steps takes a second line for the names (#81),
    # unless they go inside the steps (#90).
    rows = [
        replace(r, lines=2) if r.kind == "binned" and r.bin_names and not r.names_inside else r
        for r in rows
    ]
    binary_headers = any(r.kind in ("binary", "binary_shared") for r in rows)

    def metrics(rs: Sequence[LegendRow]) -> LegendMetrics:
        last = rs[-1] if rs else None
        return legend_metrics(
            legend_slot_count(rs),
            style,
            binary_headers=binary_headers,
            extra_lines=sum(r.lines - 1 for r in rs),
            bottom_captions=legend_caption_lines(last) if last else 0,
            bottom_lines=last.lines if last else 1,
            extra_h=sum(r.extra_h for r in rs),
        )

    # A bar taller than a swatch, for names inside its steps (#90), reaches
    # half its extra height above the swatch and half below: the row takes
    # all of it (#94), so neither the row above nor the captions under it
    # come closer than to a swatch.
    m0 = metrics(rows)
    rows = [
        replace(r, extra_h=max(0.0, binned_row_geometry(r, m0).bar_h - m0.sw_h))
        if r.kind == "binned" and r.names_inside
        else r
        for r in rows
    ]
    m = with_name_column(with_bar_width(metrics(rows), rows, style), rows)
    others = [legend_row_right(r, m) for r in rows if r.kind != "pgs"]
    right = max(others) if others else None
    max_in = legend_max_width_inches(style)
    if max_in is not None:
        # The border runs from 0 to `right + LEG_X0 + 0.1` (its rounded pad).
        # Never narrower than the other rows, which do not wrap: the panel
        # would be as wide, only taller.
        right = max(max_in / m.unit_x - LEG_X0 - 0.1, right or 0.0)
    out = []
    for row in rows:
        if row.kind == "pgs":
            if right is not None:
                row = replace(row, wrap_right=right)
                row = replace(row, lines=pgs_row_layout(row, m.icon_text_x, m).lines)
            # As tall as its icon (#159), half above its text and half below.
            above = out[-1] if out else None
            row = replace(row, extra_h=pgs_row_extra_h(row, m, above))
        out.append(row)
    if out == rows:
        return rows, m
    return out, with_name_column(with_bar_width(metrics(out), out, style), out)


BAR_WIDTHS = ("equal", "own")
BAR_ROW_KINDS = ("gradient", "binned")


def legend_gradient_width(m: LegendMetrics) -> float:
    """A gradient bar's width: the legend's common bar width (#95), or `LEG_GRAD_W`."""
    return m.bar_w or LEG_GRAD_W


def with_bar_width(
    m: LegendMetrics, rows: Sequence[LegendRow], style: PedigreeStyle
) -> LegendMetrics:
    """``m`` with one width for every bar (#95): the widest gradient or
    binned bar's, when the legend has more than one and ``legend_bar_width``
    is "equal"; the widest bar, and so the legend, keep their width."""
    bars = [r for r in rows if r.kind in BAR_ROW_KINDS]
    if style.legend_bar_width != "equal" or len(bars) < 2:
        return replace(m, bar_w=0.0)
    own = replace(m, bar_w=0.0)
    widths = [
        LEG_GRAD_W if r.kind == "gradient" else binned_row_geometry(r, own).step_w * len(r.colors)
        for r in bars
    ]
    return replace(m, bar_w=max(widths))


def with_name_column(m: LegendMetrics, rows: Sequence[LegendRow]) -> LegendMetrics:
    """``m`` with the legend's one name column (#88): `text_x`, or past the
    widest binned bar, so the names of every swatch and bar row start
    together however long a binned row's steps are."""
    ends = [binned_row_geometry(r, m).label_x for r in rows if r.kind == "binned"]
    return replace(m, name_x=max([m.text_x, *ends]))


def legend_caption_lines(row: LegendRow) -> int:
    """Lines of captions a row writes under its swatch (a gradient's ends,
    category names, a binned row's limits): what the bottom margin counts
    when it is the legend's last row (#77)."""
    if row.kind == "binned" and row.bin_names and not row.names_inside:
        return 2  # limits, then names (#81); names inside the steps take none (#90)
    return 1 if row.kind in CAPTIONED_ROW_KINDS else 0


def ptdt_legend_description() -> str:
    """Caption of the pTDT row (#52), followed by its key; the icon names
    its ±2 SD lines (#159)."""
    return f"pTDT deviation, (child{NBSP}{MINUS}{NBSP}mid-parent){NBSP}/{NBSP}SD:"


def pgs_legend_description(pgs_names: Sequence[str], scale: str = "minmax") -> str:
    """Caption of the PGS row, from its column names and the scale.

    The caption ends with a colon, and the row's key (`LegendRow.key`)
    follows it: each score's name after a line of its color (#56). The icon
    shows which way is high and names the ±2 SD lines (#159).
    """
    caption = "PGS distribution"
    if scale == "zscore":
        caption += ", z-score"
    elif scale == "percentile":
        caption += ", population percentile"
    return f"{caption}:"


def pgs_icon_tails(style: PedigreeStyle, scale: str) -> tuple[tuple[float, str], ...]:
    """The ±2 SD lines of the legend's PGS icon on ``scale``, with their
    names (#159): none on "minmax", which marks no tails."""
    names = (legend_text(style, "pgs_minus_2sd"), legend_text(style, "pgs_plus_2sd"))
    return tuple(zip(pgs_tail_positions(scale), names, strict=False))


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
        inverted = display_modes.get(name) == PhenotypeDisplayMode.GRADIENT_INVERTED
        if effective == "gradient" and style.gradient_bins.get(name):
            # Clinical thresholds (#57): the steps, each with its limits.
            cutoffs = style.gradient_bins[name]
            rows.append(
                LegendRow(
                    "binned",
                    display_name(name, style),
                    colors=tuple(
                        bin_color(to_rgb(color), i, cutoffs, inverted)
                        for i in range(len(cutoffs) + 1)
                    ),
                    bins=bin_labels(cutoffs),
                    bin_names=bin_names(name, style),
                    names_inside=bin_label_position(name, style) == "inside",
                )
            )
        elif name in style.gradient_diverging:
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
            rgb = to_rgb(color)
            # With `gradient_range` the bar's ends are the range's (#57).
            ends = style.gradient_range.get(name)
            rows.append(
                LegendRow(
                    "gradient",
                    display_name(name, style),
                    colors=tuple(blend_with_white(rgb, 1.0 - t if inverted else t) for t in steps),
                    low=f"{ends[0]:g}" if ends else legend_text(style, "low"),
                    high=f"{ends[1]:g}" if ends else legend_text(style, "high"),
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
    # is which column, so two or more binary columns share one row. Minimal
    # (#67): they get no row at all, a caption under the key says how to
    # read them, even for one column.
    minimal = style.legend == "minimal"
    compact = style.legend == "compact" or (style.legend == "auto" and bool(render.quadrants))
    if minimal and binary_in_key:
        merged = {i for i, _c in binary_in_key}
        rows = [r for i, r in enumerate(rows) if i not in merged]
    elif (compact or minimal) and len(binary_in_key) >= 2:
        first = binary_in_key[0][0]
        shared = LegendRow(
            "binary_shared",
            legend_text(style, "shared_binary"),
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
    if render.asymptomatic:
        # A white swatch with the vertical line (#104).
        rows.append(LegendRow("asymptomatic", legend_text(style, "asymptomatic")))

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

    curve_names = render.pgs_names(names)
    pgs_names = [n for n in curve_names if not is_ptdt_column(n)]
    ptdt_names = [n for n in curve_names if is_ptdt_column(n)]
    if pgs_names:
        colors = {n: style.pgs_colors.get(n, DEFAULT_PGS_COLOR) for n in pgs_names}
        rows.append(
            LegendRow(
                "pgs",
                legend_text(style, "pgs") or pgs_legend_description(pgs_names, style.pgs_scale),
                pgs=colors,
                key=tuple((pgs_trait_name(n, style), colors[n]) for n in pgs_names),
                high=legend_text(style, "pgs_high"),
                low=legend_text(style, "pgs_low"),
                tails=pgs_icon_tails(style, style.pgs_scale),
            )
        )
    if ptdt_names:
        # pTDT deviations (#52): their own row, on the z scale whatever pgs_scale.
        colors = {n: style.pgs_colors.get(n, DEFAULT_PGS_COLOR) for n in ptdt_names}
        rows.append(
            LegendRow(
                "pgs",
                legend_text(style, "ptdt") or ptdt_legend_description(),
                pgs=colors,
                key=tuple((pgs_trait_name(n, style), colors[n]) for n in ptdt_names),
                high=legend_text(style, "pgs_high"),
                low=legend_text(style, "pgs_low"),
                tails=pgs_icon_tails(style, "zscore"),
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
        rows.insert(0, LegendRow("quadrants", "", quadrants=tuple(cells), single_symbol=minimal))
        if minimal and binary_in_key:
            rows.insert(1, LegendRow("key_note", legend_text(style, "key_note")))
    return rows
