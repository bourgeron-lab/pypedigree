"""Configuration and style parameters for pedigree rendering.

All visual parameters are centralized here. No magic numbers
should exist in layout.py or renderer.py.
"""

from __future__ import annotations

import unicodedata
import warnings
from collections.abc import Sequence
from dataclasses import asdict, dataclass, field, fields, replace
from difflib import get_close_matches
from typing import Any

from pypedigree.colors import parse_rgb
from pypedigree.errors import PedigreeStyleWarning

DEFAULT_PHENOTYPE_PALETTE = [
    "#e74c3c",  # red
    "#3498db",  # blue
    "#e67e22",  # orange
    "#2ecc71",  # green
    "#9b59b6",  # purple
    "#1abc9c",  # teal
    "#f1c40f",  # yellow
    "#e91e63",  # pink
]

# More colors, used once the palette runs short: when known phenotypes and
# the colors the user chose leave too few distinct ones, or past eight
# columns. Chosen far from the palette and from the known phenotypes' colors.
EXTRA_PHENOTYPE_COLORS = [
    "#8c564b",  # brown
    "#e377c2",  # light magenta
    "#bcbd22",  # olive
]

# Marker color of a PGS column with no `pgs_colors` entry, when it is the
# figure's only PGS column (several get distinct colors, see `resolve_style`).
DEFAULT_PGS_COLOR = "#e74c3c"

# Palette pairs too alike to tell apart in a thin border or a small band,
# although their hues differ (perceptual distance alone does not catch them:
# red and orange are as far apart as green and teal).
_ALIKE_PAIRS = {
    frozenset({"#2ecc71", "#1abc9c"}),  # green, teal
    frozenset({"#f1c40f", "#bcbd22"}),  # yellow, olive
    frozenset({"#9b59b6", "#e377c2"}),  # purple, light magenta
}
# Below this CIE76 distance two colors pass for one (the palette's red is
# 6.6 from intellectual disability's red, the blue 17.9 from autism's).
_MIN_COLOR_DISTANCE = 20.0

# Phenotypes with a conventional color, recognised by name whatever the
# column prefix or case: autism blue, intellectual disability red, ADHD
# yellow, epilepsy purple. A color in `phenotype_colors` still wins.
KNOWN_PHENOTYPE_COLORS = {
    "autism": "#386eb2",
    "intellectual disability": "#d34436",
    "adhd": "#e59635",
    "epilepsy": "#5e287e",
}
_KNOWN_PHENOTYPE_ALIASES = {
    "autism": ("autism", "autisme", "asd", "tsa", "autism spectrum disorder"),
    "intellectual disability": (
        "intellectual disability",
        "id",
        "di",
        "deficience intellectuelle",
    ),
    "adhd": ("adhd", "tdah"),
    "epilepsy": ("epilepsy", "epilepsie"),
}
# Palette entries close to a known phenotype's color: skipped for the other
# columns while that phenotype is in the figure, so no two bands look alike.
_PALETTE_LOOKALIKES = {
    "autism": ("#3498db",),
    "intellectual disability": ("#e74c3c", "#e91e63"),
    "adhd": ("#e67e22", "#f1c40f"),
    "epilepsy": ("#9b59b6",),
}


def known_phenotype(name: str) -> str | None:
    """The conventional phenotype a column name stands for, or None.

    Case, accents, and ``_``/``-`` separators are ignored, so ``ASD``,
    ``autisme`` and ``Autism_spectrum_disorder`` all give ``"autism"``.

    Args:
        name: Phenotype column name, prefix already stripped.

    Returns:
        One of the keys of `KNOWN_PHENOTYPE_COLORS`, or None.
    """
    folded = unicodedata.normalize("NFKD", name).encode("ascii", "ignore").decode()
    folded = " ".join(folded.lower().replace("_", " ").replace("-", " ").split())
    for key, aliases in _KNOWN_PHENOTYPE_ALIASES.items():
        if folded in aliases:
            return key
    return None


@dataclass
class PedigreeStyle:
    """Visual style configuration for pedigree diagrams.

    All spacing values are in diagram coordinate units.
    """

    # Shape parameters
    shape_size: float = 0.8
    shape_line_width: float = 1.5
    fill_color: str = "black"
    outline_color: str = "black"
    phantom_line_style: str = "--"

    # Line parameters
    line_color: str = "black"
    line_width: float = 1.0
    double_line_gap: float = 0.06
    broken_gap_size: float = 0.15

    # Text parameters
    font_size: float = 8.0
    font_family: str = "sans-serif"
    label_offset: float = 0.15
    # Smallest size of any text, in points (None: no minimum). The secondary
    # text (notes, annotations, band and quadrant text, "PGS" captions,
    # legend) is drawn at a fraction of `font_size`; journals often ask for
    # 6 pt at print size, which those fractions miss unless the labels are
    # made much larger (#51). See `text_size`.
    min_font_size: float | None = None

    # Layout spacing
    generation_gap: float = 2.0
    sibling_gap: float = 1.5
    couple_gap: float = 1.5
    # Not used by the layout: `sibling_gap` and `couple_gap` set the spacing.
    # Kept so existing styles still load; to be removed in 1.0.
    min_node_distance: float = 1.2

    # Carrier notation: "nsgc2022" divides the symbol into patterned sections
    # (Bennett et al. 2022); "legacy" draws the `carrier` column as a
    # half-filled symbol and obligate carriers as a central dot.
    carrier_notation: str = "nsgc2022"
    carrier_dot_radius: float = 0.08

    # Multi-phenotype colors (pheno_name -> hex color)
    phenotype_colors: dict[str, str] = field(default_factory=dict)

    # Show YES/NO text labels inside binary phenotype bands (default: color only)
    show_phenotype_labels: bool = False

    # How phenotypes divide a symbol: "quadrants" (the symbol cut in four,
    # one quadrant per column, up to four), "bands" (horizontal stripes, one
    # per column), or "auto" (the default: quadrants up to four phenotype and
    # carrier sources, bands beyond, without a warning). See
    # `phenotype_quadrants`.
    phenotype_layout: str = "auto"

    # Quadrant of a column in the "quadrants" layout: column name ->
    # "top-left", "top-right", "bottom-left" or "bottom-right". Columns not
    # listed take their conventional quadrant (autism top left, ADHD top
    # right, intellectual disability bottom left, epilepsy bottom right) or
    # the next free one, in column order.
    phenotype_quadrants: dict[str, str] = field(default_factory=dict)

    # How a PGS column's value is placed on its curve: "minmax" (between the
    # pedigree's lowest and highest score), "zscore" (the value is a z-score,
    # e.g. ancestry-matched, placed on the standard normal) or "percentile"
    # (a population percentile, 0-100). With a fixed scale the curve's
    # +/-2 SD tails are marked and the tooltip gives z and percentile.
    pgs_scale: str = "minmax"

    # Write "NA" in bands with no data (-9 / NA / empty). Without it, a band
    # with no data is the same plain grey as an "unknown" (0) one.
    show_na_labels: bool = True

    # PGS curve colors (PGS column name -> hex color, e.g. {"PGS_autism": "#3498db"})
    pgs_colors: dict[str, str] = field(default_factory=dict)

    # Custom gradient normalization ranges (pheno_name -> [min, max])
    gradient_range: dict[str, list[float]] = field(default_factory=dict)

    # Diverging (tricolor) gradient config (pheno_name -> {min, mid, max, color_low, color_high})
    gradient_diverging: dict[str, dict[str, Any]] = field(default_factory=dict)

    # Zero-pad individual IDs to this width (0 = no padding)
    zero_pad_ids: int = 0

    # Circular layout parameters
    radial_gap: float = 2.5
    base_radius: float = 3.0

    # Output parameters
    dpi: int = 150
    fig_width: float = 10.0
    fig_height: float = 8.0
    margin: float = 0.5

    @classmethod
    def from_dict(cls, d: dict[str, Any], *, warn_unknown: bool = True) -> PedigreeStyle:
        """Create a style from a dictionary, ignoring unknown keys.

        Unknown keys are skipped rather than raising, so a style file stays
        forward-compatible with newer options. By default each one is
        reported as a `PedigreeStyleWarning` with a suggested spelling,
        since silently ignoring a typo such as `phenotype_color` is
        indistinguishable from the option having no effect.

        Args:
            d: Dictionary of style overrides.
            warn_unknown: Whether to warn about keys that are not style
                fields. Set False to accept a dictionary that legitimately
                carries extra data.

        Returns:
            PedigreeStyle with overrides applied over defaults.
        """
        valid_keys = {f.name for f in fields(cls)}
        unknown = [k for k in d if k not in valid_keys]
        if unknown and warn_unknown:
            for key in unknown:
                close = get_close_matches(str(key), valid_keys, n=1)
                hint = f"; did you mean {close[0]!r}?" if close else ""
                warnings.warn(
                    f"Unknown style option {key!r} ignored{hint}",
                    PedigreeStyleWarning,
                    stacklevel=2,
                )
        filtered = {k: v for k, v in d.items() if k in valid_keys}
        return cls(**filtered)

    def to_dict(self) -> dict[str, Any]:
        """Convert to a plain dictionary."""
        return asdict(self)

    def text_size(self, scale: float = 1.0) -> float:
        """Size in points of a text drawn at ``scale`` times `font_size`.

        Never below `min_font_size`. Every text of both renderers, and the
        layout's measure of the labels, take their size from here.

        Args:
            scale: Fraction of `font_size` (e.g. `SMALL_TEXT_SCALE`).

        Returns:
            The size in points.
        """
        size = self.font_size * scale
        if self.min_font_size is not None and size < self.min_font_size:
            return float(self.min_font_size)
        return size


def is_pgs_column(name: str) -> bool:
    """Check whether a phenotype column name denotes a polygenic score.

    Args:
        name: Phenotype column name.

    Returns:
        True if the name starts with "PGS" (case-insensitive).
    """
    return name.upper().startswith("PGS")


def resolve_style(
    style: PedigreeStyle | None,
    phenotype_names: Sequence[str] | None = None,
) -> PedigreeStyle:
    """Fill in the pedigree-dependent parts of a style.

    Renderers only draw phenotype bands when ``phenotype_colors`` is
    populated, which depends on the pedigree being rendered rather than on
    the style alone. Every entry point (CLI, PedigreeBuilder,
    render_pedigree) routes through this function to stay consistent.

    The input style is never mutated, and the function is idempotent:
    resolving an already-resolved style returns an equivalent style, so
    repeated renders with the same object are stable.

    Args:
        style: User-supplied style, or None for defaults.
        phenotype_names: Ordered phenotype column names of the pedigree.

    Returns:
        A new PedigreeStyle with auto-assigned phenotype colors.
    """
    resolved = PedigreeStyle() if style is None else replace(style)
    # replace() copies the dataclass but shares its mutable fields.
    resolved.phenotype_colors = dict(resolved.phenotype_colors)
    resolved.pgs_colors = dict(resolved.pgs_colors)
    resolved.phenotype_quadrants = dict(resolved.phenotype_quadrants)
    resolved.gradient_range = {k: list(v) for k, v in resolved.gradient_range.items()}
    resolved.gradient_diverging = {k: dict(v) for k, v in resolved.gradient_diverging.items()}

    all_names = list(phenotype_names or [])
    names = [n for n in all_names if not is_pgs_column(n)]

    # Conventional colors for known phenotypes (autism, ID, ADHD, epilepsy),
    # unless the user chose one.
    known = {n: k for n in names if (k := known_phenotype(n)) is not None}
    for name, key in known.items():
        resolved.phenotype_colors.setdefault(name, KNOWN_PHENOTYPE_COLORS[key])

    # Auto-assign the other columns' colors from the palette, skipping the
    # entries that look like a known phenotype's color. Each column first
    # gets the palette entry at its position among all columns, as before
    # known phenotypes existed, so a figure without one keeps its colors.
    taken = {c for key in known.values() for c in _PALETTE_LOOKALIKES[key]}
    palette = [c for c in DEFAULT_PHENOTYPE_PALETTE if c not in taken]
    spare = [c for c in DEFAULT_PHENOTYPE_PALETTE + EXTRA_PHENOTYPE_COLORS if c not in taken]
    # When some colors are already fixed (known phenotypes, or chosen by the
    # user), a column whose entry repeats one of them, or looks like one, takes
    # the first spare color that does not: two columns must not look alike.
    # Without them the palette's order is kept as it was, so such figures do
    # not change; past eight columns the spare colors come before the palette
    # starts over.
    in_use = [resolved.phenotype_colors[n] for n in names if n in resolved.phenotype_colors]
    deconflict = bool(in_use)
    for i, name in enumerate(all_names):
        if name not in names or name in resolved.phenotype_colors:
            continue
        color = palette[i % len(palette)] if palette else spare[i % len(spare)]
        if deconflict or i >= len(palette):
            clash = any(_alike(color, used) for used in in_use)
            if clash:
                # Failing a color unlike all, one not already used.
                distinct = [c for c in spare if not any(_alike(c, u) for u in in_use)]
                unused = [c for c in spare if c.lower() not in {u.lower() for u in in_use}]
                color = (distinct or unused or [color])[0]
        resolved.phenotype_colors[name] = color
        in_use.append(color)

    _assign_pgs_colors(resolved, all_names, known)
    return resolved


def _assign_pgs_colors(
    style: PedigreeStyle, all_names: Sequence[str], known: dict[str, str]
) -> None:
    """Give each PGS column of a figure with several a distinct marker color.

    A score named after a known phenotype (``PGS_autism``) takes that
    phenotype's color, the one its column has in the figure if there is one;
    the others take the first color unlike every PGS color and phenotype
    color in use, starting from `DEFAULT_PGS_COLOR`. A figure with one PGS
    column is left as it was (`DEFAULT_PGS_COLOR`), and `pgs_colors` entries
    always win.

    Args:
        style: Resolved style whose `pgs_colors` is filled in place.
        all_names: Ordered phenotype column names of the pedigree.
        known: Non-PGS columns standing for a known phenotype, with its key.
    """
    pgs = [n for n in all_names if is_pgs_column(n)]
    if len(pgs) < 2:
        return
    bands = {
        n: style.phenotype_colors[n]
        for n in all_names
        if not is_pgs_column(n) and n in style.phenotype_colors
    }
    column_color = {key: bands[n] for n, key in reversed(known.items()) if n in bands}
    used = [style.pgs_colors[n] for n in pgs if n in style.pgs_colors]
    pending = [n for n in pgs if n not in style.pgs_colors]
    for name in list(pending):
        key = known_phenotype(name[len("PGS") :])
        if key is None:
            continue
        color = column_color.get(key, KNOWN_PHENOTYPE_COLORS[key])
        # Not the color of another column, nor of a score already drawn.
        others = [c for n, c in bands.items() if known.get(n) != key]
        if not any(_alike(color, c) for c in used + others):
            style.pgs_colors[name] = color
            used.append(color)
            pending.remove(name)
    # The other scores avoid the phenotype colors too, so no marker seems to
    # belong to a band it has nothing to do with.
    candidates = [DEFAULT_PGS_COLOR, *DEFAULT_PHENOTYPE_PALETTE, *EXTRA_PHENOTYPE_COLORS]
    for name in pending:
        free = [c for c in candidates if not any(_alike(c, u) for u in used)]
        apart = [c for c in free if not any(_alike(c, b) for b in bands.values())]
        unused = [c for c in candidates if c.lower() not in {u.lower() for u in used}]
        color = (apart or free or unused or [DEFAULT_PGS_COLOR])[0]
        style.pgs_colors[name] = color
        used.append(color)


def _alike(a: str, b: str) -> bool:
    """Whether two colors would be taken for one another in a figure."""
    if a.lower() == b.lower() or frozenset({a.lower(), b.lower()}) in _ALIKE_PAIRS:
        return True
    try:
        return _lab_distance(parse_rgb(a), parse_rgb(b)) < _MIN_COLOR_DISTANCE
    except ValueError:  # an invalid user color is reported by `check_colors`
        return False


def _lab_distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    """CIE76 distance between two sRGB colors (D65 white)."""

    def lab(rgb: tuple[float, float, float]) -> tuple[float, float, float]:
        lin = [c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4 for c in rgb]
        r, g, b_ = lin
        x = (0.4124 * r + 0.3576 * g + 0.1805 * b_) / 0.95047
        y = 0.2126 * r + 0.7152 * g + 0.0722 * b_
        z = (0.0193 * r + 0.1192 * g + 0.9505 * b_) / 1.08883

        def f(t: float) -> float:
            return t ** (1 / 3) if t > 0.008856 else 7.787 * t + 16 / 116

        fx, fy, fz = f(x), f(y), f(z)
        return (116 * fy - 16, 500 * (fx - fy), 200 * (fy - fz))

    la, lb = lab(a), lab(b)
    return float(sum((p - q) ** 2 for p, q in zip(la, lb, strict=True)) ** 0.5)


def apply_pgs_spacing(
    style: PedigreeStyle,
    phenotype_names: Sequence[str] | None = None,
) -> PedigreeStyle:
    """Widen node spacing to make room for PGS distribution curves.

    PGS curves are drawn beside each node, so couples and siblings need
    extra horizontal room. Unlike `resolve_style` this is not idempotent —
    each call widens further — so exactly one caller may apply it, which
    is `render_pedigree` just before laying families out.

    Args:
        style: Style to widen.
        phenotype_names: Ordered phenotype column names of the pedigree.

    Returns:
        A new PedigreeStyle, widened if any PGS column is present, else
        `style` unchanged.
    """
    if not any(is_pgs_column(n) for n in (phenotype_names or [])):
        return style

    pgs_extra = style.shape_size * 0.55
    return replace(
        style,
        couple_gap=style.couple_gap + pgs_extra,
        sibling_gap=style.sibling_gap + pgs_extra,
    )
