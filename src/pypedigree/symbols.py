"""How standard pedigree notation is drawn, shared by both renderers.

Decisions and geometry for the notation of ``models.Notation``: pregnancies
and pregnancy losses, adoption, childlessness, assisted reproduction, gender
identity and the consultand. Both the PNG renderer (``renderer.py``) and the
interactive SVG (``svg_renderer.py``) draw from these functions, so the two
outputs cannot drift apart. No matplotlib here.

Geometry is in layout units with the y axis pointing up (as in
``FamilyLayout``); functions whose result is not symmetric top to bottom
take an ``up`` sign, so the SVG renderer, whose y axis points down, passes
``up=-1``.

Reference: Bennett RL, French KS, Resta RG, Austin J. Standardized pedigree
nomenclature update centered on sex and gender inclusivity. J Genet Couns
2022;31(6):1238-48. Figures 1-5. The 2022 revision kept the 2008 symbols for
pregnancies and assisted reproduction (Bennett et al., J Genet Couns
2008;17:424-33).
"""

from __future__ import annotations

import math
import unicodedata
from collections.abc import Sequence

from pypedigree.models import (
    Childlessness,
    Gender,
    Individual,
    Notation,
    Pregnancy,
    Sex,
)

# The unknown-sex diamond's vertices sit this far from its centre, relative
# to half the shape size, so it reads as the same size as a square.
DIAMOND_SCALE = 1.2

# A pregnancy loss is a triangle with its apex where the line from the
# parents arrives (the top of an ordinary symbol). Its base lies this far
# below the centre and it is this wide, both relative to half the shape
# size: an equilateral triangle a little shorter than a square.
TRIANGLE_BASE = 0.7
TRIANGLE_HALF_WIDTH = 0.98

# Childlessness marker: a stem below the symbol ending in a bar, two bars
# for infertility. Lengths in layout units, except the bar's half-width,
# relative to half the shape size.
CHILDLESS_STEM = 0.25
CHILDLESS_BAR_GAP = 0.07
CHILDLESS_BAR_HALF_WIDTH = 0.75
# Gap between the marker and the label below it.
CHILDLESS_LABEL_GAP = 0.05

# Adoption brackets: gap to the symbol, length of the ticks pointing
# inward, and how far they reach above and below the symbol.
BRACKET_GAP = 0.08
BRACKET_TICK = 0.1
BRACKET_OVERHANG = 0.08

# Vertical step between lines of text below the label (notes, annotations),
# in layout units: the least step, kept whatever the text's size.
TEXT_LINE_STEP = 0.25

# The step also leaves room for the text as drawn (#66): the line above, at
# `TEXT_LINE_HEIGHT_EM` of its size (DejaVu Sans' ascent plus descent), then
# `TEXT_LINE_GAP_EM` of the next line's size, so text enlarged by
# `min_font_size` keeps its gaps.
TEXT_LINE_HEIGHT_EM = 1.17
TEXT_LINE_GAP_EM = 0.3

# Documented evaluation: an asterisk beside the symbol's lower right, at this
# size relative to the style's font size. Its gap is in ems of the asterisk's
# own font, from the outside of the outline's stroke to the glyph's ink
# (#110): in layout units, to the glyph's centre, the asterisk ate the gap and
# touched circles, the more so the smaller the pedigree was drawn. The ink of
# "*" in DejaVu Sans Bold, in ems: half its width and height, and its middle
# above the baseline (it is centred on its advance).
EVALUATION_MARK = "*"
EVALUATION_MARK_FONT_SCALE = 1.3
EVALUATION_MARK_GAP_EM = 0.3
EVALUATION_MARK_INK_HALF_W_EM = 0.242
EVALUATION_MARK_INK_HALF_H_EM = 0.232
EVALUATION_MARK_INK_MIDDLE_EM = 0.51

# The "?" standing for an unknown child, relative to the style's font size.
UNKNOWN_OFFSPRING_TEXT = "?"
UNKNOWN_OFFSPRING_FONT_SCALE = 2.0
# Gap between the end of the descent line and the top of the "?".
UNKNOWN_OFFSPRING_GAP = 0.06

# Proband and consultand marker, in ems of the label font: an arrow from this
# far out (in x and in y) to the symbol's lower left corner, and the
# proband's "P" ending this far to the left of the arrow's start, level with
# it. In text units, not layout units, so that the arrow keeps its length
# and the letter its place at any scale: drawn small, as in a journal figure,
# an arrow in layout units shrank to its head and the "P" covered it (0.4.6).
# At the default font (8 pt) the arrow is 11.2 pt, as long as a large
# pedigree's was.
MARKER_ARROW_EM = 1.4
MARKER_LETTER_GAP_EM = 0.2

# Layout units per typographic point in the interactive SVG. The PNG never
# draws a pedigree smaller than this (`renderer._keep_text_scale` grows the
# figure, up to its largest size), a small one larger, so text never takes
# more room there: layout decisions that depend on the size of text convert
# it with this.
UNITS_PER_POINT = 0.02

# Text widths, in ems, for the room the proband marker needs: the bold "P"
# (0.73 in DejaVu Sans Bold), one character of a label (digits are 0.64),
# and the space kept between the "P" and the neighbour's label.
_MARKER_LETTER_EM = 0.8
_LABEL_CHAR_EM = 0.65
_MARKER_PAD_EM = 0.3

# Advance widths of DejaVu Sans, matplotlib's default font, in ems, for the
# printable ASCII characters (read from the font file; `tests/test_symbols.py`
# checks them against it). A label's width is the sum of its characters'.
# Other characters: an accented letter counts as its base letter, a wide
# (CJK) character as 1 em, anything else as `_OTHER_CHAR_EM`.
_CHAR_EM: dict[str, float] = {
    ch: em
    for chars, em in (
        ("'", 0.275),
        ("ijl", 0.278),
        ("IJ", 0.295),
        (" ,.", 0.318),
        ("/:;\\|", 0.337),
        ("f", 0.352),
        ("-", 0.361),
        ("()[]", 0.39),
        ("t", 0.392),
        ("!", 0.401),
        ("r", 0.411),
        ('"', 0.46),
        ("*_`", 0.5),
        ("s", 0.521),
        ("z", 0.525),
        ("?", 0.531),
        ("c", 0.55),
        ("L", 0.557),
        ("F", 0.575),
        ("k", 0.579),
        ("vxy", 0.592),
        ("P", 0.603),
        ("TY", 0.611),
        ("o", 0.612),
        ("a", 0.613),
        ("e", 0.615),
        ("E", 0.632),
        ("hnu", 0.634),
        ("Sbdgpq", 0.635),
        ("0123456789${}", 0.636),
        ("K", 0.656),
        ("AV", 0.684),
        ("XZ", 0.685),
        ("B", 0.686),
        ("R", 0.695),
        ("C", 0.698),
        ("U", 0.732),
        ("N", 0.748),
        ("H", 0.752),
        ("D", 0.77),
        ("G", 0.775),
        ("&", 0.78),
        ("OQ", 0.787),
        ("w", 0.818),
        ("#+<=>^~", 0.838),
        ("M", 0.863),
        ("%", 0.95),
        ("m", 0.974),
        ("W", 0.989),
        ("@", 1.0),
    )
    for ch in chars
}
# The same for DejaVu Sans Bold, in which band and quadrant text is drawn.
_CHAR_EM_BOLD = {
    ch: em
    for chars, em in (
        ("'", 0.306),
        ("ijl", 0.343),
        (" ", 0.348),
        ("/\\|", 0.365),
        ("IJ", 0.372),
        (",.", 0.38),
        (":;", 0.4),
        ("-", 0.415),
        ("f", 0.435),
        ("!", 0.456),
        ("()[]", 0.457),
        ("t", 0.478),
        ("r", 0.493),
        ("_`", 0.5),
        ('"', 0.521),
        ("*", 0.523),
        ("?", 0.58),
        ("z", 0.582),
        ("c", 0.593),
        ("s", 0.595),
        ("L", 0.637),
        ("x", 0.645),
        ("vy", 0.652),
        ("k", 0.665),
        ("a", 0.675),
        ("e", 0.678),
        ("T", 0.682),
        ("EF", 0.683),
        ("o", 0.687),
        ("$0123456789", 0.696),
        ("hnu{}", 0.712),
        ("bdgpq", 0.716),
        ("S", 0.72),
        ("Y", 0.724),
        ("Z", 0.725),
        ("P", 0.733),
        ("C", 0.734),
        ("B", 0.762),
        ("R", 0.77),
        ("X", 0.771),
        ("AV", 0.774),
        ("K", 0.775),
        ("U", 0.812),
        ("G", 0.821),
        ("D", 0.83),
        ("HN", 0.837),
        ("#+<=>^~", 0.838),
        ("OQ", 0.85),
        ("&", 0.872),
        ("w", 0.924),
        ("M", 0.995),
        ("@", 1.0),
        ("%", 1.002),
        ("m", 1.042),
        ("W", 1.103),
    )
    for ch in chars
}
# The typographic minus of the PGS values (#54), in both weights.
_CHAR_EM["\u2212"] = 0.838
_CHAR_EM_BOLD["\u2212"] = 0.838
# A no-break space, which keeps words together when a legend row wraps (#52),
# is as wide as a space.
_CHAR_EM["\u00a0"] = _CHAR_EM[" "]
_CHAR_EM_BOLD["\u00a0"] = _CHAR_EM_BOLD[" "]
_OTHER_CHAR_EM = 0.7
_WIDE_CHAR_EM = 1.0

# Space kept between the labels of two neighbours, in ems.
LABEL_GAP_EM = 0.5

# Notation notes and annotations are drawn below the ID label at this
# fraction of the label's font size.
SMALL_TEXT_SCALE = 0.7

# Other secondary text, as fractions of `font_size` (drawn at
# `PedigreeStyle.text_size(scale)`, never below `min_font_size`, #51).
PGS_CAPTION_SCALE = 0.6  # "PGS" under or above a curve
BAND_TEXT_SCALE = 0.75  # text inside a phenotype band
NA_TEXT_SCALE = 0.85  # the single "NA" of a symbol with no data at all


def text_width(text: str, font_size: float, *, bold: bool = False) -> float:
    """Width of one line of text, in layout units.

    Measured in DejaVu Sans (Bold with ``bold``) and converted with
    `UNITS_PER_POINT`, the scale of the interactive SVG, where text takes
    the most room.

    Args:
        text: The line of text.
        font_size: Font size, in points.
        bold: Measure in DejaVu Sans Bold (band and quadrant text).
    """
    table = _CHAR_EM_BOLD if bold else _CHAR_EM
    em = 0.0
    for ch in text:
        width = table.get(ch)
        if width is None:
            base = unicodedata.normalize("NFD", ch)[:1]
            if base in table:
                width = table[base]
            elif unicodedata.east_asian_width(ch) in ("W", "F"):
                width = _WIDE_CHAR_EM
            else:
                width = _OTHER_CHAR_EM
        em += width
    return em * font_size * UNITS_PER_POINT


def below_label_offsets(
    label_size: float, below_sizes: Sequence[float], units_per_point: float
) -> list[float]:
    """How far below the ID label's top each line under it starts, in layout units.

    Each step is `TEXT_LINE_STEP`, or more when the text needs it (#66): the
    height of the line above plus a gap in ems of the line's own size,
    converted with ``units_per_point`` (`UNITS_PER_POINT` in the SVG, the
    figure's scale when the PNG is drawn). At the default sizes the text
    always fits the fixed step, so nothing moves.

    Args:
        label_size: The ID label's font size, in points.
        below_sizes: The size of each line below it (notes, then annotations),
            in points.
        units_per_point: Layout units per point at the scale drawn.
    """
    offsets: list[float] = []
    offset, above = 0.0, label_size
    for size in below_sizes:
        needed = (TEXT_LINE_HEIGHT_EM * above + TEXT_LINE_GAP_EM * size) * units_per_point
        offset += max(TEXT_LINE_STEP, needed)
        offsets.append(offset)
        above = size
    return offsets


def label_width(
    label: str, below: list[str], font_size: float, below_size: float | None = None
) -> float:
    """Width of the text drawn under a symbol, in layout units.

    Args:
        label: The ID label (empty when there is none: nothing is drawn below).
        below: Notation notes and annotations, one line each, drawn below the
            label in smaller type.
        font_size: The label's font size, in points.
        below_size: The size of the lines below, in points (default:
            `SMALL_TEXT_SCALE` times ``font_size``; the renderers pass
            `PedigreeStyle.text_size(SMALL_TEXT_SCALE)`).
    """
    if not label:
        return 0.0
    if below_size is None:
        below_size = font_size * SMALL_TEXT_SCALE
    widths = [text_width(label, font_size)]
    widths += [text_width(line, below_size) for line in below]
    return max(widths)


_GENDER_SHAPE = {Gender.MAN: Sex.MALE, Gender.WOMAN: Sex.FEMALE, Gender.NONBINARY: Sex.UNKNOWN}
_ASSIGNED_AT_BIRTH = {Sex.MALE: "AMAB", Sex.FEMALE: "AFAB"}
_PREGNANCY_LOSSES = frozenset({Pregnancy.MISCARRIAGE, Pregnancy.TERMINATION, Pregnancy.ECTOPIC})


def symbol_sex(sex: Sex, notation: Notation) -> Sex:
    """The sex whose shape is drawn: square, circle or diamond.

    The 2022 revision draws gender identity, not the sex assigned at birth:
    a square for a man, a circle for a woman, a diamond for a non-binary or
    gender-diverse person.
    """
    if notation.gender is None:
        return sex
    return _GENDER_SHAPE[notation.gender]


def sex_at_birth_note(sex: Sex, notation: Notation) -> str | None:
    """ "AMAB" or "AFAB" when the symbol does not show the sex at birth.

    Nothing is written for a cisgender person, nor when the sex assigned at
    birth is unknown.
    """
    if notation.gender is None or sex not in _ASSIGNED_AT_BIRTH:
        return None
    if _GENDER_SHAPE[notation.gender] == sex:
        return None
    return _ASSIGNED_AT_BIRTH[sex]


def notes(ind: Individual) -> list[str]:
    """Notation text written below the label, top to bottom.

    "SB" for a stillbirth and "ECT" for an ectopic pregnancy, then the sex
    assigned at birth when the symbol does not show it, then the result of
    a documented evaluation (written as given, without the "E" the 2022
    revision dropped).
    """
    out = []
    if ind.notation.pregnancy in (Pregnancy.STILLBIRTH, Pregnancy.ECTOPIC):
        out.append(ind.notation.pregnancy.value)
    note = sex_at_birth_note(ind.sex, ind.notation)
    if note:
        out.append(note)
    if ind.notation.evaluation:
        out.append(ind.notation.evaluation)
    return out


def is_pregnancy_loss(notation: Notation) -> bool:
    """True for a miscarriage, an ectopic pregnancy or a termination: a triangle."""
    return notation.pregnancy in _PREGNANCY_LOSSES


def is_struck_through(ind: Individual) -> bool:
    """Whether the symbol gets the deceased line; a stillbirth always does."""
    return ind.is_deceased or ind.notation.pregnancy == Pregnancy.STILLBIRTH


def inner_letter(notation: Notation) -> str | None:
    """Text written inside the symbol: the number of individuals it stands
    for ("5", "n", #104), else "P" (pregnancy), "D" or "S"."""
    if notation.count:
        return notation.count
    if notation.pregnancy == Pregnancy.ONGOING:
        return "P"
    if notation.reproductive_role is not None:
        return notation.reproductive_role.value
    return None


# The inner letter ("P", "D", "S") is written at the centre of the symbol,
# as NSGC draws it, whatever the phenotype columns (#105). On a symbol divided
# into bands, quadrants or carrier sections it has a white backing; band or
# quadrant text keeps this far from the letter's middle, in ems of the
# letter's font (half the backing's width and height, and a gap).
LETTER_CLEAR_HALF_W_EM = 0.55
LETTER_CLEAR_HALF_H_EM = 0.75
LETTER_CLEAR_GAP_EM = 0.15
# Each character after the first ("12" individuals in one symbol, #104)
# widens it by this much on each side.
LETTER_CLEAR_CHAR_EM = 0.32


def triangle_points(x: float, y: float, half: float, up: float = 1.0) -> list[tuple[float, float]]:
    """Vertices of a pregnancy-loss triangle: apex, base left, base right."""
    base_y = y - up * half * TRIANGLE_BASE
    w = half * TRIANGLE_HALF_WIDTH
    return [(x, y + up * half), (x - w, base_y), (x + w, base_y)]


def termination_stroke(
    x: float, y: float, half: float, up: float = 1.0
) -> tuple[tuple[float, float], tuple[float, float]]:
    """The line struck through a termination's triangle, lower left to upper right."""
    w = half * TRIANGLE_HALF_WIDTH
    return (x - w * 1.1, y - up * half * 0.95), (x + w * 0.8, y + up * half * 0.75)


def marker_points(tip_x: float, tip_y: float, em: float, up: float = 1.0) -> tuple[float, float]:
    """Start of the proband or consultand arrow, in the units of ``em``.

    The SVG passes the font size in layout units (`UNITS_PER_POINT` per
    point); the PNG draws the same offsets in points from the tip. The
    proband's "P" is level with this point, its right edge
    `MARKER_LETTER_GAP_EM` ems further left.

    Args:
        tip_x: X of the point the arrow ends at, the symbol's lower left
            corner (or the outline's point in that direction).
        tip_y: Y of that point.
        em: The label font size, in the units of ``tip_x`` and ``tip_y``.
        up: +1 with the y axis pointing up (PNG), -1 pointing down (SVG).

    Returns:
        (x, y) where the arrow starts.
    """
    return tip_x - em * MARKER_ARROW_EM, tip_y - up * em * MARKER_ARROW_EM


def marker_clearance(half: float, font_size: float, label_chars: int, letter: bool) -> float:
    """Distance to the left neighbour that keeps the marker off their label.

    The marker sits at the lower left of the symbol, level with the labels,
    so with tight spacing or large text the "P" would run into the label of
    the individual on the left. The layout keeps a proband's or a
    consultand's left neighbour at least this far, centre to centre.

    Args:
        half: Half the shape size.
        font_size: Label and marker font size, in points.
        label_chars: Length of the longest label the neighbour may have.
        letter: Whether the marker has a "P" (a proband) or is an arrow
            alone (a consultand).
    """
    em = font_size * UNITS_PER_POINT
    reach = half + em * MARKER_ARROW_EM
    if letter:
        reach += em * (MARKER_LETTER_GAP_EM + _MARKER_LETTER_EM)
    return reach + em * (_MARKER_PAD_EM + _LABEL_CHAR_EM * label_chars / 2)


def shape_half_width_at(sex: Sex, half: float, y: float) -> float:
    """Half the width of an upright symbol at height ``y`` from its centre.

    A square keeps its width, a circle narrows as a chord, the unknown-sex
    diamond linearly to its vertices. Used to fit text inside bands and
    quadrants (#55).

    Args:
        sex: The symbol's shape (`Sex.MALE` square, `FEMALE` circle, else
            the diamond).
        half: Half the shape size.
        y: Vertical distance from the centre, either sign.
    """
    y = abs(y)
    if sex == Sex.MALE:
        return half if y <= half else 0.0
    if sex == Sex.FEMALE:
        return math.sqrt(max(0.0, half * half - y * y))
    return max(0.0, half * DIAMOND_SCALE - y)


def shape_half_width(sex: Sex, half: float) -> float:
    """Horizontal half-extent of a shape (the diamond reaches further)."""
    return half * DIAMOND_SCALE if sex not in (Sex.MALE, Sex.FEMALE) else half


def grown_half(sex: Sex, half: float, distance: float) -> float:
    """Half size of the shape whose outline runs `distance` outside this one's.

    Used for the colored borders drawn around a symbol (#46). The square and
    the circle grow by `distance`; the diamond's sides are at 45°, so its
    vertices move √2 times further.
    """
    if sex in (Sex.MALE, Sex.FEMALE):
        return half + distance
    return half + distance * math.sqrt(2) / DIAMOND_SCALE


def shape_reach(sex: Sex, half: float, direction: float, rotation: float = 0.0) -> float:
    """Distance from a shape's centre to its outline in a given direction.

    Args:
        sex: Decides the shape: square (male), circle (female) or diamond.
        half: Half the shape size.
        direction: Direction to measure along, in radians.
        rotation: Rotation applied to the shape, in radians.
    """
    theta = direction - rotation  # the direction in the shape's own frame
    c, s = abs(math.cos(theta)), abs(math.sin(theta))
    if sex == Sex.MALE:
        return half / max(c, s)
    if sex == Sex.FEMALE:
        return half
    # The diamond is a square standing on a vertex: its outline is |x|+|y| = R.
    return half * DIAMOND_SCALE / (c + s)


def bracket_lines(
    x: float, y: float, half: float, sex: Sex, triangle: bool = False
) -> list[list[tuple[float, float]]]:
    """The two adoption brackets around a symbol, as polylines.

    Symmetric top to bottom, so the same points serve both y directions.
    """
    if triangle:
        reach_x, reach_y = half * TRIANGLE_HALF_WIDTH, half
    else:
        # Square and circle reach `half` both ways; the diamond further.
        reach_x = reach_y = shape_half_width(sex, half)
    left = x - reach_x - BRACKET_GAP
    right = x + reach_x + BRACKET_GAP
    top = y + reach_y + BRACKET_OVERHANG
    bottom = y - reach_y - BRACKET_OVERHANG
    lt, rt = left + BRACKET_TICK, right - BRACKET_TICK
    return [
        [(lt, top), (left, top), (left, bottom), (lt, bottom)],
        [(rt, top), (right, top), (right, bottom), (rt, bottom)],
    ]


def childless_depth(notation: Notation) -> float:
    """Room the childlessness marker takes below the symbol, label gap included.

    None for a couple's childlessness, which hangs from the relationship line
    (#104).
    """
    if notation.childlessness is None or notation.childless_partner is not None:
        return 0.0
    depth = CHILDLESS_STEM + CHILDLESS_LABEL_GAP
    if notation.childlessness == Childlessness.INFERTILITY:
        depth += CHILDLESS_BAR_GAP
    return depth


def childless_segments(
    x: float,
    y: float,
    sex: Sex,
    half: float,
    childlessness: Childlessness,
    direction: float = -math.pi / 2,
) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    """Segments of a childlessness marker: the stem, then one or two bars.

    The stem leaves the symbol in ``direction`` (straight down by default;
    outward in the circular layout) and the bars cross its end.
    """
    c, s = math.cos(direction), math.sin(direction)
    start = shape_reach(sex, half, direction)
    end = start + CHILDLESS_STEM
    w = half * CHILDLESS_BAR_HALF_WIDTH

    def at(d: float, across: float) -> tuple[float, float]:
        # `across` runs perpendicular to the stem.
        return (x + d * c - across * s, y + d * s + across * c)

    segments = [(at(start, 0.0), at(end, 0.0)), (at(end, -w), at(end, w))]
    if childlessness == Childlessness.INFERTILITY:
        d = end + CHILDLESS_BAR_GAP
        segments.append((at(d, -w), at(d, w)))
    return segments


def couple_childless_segments(
    x: float,
    y: float,
    half: float,
    childlessness: Childlessness,
    direction: float = -math.pi / 2,
) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    """A couple's childlessness marker (#104): from the middle of their
    relationship line at (``x``, ``y``), where a descent line would start, a
    stem ending in one or two bars, as `childless_segments` draws it under
    one symbol.

    The bars hang as low as an individual's marker under a square, below the
    partners' symbols (`half` + `CHILDLESS_STEM` from the line): just under
    the line they filled the gap between the partners and read as a second
    relationship line, a consanguineous couple (#124)."""
    c, s = math.cos(direction), math.sin(direction)
    end = half + CHILDLESS_STEM
    w = half * CHILDLESS_BAR_HALF_WIDTH

    def at(d: float, across: float) -> tuple[float, float]:
        return (x + d * c - across * s, y + d * s + across * c)

    segments = [(at(0.0, 0.0), at(end, 0.0)), (at(end, -w), at(end, w))]
    if childlessness == Childlessness.INFERTILITY:
        d = end + CHILDLESS_BAR_GAP
        segments.append((at(d, -w), at(d, w)))
    return segments


# The asymptomatic or presymptomatic carrier's line (#104): vertical, through
# the symbol's centre, and past its top and bottom by this many half-sizes,
# so it is not taken for the line between two quadrants or carrier sections.
ASYMPTOMATIC_OVERSHOOT = 0.2


def asymptomatic_line(
    x: float, y: float, half: float, sex: Sex
) -> tuple[tuple[float, float], tuple[float, float]]:
    """The vertical line of an asymptomatic or presymptomatic carrier (#104).

    Symmetric top to bottom, so the same points serve both y directions.
    """
    reach = shape_half_width(sex, half) + half * ASYMPTOMATIC_OVERSHOOT
    return (x, y - reach), (x, y + reach)


def unknown_offspring_anchor(
    x: float, y: float, half: float, up: float = 1.0
) -> tuple[float, float]:
    """Top centre of the "?" of an unknown child: just below its descent line.

    The descent line stops where the top of a symbol would be; the "?" hangs
    from there, anchored by its top edge.
    """
    return x, y + up * (half - UNKNOWN_OFFSPRING_GAP)


# The PGS curve's box beside a symbol: its width and its gap from the symbol,
# in symbol half-sizes. It goes on the left, or on the right for the
# right-hand partner of a couple (`NodeLayout.pgs_side`, #78), off the line
# between the partners.
PGS_BOX_WIDTH = 0.8
PGS_BOX_GAP = 0.15

# The deceased line runs corner to corner and this far past each corner, in
# half-sizes; a PGS box beside a struck-through symbol keeps this much more
# from it, so the line's end does not reach the box (#80).
DECEASED_OVERSHOOT = 0.2
PGS_DECEASED_MARGIN = 0.12


def pgs_box_near(
    half: float,
    sex: Sex,
    evaluated: bool,
    mark_size: float,
    side: str,
    *,
    struck: bool = False,
    ring: float = 0.0,
    line_width: float = 0.0,
) -> float:
    """From a symbol's centre to the near edge of its PGS box, in layout units.

    The gap is `PGS_BOX_GAP` half-sizes beyond the shape's reach on that
    side: a diamond reaches `DIAMOND_SCALE` half-sizes (#111), and outer
    borders (``ring``, #46) grow the shape with `grown_half`. Beside a
    struck-through (deceased) symbol, the box also clears the line's end,
    `DECEASED_OVERSHOOT` past the corner, by `PGS_DECEASED_MARGIN` (#80). A
    box on the right also clears the documented-evaluation asterisk at the
    symbol's lower right (#78), measured at the SVG's scale.

    Args:
        half: Half the symbol's size.
        sex: The shape drawn.
        evaluated: Whether the individual has the asterisk.
        mark_size: The asterisk's font size, in points.
        side: "left" or "right" (`NodeLayout.pgs_side`).
        struck: Whether the symbol is struck through (deceased).
        ring: Room the borders after the first take, in layout units.
        line_width: Width of the outermost outline, in points.
    """
    grown = grown_half(sex, half, ring) if ring else half
    near = shape_half_width(sex, grown) + half * PGS_BOX_GAP
    if struck:
        near = max(near, ring + half * (1 + DECEASED_OVERSHOOT + PGS_DECEASED_MARGIN))
    if side == "right" and evaluated:
        mark_x, _ = evaluation_mark_position(
            0.0, 0.0, grown, sex, mark_size, line_width, UNITS_PER_POINT
        )
        mark_right = mark_x + EVALUATION_MARK_INK_HALF_W_EM * mark_size * UNITS_PER_POINT
        near = max(near, mark_right + half * PGS_BOX_GAP)
    return near


def evaluation_mark_placement(
    x: float,
    y: float,
    half: float,
    sex: Sex,
    size: float,
    line_width: float,
    triangle: bool = False,
    up: float = 1.0,
) -> tuple[float, float, float, float]:
    """Where the documented-evaluation asterisk goes, beside the lower right.

    Bennett et al. 2022, Figure 2 ("Documented evaluation, records
    reviewed"): just outside the symbol's lower-right corner, on the
    diagonal for a circle or the lower-right edge of a diamond; for a
    pregnancy loss, beside the right end of the triangle's base. The point
    is on the outline, in layout units; the offset, in points, goes from
    there to the middle of the asterisk's ink: half the outline's stroke,
    the ink's half-extent in that direction (more on the diagonal) and
    `EVALUATION_MARK_GAP_EM` (#110). The PNG draws it as an offset in points,
    the SVG converts it with its scale (`evaluation_mark_position`).

    Args:
        x: X of the symbol's centre.
        y: Y of the symbol's centre.
        half: Half the shape size (grown round outer borders, if any).
        sex: The shape drawn.
        size: The asterisk's font size, in points.
        line_width: Width of the outline it sits beside, in points.
        triangle: Whether the symbol is a pregnancy-loss triangle.
        up: +1 with the y axis pointing up (PNG), -1 pointing down (SVG).

    Returns:
        (x, y) on the outline, then (dx, dy) in points, in the same y
        direction.
    """
    if triangle or sex == Sex.MALE:
        if triangle:
            px, py = x + half * TRIANGLE_HALF_WIDTH, y - up * half * TRIANGLE_BASE
        else:
            # Beside the right side, near the bottom, as in the figure.
            px, py = x + half, y - up * half * 0.8
        dx, dy, extent = 1.0, 0.0, EVALUATION_MARK_INK_HALF_W_EM
    else:
        d = shape_reach(sex, half, -math.pi / 4) / math.sqrt(2)
        px, py = x + d, y - up * d
        dx = dy = 1 / math.sqrt(2)
        dy = -up * dy
        extent = (EVALUATION_MARK_INK_HALF_W_EM + EVALUATION_MARK_INK_HALF_H_EM) / math.sqrt(2)
    distance = line_width / 2 + (extent + EVALUATION_MARK_GAP_EM) * size
    return px, py, dx * distance, dy * distance


def evaluation_mark_position(
    x: float,
    y: float,
    half: float,
    sex: Sex,
    size: float,
    line_width: float,
    units_per_point: float,
    triangle: bool = False,
    up: float = 1.0,
) -> tuple[float, float]:
    """Middle of the asterisk's ink at a given scale (`evaluation_mark_placement`).

    Args:
        units_per_point: Layout units per point at the drawing's scale (the
            SVG's `UNITS_PER_POINT`).
    """
    px, py, dx, dy = evaluation_mark_placement(x, y, half, sex, size, line_width, triangle, up)
    return px + dx * units_per_point, py + dy * units_per_point


def text_lines_below(ind: Individual) -> int:
    """Lines of text below the label: notes, then annotations."""
    return len(notes(ind)) + sum(1 for v in ind.annotations.values() if v)
