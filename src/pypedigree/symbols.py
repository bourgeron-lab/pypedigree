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

# Vertical step between lines of text below the label (notes, annotations).
TEXT_LINE_STEP = 0.25

# Documented evaluation: an asterisk this far outside the symbol's lower
# right, at this size relative to the style's font size.
EVALUATION_MARK = "*"
EVALUATION_MARK_GAP = 0.08
EVALUATION_MARK_FONT_SCALE = 1.3

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

# Layout units per typographic point in the interactive SVG. The PNG draws a
# large pedigree at 0.8 inch (57.6 pt) per unit and a small one larger, so
# text never takes more room there: layout decisions that depend on the size
# of text convert it with this.
UNITS_PER_POINT = 0.02

# Text widths, in ems, for the room the proband marker needs: the bold "P"
# (0.73 in DejaVu Sans Bold), one character of a label (digits are 0.64),
# and the space kept between the "P" and the neighbour's label.
_MARKER_LETTER_EM = 0.8
_LABEL_CHAR_EM = 0.65
_MARKER_PAD_EM = 0.3

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
    """Letter written inside the symbol: "P" (pregnancy), "D" or "S"."""
    if notation.pregnancy == Pregnancy.ONGOING:
        return "P"
    if notation.reproductive_role is not None:
        return notation.reproductive_role.value
    return None


def letter_position(
    x: float, y: float, half: float, sex: Sex, outside: bool, up: float = 1.0
) -> tuple[float, float]:
    """Where the inner letter goes: the centre, or beside the upper right.

    A symbol divided into phenotype bands has text of its own in the
    middle, so the letter moves outside, next to its upper-right corner;
    the renderer anchors it there by its lower-left corner.
    """
    if not outside:
        return (x, y)
    return (x + shape_half_width(sex, half) + LETTER_OFFSET, y + up * half * 0.55)


# Gap between a symbol and a letter written beside it.
LETTER_OFFSET = 0.04


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


def shape_half_width(sex: Sex, half: float) -> float:
    """Horizontal half-extent of a shape (the diamond reaches further)."""
    return half * DIAMOND_SCALE if sex not in (Sex.MALE, Sex.FEMALE) else half


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
    """Room the childlessness marker takes below the symbol, label gap included."""
    if notation.childlessness is None:
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


def unknown_offspring_anchor(
    x: float, y: float, half: float, up: float = 1.0
) -> tuple[float, float]:
    """Top centre of the "?" of an unknown child: just below its descent line.

    The descent line stops where the top of a symbol would be; the "?" hangs
    from there, anchored by its top edge.
    """
    return x, y + up * (half - UNKNOWN_OFFSPRING_GAP)


def evaluation_mark_position(
    x: float, y: float, half: float, sex: Sex, triangle: bool = False, up: float = 1.0
) -> tuple[float, float]:
    """Centre of the documented-evaluation asterisk: beside the lower right.

    Bennett et al. 2022, Figure 2 ("Documented evaluation, records
    reviewed"): just outside the symbol's lower-right corner, or the
    lower-right edge of a diamond; for a pregnancy loss, beside the right
    end of the triangle's base.
    """
    if triangle:
        return (
            x + half * TRIANGLE_HALF_WIDTH + EVALUATION_MARK_GAP,
            y - up * half * TRIANGLE_BASE,
        )
    if sex == Sex.MALE:
        # Beside the right side, near the bottom, as in the figure.
        return x + half + EVALUATION_MARK_GAP, y - up * half * 0.8
    reach = shape_reach(sex, half, -math.pi / 4) + EVALUATION_MARK_GAP
    d = reach / math.sqrt(2)
    return x + d, y - up * d


def text_lines_below(ind: Individual) -> int:
    """Lines of text below the label: notes, then annotations."""
    return len(notes(ind)) + sum(1 for v in ind.annotations.values() if v)
