"""Color parsing without matplotlib.

The SVG writer, the layout and the parser must not import matplotlib (its
import is slow and ``pyplot`` holds global state, see
bourgeron-lab/pypedigree#8). This module reads the colors a style can hold
the way ``matplotlib.colors.to_rgb`` does, with the same values:

- hex: ``#rrggbb``, ``#rgb``, ``#rrggbbaa``, ``#rgba`` (alpha ignored);
- names: the 148 CSS4 names (``grey`` and ``gray`` spellings), the
  single letters ``b g r c m y k w`` and the ``tab:`` palette, case
  ignored except for single letters;
- a grey level as a string, ``"0"`` to ``"1"``;
- an RGB or RGBA tuple of floats in [0, 1].

Anything else (``xkcd:`` names, ``C0`` cycle colors) is handed to
matplotlib when it is installed, so a style accepts the same colors as
before; without matplotlib it is an invalid color.
"""

from __future__ import annotations

import re
from collections.abc import Sequence

RGB = tuple[float, float, float]

# matplotlib.colors.BASE_COLORS: RGB tuples, not hex (green is 0.5, not 128/255).
_BASE_COLORS: dict[str, RGB] = {
    "b": (0.0, 0.0, 1.0),
    "g": (0.0, 0.5, 0.0),
    "r": (1.0, 0.0, 0.0),
    "c": (0.0, 0.75, 0.75),
    "m": (0.75, 0.0, 0.75),
    "y": (0.75, 0.75, 0.0),
    "k": (0.0, 0.0, 0.0),
    "w": (1.0, 1.0, 1.0),
}

# matplotlib.colors.TABLEAU_COLORS and CSS4_COLORS.
_NAMED_COLORS: dict[str, str] = {
    "tab:blue": "#1f77b4",
    "tab:orange": "#ff7f0e",
    "tab:green": "#2ca02c",
    "tab:red": "#d62728",
    "tab:purple": "#9467bd",
    "tab:brown": "#8c564b",
    "tab:pink": "#e377c2",
    "tab:gray": "#7f7f7f",
    "tab:olive": "#bcbd22",
    "tab:cyan": "#17becf",
    "aliceblue": "#f0f8ff",
    "antiquewhite": "#faebd7",
    "aqua": "#00ffff",
    "aquamarine": "#7fffd4",
    "azure": "#f0ffff",
    "beige": "#f5f5dc",
    "bisque": "#ffe4c4",
    "black": "#000000",
    "blanchedalmond": "#ffebcd",
    "blue": "#0000ff",
    "blueviolet": "#8a2be2",
    "brown": "#a52a2a",
    "burlywood": "#deb887",
    "cadetblue": "#5f9ea0",
    "chartreuse": "#7fff00",
    "chocolate": "#d2691e",
    "coral": "#ff7f50",
    "cornflowerblue": "#6495ed",
    "cornsilk": "#fff8dc",
    "crimson": "#dc143c",
    "cyan": "#00ffff",
    "darkblue": "#00008b",
    "darkcyan": "#008b8b",
    "darkgoldenrod": "#b8860b",
    "darkgray": "#a9a9a9",
    "darkgreen": "#006400",
    "darkgrey": "#a9a9a9",
    "darkkhaki": "#bdb76b",
    "darkmagenta": "#8b008b",
    "darkolivegreen": "#556b2f",
    "darkorange": "#ff8c00",
    "darkorchid": "#9932cc",
    "darkred": "#8b0000",
    "darksalmon": "#e9967a",
    "darkseagreen": "#8fbc8f",
    "darkslateblue": "#483d8b",
    "darkslategray": "#2f4f4f",
    "darkslategrey": "#2f4f4f",
    "darkturquoise": "#00ced1",
    "darkviolet": "#9400d3",
    "deeppink": "#ff1493",
    "deepskyblue": "#00bfff",
    "dimgray": "#696969",
    "dimgrey": "#696969",
    "dodgerblue": "#1e90ff",
    "firebrick": "#b22222",
    "floralwhite": "#fffaf0",
    "forestgreen": "#228b22",
    "fuchsia": "#ff00ff",
    "gainsboro": "#dcdcdc",
    "ghostwhite": "#f8f8ff",
    "gold": "#ffd700",
    "goldenrod": "#daa520",
    "gray": "#808080",
    "green": "#008000",
    "greenyellow": "#adff2f",
    "grey": "#808080",
    "honeydew": "#f0fff0",
    "hotpink": "#ff69b4",
    "indianred": "#cd5c5c",
    "indigo": "#4b0082",
    "ivory": "#fffff0",
    "khaki": "#f0e68c",
    "lavender": "#e6e6fa",
    "lavenderblush": "#fff0f5",
    "lawngreen": "#7cfc00",
    "lemonchiffon": "#fffacd",
    "lightblue": "#add8e6",
    "lightcoral": "#f08080",
    "lightcyan": "#e0ffff",
    "lightgoldenrodyellow": "#fafad2",
    "lightgray": "#d3d3d3",
    "lightgreen": "#90ee90",
    "lightgrey": "#d3d3d3",
    "lightpink": "#ffb6c1",
    "lightsalmon": "#ffa07a",
    "lightseagreen": "#20b2aa",
    "lightskyblue": "#87cefa",
    "lightslategray": "#778899",
    "lightslategrey": "#778899",
    "lightsteelblue": "#b0c4de",
    "lightyellow": "#ffffe0",
    "lime": "#00ff00",
    "limegreen": "#32cd32",
    "linen": "#faf0e6",
    "magenta": "#ff00ff",
    "maroon": "#800000",
    "mediumaquamarine": "#66cdaa",
    "mediumblue": "#0000cd",
    "mediumorchid": "#ba55d3",
    "mediumpurple": "#9370db",
    "mediumseagreen": "#3cb371",
    "mediumslateblue": "#7b68ee",
    "mediumspringgreen": "#00fa9a",
    "mediumturquoise": "#48d1cc",
    "mediumvioletred": "#c71585",
    "midnightblue": "#191970",
    "mintcream": "#f5fffa",
    "mistyrose": "#ffe4e1",
    "moccasin": "#ffe4b5",
    "navajowhite": "#ffdead",
    "navy": "#000080",
    "oldlace": "#fdf5e6",
    "olive": "#808000",
    "olivedrab": "#6b8e23",
    "orange": "#ffa500",
    "orangered": "#ff4500",
    "orchid": "#da70d6",
    "palegoldenrod": "#eee8aa",
    "palegreen": "#98fb98",
    "paleturquoise": "#afeeee",
    "palevioletred": "#db7093",
    "papayawhip": "#ffefd5",
    "peachpuff": "#ffdab9",
    "peru": "#cd853f",
    "pink": "#ffc0cb",
    "plum": "#dda0dd",
    "powderblue": "#b0e0e6",
    "purple": "#800080",
    "rebeccapurple": "#663399",
    "red": "#ff0000",
    "rosybrown": "#bc8f8f",
    "royalblue": "#4169e1",
    "saddlebrown": "#8b4513",
    "salmon": "#fa8072",
    "sandybrown": "#f4a460",
    "seagreen": "#2e8b57",
    "seashell": "#fff5ee",
    "sienna": "#a0522d",
    "silver": "#c0c0c0",
    "skyblue": "#87ceeb",
    "slateblue": "#6a5acd",
    "slategray": "#708090",
    "slategrey": "#708090",
    "snow": "#fffafa",
    "springgreen": "#00ff7f",
    "steelblue": "#4682b4",
    "tan": "#d2b48c",
    "teal": "#008080",
    "thistle": "#d8bfd8",
    "tomato": "#ff6347",
    "turquoise": "#40e0d0",
    "violet": "#ee82ee",
    "wheat": "#f5deb3",
    "white": "#ffffff",
    "whitesmoke": "#f5f5f5",
    "yellow": "#ffff00",
    "yellowgreen": "#9acd32",
}

_HEX_RE = re.compile(r"#(?:[0-9a-fA-F]{3,4}|[0-9a-fA-F]{6}|[0-9a-fA-F]{8})")


def _from_hex(value: str) -> RGB:
    digits = value[1:]
    if len(digits) in (3, 4):
        digits = "".join(d * 2 for d in digits)
    r, g, b = (int(digits[i : i + 2], 16) / 255 for i in (0, 2, 4))
    return (r, g, b)


def _from_string(value: str) -> RGB | None:
    if value.lower() == "none":
        return (0.0, 0.0, 0.0)
    if value in _BASE_COLORS:
        return _BASE_COLORS[value]
    named = _NAMED_COLORS.get(value) or (
        _NAMED_COLORS.get(value.lower()) if len(value) != 1 else None
    )
    if named is not None:
        return _from_hex(named)
    if _HEX_RE.fullmatch(value):
        return _from_hex(value)
    try:
        grey = float(value)
    except ValueError:
        return None
    if not 0 <= grey <= 1:
        raise ValueError(f"Invalid string grayscale value {value!r}: must be within 0-1")
    return (grey, grey, grey)


def _from_tuple(value: Sequence[object]) -> RGB:
    if len(value) not in (3, 4):
        raise ValueError(f"RGBA sequence should have length 3 or 4, not {len(value)}")
    channels: list[float] = []
    for channel in value:
        if isinstance(channel, bool) or not isinstance(channel, (int, float)):
            raise ValueError(f"Invalid RGBA argument: {tuple(value)!r}")
        channels.append(float(channel))
    if any(not 0 <= c <= 1 for c in channels):
        raise ValueError("RGBA values should be within 0-1 range")
    return (channels[0], channels[1], channels[2])


def parse_rgb(value: object) -> RGB:
    """(r, g, b) floats in [0, 1] for a color, as ``matplotlib.colors.to_rgb``.

    Raises:
        ValueError: If the color is not recognised.
    """
    if isinstance(value, str):
        rgb = _from_string(value)
        if rgb is not None:
            return rgb
        return _matplotlib_rgb(value)
    if isinstance(value, (tuple, list)):
        return _from_tuple(value)
    raise ValueError(f"Invalid RGBA argument: {value!r}")


def _matplotlib_rgb(value: str) -> RGB:
    """Colors only matplotlib knows (``xkcd:``, ``C0``), if it is installed."""
    try:
        import matplotlib.colors as mcolors
    except ImportError:
        raise ValueError(f"Invalid RGBA argument: {value!r}") from None
    r, g, b = mcolors.to_rgb(value)
    return (float(r), float(g), float(b))


def rgb_to_hex(rgb: Sequence[float]) -> str:
    """``#rrggbb`` for an RGB(A) tuple, rounded as ``matplotlib.colors.to_hex``."""
    return "#" + "".join(format(round(c * 255), "02x") for c in rgb[:3])
