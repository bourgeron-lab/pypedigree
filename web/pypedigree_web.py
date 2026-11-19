"""Glue between the web page and pypedigree, run by Pyodide in the browser.

The page calls `render` with the text of a pedigree file and the options
chosen in the form, and gets back JSON: the interactive SVG to show, the
static SVG to save or turn into a PNG, the families found and any warnings.
Nothing here touches the network or the page: the browser does that.

Only the parser, the layout and the SVG writer are used, so matplotlib,
numpy and click are never imported (pypedigree is installed without its
dependencies). The module also runs under CPython, where the tests use it.
"""

from __future__ import annotations

import json
import warnings
from typing import Any

import pypedigree
from pypedigree import PedigreeStyle, parse_pedigree_from_string
from pypedigree.errors import PedigreeError
from pypedigree.svg_renderer import render_interactive_svg, render_svg

#: Form options that map one to one onto `PedigreeStyle` fields, with the
#: values the page may send.
STYLE_CHOICES: dict[str, tuple[str, ...]] = {
    "phenotype_layout": ("auto", "quadrants", "bands"),
    "pgs_scale": ("minmax", "zscore", "percentile"),
    "carrier_notation": ("nsgc2022", "legacy"),
}
#: Prefix of every id in the interactive SVG, so it cannot clash with the page.
ID_PREFIX = "pv-"


def version() -> str:
    """The version of pypedigree that draws the figures."""
    return pypedigree.__version__


def render(text: str, options_json: str = "{}") -> str:
    """Draw a pedigree file.

    Args:
        text: Content of a pedigree file (.ped, .fam, TSV, CSV, or space
            separated, with or without a header).
        options_json: JSON object with any of: ``family`` (an ID, or empty
            for all), ``phenotype_layout``, ``pgs_scale``,
            ``carrier_notation``, ``min_font_size`` (points, or empty for
            none), ``show_labels``, ``show_legend``, ``show_na_labels`` and
            ``style`` (a JSON object of
            `PedigreeStyle` fields, as in a ``--style`` file).

    Returns:
        A JSON object. On success: ``ok`` true, ``families`` (all family IDs
        in the file), ``family`` (the one drawn, or empty for all),
        ``individuals`` (how many rows were drawn: phantom parents and
        duplicates are not counted), ``svg`` (interactive),
        ``static_svg`` (no script) and ``warnings``. On failure: ``ok``
        false, ``error`` and ``warnings``.
    """
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            result = _render(text, json.loads(options_json or "{}"))
        except (PedigreeError, ValueError, TypeError, KeyError) as exc:
            result = {"ok": False, "error": _message(exc)}
    result["warnings"] = list(dict.fromkeys(str(w.message) for w in caught))
    return json.dumps(result)


def _render(text: str, options: dict[str, Any]) -> dict[str, Any]:
    if not text.strip():
        raise ValueError("The file is empty.")
    pedigree = parse_pedigree_from_string(text)
    families = sorted(pedigree.families)
    if not families:
        raise ValueError("No individual was found in the file.")
    family = str(options.get("family") or "")
    if family:
        if family not in pedigree.families:
            raise ValueError(f"Family '{family}' not found. Available: {', '.join(families)}")
        pedigree = pedigree.subset([family])

    style = _style(options)
    kwargs = {
        "style": style,
        "show_labels": bool(options.get("show_labels", True)),
        "show_legend": bool(options.get("show_legend", True)),
    }
    return {
        "ok": True,
        "families": families,
        "family": family,
        "individuals": sum(
            1
            for f in pedigree.families.values()
            for i in f.individuals.values()
            if not (i.is_phantom or i.is_duplicate or i.is_unknown_offspring)
        ),
        "svg": render_interactive_svg(pedigree, id_prefix=ID_PREFIX, responsive=True, **kwargs),
        "static_svg": render_svg(pedigree, **kwargs),
    }


def _style(options: dict[str, Any]) -> PedigreeStyle:
    """The style from the advanced JSON, then the form's own choices."""
    overrides = options.get("style") or {}
    if isinstance(overrides, str):
        overrides = json.loads(overrides) if overrides.strip() else {}
    if not isinstance(overrides, dict):
        raise ValueError("The style must be a JSON object, as in a --style file.")
    style = PedigreeStyle.from_dict(overrides)
    for name, allowed in STYLE_CHOICES.items():
        value = options.get(name)
        if value is None:
            continue
        if value not in allowed:
            raise ValueError(f"Unknown {name} '{value}'. Choose one of: {', '.join(allowed)}")
        setattr(style, name, value)
    if "show_na_labels" in options:
        style.show_na_labels = bool(options["show_na_labels"])
    size = options.get("min_font_size")
    if size not in (None, ""):
        try:
            style.min_font_size = float(size)
        except (TypeError, ValueError):
            raise ValueError(
                f"The smallest text size must be a number of points, not '{size}'."
            ) from None
    return style


def _message(exc: Exception) -> str:
    if isinstance(exc, json.JSONDecodeError):
        return f"The style is not valid JSON: {exc.msg} (line {exc.lineno}, column {exc.colno})."
    return str(exc) or type(exc).__name__
