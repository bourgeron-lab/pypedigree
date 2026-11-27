"""One file per family, for a cohort (#49).

`render_families` draws each family of a pedigree into its own file, with
the colors, column modes and legend of the whole file
(`Pedigree.subset` keeps it as `Pedigree.column_scope`), so the figures of
one cohort can be compared. A family that fails is reported and the others
are still drawn. matplotlib is imported only for PNG, PDF and matplotlib SVG.
"""

from __future__ import annotations

import warnings
from collections.abc import Iterable
from pathlib import Path

from pypedigree.carriers import check_notation
from pypedigree.config import PedigreeStyle, check_text_scale, resolve_style
from pypedigree.errors import PedigreeBatchError, PedigreeDataWarning, PedigreeStyleWarning
from pypedigree.models import Pedigree
from pypedigree.phenotypes import (
    check_band_text,
    check_curve_options,
    check_gradient_bins,
    check_legend,
    check_min_font_size,
    check_pgs_scale,
    check_phenotype_layout,
    validate_style_colors,
)
from pypedigree.transmission import check_transmission

FORMATS = ("png", "pdf", "svg")

# Device names Windows reserves whatever the extension ("CON.png" cannot be
# created there).
_RESERVED = {"con", "prn", "aux", "nul"} | {f"{p}{i}" for p in ("com", "lpt") for i in range(1, 10)}
_MAX_NAME = 100


def family_file_names(family_ids: Iterable[str]) -> dict[str, str]:
    """A distinct, portable file name (without extension) for each family.

    Letters, digits, ``-``, ``_`` and ``.`` are kept and every other
    character becomes ``_``; leading and trailing dots and spaces are
    dropped, Windows device names (``CON``, ``COM1``...) get a trailing
    ``_``, and names are cut to 100 characters. Two IDs giving the same
    name, ignoring case (``F/1`` and ``F_1``, or ``f1`` and ``F1`` on a
    file system that ignores case), are told apart by ``_2``, ``_3``... in
    file order, with a `PedigreeDataWarning`.

    Args:
        family_ids: Family IDs, in the order of the files.

    Returns:
        Family ID → file name.
    """
    names: dict[str, str] = {}
    taken: set[str] = set()
    for fid in dict.fromkeys(family_ids):
        base = "".join(c if c.isalnum() or c in "-_." else "_" for c in fid)
        base = base.strip(". ")[:_MAX_NAME] or "family"
        if base.split(".")[0].lower() in _RESERVED:
            base += "_"
        name, n = base, 1
        while name.lower() in taken:
            n += 1
            name = f"{base}_{n}"
        if name != base:
            warnings.warn(
                f"Family '{fid}' is written as '{name}': '{base}' is the file name "
                "of an earlier family.",
                PedigreeDataWarning,
                stacklevel=2,
            )
        taken.add(name.lower())
        names[fid] = name
    return names


def render_families(
    pedigree: Pedigree,
    output_dir: str | Path,
    *,
    fmt: str = "png",
    family_ids: Iterable[str] | None = None,
    style: PedigreeStyle | None = None,
    show_labels: bool = True,
    show_legend: bool = True,
    circular: bool = False,
    interactive: bool = False,
    script: bool = True,
    id_prefix: str = "",
    responsive: bool = False,
) -> dict[str, Path]:
    """Draw each family into its own file in ``output_dir``.

    Each file is named after its family (`family_file_names`). Every family
    is drawn with the colors, column modes and legend of the whole
    pedigree, so the figures of one cohort can be compared. The directory
    is created if needed; existing files are overwritten.

    Args:
        pedigree: The parsed pedigree.
        output_dir: Directory to write into.
        fmt: ``"png"``, ``"pdf"`` or ``"svg"`` (matplotlib's SVG). Ignored
            with ``interactive``.
        family_ids: Families to draw, in order (default: all).
        style: Visual style parameters. Uses defaults if None.
        show_labels: Whether to draw sample ID labels.
        show_legend: Whether to draw the legend.
        circular: Circular layout (PNG, PDF and matplotlib SVG only).
        interactive: Write pypedigree's SVG (`render_interactive_svg`)
            instead of a matplotlib figure.
        script: With ``interactive``: False writes the same SVG without
            script or tooltip (`render_svg`).
        id_prefix: With ``interactive``: prefix of every id in each SVG.
        responsive: With ``interactive``: ``width="100%"`` and no height.

    Returns:
        Family ID → path of its file, in drawing order.

    Raises:
        PedigreeBatchError: If some families could not be drawn, after the
            others were; it lists the files written and the failures.
        PedigreeError: If the style is invalid (nothing is drawn).
        KeyError: If a family ID is not in the pedigree.
        ValueError: For an unknown format or an invalid ``id_prefix``.
    """
    ids = list(dict.fromkeys(pedigree.family_ids if family_ids is None else family_ids))
    missing = [f for f in ids if f not in pedigree.families]
    if missing:
        raise KeyError(f"No family {', '.join(map(repr, missing))} in this pedigree")
    ext = "svg" if interactive else fmt.lower()
    if ext not in FORMATS:
        raise ValueError(f"Unknown format {fmt!r}: expected one of {', '.join(FORMATS)}")

    # Checked once, not once per family: an invalid style fails them all.
    style = resolve_style(style, pedigree.phenotype_names, pedigree.phenotype_display_modes)
    validate_style_colors(style)
    check_notation(style.carrier_notation)
    check_phenotype_layout(style)
    check_pgs_scale(style)
    check_curve_options(style)
    check_min_font_size(style)
    check_text_scale(style)
    check_band_text(style)
    check_legend(style)
    check_gradient_bins(style)
    check_transmission(style)
    if interactive:
        from pypedigree.svg_renderer import check_id_prefix, render_interactive_svg, render_svg

        check_id_prefix(id_prefix)
        if circular:
            warnings.warn(
                "The circular layout is not available in the interactive SVG; "
                "drawing the hierarchical layout.",
                PedigreeStyleWarning,
                stacklevel=2,
            )
        svg = render_interactive_svg if script else render_svg

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    names = family_file_names(ids)
    written: dict[str, Path] = {}
    failures: dict[str, str] = {}
    for fid in ids:
        path = out / f"{names[fid]}.{ext}"
        one = pedigree.subset([fid])
        try:
            if interactive:
                svg(
                    one,
                    style=style,
                    output=path,
                    show_labels=show_labels,
                    show_legend=show_legend,
                    id_prefix=id_prefix,
                    responsive=responsive,
                )
            else:
                _save_figure(one, path, style, show_labels, show_legend, circular)
        except Exception as exc:  # one family's failure must not stop the others
            failures[fid] = str(exc) or type(exc).__name__
        else:
            written[fid] = path
    if failures:
        raise PedigreeBatchError(written, failures)
    return written


def _save_figure(
    pedigree: Pedigree,
    path: Path,
    style: PedigreeStyle,
    show_labels: bool,
    show_legend: bool,
    circular: bool,
) -> None:
    """Draw one pedigree with matplotlib, save it and close its figures."""
    import matplotlib.pyplot as plt

    from pypedigree.renderer import render_pedigree

    before = set(plt.get_fignums())
    try:
        render_pedigree(
            pedigree,
            style=style,
            output=path,
            show_labels=show_labels,
            show_legend=show_legend,
            circular=circular,
        )
    finally:
        # Hundreds of open figures would exhaust memory; the caller's own
        # figures are left alone.
        for num in set(plt.get_fignums()) - before:
            plt.close(num)
