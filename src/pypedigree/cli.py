"""Command-line interface for pypedigree."""

from __future__ import annotations

import json
import warnings
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

import click

from pypedigree import __version__
from pypedigree.config import PedigreeStyle, resolve_style
from pypedigree.errors import PedigreeError
from pypedigree.parser import (
    apply_phenotypes_to_pedigree,
    parse_fam_file,
    parse_pedigree,
    parse_phenotype_file,
    parse_plink_ped,
)


@contextmanager
def _reported_warnings() -> Iterator[None]:
    """Report the warnings raised inside the block as ``Warning: <message>``.

    Parsing, style resolution and rendering all warn through Python's
    `warnings` (`PedigreeDataWarning`, `PedigreeStyleWarning`); left alone,
    they would print as raw Python warnings carrying an internal path and an
    unrelated source line. Each message is written once, on stderr, when the
    block ends, also when it ends with an error.
    """
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            yield
        finally:
            for message in dict.fromkeys(str(entry.message) for entry in caught):
                click.echo(f"Warning: {message}", err=True)


@click.command()
@click.argument("pedigree_file", type=click.Path(exists=True))
@click.option("-o", "--output", type=click.Path(), default=None, help="Output file path.")
@click.option(
    "--format",
    "fmt",
    type=click.Choice(["png", "svg", "pdf"], case_sensitive=False),
    default=None,
    help="Output format (auto-detected from extension).",
)
@click.option("--family", default=None, help="Render only this family ID.")
@click.option("--dpi", type=int, default=None, help="Output resolution (default: 150).")
@click.option("--show", is_flag=True, help="Display interactively instead of saving.")
@click.option(
    "--style",
    "style_file",
    type=click.Path(exists=True),
    default=None,
    help="JSON file with style overrides.",
)
@click.option(
    "--separator",
    type=click.Choice(["tab", "comma", "space"]),
    default=None,
    help="Force field separator.",
)
@click.option(
    "--pheno",
    "pheno_file",
    type=click.Path(exists=True),
    default=None,
    help="Separate .pheno TSV file with extra phenotype data.",
)
@click.option(
    "--annotations",
    "annotations_file",
    type=click.Path(exists=True),
    default=None,
    help="Separate TSV file with text annotations to display below individuals.",
)
@click.option(
    "--outline-color",
    default=None,
    help="Shape and deceased marker outline color (e.g. 'grey', '#888888'). Default: black.",
)
@click.option(
    "--carrier-notation",
    type=click.Choice(["nsgc2022", "legacy"]),
    default=None,
    help=(
        "Carrier status notation. nsgc2022 (default): the symbol is divided into "
        "patterned sections, one per carrier result. legacy: half-filled symbol "
        "for carriers, central dot for obligate carriers."
    ),
)
@click.option("--no-labels", is_flag=True, help="Hide sample ID labels.")
@click.option("--no-legend", is_flag=True, help="Hide the phenotype legend.")
@click.option(
    "--show-phenotype-labels",
    is_flag=True,
    help="Show YES/NO text labels inside binary phenotype bands.",
)
@click.option(
    "--phenotype-layout",
    type=click.Choice(["auto", "quadrants", "bands"]),
    default=None,
    help="How phenotypes divide a symbol: four quadrants, one per column (up to "
    "four), or horizontal bands. Default 'auto': quadrants up to four columns, "
    "bands beyond.",
)
@click.option(
    "--pgs-scale",
    type=click.Choice(["minmax", "zscore", "percentile"]),
    default=None,
    help="How PGS values are placed on their curve: between the pedigree's lowest and "
    "highest score (minmax, default), as z-scores, or as population percentiles.",
)
@click.option(
    "--no-na-labels",
    is_flag=True,
    help='Do not write "NA" in bands with no data (they stay grey, like unknown).',
)
@click.option(
    "--interactive",
    is_flag=True,
    help="Generate interactive SVG with tooltips and click highlighting.",
)
@click.option(
    "--no-script",
    is_flag=True,
    help="With --interactive: the same SVG without script or tooltip, "
    "for print or for a page that adds its own behaviour.",
)
@click.option(
    "--id-prefix",
    default="",
    metavar="PREFIX",
    help="With --interactive: prepend PREFIX to every id in the SVG, "
    "so several figures can share one page.",
)
@click.option(
    "--responsive",
    is_flag=True,
    help='With --interactive: write width="100%" and no height, '
    "so the figure fills the width of its container.",
)
@click.option(
    "--plink",
    is_flag=True,
    help="Parse as PLINK format (only first 6 columns, ignore genotype data).",
)
@click.option(
    "--circular",
    is_flag=True,
    help="Use circular (radial) layout with founders at the centre.",
)
@click.version_option(version=__version__, prog_name="pypedigree")
def main(
    pedigree_file: str,
    output: str | None,
    fmt: str | None,
    family: str | None,
    dpi: int | None,
    show: bool,
    style_file: str | None,
    separator: str | None,
    pheno_file: str | None,
    annotations_file: str | None,
    outline_color: str | None,
    carrier_notation: str | None,
    no_labels: bool,
    no_legend: bool,
    show_phenotype_labels: bool,
    no_na_labels: bool,
    pgs_scale: str | None,
    phenotype_layout: str | None,
    interactive: bool,
    no_script: bool,
    id_prefix: str,
    responsive: bool,
    plink: bool,
    circular: bool,
) -> None:
    """Generate pedigree diagrams from pedigree files.

    PEDIGREE_FILE is a tab/comma/space-separated pedigree file.
    Supports .ped, .fam (PLINK), and extended formats.
    """
    try:
        input_path = Path(pedigree_file)
        with _reported_warnings():
            # Build style
            style = PedigreeStyle()
            if style_file:
                with open(style_file) as f:
                    overrides = json.load(f)
                style = PedigreeStyle.from_dict(overrides)
            if dpi is not None:
                style.dpi = dpi
            if outline_color is not None:
                style.outline_color = outline_color
            if carrier_notation is not None:
                style.carrier_notation = carrier_notation
            if show_phenotype_labels:
                style.show_phenotype_labels = True
            if no_na_labels:
                style.show_na_labels = False
            if pgs_scale:
                style.pgs_scale = pgs_scale
            if phenotype_layout:
                style.phenotype_layout = phenotype_layout

            # Parse — auto-detect .fam or use --plink flag
            if input_path.suffix.lower() == ".fam":
                pedigree = parse_fam_file(pedigree_file)
            elif plink:
                pedigree = parse_plink_ped(pedigree_file)
            else:
                pedigree = parse_pedigree(pedigree_file)

            # Apply external phenotype file if provided
            if pheno_file:
                pheno_data, pheno_names = parse_phenotype_file(pheno_file)
                apply_phenotypes_to_pedigree(pedigree, pheno_data, pheno_names)

            # Apply external annotation file if provided
            if annotations_file:
                from pypedigree.parser import (
                    apply_annotations_to_pedigree,
                    parse_annotation_file,
                )

                ann_data, ann_names = parse_annotation_file(annotations_file)
                apply_annotations_to_pedigree(pedigree, ann_data, ann_names)

            # Fill in palette colors and PGS spacing from the pedigree's phenotypes
            style = resolve_style(style, pedigree.phenotype_names)

        # Filter to specific family if requested
        if family:
            fam = pedigree.get_family(family)
            if fam is None:
                available = ", ".join(pedigree.family_ids)
                raise click.ClickException(f"Family '{family}' not found. Available: {available}")
            from pypedigree.models import Pedigree

            pedigree = Pedigree(families={family: fam})

        # Interactive SVG mode
        if (no_script or id_prefix or responsive) and not interactive:
            raise click.UsageError("--no-script, --id-prefix and --responsive need --interactive.")
        if interactive:
            from pypedigree.svg_renderer import check_id_prefix, render_interactive_svg, render_svg

            try:
                check_id_prefix(id_prefix)
            except ValueError as e:
                raise click.BadParameter(str(e), param_hint="--id-prefix") from e

            if circular:
                click.echo(
                    "Warning: --circular is not supported with --interactive; "
                    "drawing the hierarchical layout.",
                    err=True,
                )

            if output is None:
                output = str(input_path.with_suffix(".svg"))
            elif not output.endswith(".svg"):
                output = str(Path(output).with_suffix(".svg"))
            render = render_svg if no_script else render_interactive_svg
            with _reported_warnings():
                render(
                    pedigree,
                    style=style,
                    output=output,
                    show_labels=not no_labels,
                    show_legend=not no_legend,
                    id_prefix=id_prefix,
                    responsive=responsive,
                )
            kind = "static" if no_script else "interactive"
            click.echo(f"Saved {kind} SVG: {output}")
            return

        # Determine output path
        if output is None and not show:
            ext = f".{fmt}" if fmt else ".png"
            output = str(input_path.with_suffix(ext))

        if fmt and output and not output.endswith(f".{fmt}"):
            output = str(Path(output).with_suffix(f".{fmt}"))

        # Render. matplotlib is only imported here, so --interactive never
        # imports it (#8).
        import matplotlib

        matplotlib.use("TkAgg" if show else "Agg")
        import matplotlib.pyplot as plt

        from pypedigree.renderer import render_pedigree

        with _reported_warnings():
            fig = render_pedigree(
                pedigree,
                style=style,
                output=None if show else output,
                show_labels=not no_labels,
                show_legend=not no_legend,
                circular=circular,
            )
        if show:
            plt.show()
        else:
            click.echo(f"Saved: {output}")
            plt.close(fig)

    except PedigreeError as e:
        raise click.ClickException(str(e)) from e
