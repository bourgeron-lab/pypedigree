"""Draw the notation images of `docs/pedigree_standard.md` with pypedigree.

    uv run python docs/imgs/make_notation_images.py

Every image in `docs/imgs/` (not `docs/imgs/examples/`) is written by this
script from invented individuals, so the reference images show exactly what
pypedigree draws. Labels and the legend are hidden; the drawing is cropped to
its content.
"""

from __future__ import annotations

import io
from collections.abc import Callable
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from PIL import Image, ImageChops  # noqa: E402  (a matplotlib dependency)

from pypedigree import PedigreeBuilder, PedigreeStyle, compute_layout, render_pedigree  # noqa: E402
from pypedigree.config import resolve_style  # noqa: E402
from pypedigree.models import EdgeType  # noqa: E402

OUT = Path(__file__).resolve().parent
DPI = 150
PAD = 6  # pixels of white kept around the drawing
UNIT = 0.8  # inches per layout unit, as `renderer._auto_figsize` for large pedigrees
STYLE = PedigreeStyle(dpi=DPI)


def _save(
    builder: PedigreeBuilder,
    name: str,
    annotate: Callable | None = None,
    extra: tuple[float, float, float, float] = (0, 0, 0, 0),
) -> None:
    """Draw at a fixed scale (`UNIT` inches per layout unit), crop to the content.

    `extra` widens the view (left, right, bottom, top, in layout units) for
    text that `annotate` puts outside the pedigree.
    """
    pedigree = builder.build()
    fig, ax = plt.subplots()
    render_pedigree(pedigree, STYLE, ax=ax, show_labels=False, show_legend=False)
    if annotate is not None:
        annotate(ax, pedigree)
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    left, right, bottom, top = extra
    ax.set_xlim(x0 - left, x1 + right)
    ax.set_ylim(y0 - bottom, y1 + top)
    ax.set_position((0, 0, 1, 1))
    fig.set_size_inches((x1 - x0 + left + right) * UNIT, (y1 - y0 + bottom + top) * UNIT)
    buf = io.BytesIO()
    fig.savefig(buf, dpi=DPI, facecolor="white")
    plt.close(fig)
    image = Image.open(buf).convert("RGB")
    box = ImageChops.difference(image, Image.new("RGB", image.size, "white")).getbbox()
    left_px, top_px, right_px, bottom_px = box
    image.crop((left_px - PAD, top_px - PAD, right_px + PAD, bottom_px + PAD)).save(
        OUT / f"{name}.png"
    )
    print(f"wrote {name}.png")


def single(
    name: str, extra: tuple[float, float, float, float] = (0, 0, 0, 0), **kw: object
) -> None:
    _save(PedigreeBuilder().add_individual("A", **kw), name, extra=extra)


def couple(name: str, **union: bool) -> None:
    b = (
        PedigreeBuilder()
        .add_individual("F", sex="male", phenotype="unaffected")
        .add_individual("M", sex="female", phenotype="unaffected")
        .add_individual("C", sex="female", phenotype="unaffected", father="F", mother="M")
        .add_union("F", "M", **union)
    )
    _save(b, name)


def sibship(name: str, sexes: list[str], twin_type: str | None = None) -> None:
    b = (
        PedigreeBuilder()
        .add_individual("F", sex="male", phenotype="unaffected")
        .add_individual("M", sex="female", phenotype="unaffected")
    )
    kids = [f"C{i}" for i in range(len(sexes))]
    for kid, sex in zip(kids, sexes, strict=True):
        b.add_individual(kid, sex=sex, phenotype="unaffected", father="F", mother="M")
    if twin_type:
        b.add_twins(kids, twin_type)
    _save(b, name)


def diagram() -> None:
    """The four kinds of line, labelled (relationship, descent, sibship, individual)."""
    b = (
        PedigreeBuilder()
        .add_individual("F", sex="male", phenotype="unaffected")
        .add_individual("M", sex="female", phenotype="unaffected")
    )
    for kid, sex in (("C1", "male"), ("C2", "female"), ("C3", "male")):
        b.add_individual(kid, sex=sex, phenotype="unaffected", father="F", mother="M")
    b.add_individual("C4", phenotype="unaffected", father="F", mother="M", pregnancy="SAB")

    def annotate(ax, pedigree) -> None:
        family = next(iter(pedigree.families.values()))
        layout = compute_layout(family, resolve_style(STYLE))
        edges = {t: [e for e in layout.edges if e.edge_type == t] for t in EdgeType}
        rel = edges[EdgeType.RELATIONSHIP][0]
        desc = edges[EdgeType.DESCENT][0]
        sib = edges[EdgeType.SIBSHIP][0]
        ind = sorted(edges[EdgeType.INDIVIDUAL], key=lambda e: e.start[0])[-1]
        size = STYLE.shape_size
        arrow = dict(arrowstyle="-|>", lw=0.9, color="black", shrinkA=0, shrinkB=1)
        text = dict(fontsize=8, va="center")
        xm = (rel.start[0] + rel.end[0]) / 2
        ym = (desc.start[1] + desc.end[1]) / 2
        xs = min(sib.start[0], sib.end[0])
        xr = max(sib.start[0], sib.end[0])
        yi = (ind.start[1] + ind.end[1]) / 2
        ax.annotate(
            "1. relationship line",
            xy=(xm - 0.25 * size, rel.start[1]),
            xytext=(xm - 0.25 * size, rel.start[1] + 0.9 * size),
            arrowprops={**arrow, "relpos": (1, 0)},
            ha="right",
            **text,
        )
        ax.annotate(
            "2. line of descent",
            xy=(desc.start[0], ym),
            xytext=(xr + 0.6 * size, ym),
            arrowprops=arrow,
            ha="left",
            **text,
        )
        ax.annotate(
            "3. sibship line",
            xy=(xs + 0.4 * size, sib.start[1]),
            xytext=(xs - 0.2 * size, sib.start[1] + 0.55 * size),
            arrowprops={**arrow, "relpos": (1, 0)},
            ha="right",
            **text,
        )
        ax.annotate(
            "4. individual's line",
            xy=(ind.start[0], yi),
            xytext=(xr + 0.6 * size, yi),
            arrowprops=arrow,
            ha="left",
            **text,
        )

    _save(b, "diagram_informations", annotate, extra=(2.6, 3.2, 0, 0.8))


def main() -> None:
    for sex in ("male", "female", "unknown"):
        tag = "unspecified_sex" if sex == "unknown" else sex
        single(f"individual_unaffected_{tag}", sex=sex, phenotype="unaffected")
        single(f"individual_affected_{tag}", sex=sex, phenotype="affected")
        single(f"deceased_{tag}", sex=sex, phenotype="unaffected", deceased=True)
    # The proband's arrow and "P" sit below and left of the symbol.
    marker = (0.8, 0, 0.5, 0)
    single("male_proband", marker, sex="male", phenotype="affected", proband=True)
    single("female_proband", marker, sex="female", phenotype="affected", proband=True)
    couple("child_with_parents")
    couple("consanguinous_relationship", consanguineous=True)
    couple("broken_relationship", broken=True)
    sibship("monozygotic_twins", ["male", "male"], "MZ")
    sibship("dizygotic_twins", ["male", "female"], "DZ")
    sibship("unknown_zygocity_twins", ["male", "male"], "TW")
    sibship("trizygotic_twins", ["male", "male", "female"], "DZ")
    diagram()


if __name__ == "__main__":
    main()
