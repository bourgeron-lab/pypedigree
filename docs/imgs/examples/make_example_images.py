"""Draw the example images of the documentation with pypedigree.

    uv run python docs/imgs/examples/make_example_images.py

Every image in `docs/imgs/examples/` is written by this script, from the
synthetic fixtures in `tests/data/` and the options the examples gallery
(`docs/examples_gallery.md`) gives for them, so the documentation shows what
the current version draws. Rerun it before a release. `render_pedigree` draws
each figure at its own size (fitted to the pedigree and the legend, #123);
`multigenerational_detail.png` is a crop of the large example to the branch
holding its first-cousin loop.
"""

from __future__ import annotations

import io
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from PIL import Image  # noqa: E402  (a matplotlib dependency)

from pypedigree import (  # noqa: E402
    PedigreeBuilder,
    PedigreeStyle,
    compute_layout,
    parse_pedigree,
    render_pedigree,
)

OUT = Path(__file__).resolve().parent
DATA = OUT.parents[2] / "tests" / "data"
DPI = 150

# name -> (input under tests/data, style JSON under tests/data or None, style overrides)
EXAMPLES: dict[str, tuple[str, str | None, dict[str, Any]]] = {
    "simple_trio": ("basic/simple_trio.ped", None, {}),
    "three_gen": ("basic/three_gen.ped", None, {}),
    "twins_mz": ("twins/twins_mz.ped", None, {}),
    "consanguinity": ("consanguinity/consanguinity.ped", None, {}),
    "deceased": ("extended/extended_deceased.ped", None, {}),
    "carrier": ("phenotype/carrier.ped", None, {}),
    "carrier_genes": ("phenotype/carrier_genes.ped", None, {}),
    "multi_phenotype": (
        "phenotype/multi_phenotype_inline.ped",
        None,
        {"phenotype_layout": "bands"},
    ),
    "multi_phenotype_pgs": (
        "complex/multi_phenotype_pgs.ped",
        "complex/multi_phenotype_pgs_style.json",
        {},
    ),
    "annotations": ("phenotype/test_all.ped", "phenotype/test_all_style.json", {}),
    "gradient_bins": ("phenotype/gradient_bins.ped", "phenotype/gradient_bins.json", {}),
    "legend_minimal": ("phenotype/legend_minimal.ped", "phenotype/legend_minimal.json", {}),
    "ptdt": ("phenotype/ptdt.ped", "phenotype/ptdt.json", {}),
    "curve_traits": ("phenotype/curve_traits.ped", "phenotype/curve_traits.json", {}),
    "canrisk": ("format/canrisk_v4.txt", None, {}),
    "quadrants": ("phenotype/quadrants.ped", None, {"phenotype_layout": "quadrants"}),
    "outline_rings": ("phenotype/outline_rings.ped", None, {}),
    "nomenclature": ("nomenclature/nomenclature.ped", None, {}),
    "multigenerational": ("complex/multigenerational.ped", None, {"phenotype_layout": "bands"}),
}

# The detail of the large example: the branch of 201, whose grandchildren 401
# and 406 marry (the figure caption names them), with the partners.
DETAIL_ROOT = "201"
DETAIL_PAD = 0.9  # layout units around those individuals


def _style(style_file: str | None, overrides: dict[str, Any]) -> PedigreeStyle:
    settings: dict[str, Any] = {}
    if style_file:
        settings = json.loads((DATA / style_file).read_text(encoding="utf-8"))
    settings.update(overrides)
    settings["dpi"] = DPI
    return PedigreeStyle.from_dict(settings)


def _write(fig: plt.Figure, name: str, **savefig: Any) -> None:
    buf = io.BytesIO()
    fig.savefig(buf, dpi=DPI, facecolor="white", **savefig)
    plt.close(fig)
    Image.open(buf).convert("RGB").save(OUT / f"{name}.png", optimize=True)
    print(f"wrote {name}.png")


def example(name: str) -> None:
    source, style_file, overrides = EXAMPLES[name]
    fig = render_pedigree(parse_pedigree(DATA / source), style=_style(style_file, overrides))
    _write(fig, name)


def api_example() -> None:
    """The builder example of the documentation (`PedigreeBuilder`)."""
    builder = PedigreeBuilder("FAM1")
    builder.add_individual("father", sex="male")
    builder.add_individual("mother", sex="female")
    builder.add_individual(
        "child", sex="male", father="father", mother="mother", phenotype="affected"
    )
    fig = builder.render(style=PedigreeStyle(dpi=DPI))
    _write(fig, "api_example")


def multigenerational_detail() -> None:
    """The large example, cropped to the branch of `DETAIL_ROOT`."""
    source, style_file, overrides = EXAMPLES["multigenerational"]
    style = _style(style_file, overrides)
    pedigree = parse_pedigree(DATA / source)
    fig = render_pedigree(pedigree, style=style, show_legend=False)
    ax = fig.axes[0]
    family = next(iter(pedigree.families.values()))
    branch = {DETAIL_ROOT}
    grew = True
    while grew:  # descendants, then everyone partnered with one of them
        grew = False
        for union in family.unions:
            if branch & set(union.partner_ids):
                new = (set(union.partner_ids) | set(union.children_ids)) - branch
                if new:
                    branch |= new
                    grew = True
    # Not the parents of a partner who married in (they are outside the branch).
    nodes = [n for n in compute_layout(family, style).nodes if n.individual_id in branch]
    xs, ys = [n.x for n in nodes], [n.y for n in nodes]
    fig.canvas.draw()
    corners = ax.transData.transform(
        [(min(xs) - DETAIL_PAD, min(ys) - DETAIL_PAD), (max(xs) + DETAIL_PAD, max(ys) + DETAIL_PAD)]
    )
    from matplotlib.transforms import Bbox

    box = Bbox(corners).transformed(fig.dpi_scale_trans.inverted())
    _write(fig, "multigenerational_detail", bbox_inches=box)


def main() -> None:
    for name in EXAMPLES:
        example(name)
    api_example()
    multigenerational_detail()


if __name__ == "__main__":
    main()
