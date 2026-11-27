"""Assemble the web app into a static site (default ``web/_site``).

    uv run python web/build.py              # build the wheel with uv, then the site
    uv run python web/build.py --wheel dist/pypedigree-0.4.6-py3-none-any.whl
    python3 -m http.server -d web/_site 8000

The site holds the page, the wheel of this checkout (installed in the browser
by Pyodide, without its dependencies), the synthetic examples from
``tests/data`` and the HTML documentation. GitHub Pages serves it
(``.github/workflows/pages.yml``). Nothing in it is fetched at run time except
Pyodide itself, from its CDN.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

WEB = Path(__file__).resolve().parent
ROOT = WEB.parent
DATA = ROOT / "tests" / "data"

#: Pyodide release loaded from the jsDelivr CDN.
PYODIDE_VERSION = "314.0.8"

#: Files of the page itself, copied as they are.
PAGE_FILES = ["index.html", "help.html", "app.js", "app.css", "pypedigree_web.py"]

#: Examples offered on the page: all synthetic fixtures of the test suite.
#: ``options`` presets the form; ``style`` is the fixture's style file, if any.
EXAMPLES: list[dict[str, object]] = [
    {
        "file": "basic/three_gen.ped",
        "title": "Three generations",
        "description": "Grandparents, parents and an affected child.",
    },
    {
        "file": "complex/multi_phenotype_pgs.ped",
        "title": "Phenotypes and polygenic score",
        "description": "Two binary and two continuous phenotypes, and a polygenic score.",
        "style": "complex/multi_phenotype_pgs_style.json",
    },
    {
        "file": "phenotype/curve_traits.ped",
        "title": "IQ and polygenic scores on curves",
        "description": "IQ on its population curve (mean 100, SD 15) beside two polygenic "
        "scores: one box per kind of curve.",
        "style": "phenotype/curve_traits.json",
    },
    {
        "file": "phenotype/quadrants.ped",
        "title": "Four phenotypes in quadrants",
        "description": "Autism, ADHD, intellectual disability and epilepsy, one per quadrant, "
        "and an outline column.",
    },
    {
        "file": "phenotype/carrier_genes.ped",
        "title": "Carriers of two genes",
        "description": "Carrier status per gene, in the NSGC 2022 notation.",
    },
    {
        "file": "consanguinity/consanguinity_first_cousins.ped",
        "title": "First-cousin marriage",
        "description": "A consanguineous couple (double line) and their affected child.",
    },
    {
        "file": "twins/twins.ped",
        "title": "Twins",
        "description": "Monozygotic twins and their older sibling.",
    },
    {
        "file": "nomenclature/nomenclature.ped",
        "title": "Standard notation",
        "description": "Pregnancies, adoption, gamete donor, gender identity, consultand.",
    },
    {
        "file": "complex/multigenerational.ped",
        "title": "Large pedigree with loops",
        "description": "81 people over five generations, with two loops.",
    },
    {
        "file": "format/canrisk_v4.txt",
        "title": "CanRisk file",
        "description": "Two families from the CanRisk tool: cancers, gene tests and PRS.",
    },
    {
        "file": "format/phenopacket_consanguineous.json",
        "title": "Phenopackets family",
        "description": "A GA4GH Phenopackets v2 Family: first-cousin parents, the proband's "
        "disease and the HPO terms of every relative.",
        "options": {"hpo_features": True},
    },
]


def build_wheel(dest: Path) -> Path:
    """Build the wheel of this checkout with uv, into ``dest``."""
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run(["uv", "build", "--wheel", "-o", tmp, str(ROOT)], check=True)
        (wheel,) = Path(tmp).glob("pypedigree-*.whl")
        return Path(shutil.copy2(wheel, dest / wheel.name))


def wheel_version(wheel: Path) -> str:
    """The version in a wheel's file name (``pypedigree-0.4.6-py3-none-any.whl``)."""
    return wheel.name.split("-")[1]


def build(out: Path, wheel: Path | None = None) -> Path:
    """Write the site into ``out``, replacing what was there."""
    if out.exists():
        shutil.rmtree(out)
    (out / "wheels").mkdir(parents=True)
    (out / "examples").mkdir()

    for name in PAGE_FILES:
        shutil.copy2(WEB / name, out / name)

    if wheel is None:
        wheel = build_wheel(out / "wheels")
    else:
        wheel = Path(shutil.copy2(wheel, out / "wheels" / Path(wheel).name))

    examples = []
    for example in EXAMPLES:
        source = DATA / str(example["file"])
        target = f"examples/{source.name}"
        shutil.copy2(source, out / target)
        entry = {k: v for k, v in example.items() if k not in ("file", "style")}
        entry["file"] = target
        if "style" in example:
            entry["style"] = json.loads((DATA / str(example["style"])).read_text(encoding="utf-8"))
        examples.append(entry)

    config = {
        "version": wheel_version(wheel),
        "wheel": f"wheels/{wheel.name}",
        "pyodide": PYODIDE_VERSION,
        "examples": examples,
    }
    (out / "config.json").write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")

    # The full documentation, with its images, next to the app.
    docs = out / "docs"
    docs.mkdir()
    shutil.copy2(ROOT / "docs" / "pypedigree_documentation.html", docs / "index.html")
    shutil.copytree(ROOT / "docs" / "imgs", docs / "imgs")
    (out / ".nojekyll").write_text("")
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("-o", "--output", type=Path, default=WEB / "_site")
    parser.add_argument("--wheel", type=Path, help="Use this wheel instead of building one.")
    args = parser.parse_args(argv)
    out = build(args.output, args.wheel)
    config = json.loads((out / "config.json").read_text(encoding="utf-8"))
    print(f"Site written to {out} (pypedigree {config['version']}, Pyodide {config['pyodide']})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
