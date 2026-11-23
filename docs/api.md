# API Reference

## Parsing

### `parse_pedigree(file_path, *, canrisk_names=False)`

Parse a pedigree from a file. Auto-detects separator (tab, comma, space) and header rows. When a header is present, columns are matched by name (not position), so phenotype columns can appear directly after Sex without needing the full 10-column format. The file is read as UTF-8 whatever the system's default encoding (cp1252 on most Windows systems); a byte order mark is ignored, and a file that is not UTF-8 raises `PedigreeParseError`. `parse_fam_file`, `parse_plink_ped`, `--pheno` and `--annotations` files are read the same way, and `export_fam` writes UTF-8 with `\n` line endings.

A CanRisk file (first line `##CanRisk <version>`, formats 1 to 4, #109) is read too, by `parse_pedigree` and `parse_pedigree_from_string`: `pypedigree.canrisk.canrisk_to_text(text)` rewrites it as a pedigree file with a header (cancers as binary columns, "not applicable" for ovarian and prostate cancer of the other sex assigned at birth, ages, gene tests and pathology as annotations, the target as consultand, the `##PRS_*` z-scores as `PGS_` columns) and `pypedigree.canrisk.is_canrisk(text)` detects it. See `docs/pedigree_standard.md`, "CanRisk files".

A GA4GH Phenopackets v2 JSON `Family` or `Pedigree` (#108) is read the same way (`pypedigree.phenopackets`): persons, sex and affected status (camelCase or snake_case, enums by name or by the Phenopackets numbers), the proband, deaths, and `consanguinousParents` (a double line, with a warning when the structure does not show it). See `docs/pedigree_standard.md`, "GA4GH Phenopackets v2".

### `export_phenopacket_pedigree(pedigree, file_path)`

Write a Phenopackets v2 `Pedigree` message (`{"persons": [...]}`, camelCase, enums by name, `UNKNOWN_SEX` and `MISSING` written out) without duplicates, as `export_fam` does. A parent named in the input but with no row of their own is written as a person (sex from their role, `MISSING`), so a strict validator finds every `paternalId` and `maternalId`; the unnamed `_phantom_*` parents are not written. Returns the `Path` written.

```python
from pypedigree import export_phenopacket_pedigree, parse_pedigree

export_phenopacket_pedigree(parse_pedigree("family.ped"), "pedigree.json")
```

```python
from pypedigree import parse_pedigree

ped = parse_pedigree("family.ped")
```

**Args**: `file_path` — Path to the pedigree file (str or Path). `canrisk_names` — for a CanRisk file, label individuals with its `Name` column instead of `IndivID`, where given and unique in the family (#120); ignored for other files.

**Returns**: `Pedigree` object containing all families.

**Raises**: `PedigreeParseError`, `PedigreeValidationError`, `FileNotFoundError`.

---

### `parse_pedigree_from_string(text, *, canrisk_names=False)`

Parse a pedigree from a string (same format as file contents).

```python
from pypedigree import parse_pedigree_from_string

text = "FAM1\tfather\t0\t0\t1\t1\nFAM1\tmother\t0\t0\t2\t1\nFAM1\tchild\tfather\tmother\t2\t2\n"
ped = parse_pedigree_from_string(text)
```

**Args**: `text` — Pedigree data as a string. `canrisk_names` — as for `parse_pedigree`.

**Returns**: `Pedigree` object.

---

### `parse_pedigree_from_dataframe(df)`

Parse a pedigree from a pandas DataFrame. Requires `pip install pypedigree[pandas]`.

The DataFrame is read as a pedigree file would be:
- **By name**, when its first column is named like a family ID column (`FID`, `family_id`,
  `family`) or another known header name (`IID`, `sex`...): the column names are the file's
  header. Columns are matched by name in any order, and extra phenotype (`b_`, `g_`, `c_`,
  `o_`, `cr_`, `PGS`, `pTDT`, no prefix), annotation (`a_`), notation, `css_class` and `data_`
  columns are read, as from a file. This is what `pd.read_csv("family.ped", sep="\t")` of a
  file with a header gives.
- **By position** otherwise, e.g. `pd.read_csv(..., header=None)` (columns 0, 1, 2...):
  1. family_id, 2. sample_id, 3. paternal_id, 4. maternal_id, 5. sex, 6. phenotype,
  7. date_of_birth (opt), 8. twins (opt), 9. deceased (opt), 10. proband (opt). Further
  columns are ignored with a `PedigreeDataWarning`.

Missing values (NaN, `None`, `pd.NA`) are missing: an unknown parent, sex or phenotype, no
data for an extra phenotype (never a value of 0). `pd.read_csv` reads `N/A` and `n/a` as NaN by
default, so "not applicable" (#107) becomes no data: to keep it, read the file with
`keep_default_na=False, na_values=["NA", ""]` (#121). A float with no fractional part (`2.0`, as
pandas stores an integer column holding NaN) is read as the integer, `True`/`False` as 1/0.
A tab or a line break in a cell raises `PedigreeParseError`.

```python
import pandas as pd
from pypedigree import parse_pedigree_from_dataframe

df = pd.DataFrame({
    "fid": ["FAM1", "FAM1", "FAM1"],
    "iid": ["father", "mother", "child"],
    "pat": ["0", "0", "father"],
    "mat": ["0", "0", "mother"],
    "sex": [1, 2, 2],
    "pheno": [1, 1, 2],
    "b_autism": [1, 1, 2],     # read by name: the first column is "fid"
    "PGS_autism": [0.2, -0.4, 1.8],
})
ped = parse_pedigree_from_dataframe(df)
```

**Args**: `df` — pandas DataFrame with at least 6 columns.

**Returns**: `Pedigree` object.

---

### `parse_fam_file(file_path)`

Parse a PLINK `.fam` file. The `.fam` format is whitespace-delimited with exactly 6 columns (FID, IID, PAT, MAT, SEX, PHENO) and no header row.

```python
from pypedigree import parse_fam_file

ped = parse_fam_file("family.fam")
```

**Args**: `file_path` — Path to the `.fam` file (str or Path).

**Returns**: `Pedigree` object.

---

### `parse_plink_ped(file_path)`

Parse a PLINK `.ped` file, ignoring genotype columns after column 6. Only the standard 6 pedigree columns are read; no header row is expected.

```python
from pypedigree import parse_plink_ped

ped = parse_plink_ped("family.ped")
```

**Args**: `file_path` — Path to the PLINK `.ped` file (str or Path).

**Returns**: `Pedigree` object.

---

### `export_fam(pedigree, file_path)`

Export a `Pedigree` to PLINK `.fam` format. Writes tab-delimited, 6 columns, no header. Phantom parents and duplicate nodes are excluded; internal references (`_dup_`, `_phantom_`) are resolved back to original IDs.

```python
from pypedigree import parse_pedigree, export_fam

ped = parse_pedigree("family.ped")
export_fam(ped, "output.fam")
```

**Args**:
- `pedigree` — `Pedigree` to export.
- `file_path` — Output file path (str or Path).

**Returns**: `Path` of the written file.

---

### `parse_phenotype_file(file_path)`

Parse a separate phenotype file (.pheno TSV). First row is a header (`sample_id`, `pheno1`, `pheno2`, ...), subsequent rows are data.

```python
from pypedigree import parse_phenotype_file

pheno_data, pheno_names = parse_phenotype_file("data.pheno")
# pheno_data: {"1": {"epilepsy": 1.0, "IQ": 105.0}, ...}
# pheno_names: ["epilepsy", "IQ", "SRS"]
```

**Args**: `file_path` — Path to the .pheno file (str or Path).

**Returns**: Tuple of `(pheno_data, pheno_names)` where `pheno_data` maps `sample_id -> {pheno_name: raw_value}`.

---

### `apply_phenotypes_to_pedigree(pedigree, pheno_data, pheno_names)`

Apply external phenotype data to individuals in a pedigree and normalize values.

```python
from pypedigree import parse_pedigree, parse_phenotype_file, apply_phenotypes_to_pedigree

ped = parse_pedigree("family.ped")
pheno_data, pheno_names = parse_phenotype_file("data.pheno")
apply_phenotypes_to_pedigree(ped, pheno_data, pheno_names)
```

**Args**:
- `pedigree` — `Pedigree` to update (modified in place).
- `pheno_data` — Map of `sample_id -> {pheno_name: raw_value}`.
- `pheno_names` — Ordered list of phenotype names.

**Warns**: `PedigreeDataWarning` when sample IDs in `pheno_data` match no
individual in the pedigree, giving their number and the first five. Their
values are ignored, and an individual with no phenotype data is drawn with
"NA" bands, so a file whose IDs match nobody gives an all-"NA" figure.

---

### `parse_annotation_file(file_path)`

Parse an annotation file (TSV). First row is a header (`sample_id`, followed by annotation column names), subsequent rows are data. Annotation columns provide text displayed below the individual's ID label (not phenotypes).

```python
from pypedigree import parse_annotation_file

anno_data, anno_names = parse_annotation_file("annotations.tsv")
# anno_data: {"1": {"Role": "Proband", "Clinic": "Paris"}, ...}
# anno_names: ["Role", "Clinic"]
```

**Args**: `file_path` -- Path to the annotation file (str or Path).

**Returns**: Tuple of `(anno_data, anno_names)` where `anno_data` maps `sample_id -> {anno_name: value}` and `anno_names` is the ordered list of annotation column names.

---

### `apply_annotations_to_pedigree(pedigree, anno_data, anno_names)`

Apply external annotation data to individuals in a pedigree.

```python
from pypedigree import parse_pedigree, parse_annotation_file, apply_annotations_to_pedigree

ped = parse_pedigree("family.ped")
anno_data, anno_names = parse_annotation_file("annotations.tsv")
apply_annotations_to_pedigree(ped, anno_data, anno_names)
```

**Args**:
- `pedigree` -- `Pedigree` to update (modified in place).
- `anno_data` -- Map of `sample_id -> {anno_name: value}`.
- `anno_names` -- Ordered list of annotation column names.

### `build_notation(*, pregnancy=None, adopted=None, childless=None, art=None, gender=None, consultand=False, evaluated=False)`

Build a `Notation` from enum values or the codes of the `.ped` notation columns
(case-insensitive): `pregnancy` `P`/`SAB`/`TOP`/`ECT`/`SB`, `adopted` `in`/`out`,
`childless` `choice`/`infertility`, `art` `D`/`S`, `gender` `man`/`woman`/`nonbinary`,
`evaluated` `True`/`1`/`yes` (documented evaluation) or a result text written below the symbol.
Raises `PedigreeParseError` on an unknown code. `PedigreeBuilder.add_individual()`
uses it.

```python
from pypedigree import build_notation

build_notation(pregnancy="SAB", adopted="in")
# Notation(pregnancy=<Pregnancy.MISCARRIAGE: 'SAB'>, adoption=<Adoption.IN: 'in'>, ...)
```

---

## Layout

### `compute_layout(family, style=None, *, labels=True)`

Compute (x, y) coordinates for all individuals and connection geometry. No matplotlib dependency.

```python
from pypedigree import compute_layout, PedigreeStyle

layout = compute_layout(family, style=PedigreeStyle())
```

**Args**: `family` — `Family` object. `style` — `PedigreeStyle` (optional, defaults used if None). `labels` — whether labels will be drawn: neighbours are then kept far enough apart for their labels (the renderers pass `show_labels`).

**Returns**: `FamilyLayout` with `.nodes`, `.edges`, `.twin_connectors`, `.duplicate_arcs`.

#### Layout pipeline

0. **Couple alone in its family** — A family holding only two parents, a man and a woman, both founders, with no child in the file (a trio whose child is missing, say) is drawn as their couple, father on the left, with a descent line to a "?" instead of two unrelated symbols. The pedigree itself is not changed. `with_unknown_offspring(family)` (in `pypedigree.layout`) returns a copy with a placeholder child (`is_unknown_offspring=True`) and the couple's union, or the family itself.
1. **Assign generations** — BFS from founders using a parent→child index, then partners are aligned and children propagated until every child sits below both parents. A duplicate node stands on the generation of the individual it repeats.
2. **Compute horizontal positions** — Each union's descendants form a block. Sibling blocks are packed left to right against each other's contours, generation by generation, and every couple is centred over its children. Nodes on a row stay at least `sibling_gap` apart (`couple_gap` between partners), and the children of one union occupy a stretch of their row that no other sibship enters. Where two neighbours' labels would overlap at that distance (long IDs, a large `font_size`, long annotations), those two are moved apart by the room their labels need: each label is measured with the character widths of DejaVu Sans at the SVG's scale, the widest of its lines counting (ID, then notes and annotations at 70 % of the size), with half an em between labels.
3. **Break tangled loops** — A loop partner (both partners descend from people in the pedigree) is drawn next to their partner, with their parents' sibship line stretched to reach them. Step 2 orders siblings so the two branches of a loop meet: each ancestor on one side is placed last among their siblings, each on the other side first. When both branches descend from the same sibship (cousins marrying), those two siblings are kept together, before any sibling who marries into a family placed later. Where the stretched line would still overlap another sibship, the partner is drawn twice instead: among their siblings, and as a duplicate next to their partner, joined by a dotted arc. One loop is broken per group of overlapping lines, then steps 1-2 are repeated, so no loop is broken that another break made unnecessary. The input family is never modified; the duplicates appear as extra nodes with `is_duplicate=True`.
4. **Build geometry** — Nodes, edges, twin connectors, duplicate arcs. Every edge is horizontal or vertical: a single child that is not centred under its parents (possible only in a loop) hangs from a short sibship-level segment rather than a diagonal. Sibship lines of different unions never overlap. Duplicate arcs bow by at most 0.4 generation gaps.

---

### `compute_layout_circular(family, style=None)`

Compute a circular/radial layout. Founders sit at the centre and each later generation on its own concentric ring, in the left-to-right order of the hierarchical layout wrapped around the circle (counter-clockwise, centred at the bottom). It therefore inherits that layout's guarantees: siblings stay together, sibship arcs of different couples never overlap, partners sit side by side, and tangled loops are drawn with duplicated individuals. A small family occupies only the arc it needs; a large one the whole circle.

```python
from pypedigree import compute_layout_circular, PedigreeStyle

layout = compute_layout_circular(family, style=PedigreeStyle())
```

**Args**: `family` — `Family` object. `style` — `PedigreeStyle` (optional; `base_radius` and `radial_gap` control ring sizing).

`base_radius` is expanded automatically when a ring would be too crowded for
the shapes to fit, so nodes never overlap. Shapes stay upright (a square turned
to face outward would read as the unknown-sex diamond), so phenotype bands are
horizontal as in the hierarchical layout. ID labels are set outside each shape,
away from the lines arriving from the centre, on a white backing: radially
outward (`NodeLayout.label_angle`), or turned aside to the first direction where
no line, shape or label is in the way. `render_pedigree` places them when the
figure is drawn, at its final size.
Twins hang from one point of the sibship arc, as a V (with a bar for MZ
twins), and duplicates are joined to their original by a dotted arc.

**Returns**: `FamilyLayout` with `.nodes`, `.edges` (edges may contain `waypoints`), `.twin_connectors`, `.duplicate_arcs`.

---

## PedigreeBuilder

### `PedigreeBuilder(family_id="FAM1")`

Fluent builder for constructing pedigrees programmatically. All mutation methods return `self` for chaining.

```python
from pypedigree import PedigreeBuilder

builder = PedigreeBuilder("FAM1")
builder.add_individual("father", sex="male")
builder.add_individual("mother", sex="female")
builder.add_individual("child", sex="male", father="father", mother="mother",
                       phenotype="affected", phenotypes={"IQ": 85})
fig = builder.render(output="family.png")
```

### `add_individual(sample_id, ...)`

Add an individual to the pedigree.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `sample_id` | `str` | required | Unique identifier |
| `sex` | `str \| Sex` | `"unknown"` | `"male"`, `"female"`, `"unknown"` (or `"m"`, `"f"`, `"1"`, `"2"`) |
| `phenotype` | `str \| Phenotype` | `"unknown"` | `"affected"`, `"unaffected"`, `"unknown"` |
| `father` | `str \| None` | `None` | Father's sample_id |
| `mother` | `str \| None` | `None` | Mother's sample_id |
| `family_id` | `str \| None` | `None` | Override default family ID |
| `deceased` | `bool` | `False` | Deceased marker |
| `proband` | `bool` | `False` | Proband marker |
| `inferred` | `bool` | `False` | Never sampled or sequenced, only links real members ("dummy"): dashed outline like a phantom, but label, parents, phenotypes and `.fam` export kept |
| `carrier` | `bool` | `False` | Carrier: patterned section (NSGC 2022 notation) or half-filled shape (legacy) |
| `obligate_carrier` | `bool` | `False` | Obligate carrier: sparser pattern (NSGC 2022) or center dot (legacy) |
| `date_of_birth` | `str \| None` | `None` | Date of birth |
| `phenotypes` | `dict \| None` | `None` | Extra phenotypes as `{name: value}`; `"N/A"` for a phenotype that does not apply (#107) |
| `annotations` | `list[str] \| None` | `None` | Annotation lines displayed below the ID label (italic grey text) |
| `pregnancy` | `str \| Pregnancy \| None` | `None` | `"P"` ongoing pregnancy, `"SAB"` miscarriage, `"TOP"` termination, `"ECT"` ectopic, `"SB"` stillbirth |
| `adopted` | `str \| Adoption \| None` | `None` | `"in"` (`father`/`mother` are the adoptive parents; dashed line) or `"out"` (they are the biological parents); brackets around the symbol |
| `childless` | `str \| Childlessness \| None` | `None` | `"choice"` (one bar) or `"infertility"` (two bars) below the symbol; `"choice:<ID>"` for a couple without children, the marker hanging from their line (#104) |
| `art` | `str \| ReproductiveRole \| None` | `None` | `"D"` gamete donor or `"S"` surrogate, written inside the symbol |
| `gender` | `str \| Gender \| None` | `None` | Gender identity `"man"`, `"woman"`, `"nonbinary"`: the symbol follows it and `sex` is the sex assigned at birth (AMAB/AFAB below when they differ) |
| `consultand` | `bool` | `False` | Arrow without "P" |
| `count` | `int \| str \| None` | `None` | Several individuals in one symbol (2 or more, or `"n"`), written inside; cannot be a parent (#104). Unlike a file's `count` column, any other value is an error |
| `asymptomatic` | `bool` | `False` | Asymptomatic or presymptomatic carrier: a vertical line through the symbol (#104) |
| `evaluated` | `bool \| str` | `False` | Documented evaluation: `True` draws an asterisk at the lower right; text (e.g. `"WGS: negative"`) is the result, written below |
| `twin_of` | `str \| list[str] \| None` | `None` | Sample ID of this individual's twin (a list for triplets and more); the twin may be added before or after |
| `css_class` | `str \| list[str] \| None` | `None` | Extra classes on the individual's node in the interactive SVG (`"flagged"`, several separated by spaces, or a list) |
| `data` | `Mapping[str, object] \| None` | `None` | Extra `data-*` attributes on that node: `{"qc": "high"}` writes `data-qc="high"`; `_` in a name becomes `-`, `None` values are left out |
| `twin_type` | `str \| TwinType` | `"TW"` | `"MZ"` (monozygotic), `"DZ"` (dizygotic) or `"TW"` (zygosity unknown); also `"monozygotic"`, `"dizygotic"`, `"unknown"` |

An unknown notation code or twin type raises `PedigreeError`.

### `add_twins(sample_ids, twin_type="TW")`

Declare a group of twins (or triplets, ...) of one zygosity: every member is
recorded as the twin of every other, as the `twins` column does
(`MZ:S8` on S7's row, `MZ:S7` on S8's). The individuals may be added before
or after. Fewer than two distinct IDs, or an unknown type, raise
`PedigreeError`. At `build()`, a twin that was never added, twins in
different families, or twins that do not share their parents give a
`PedigreeDataWarning`.

```python
builder.add_twins(["S7", "S8"], "MZ")
# or, on one of them:
builder.add_individual("S8", father="F", mother="M", twin_of="S7", twin_type="MZ")
```

### `add_union(partner1, partner2, ...)`

Explicitly add a union. Only needed for childless couples or consanguineous/broken flags.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `partner1` | `str` | required | First partner's sample_id |
| `partner2` | `str` | required | Second partner's sample_id |
| `consanguineous` | `bool` | `False` | Double-line (shared ancestors) |
| `broken` | `bool` | `False` | Broken-line (dissolved) |
| `childless` | `str \| Childlessness \| None` | `None` | A couple without children: `"choice"` (one bar) or `"infertility"` (two), hanging from their relationship line (#104) |

### `set_style(style)` / `build()` / `render(output=None, ...)` / `render_interactive_svg(output=None, ...)` / `render_svg(output=None, ...)`

- `set_style(style)` — Set custom `PedigreeStyle`.
- `build()` — Build a `Pedigree` object. Auto-infers unions and phantom parents, as the parser does: a parent named but not added becomes a phantom with that ID, and a child with only one parent gets an unlabelled phantom for the other.
- `render(output, *, style, show_labels, show_legend, circular=False)` — Build and render as matplotlib Figure. Auto-assigns phenotype colors. Pass `circular=True` for radial layout.
- `render_interactive_svg(output, *, style, show_labels, show_legend, id_prefix, responsive)` — Build and render as interactive SVG string with tooltips and click highlighting.
- `render_svg(output, *, style, show_labels, show_legend, id_prefix, responsive)` — Build and render as a static SVG string (no script, no tooltip), for print or host pages.

### Jupyter integration

The builder supports inline display in Jupyter notebooks:

```python
builder  # displays the pedigree inline
```

---

## Rendering

### `render_pedigree(pedigree, style=None, output=None, ax=None, show_labels=True, show_legend=True, circular=False)`

Render all families in a pedigree to a matplotlib figure.

```python
from pypedigree import parse_pedigree, render_pedigree

ped = parse_pedigree("family.ped")
fig = render_pedigree(ped, output="family.png")
```

**Args**:
- `pedigree` — `Pedigree` object.
- `style` — `PedigreeStyle` (optional).
- `output` — File path to save (optional). Format inferred from extension.
- `ax` — Existing matplotlib `Axes` to render onto (optional).
- `show_labels` — Whether to draw sample ID labels (default: True).
- `show_legend` — Whether to draw phenotype legend (default: True).
- `circular` — Use circular/radial layout instead of the default top-down layout (default: False).

**Returns**: matplotlib `Figure`.

**Note**: For large pedigrees, figure dimensions are auto-scaled based on layout bounds.

`render_pedigree` and `render_family` are imported on first use: `import pypedigree` does not import matplotlib, and parsing, layout and SVG output never do.

---

### `render_family(family, style=None, ax=None, show_labels=True)`

Render a single family, returning both the figure and the computed layout.

```python
from pypedigree import parse_pedigree, render_family

ped = parse_pedigree("family.ped")
fam = ped.get_family("FAM1")
fig, layout = render_family(fam)
```

`render_family` draws the family's structure only, without phenotype bands, PGS curves or legend, which are resolved for the whole pedigree. To draw one family with them, use `render_pedigree(ped.subset(["FAM1"]))`.

**Returns**: Tuple of (`Figure`, `FamilyLayout`).

---

### `render_interactive_svg(pedigree, style=None, output=None, show_labels=True, show_legend=True, *, id_prefix="", responsive=False)`

Render all families as an interactive SVG with embedded JavaScript for tooltips, click highlighting, zoom, and pan. The SVG is written directly, without matplotlib (colors are read by `pypedigree.colors`), and draws the same figure as `render_pedigree()`: both renderers take their decisions from `pypedigree.phenotypes` and `pypedigree.carriers`.

```python
from pypedigree import parse_pedigree, render_interactive_svg

ped = parse_pedigree("family.ped")
svg_str = render_interactive_svg(ped, output="family.svg")
```

**Args**:
- `pedigree` — `Pedigree` object.
- `style` — `PedigreeStyle` (optional).
- `output` — File path to save (optional).
- `show_labels` — Whether to draw sample ID labels (default: True).
- `show_legend` — Whether to draw the phenotype and carrier legend (default: True).
- `id_prefix` — Prepended to every `id` in the SVG and to every `url(#...)` / `href="#..."` reference to one, so several figures can share a page (default: `""`, ids `main-group`, `tooltip`, ...). Letters, digits, `_`, `-`, `.`, starting with a letter or `_`; anything else raises `ValueError`.
- `responsive` — Write `width="100%"` and no height instead of a size in pixels, so the figure fills the width of its container with the proportions of its `viewBox` (default: False).

**Returns**: SVG content as a string.

**Features**:
- Hover tooltips (sample_id, sex, phenotype, DOB, extra phenotypes, carrier status)
- Host hooks: each node is a `<g class="node …">` with `data-id`, `data-parents`, `data-children`, `data-siblings`, plus the individual's own classes and `data-*` attributes (`css_class` / `data_` columns, `css_class=` / `data=` builder keywords), so a page can style them (`.node.flagged`) or select them (`[data-qc="high"]`) without re-rendering. The click highlighting adds and removes only the `dimmed` class, so a host's classes stay
- Click node to highlight parents, children, siblings (click again to reset)
- Mouse wheel zoom when the SVG is opened on its own; inlined in a page, Ctrl/Cmd + wheel (or a trackpad pinch) zooms and the plain wheel scrolls the page. Ctrl+drag or the middle button pans
- Carrier notation: patterned sections defined once as SVG `<pattern>`s (NSGC 2022), or half-fill and center dot with `carrier_notation="legacy"`
- Phenotype bands (binary, gradient, diverging, categorical), `o_` outline colors, PGS distribution curves and `a_` annotations, from the same `style.phenotype_colors`, `gradient_range`, `gradient_diverging` and `pgs_colors` as the PNG
- The legend below the pedigree, with the same rows and proportions as the PNG legend
- Several families side by side, each titled

**Node attributes, for host pages.** Each individual is a `<g class="node">`
with `data-id` (sample ID), `data-parents`, `data-children`, `data-siblings`
(comma-separated IDs), and `data-x`, `data-y` (the centre of the symbol) and
`data-size` (its size). Positions are in the coordinates of `#main-group`, the
drawing's user units as set by the `viewBox`; zoom and pan transform that
group, so an element a host appends inside it, at those coordinates, stays on
the node:

```js
const node = svg.querySelector('.node[data-id="S7"]');
const x = +node.dataset.x, y = +node.dataset.y, r = +node.dataset.size / 2;
const ring = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
ring.setAttribute('cx', x); ring.setAttribute('cy', y); ring.setAttribute('r', r * 1.4);
svg.getElementById('main-group').appendChild(ring);
```

The circular layout is not available in the interactive SVG. Text is placed
with explicit offsets rather than `dominant-baseline`, so the file renders the
same in browsers, Inkscape and librsvg.

**Inlining in a page.** The SVG can be put in an HTML page, next to other
SVGs, as `innerHTML` or in a template. Its root has the class `pypedigree`,
its CSS only applies under that class, and its script only acts on its own
root (`document.currentScript.closest('svg')`) with no listener on the
window. Give each figure its own `id_prefix` so their ids do not collide:

```python
from pypedigree import render_interactive_svg

parts = [
    render_interactive_svg(ped_a, id_prefix="fam-a-", responsive=True),
    render_interactive_svg(ped_b, id_prefix="fam-b-", responsive=True),
]
```

A host then finds a figure's group with `svg.querySelector('[id="fam-a-main-group"]')`.
Scripts inserted through `innerHTML` do not run: such a page gets the drawing,
the `<title>` tooltips and the `data-*` attributes, and adds its own behaviour,
or uses `render_svg()`.

---

### `render_svg(pedigree, style=None, output=None, show_labels=True, show_legend=True, *, id_prefix="", responsive=False)`

The same drawing as `render_interactive_svg()`, without script and without the
tooltip group: for a printed report, or for a page that adds its own
behaviour. Each individual keeps its `<title>` (shown by browsers as a native
tooltip), its classes and its `data-*` attributes, and the `dimmed` class still
fades a node or an edge. Same arguments as `render_interactive_svg()`. CLI:
`--interactive --no-script`.

```python
from pypedigree import parse_pedigree, render_svg

ped = parse_pedigree("family.ped")
svg_str = render_svg(ped, id_prefix="fam1-", responsive=True)
```

### `render_families(pedigree, output_dir, *, fmt="png", family_ids=None, style=None, show_labels=True, show_legend=True, circular=False, interactive=False, script=True, id_prefix="", responsive=False)`

Draw each family into its own file in `output_dir` (created if needed; existing files are overwritten), and return a `dict` of family ID → `Path`, in drawing order. CLI: `--output-dir`.

Every family is drawn with the colors, column modes and legend of the whole pedigree (see `Pedigree.subset`), so the figures of one cohort can be compared. The style is checked once, before anything is drawn. A family that fails does not stop the others: once all are drawn, `PedigreeBatchError` is raised, with `.written` (family ID → path) and `.failures` (family ID → message). Figures are closed as they are saved, so hundreds of families do not fill the memory; the caller's own figures are left open.

| Argument | Description |
|----------|-------------|
| `fmt` | `"png"`, `"pdf"` or `"svg"` (matplotlib's SVG); ignored with `interactive` |
| `family_ids` | Families to draw, in order (default: all). `KeyError` for an unknown ID |
| `circular` | Circular layout (PNG, PDF, matplotlib SVG); with `interactive` it warns and draws the hierarchical layout |
| `interactive` | pypedigree's SVG (`render_interactive_svg`); matplotlib is not imported |
| `script` | With `interactive`: `False` writes `render_svg`'s static SVG |
| `id_prefix`, `responsive` | With `interactive`: as in `render_interactive_svg` |

```python
from pypedigree import parse_pedigree, render_families

ped = parse_pedigree("cohort.ped")
files = render_families(ped, "figures", fmt="svg")   # {'FAM1': Path('figures/FAM1.svg'), ...}
```

### `family_file_names(family_ids)`

The file name (without extension) `render_families` gives each family. Letters (accented ones too), digits, `-`, `_` and `.` are kept, every other character becomes `_`; leading and trailing dots and spaces are dropped (an empty name becomes `family`), Windows device names (`CON`, `NUL`, `COM1`...) get a trailing `_`, and names are cut to 100 characters. Two IDs giving the same name, ignoring case (`F/1` and `F_1`, or `f1` and `F1`), are told apart by `_2`, `_3`... in order, with a `PedigreeDataWarning`.

---

## Data Model

### Enums

| Enum | Values | Description |
|------|--------|-------------|
| `Sex` | `MALE (1)`, `FEMALE (2)`, `UNKNOWN (0)` | Biological sex |
| `Phenotype` | `UNAFFECTED (1)`, `AFFECTED (2)`, `UNKNOWN (0)` | Phenotypic status |
| `TwinType` | `MONOZYGOTIC`, `DIZYGOTIC`, `UNKNOWN` | Twin zygosity |
| `EdgeType` | `RELATIONSHIP`, `DESCENT`, `SIBSHIP`, `INDIVIDUAL`, `CHILDLESS` | Connection line type (`CHILDLESS`: the stem and bars of a childlessness marker) |
| `Pregnancy` | `ONGOING ("P")`, `MISCARRIAGE ("SAB")`, `TERMINATION ("TOP")`, `ECTOPIC ("ECT")`, `STILLBIRTH ("SB")` | A pregnancy or how it ended |
| `Adoption` | `IN`, `OUT` | Adopted in (parents given are adoptive) or out (they are biological) |
| `Childlessness` | `BY_CHOICE`, `INFERTILITY` | No children by choice or reason unknown, or infertility |
| `ReproductiveRole` | `DONOR ("D")`, `SURROGATE ("S")` | Role in assisted reproduction |
| `Gender` | `MAN`, `WOMAN`, `NONBINARY` | Gender identity; decides the symbol |
| `LineStyle` | `SOLID`, `DOUBLE`, `BROKEN`, `DASHED` | Line rendering style |

### `Individual`

Represents one person in the pedigree.

| Field | Type | Description |
|-------|------|-------------|
| `sample_id` | `str` | Unique identifier within family |
| `family_id` | `str` | Family this individual belongs to |
| `sex` | `Sex` | Biological sex |
| `phenotype` | `Phenotype` | Affected/unaffected status |
| `father_id` | `str \| None` | Father's sample_id |
| `mother_id` | `str \| None` | Mother's sample_id |
| `date_of_birth` | `str \| None` | DOB string (YYYYMMDD/YYYYMM/YYYY) |
| `twin_info` | `list[TwinInfo]` | Twin relationship declarations |
| `is_deceased` | `bool` | Whether the individual is deceased |
| `is_proband` | `bool` | Whether this is the proband |
| `is_phantom` | `bool` | Whether this is a generated phantom parent |
| `is_inferred` | `bool` | A row never sampled (`not_sampled` column, `inferred=` keyword): dashed outline, label and phenotypes kept |
| `is_unknown_offspring` | `bool` | The placeholder child the layout adds under a couple alone in its family (drawn as "?") |
| `css_classes` | `list[str]` | Extra classes on the individual's SVG node (`css_class` column or keyword) |
| `data_attributes` | `dict[str, str]` | Extra `data-*` attributes on the SVG node, name without `data-` (`data_` columns, `data` keyword) |
| `is_duplicate` | `bool` | Whether this is a duplicate node |
| `is_carrier` | `bool` | Whether this is a carrier (`carrier` column = 1) |
| `is_obligate_carrier` | `bool` | Whether this is an obligate carrier (`carrier` column = 2) |
| `duplicate_of` | `str \| None` | Original sample_id if duplicate |
| `extra_phenotypes` | `dict[str, float \| None]` | Normalized phenotype values [0, 1] |
| `extra_phenotypes_raw` | `dict[str, float \| None]` | Original raw phenotype values |
| `not_applicable` | `set[str]` | Phenotype and carrier columns that do not apply to the individual (`N/A` cells, #107, #121); their value is stored as no data and drawn white, as unaffected (#120) |
| `annotations` | `dict[str, str]` | Annotation key-value pairs (displayed as italic grey text below ID label) |
| `notation` | `Notation` | Pregnancy, adoption, childlessness, assisted reproduction, gender identity and consultand (see below) |

### `Notation`

Frozen dataclass of the standard notation beyond sex, affection status and death
(NSGC nomenclature, Bennett et al. 2008 and 2022). Every field defaults to "nothing
to note". How each is drawn is decided in `pypedigree.symbols`.

| Field | Type | Drawn as |
|-------|------|----------|
| `pregnancy` | `Pregnancy \| None` | "P" inside the symbol; a triangle for SAB, ECT ("ECT" below) and TOP (struck through); the symbol struck through with "SB" below for a stillbirth |
| `adoption` | `Adoption \| None` | Brackets; the line to adoptive parents is dashed |
| `childlessness` | `Childlessness \| None` | A line below the symbol ending in one bar, or two |
| `reproductive_role` | `ReproductiveRole \| None` | "D" or "S" inside the symbol |
| `gender` | `Gender \| None` | Symbol shape; AMAB/AFAB below when it differs from `Individual.sex` |
| `is_consultand` | `bool` | Arrow without "P" |
| `is_evaluated` | `bool` | Documented evaluation, records reviewed: an asterisk at the lower right of the symbol |
| `evaluation` | `str \| None` | The evaluation's result, written below the symbol (no "E" prefix, as in the 2022 revision) |
| `count` | `str \| None` | Several individuals in one symbol: their number, or `"n"` (#104) |
| `is_asymptomatic` | `bool` | Asymptomatic or presymptomatic carrier: a vertical line through the symbol (#104) |
| `childless_partner` | `str \| None` | With `childlessness`, the partner of a couple without children (#104) |

### `Union`

Represents a mating pair and their children.

| Field | Type | Description |
|-------|------|-------------|
| `partner_ids` | `tuple[str, str]` | (male_id, female_id) |
| `children_ids` | `list[str]` | Children of this union |
| `is_consanguinous` | `bool` | Auto-detected shared ancestry |
| `is_broken` | `bool` | Separated/divorced (API-only) |
| `childlessness` | `Childlessness \| None` | A couple declared without children: the marker hangs from their relationship line, its bars below the partners' symbols (#104, #124) |

### `Family`

Collection of individuals and unions sharing a family_id.

| Method | Returns | Description |
|--------|---------|-------------|
| `get_founders()` | `list[Individual]` | Individuals with no parents |
| `get_children(sample_id)` | `list[Individual]` | Children of a given parent |
| `get_union_for_parents(id1, id2)` | `Union \| None` | Find union by partner IDs |

### `Pedigree`

Top-level container for multiple families.

| Property / Method | Returns | Description |
|-------------------|---------|-------------|
| `family_ids` | `list[str]` | List of family identifiers |
| `get_family(family_id)` | `Family \| None` | Look up a family by ID |
| `subset(family_ids)` | `Pedigree` | Only these families, with every phenotype and annotation column kept (to draw one family with `render_pedigree`; `Pedigree(families=...)` would drop the columns). `KeyError` for an unknown ID |
| `source` | `Pedigree \| None` | For a subset, the whole pedigree it was taken from (not compared, not in `repr`) |
| `column_scope` | `Pedigree` | `source`, or the pedigree itself: the renderers decide each column's mode, category colors and carrier legend over its individuals, as the parser normalised the values over the whole file. So a subset is drawn as in the figure of the whole file |
| `individuals()` | `Iterator[Individual]` | Every individual of every family, in order |
| `phenotype_names` | `list[str]` | Names of extra phenotype columns |
| `annotation_names` | `list[str]` | Names of annotation columns (from `a_` prefixed headers) |

### Layout Types

| Type | Key Fields | Description |
|------|------------|-------------|
| `NodeLayout` | `individual_id, x, y, sex, phenotype, is_deceased, is_proband, is_phantom, is_duplicate, is_carrier, is_obligate_carrier, label, extra_phenotypes, extra_phenotypes_raw, not_applicable, annotations, pgs_side, notation, notes, marker_depth, is_unknown_offspring` | Position and attributes of one node. `is_unknown_offspring` marks the placeholder child of a couple alone in its family: drawn as a "?", with no symbol, label or tooltip. `sex` is the shape drawn (from the gender identity when given); `is_deceased` includes stillbirths. `annotations: list[str]` contains annotation lines to display below the ID label, after `notes` ("SB", "ECT", "AMAB", "AFAB"). `marker_depth` is the room a childlessness marker takes below the symbol; the label moves down by it. `pgs_side: str` indicates PGS curve placement ("left" or "right"). |
| `EdgeLayout` | `start, end, edge_type, line_style, waypoints` | A connection line segment (waypoints for polyline/curved edges) |
| `TwinConnector` | `convergence_point, child_positions, twin_type` | V-lines for twin groups |
| `DuplicateArc` | `from_pos, to_pos` | Dotted arc between duplicate nodes |
| `FamilyLayout` | `nodes, edges, twin_connectors, duplicate_arcs` | Complete layout for one family |

---

## Configuration

### `PedigreeStyle`

All visual parameters in one dataclass. Pass to `compute_layout()` or `render_pedigree()`.

```python
from pypedigree import PedigreeStyle

style = PedigreeStyle(shape_size=1.0, dpi=300)
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `shape_size` | `0.8` | Size of individual shapes |
| `shape_line_width` | `1.5` | Shape border width |
| `fill_color` | `"black"` | Fill color for affected individuals |
| `outline_color` | `"black"` | Shape outline and deceased marker color |
| `phantom_line_style` | `"--"` | Line style for phantom parents |
| `line_color` | `"black"` | Connection line color |
| `line_width` | `1.0` | Connection line width |
| `double_line_gap` | `0.06` | Gap between double lines (consanguinity) |
| `broken_gap_size` | `0.15` | Gap size in broken lines |
| `carrier_notation` | `"nsgc2022"` | `"nsgc2022"`: the symbol is divided into patterned sections, one per carrier result, defined in the legend (Bennett et al., J Genet Couns 2022). `"legacy"`: half-filled shape for carriers, center dot for obligate carriers. Any other value raises `PedigreeError` |
| `carrier_dot_radius` | `0.08` | Radius of the obligate carrier center dot (legacy notation) |
| `font_size` | `8.0` | Label font size, in points. The legend panel is sized from it (text, swatches and row height), so raise it to keep a figure legible once scaled down for print |
| `min_font_size` | `None` | Smallest size of any text, in points (e.g. `6`, a common journal rule). The secondary text is drawn at fractions of `font_size` (notes and annotations 0.7, band text 0.75, quadrant text 0.65, "PGS" caption 0.6, legend 0.85 and less); each is raised to this minimum. The layout measures the lines below labels at that size and the legend grows with its text. `None`: no minimum. CLI: `--min-font-size` |
| `font_family` | `"sans-serif"` | Label font family |
| `label_offset` | `0.15` | Vertical offset for labels below shapes |
| `generation_gap` | `2.0` | Vertical distance between generations |
| `sibling_gap` | `1.5` | Horizontal distance between siblings (more where long labels need it) |
| `couple_gap` | `1.5` | Horizontal distance between partners (more where long labels need it) |
| `min_node_distance` | `1.2` | Deprecated, no effect: `sibling_gap` and `couple_gap` set the spacing. To be removed in 1.0 |
| `radial_gap` | `2.5` | Distance between concentric rings in circular layout |
| `base_radius` | `3.0` | Radius of the innermost ring in circular layout |
| `phenotype_colors` | `{}` | Map of phenotype name to hex color (e.g. `{"IQ": "#3498db"}`) |
| `show_phenotype_labels` | `False` | Show YES/NO text labels inside binary phenotype bands |
| `phenotype_layout` | `"auto"` | `"quadrants"` (the symbol cut in four, one quadrant per phenotype or `cr_` carrier source, at most four; more fall back to bands with a `PedigreeStyleWarning`), `"bands"` (one horizontal band per phenotype; the default before 0.5.0) or `"auto"` (quadrants up to four sources, bands beyond, without a warning). With carrier sources and no phenotype column nothing is cut: the symbol and the legend are drawn as with bands. Any other value raises `PedigreeError` |
| `phenotype_quadrants` | `{}` | Column name → `"top-left"`, `"top-right"`, `"bottom-left"` or `"bottom-right"`, for the quadrant layout. Unlisted columns take their conventional quadrant (autism top left, ADHD top right, intellectual disability bottom left, epilepsy bottom right) or the next free one in column order. Two columns in one quadrant raise `PedigreeError` |
| `show_na_labels` | `True` | Write "NA" in bands with no data (`-9`, `NA`, empty). When `False` they stay grey with no text, so they look like an "unknown" (`0`) band |
| `band_text` | `"auto"` | Text in bands and quadrants. `"auto"`: "name: value", or the value alone, or the value smaller (down to `min_font_size`, or half `font_size`), whichever fits; else none. `"name_value"`: always "name: value" (before 0.7.0). `"value"`: the value alone. `"none"`: no value ("NA" still follows `show_na_labels`) |
| `pgs_scale` | `"zscore"` | Position of a PGS value on its curve: `"zscore"` (a z-score on the drawn standard normal, clipped at about ±2.8 SD; the default since 0.10.0), `"minmax"` (between the pedigree's lowest and highest score, for raw scores) or `"percentile"` (a population percentile, 0-100). On the fixed scales the ±2 SD tails are marked and the SVG tooltip gives z and percentile. On `"zscore"`, a PGS column whose values do not look like z-scores draws a `PedigreeDataWarning` suggesting `"minmax"`. Any other value raises `PedigreeError` |
| `phenotype_labels` | `{}` | Name shown for a column (after its prefix), e.g. `{"SRS": "SRS T-score", "PGS_IQ": "IQ"}`: legend, quadrant key, band text, SVG tooltips, carrier rows. Others show their name with underscores as spaces; a PGS column in the legend's key, without "PGS_". Other style options keep using the column names |
| `legend` | `"auto"` | `"compact"`: with the quadrant layout, the binary columns of the key share one Affected / Unaffected / Unknown row (the key gives their colors). `"full"`: one row per column. `"auto"`: compact with quadrants, full with bands. `"minimal"` (journal figures): compact, with no title, no sex symbols and `legend_spacing` 0.7 unless those are set, a frame that hugs its content (no minimum width), a quadrant key of one symbol, and the key's binary columns read from a caption under it ("Colored: affected · white: unaffected · grey: unknown") instead of a row of swatches |
| `legend_title` | `"Legend"` | The legend's title; `None` or `""` for none |
| `legend_sex_symbols` | `True` | The "□ = Male ○ = Female" row; `False` drops it |
| `legend_max_width` | `None` | Widest the legend may be: a number in centimetres (`8.5`) or text with its unit (`"85mm"`, `"8.5cm"`, `"3.3in"`). The PGS row may grow up to it before wrapping (never wrapping before the widest other row ends, since the others do not wrap). `None`: the PGS row wraps where the widest other row ends: the caption on its first line, the score key on the lines below, broken between scores |
| `legend_text` | `{}` | The legend's fixed wording, by role: `key_note`, `shared_binary`, `affected`, `unaffected`, `unknown`, `pgs`, `ptdt`, `male`, `female`, `low`, `high` (`not_applicable` is still accepted, unused since 0.10.1: N/A has no legend row, #120) (e.g. `{"key_note": "Coloured: affected · white: unaffected · grey: unknown"}` for British spelling, or another language). `pgs` and `ptdt` replace the whole caption before the scores' key; a longer `male` moves the female symbol along. The defaults follow the docs' spelling ("color", "grey"). An unknown role raises `PedigreeError` with a suggested spelling. The title is `legend_title` |
| `legend_bar_width` | `"equal"` | With more than one gradient or binned row, `"equal"` draws every bar as wide as the widest (a gradient stretched, its end captions moving with its ends; a binned bar's steps sharing it equally, never narrower than they need); `"own"` keeps each bar at its own width, as before. The legend keeps its width. Anything else raises `PedigreeError` |
| `legend_spacing` | `1.0` | Factor on the legend's row height and gaps (header, captions, padding), 0.6 to 3, e.g. 0.7 for print. A row never gets shorter than its swatch and the captions under it. The figure is as much shorter as the legend |
| `show_pgs_values` | `False` | Write each PGS value beside its marker, in its color, left of the curve: "+2.1" / "−0.4" on `"zscore"`, the rounded percentile on `"percentile"`, three significant digits on `"minmax"`. Close values are spread apart; the layout widens the gaps by the widest value. CLI: `--show-pgs-values` |
| `pgs_colors` | `{}` | Map of PGS column name to hex color for the marker line. Unset: `#e74c3c` for a single PGS column; with several, distinct colors (known phenotypes' for `PGS_autism` and the like) |
| `gradient_range` | `{}` | Custom min/max for gradient normalization, e.g. `{"SRS": [40, 90]}`; the legend's bar then writes these ends instead of "Low" / "High" |
| `gradient_bins` | `{}` | Clinical thresholds: column → ascending cut-offs, e.g. `{"SRS": [60, 66, 76], "IQ": [70]}`. The band takes one of `len(cut-offs) + 1` steps of the column's color, white to full (a `gi_` column the other way), a value at a cut-off going to the step above; the value is still written. The legend shows the steps with their limits ("< 60", "60–65", "66–75", "≥ 76"; with fractional cut-offs a step is written up to the next one, which it excludes) and the SVG tooltip adds the step ("SRS: 70 (66–75)"). Takes precedence over `gradient_diverging` and `gradient_range` for the color. Anything but ascending numbers raises `PedigreeError` |
| `gradient_bin_label_position` | `"below"` | Where the legend writes the step names: `"below"` (under the limits) or `"inside"` (centred in each colored step, in a color that reads on it, `text_color_on`; the limits stay under the bar on one line, so the row is shorter; the bar grows to a line of names when needed and each step to its name). For all binned columns, or by column, e.g. `{"SRS": "inside"}` (others `"below"`). Anything else raises `PedigreeError` |
| `gradient_bin_labels` | `{}` | Names of the steps of `gradient_bins`: column → one name per step, from the lowest up (`len(cut-offs) + 1` entries, whatever direction a `gi_` column's color runs), e.g. `{"SRS": [None, "mild", "moderate", "severe"]}`; `None` or `""` leaves a step unnamed. The legend writes each name under its step's limits (the row takes a second caption line) and the SVG tooltip gives it with the step ("SRS: 81 (≥ 76, severe)"). A column without `gradient_bins`, or a wrong number of names, raises `PedigreeError` |
| `gradient_diverging` | `{}` | Tricolor diverging gradient config, e.g. `{"IQ": {"min": 40, "mid": 100, "max": 160, "color_low": "#e74c3c", "color_high": "#2ecc71"}}` |
| `zero_pad_ids` | `0` | Zero-pad sample IDs to this width in labels (0 = no padding) |
| `labels` | `"id"` | What labels a symbol (#125): `"id"` its sample ID, `"generation"` its NSGC number from the drawing (generation in Roman numerals from the top, position from the left: "II-3"; a duplicate has its original's), `"both"` the number with the ID right under it. `NodeLayout.number`; SVG tooltip "Number: II-3" and `data-number`. Another value raises `PedigreeError` |
| `dpi` | `150` | Output resolution |
| `fig_width` | `10.0` | Width in inches a small pedigree is drawn to fill (sets its scale; the figure is then fitted to the drawing) |
| `fig_height` | `8.0` | Height in inches a small pedigree is drawn to fill (sets its scale; the figure is then fitted to the drawing) |
| `margin` | `0.5` | Margin around the diagram |

**Methods**: `PedigreeStyle.from_dict(d)` — create from dictionary. `style.to_dict()` — convert to dictionary.

`PedigreeStyle.text_size(scale=1.0)` returns the size in points of a text drawn at `scale` times `font_size`, never below `min_font_size`: `PedigreeStyle(font_size=6, min_font_size=6).text_size(0.6)` is `6.0`. Both renderers and the layout take every text size from it.

### `DEFAULT_PHENOTYPE_PALETTE` and `EXTRA_PHENOTYPE_COLORS`

In `pypedigree.config`. `DEFAULT_PHENOTYPE_PALETTE`: the 8 hex colors given to columns without one in `phenotype_colors`, by their position among the columns: red, blue, orange, green, purple, teal, yellow, pink. `EXTRA_PHENOTYPE_COLORS`: brown, light magenta and olive, used when the palette runs short.

When some colors are fixed (a known phenotype's, or one in `phenotype_colors`), a column whose palette color repeats one of them, or looks like one, takes the next color that does not, from the palette then the extra colors. Colors look alike when their CIE76 distance is below 20 or they are green and teal, yellow and olive, or purple and light magenta. A figure with no known phenotype and no chosen color keeps the palette's order; past eight columns the extra colors come before it repeats.

### `KNOWN_PHENOTYPE_COLORS` and `known_phenotype(name)`

In `pypedigree.config`. Phenotypes with a conventional color, applied to any column without one in `phenotype_colors`, in both layouts: autism `#386eb2` (blue), intellectual disability `#d34436` (red), ADHD `#e59635` (yellow), epilepsy `#5e287e` (purple). `known_phenotype(name)` returns `"autism"`, `"intellectual disability"`, `"adhd"`, `"epilepsy"` or `None`; case, accents, `_` and `-` are ignored, and the aliases are `ASD`, `TSA`, `autisme`, `autism spectrum disorder`; `ID`, `DI`, `déficience intellectuelle`; `TDAH`; `épilepsie`. When such a column is present, the palette skips its own lookalike colors for the other columns.

### `PhenotypeDisplayMode`

Enum controlling how extra phenotype columns are rendered. Set via column name prefix in the .ped header:

| Prefix | Mode | Rendering |
|--------|------|-----------|
| `b_` | `BINARY` | Color (affected) / white (unaffected) / grey (unknown) / grey with "NA" (no data) bands. When all bands are NA, a single centered "NA" is shown. |
| `g_` | `GRADIENT` | White to color continuous gradient with value labels |
| `gi_` | `GRADIENT_INVERTED` | Color to white gradient (high value = white, low value = full color). Useful for metrics where low = bad (e.g., IQ: low IQ = more red). |
| `o_` | `OUTLINE` | Colored shape border (not a band). An individual affected in several `o_` columns gets one border per column, inner first, in column order (#46) |
| `c_` | `CATEGORICAL` | Distinct color per unique string value |
| `cr_` | `CARRIER` | Carrier status of one gene or condition (0 = not, 1 = carrier, 2 = obligate carrier), drawn as a patterned section, not a phenotype band |
| *(none)* | `AUTO` | Auto-detect from data (default, backward compatible) |

**Example header:**
```
FID  ID  Father  Mother  Sex  b_ASD  g_IQ  o_Proband  c_Diagnosis
```

The prefix is stripped from the display name (e.g. `b_ASD` is shown as `ASD`).

**Multiple PGS columns**: Any column whose name starts with `PGS` (case-insensitive) is treated as a PGS column -- e.g. `PGS_autism`, `PGS_IQ`, `PGS_ADHD`. All PGS values are rendered on a single gaussian curve with multiple colored markers. Marker colors are set via the `pgs_colors` style option. With several PGS columns and no `pgs_colors` entry, each gets its own color: a score named after a known phenotype (`PGS_autism`, `PGS_ADHD`, `PGS_ID`, `PGS_epilepsy`) takes that phenotype's color, the one its column has in the figure if there is one; the others take the first color unlike the other scores and the phenotype bands, starting from `#e74c3c`. A single PGS column keeps `#e74c3c`. `pgs_colors` entries always win. When PGS columns are present, `couple_gap` and `sibling_gap` are auto-increased to accommodate the curve. The curve is drawn left of the symbol, or right of it for the right-hand partner of a couple (`NodeLayout.pgs_side`), so that no box or value sits on the line between partners; its values go on the far side of the box from the symbol, and the layout keeps the box's room on that side.

**pTDT deviations** (#52): a column whose name starts with `pTDT` (case-insensitive, e.g. `pTDT_autism`) holds each child's polygenic transmission deviation, (PGS_child − PGS_mid-parent) / SD(PGS_mid-parent) (Weiner et al. 2017), N(0, 1) without over-transmission. It is drawn with the PGS curve, always on the z scale whatever `pgs_scale` (±2 SD marked), captioned "pTDT" under the curve; the legend has its own row, "pTDT deviation, (child − mid-parent) / SD, ±2 SD marked:" then its key, and the SVG tooltip reads "pTDT = +2.12" (no percentile). An individual's curve shows one kind: its pTDT deviations if it has any, otherwise its PGS, so the parents' PGS and the children's pTDT share a figure; an individual with both draws a `PedigreeDataWarning`. `pTDT_autism` may take the color `PGS_autism` has. pypedigree does not compute the deviations: give them in the file.

**Gradient value display**: Gradient band values (e.g. "IQ: 90", "SRS: 76") are always written, whatever the band's color intensity. Text in bands and quadrants fits its place (`band_text: "auto"`, the default since 0.7.0): a gradient band writes "IQ: 90", or "90" alone, or the value smaller (down to `min_font_size`, or half `font_size`), whichever first fits the band's width at that height (a circle's chord), and nothing when even that does not fit; quadrant text likewise. `band_text: "name_value"` always writes "IQ: 90" (as before 0.7.0), `"value"` the value alone, `"none"` no value.

**AUTO mode resolution**: A column with no prefix is resolved **once per column**, from all of its values across the pedigree:

| Values found in the column | Effective mode |
|---------------------------|----------------|
| Any string value | categorical |
| Any numeric value outside `0`, `1`, `2` | gradient |
| Otherwise (only `0`/`1`/`2`, or empty) | binary |

The decision is never taken per individual, so a column cannot mix binary and
gradient bands within the same diagram, and its legend row always matches its
bands. An explicit `gradient_range`, `gradient_diverging` or `gradient_bins`
entry for the column also forces gradient.

---

### `resolve_style(style, phenotype_names=None)`

Fill in the pedigree-dependent parts of a style: phenotype columns with no
explicit color get one from `DEFAULT_PHENOTYPE_PALETTE`, even when other
columns have one (they were grey before 0.6.1), never one that looks like a
color already in the figure. Bands are only drawn
when `phenotype_colors` is populated, so this is what makes a multi-phenotype
figure render. `render_pedigree()`, the CLI and `PedigreeBuilder` all call it,
which is why they produce identical figures.

```python
from pypedigree import parse_pedigree, resolve_style

ped = parse_pedigree("multi_pheno.ped")
style = resolve_style(None, ped.phenotype_names)
print(style.phenotype_colors)   # {'epilepsy': '#e74c3c', 'IQ': '#3498db', ...}
```

**Args**: `style` — `PedigreeStyle` or None for defaults. `phenotype_names` — ordered phenotype column names (e.g. `pedigree.phenotype_names`).

**Returns**: a new `PedigreeStyle`. The input is never mutated, and the
function is idempotent.

---

### `apply_pgs_spacing(style, phenotype_names=None, values_width=0.0)`

Widen `couple_gap` and `sibling_gap` to leave room for PGS distribution
curves drawn beside each node, and by `values_width` (layout units) for the
values written left of them with `show_pgs_values` (the renderers pass the
widest value's width). Returns `style` unchanged when no PGS column is present.

Unlike `resolve_style()` this is **not** idempotent — each call widens
further. `render_pedigree()` is its only caller; do not apply it yourself
before calling `render_pedigree()`.

**Args**: `style` — `PedigreeStyle`. `phenotype_names` — ordered phenotype column names. `values_width` — room for the PGS values, in layout units.

**Returns**: a new `PedigreeStyle`, or the input unchanged.

---

### `is_pgs_column(name)`

Whether a phenotype column name denotes a polygenic score, i.e. starts with
`PGS` (case-insensitive).

### `is_ptdt_column(name)`

Whether a phenotype column name holds pTDT deviations (#52), i.e. starts with
`pTDT` (case-insensitive). Such a column is drawn as a curve too, but
`is_pgs_column` is false for it.

**Returns**: `bool`.

---

## Exceptions

| Exception | Parent | When raised |
|-----------|--------|-------------|
| `PedigreeError` | `Exception` | Base class for all pypedigree errors; also raised for an invalid color in a style |
| `PedigreeParseError` | `PedigreeError` | File format or parsing issues |
| `PedigreeValidationError` | `PedigreeError` | Data validation failures (cycles, invalid references) |
| `PedigreeLayoutError` | `PedigreeError` | Layout computation failures |
| `PedigreeBatchError` | `PedigreeError` | `render_families`: some families could not be drawn (the others were). `.written`: family ID → path; `.failures`: family ID → message |

Every color in a style is validated before drawing starts, and all offending
values are reported in one error:

```
PedigreeError: Invalid color(s) in style: outline_color = '#12345';
phenotype_colors['Epilepsy'] = '#zzzzzz'. Expected hex values like
'#e74c3c' or named colors like 'grey'.
```

## Warnings

Non-fatal diagnostics use dedicated categories, so they can be filtered or
escalated independently of other libraries' warnings.

| Warning | Parent | When emitted |
|---------|--------|--------------|
| `PedigreeWarning` | `UserWarning` | Base class for all pypedigree warnings |
| `PedigreeStyleWarning` | `PedigreeWarning` | A style definition contains something unusable, e.g. an unknown option |
| `PedigreeDataWarning` | `PedigreeWarning` | Pedigree data is questionable but still renderable, e.g. a father recorded as female |

```python
import warnings
from pypedigree import PedigreeDataWarning, parse_pedigree

# Treat questionable data as a hard failure
warnings.simplefilter("error", PedigreeDataWarning)
parse_pedigree("family.ped")

# Or silence only pypedigree's data warnings
warnings.simplefilter("ignore", PedigreeDataWarning)
```

An unknown style option is ignored but reported, with a suggestion when one is
close:

```python
PedigreeStyle.from_dict({"phenotype_color": {...}})
# PedigreeStyleWarning: Unknown style option 'phenotype_color' ignored;
# did you mean 'phenotype_colors'?
```

Pass `from_dict(d, warn_unknown=False)` to accept a dictionary that
legitimately carries extra data.

The CLI captures these and prints them as `Warning: <message>` on stderr.
