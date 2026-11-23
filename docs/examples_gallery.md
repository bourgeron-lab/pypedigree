# Examples Gallery

Each example file in `tests/data/` demonstrates specific pypedigree features, organized by subdirectory. Render any example with:

```bash
pypedigree tests/data/<subdir>/<file>.ped -o output.png
```

---

## simple_trio.ped

**Description**: Father, mother, and one affected child.

**Features**: Basic trio layout, affected (filled) shape, label rendering.

```bash
pypedigree tests/data/basic/simple_trio.ped -o simple_trio.png
```

---

## two_gen_siblings.ped

**Description**: Two parents and three children.

**Features**: Sibship line connecting multiple siblings, horizontal ordering.

```bash
pypedigree tests/data/basic/two_gen_siblings.ped -o two_gen_siblings.png
```

---

## three_gen.ped

**Description**: Grandparents, parents, and two grandchildren.

**Features**: Multi-generation layout with three vertical levels, descent lines.

```bash
pypedigree tests/data/basic/three_gen.ped -o three_gen.png
```

---

## extended_markers.ped

**Description**: Three generations with a deceased grandfather and a proband child. Uses the extended 10-column format.

**Features**: Deceased marker (diagonal slash), proband marker ("P" annotation), extended column parsing (columns 7-10: DOB, twins, deceased, proband).

```bash
pypedigree tests/data/extended/extended_markers.ped -o extended_markers.png
```

---

## consanguinity.ped

**Description**: First-cousin marriage across three generations. Grandparents produce two children who marry unrelated spouses. Their children (first cousins) marry each other.

**Features**: Automatic consanguinity detection via shared ancestor analysis, double relationship line for consanguinous unions.

```bash
pypedigree tests/data/consanguinity/consanguinity.ped -o consanguinity.png
```

---

## twins.ped

**Description**: Parents with MZ (monozygotic) twins and additional siblings, with dates of birth for ordering.

**Features**: Twin connectors (V-lines from sibship to children), MZ horizontal bar, DOB-based sibling ordering (older siblings to the left).

```bash
pypedigree tests/data/twins/twins.ped -o twins.png
```

---

## multiple_relationships.ped

**Description**: One individual (father) with two different partners and children from each union.

**Features**: Duplicate node creation (same individual appears twice), dotted arc connecting duplicate instances, separate sub-tree layout for each union.

```bash
pypedigree tests/data/multi_partner/multiple_relationships.ped -o multiple_relationships.png
```

---

## Style customization

**Description**: Any pedigree can be rendered with custom visual styles using the `--style` JSON file or direct CLI options.

**Example — grey outlines**: Render shapes and deceased markers with grey outlines instead of black.

```bash
# Via CLI option
pypedigree tests/data/extended/extended_markers.ped --outline-color grey -o grey_outlines.png

# Via JSON style file
echo '{"outline_color": "grey"}' > style.json
pypedigree tests/data/extended/extended_markers.ped --style style.json -o grey_outlines.png
```

See the [PedigreeStyle reference](api.md#pedigreestyle) for all available parameters.

---

## Multi-phenotype rendering

**Description**: Display multiple phenotypes inside each individual's shape: one quadrant each up to four phenotypes (the default since 0.5.0), or horizontal color bands (`--phenotype-layout bands`, and automatically beyond four). Binary phenotypes show color (affected), white (unaffected), grey (unknown/0), or grey with "NA" (no data/-9). When all bands are NA, a single centered "NA" is shown. Continuous phenotypes show a color gradient proportional to the normalized value. The legend is sized from `font_size`: its text, swatches and rows keep the same physical size whatever the size of the pedigree, and the figure grows downward to make room for it.

**With inline phenotypes** (columns 11+ in the .ped file, requires header):

```bash
pypedigree mockup_v2.ped -o mockup_v2.png
```

**With a separate .pheno file** + custom colors:

```bash
pypedigree family.ped --pheno data.pheno --style style.json -o output.png
```

Where `style.json` maps phenotype names to colors:

```json
{"phenotype_colors": {"epilepsy": "#e74c3c", "IQ": "#3498db", "SRS": "#e67e22"}}
```

Colors are auto-assigned from a default palette if no `phenotype_colors` are provided. Autism (`ASD`, `TSA`), intellectual disability (`ID`, `DI`), ADHD (`TDAH`) and epilepsy columns get conventional colors instead: blue, red, yellow and purple; the palette then skips its own blue, red, orange/yellow and purple for the other columns, and no two columns get colors that look alike (green and teal, say): the next distinct color is taken, brown, light magenta or olive when the palette runs short.

---

## couple_only.ped (parents without children)

**Description**: A family with only the two parents in the file, as when a trio's child is missing. Instead of two unrelated symbols, pypedigree draws their couple, father on the left even though the mother comes first in the file, with a descent line to a "?". The rule applies to a family of exactly two founders, a man and a woman; the pedigree data are not changed (an exported `.fam` holds the two parents only).

```bash
pypedigree tests/data/basic/couple_only.ped --phenotype-layout quadrants -o couple_only.png
```

---

## quadrants.ped (four phenotypes in quadrants)

![Four phenotypes, one quadrant each](imgs/examples/quadrants.png)

**Description**: Synthetic three-generation pedigree, written by `tests/data/generate_fixtures.py` from a fixed seed, with four binary phenotype columns (`b_Autism`, `b_ADHD`, `b_ID`, `b_Epilepsy`) and an `o_WGS` outline for the individuals with sequencing data, green: the automatic colors keep it unlike the four phenotypes' (a red outline would look like the intellectual disability quadrant).

**Features**: the quadrant layout (the default since 0.5.0) cuts each symbol in four. Autism is top left in blue, ADHD top right in yellow, intellectual disability bottom left in red, epilepsy bottom right in purple: these four phenotypes are recognised by name (whatever the column order) and get their conventional quadrant and color. Other columns take the free quadrants in column order, and `phenotype_quadrants` in a style file places any column explicitly. A `cr_` carrier source takes a quadrant too, patterned. The legend starts with a key, a square and a circle cut in four with the column names around them. Since the key gives each column's color, the binary columns then share one Affected / Unaffected / Unknown row, its first swatch striped with their colors (`"legend": "full"` in a style file gives one row per column). A symbol with no data at all (316) is grey with a single "NA". `"phenotype_labels": {"ID": "Intellectual disability"}` shows a column under another name.

```bash
pypedigree tests/data/phenotype/quadrants.ped --phenotype-layout quadrants -o quadrants.png
```

---

## outline_rings.ped (several outline columns)

![One border per outline column](imgs/examples/outline_rings.png)

**Description**: Synthetic three-generation family with one binary phenotype (`b_ASD`) and three outline columns saying which data each individual has: `o_WGS`, `o_WES` and `o_Array`.

**Features**: an individual affected in several `o_` columns gets one colored border per column, the first on the outline, the others around it, outwards, in column order (III-1 has the three, inner orange for WGS). One affected in a single column keeps the usual border (II-3, III-2). All symbols are spaced and labelled for the most borders anyone has, so the labels of a generation stay level, and the proband's arrow points at the outer border. The legend adds "Several borders: one per column, inner first"; in the interactive SVG the tooltip lists them ("Borders: WGS, WES, Array"). More than three borders draw a warning: the symbols then take much more room.

```bash
pypedigree tests/data/phenotype/outline_rings.ped -o outline_rings.png
```

---

## gradient_bins.ped (clinical thresholds)

![Scores read against their cut-offs](imgs/examples/gradient_bins.png)

**Description**: Synthetic family of five with an autism diagnosis (`b_autism`) and three scores: an SRS T-score (`g_SRS`), an IQ (`gi_IQ`, low is atypical) and an ADHD rating scale (`g_ADHD_RS`). The style (`gradient_bins.json`) gives the SRS its cut-offs (60 mild, 66 moderate, 76 severe) and their names, the IQ one at 70, and the ADHD-RS its range (0 to 54).

**Features**: with `gradient_bins`, a score is drawn in one of `len(cut-offs) + 1` steps of its color instead of a continuous shade, so the reader sees which side of a cut-off it falls (a value at a cut-off goes to the step above). The legend shows each step with its limits ("< 60", "60–65", "66–75", "≥ 76"); the inverted IQ is colored below 70 only. The value is still written in the band, and the SVG tooltip adds the step ("SRS: 70 (66–75, moderate)"). `gradient_bin_labels` names the steps, from the lowest up: "mild", "moderate", "severe" under the SRS limits (none below 60), "intellectual disability" under the IQ's first step. With `gradient_range`, the legend's bar writes the range's ends ("0", "54") instead of "Low" / "High".

```bash
pypedigree tests/data/phenotype/gradient_bins.ped --style tests/data/phenotype/gradient_bins.json -o gradient_bins.png
```

---

## legend_minimal.ped (a compact legend for journal figures)

![The minimal legend](imgs/examples/legend_minimal.png)

**Description**: Synthetic quartet with two binary phenotypes (`b_autism`, `b_ADHD`), a gradient (`g_SRS`) and two polygenic scores, drawn at 6 pt as in a journal figure (`legend_minimal.json`: `"legend": "minimal"`, `font_size` and `min_font_size` 6).

**Features**: `legend: "minimal"` drops the "Legend" title and the "= Male / = Female" row, tightens the rows (`legend_spacing` 0.7), draws the quadrant key with one symbol, and replaces the Affected / Unaffected / Unknown swatches of the key's binary columns by one caption under the key. Here the legend is 2.7 cm tall instead of 4.5 cm, and the figure 1.8 cm shorter. Each part can also be set alone: `legend_title` (`null` for none), `legend_sex_symbols` (`false`) and `legend_spacing`.

```bash
pypedigree tests/data/phenotype/legend_minimal.ped --style tests/data/phenotype/legend_minimal.json -o legend_minimal.png
```

---

## PGS distribution

**Description**: If a column named `PGS` (case-insensitive) is present in the data, an upright gaussian distribution curve, its peak to the left (turned clockwise, the usual bell with positive values on the right), is drawn to the left of each individual with a red marker indicating their score.

```bash
# PGS as column 14 (inline)
pypedigree mockup_with_pgs.ped -o output.png

# PGS in a separate .pheno file
pypedigree family.ped --pheno data_with_pgs.pheno -o output.png

# Ancestry-matched z-scores on a fixed scale, +/-2 SD marked
pypedigree family.ped --pgs-scale zscore -o output.png

# The same, with each value written beside its marker (+2.1, −0.4)
pypedigree family.ped --pgs-scale zscore --show-pgs-values -o output.png
```

By default the marker sits between the pedigree's lowest and highest score. With `--pgs-scale zscore` or `percentile` it is placed on the population's standard normal, so the same score is drawn at the same place in every pedigree.

---

## Annotations

**Description**: Columns with the `a_` prefix in the .ped header are annotation columns. Their values are displayed as italic grey text below the individual's ID label. Missing values (0, -9, NA, empty) produce no annotation line. Annotations can also be loaded from a separate file via `--annotations`.

**With inline annotations** (header-based `a_` prefix columns):

```bash
pypedigree tests/data/phenotype/test_all.ped -o annotations.png
```

**With a separate annotation file**:

```bash
pypedigree family.ped --annotations annotations.tsv -o annotations.png
```

**Python API**:

```python
from pypedigree import PedigreeBuilder

builder = PedigreeBuilder("FAM1")
builder.add_individual("father", sex="male", annotations=["Clinic: Paris", "Ref: 42"])
builder.add_individual("mother", sex="female")
builder.add_individual("child", sex="male", father="father", mother="mother")
fig = builder.render(output="annotations.png")
```

---

## Diverging gradient

**Description**: Use `gradient_diverging` in the style JSON to define a tricolor gradient for continuous phenotypes. Values below the midpoint blend from `color_low` to white, values above blend from white to `color_high`.

```bash
echo '{"gradient_diverging": {"IQ": {"min": 40, "mid": 100, "max": 160, "color_low": "#e74c3c", "color_high": "#2ecc71"}}}' > style.json
pypedigree family.ped --style style.json -o diverging.png
```

---

## ptdt.ped (pTDT deviations)

![The parents' PGS, the children's pTDT](imgs/examples/ptdt.png)

**Description**: Synthetic family with two autistic children and two unaffected siblings. The parents have an autism PGS (`PGS_autism`), the children their pTDT deviation (`pTDT_autism`): (PGS_child − PGS_mid-parent) / SD(PGS_mid-parent), +2.1 and +1.4 for the autistic children, about 0 for their siblings (`ptdt.json`: `pgs_scale` "zscore", `show_pgs_values`).

**Features**: a `pTDT_` column is drawn with the PGS curve, on the z scale whatever `pgs_scale`, and captioned "pTDT", so the reader sees what each curve shows. The legend has a row for each kind; the SVG tooltip reads "pTDT = +2.10". An individual's curve shows its pTDT deviations if it has any, otherwise its PGS.

```bash
pypedigree tests/data/phenotype/ptdt.ped --style tests/data/phenotype/ptdt.json -o ptdt.png
```

---

## Multiple PGS columns

**Description**: Any column starting with `PGS` (case-insensitive) is treated as a PGS column. Multiple PGS values are displayed on a single gaussian curve with distinct colored markers.

```
FID	ID	Father	Mother	Sex	PGS_autism	PGS_IQ	PGS_ADHD
FAM1	1	0	0	1	0.72	-0.35	0.10
```

```bash
echo '{"pgs_colors": {"PGS_autism": "#e74c3c", "PGS_IQ": "#3498db", "PGS_ADHD": "#e67e22"}}' > style.json
pypedigree family.ped --style style.json -o multi_pgs.png
```

Without `pgs_colors`, the scores still get distinct colors: `PGS_autism` autism's blue, `PGS_ADHD` ADHD's orange (unless a band of the figure already uses an orange), and `PGS_IQ` the red.

---

## One file per family

**Description**: `--output-dir` draws each family of a file into its own file, named after the family ID, in one run. Every file has the colors, column modes and legend of the whole file, so the figures of a cohort can be compared side by side.

```bash
pypedigree tests/data/basic/multi_family.ped --output-dir figures/
# Saved: figures/FAM1.png
# Saved: figures/FAM2.png
pypedigree tests/data/basic/multi_family.ped --output-dir figures/ --interactive
```

From Python: `render_families(pedigree, "figures", fmt="svg")`.

---

## Minimum text size (journal figures)

**Description**: Journals often require every text of a figure to be at least 6 pt at print size. Notes, annotations, values in the symbols, "PGS" captions and the legend are drawn smaller than the ID labels; `--min-font-size` raises them to a minimum without enlarging the labels. The layout and the legend make room for the larger text.

```bash
echo '{"font_size": 6}' > style.json
pypedigree tests/data/phenotype/test_all.ped --style style.json --min-font-size 6 -o journal.pdf
```

---

## Zero-padded IDs

**Description**: Use the `zero_pad_ids` style option to zero-pad sample IDs in labels. For example, with `zero_pad_ids: 3`, ID `1` displays as `001`.

```bash
echo '{"zero_pad_ids": 3}' > style.json
pypedigree family.ped --style style.json -o padded_ids.png
```

---

## Gradient inverted (gi_ prefix)

**Description**: The `gi_` prefix activates `GRADIENT_INVERTED` mode, producing a color-to-white scale where high values are white and low values are the full color. This is useful for metrics where low values indicate concern (e.g., IQ: low IQ = more red).

```
FID	ID	Father	Mother	Sex	gi_IQ
FAM1	1	0	0	1	70
FAM1	2	0	0	2	120
```

```bash
pypedigree family_with_gi.ped -o gradient_inverted.png
```

---

## Rendering all examples at once

```bash
for f in tests/data/**/*.ped; do
    pypedigree "$f" -o "/tmp/$(basename "${f%.ped}.png")"
done
```

## carrier.ped

**Description**: Three-generation family demonstrating carrier notation. Grandmother, father and child2 are carriers, aunt and cousin1 are obligate carriers.

**Features**: NSGC 2022 carrier notation: carriers have a patterned symbol, obligate carriers a sparser version of the same pattern, both defined in the legend. `--carrier-notation legacy` draws the pre-2022 half-filled shapes and center dots instead.

```bash
pypedigree tests/data/phenotype/carrier.ped -o carrier.png
pypedigree tests/data/phenotype/carrier.ped --carrier-notation legacy -o carrier_legacy.png
```

---

## carrier_genes.ped

**Description**: Carrier status for two genes, one `cr_` column each (`cr_CFTR`, `cr_HEXA`). The father carries both and his symbol is split in two patterns; child1 is affected and an obligate CFTR carrier, so his symbol is half filled, half patterned. A missing value (`-9`) is not drawn.

**Features**: `cr_` carrier columns, one pattern per gene, symbols divided as NSGC 2022 recommends.

```bash
pypedigree tests/data/phenotype/carrier_genes.ped -o carrier_genes.png
```

---

## nomenclature.ped

**Description**: Every standard notation symbol beyond sex, affection status and death (NSGC nomenclature, 2008 and its 2022 revision). FAM1: II-1 is infertile (two bars, reason in an `a_` column), II-2 the consultand (arrow without "P"), II-4 a transgender woman (circle, AMAB below), II-5 non-binary and childless by choice (diamond, one bar, AFAB below), II-6 a miscarriage (triangle). In the third generation: the proband, a child adopted in (brackets, dashed line), a child adopted out (brackets), a stillbirth ("SB"), an affected termination (filled, struck-through triangle), an ectopic pregnancy ("ECT") and an ongoing pregnancy ("P"). FAM2: a sperm donor ("D") and an affected pregnancy.

**Features**: `pregnancy`, `adopted`, `childless`, `art`, `gender_identity` and `consultand` columns; gestational ages in an `a_` column.

```bash
pypedigree tests/data/nomenclature/nomenclature.ped -o nomenclature.png
pypedigree tests/data/nomenclature/nomenclature.ped --interactive -o nomenclature.svg
```

---

## canrisk_v4.txt

**Description**: Two invented families saved as a CanRisk v4 file (#109). FAM1: the target (the consultand, also the proband since she has breast cancer at 35 and 37, BRCA1 positive, triple negative), her mother (ovarian cancer, deceased), her father (prostate cancer, HOXB13 heterozygous), identical twin sisters, a brother who tested negative, and grandparents with pancreatic and breast cancer. Breast and ovarian PRS z-scores for the target. FAM2: an unaffected target with a prostate PRS, his father with pancreatic cancer.

**Features**: the CanRisk reader: one column per cancer, ovarian and prostate cancer "not applicable" for the other sex assigned at birth (white, as unaffected), ages, gene tests and pathology as annotations in NSGC order, PRS on the z-score scale.

![A CanRisk file](imgs/examples/canrisk.png)

```bash
pypedigree tests/data/format/canrisk_v4.txt -o canrisk.png
pypedigree tests/data/format/canrisk_v4.txt --interactive -o canrisk.svg
```

---

## Interactive SVG

**Description**: Any pedigree can be exported as an interactive SVG with tooltips and click-to-highlight. It carries the same phenotype bands, legend, PGS curves and annotations as the PNG.

```bash
# Generate interactive SVG
pypedigree tests/data/basic/three_gen.ped --interactive -o interactive.svg

# Multi-phenotype pedigree with its style file
pypedigree tests/data/phenotype/test_all.ped --style tests/data/phenotype/test_all_style.json \
    --interactive -o test_all.svg

# Open in browser for tooltips and click highlighting
open interactive.svg
```

**Python API**:

```python
from pypedigree import parse_pedigree, render_interactive_svg

ped = parse_pedigree("tests/data/basic/three_gen.ped")
render_interactive_svg(ped, output="interactive.svg")
```

---

## Programmatic Builder API

Build pedigrees from code without writing .ped files:

```python
from pypedigree import PedigreeBuilder

builder = PedigreeBuilder("FAM1")
builder.add_individual("gf", sex="male")
builder.add_individual("gm", sex="female")
builder.add_individual("father", sex="male", father="gf", mother="gm")
builder.add_individual("mother", sex="female")
builder.add_individual("child", sex="male", father="father", mother="mother",
                       phenotype="affected", phenotypes={"IQ": 85, "epilepsy": 1})
fig = builder.render(output="builder_example.png")
```

---

## large_multi_union.ped

**Description**: Synthetic pedigree, written by `tests/data/generate_fixtures.py` from a fixed seed — 55 individuals across 3 generations (a founder couple, their thirteen children and twelve spouses, twenty-eight grandchildren) with 3 phenotype columns.

**Features**: Smart header-based column mapping, prefix-based display modes (`b_ASD`, `b_Dyslexia`, `o_Sleep disorder`), Deceased and Proband columns, 0 vs -9 distinction (unknown vs no data), legend sized from the font, large family layout.

```bash
# --style takes a path to a JSON file:
cat > style.json <<'EOF'
{"phenotype_colors": {"ASD": "#e74c3c", "Dyslexia": "#3498db", "Sleep disorder": "#e67e22"}}
EOF
pypedigree tests/data/complex/large_multi_union.ped --style style.json -o large_multi_union.png
```

---

## multigenerational.ped

![Five-generation pedigree with two loops](imgs/examples/multigenerational.png)

Detail of the first branch: 401 and 406 are first cousins, grandchildren of 201. Their parents' sibships are ordered so they meet: 401 is the last of his sibship, 406 the first of hers.

![Detail: the first-cousin loop](imgs/examples/multigenerational_detail.png)

**Description**: Synthetic pedigree, written by `tests/data/generate_fixtures.py` from a fixed seed. 81 individuals over five generations descending from one founder couple. Two loops: a first-cousin marriage (401 and 406, double line) and a marriage between the second and fourth branches (412 and 429, second cousins). Both are drawn with the partners side by side and no duplicated individual: each sibship on the way down from the common ancestors is ordered so the two branches meet. The affected children below both loops (501, the proband, 504, 505 and 506) have unaffected obligate-carrier parents; their brother 502 is a carrier.

**Features**: loops, `b_ASD` binary band, `g_IQ` gradient band with values, `cr_HEXA` carrier band (carrier and obligate carrier), dizygotic (318-319) and monozygotic (410-411) twins, deceased individuals with their age at death in an `a_` column, a miscarriage (503) and a child adopted in (436), `-9` no-data values for the oldest generations. Identifiers are cohort-style (generation, then rank in the file), not the clinical left-to-right numbering, since loop ordering moves siblings.

```bash
pypedigree tests/data/complex/multigenerational.ped --phenotype-layout bands -o multigenerational.png
# Interactive version, to zoom in and hover over individuals
pypedigree tests/data/complex/multigenerational.ped --interactive -o multigenerational.svg
```

---

## Circular / Radial Layout

**Description**: Any pedigree can be rendered in a circular/radial layout where founders sit at the center and each subsequent generation occupies a concentric ring outward, in the left-to-right order of the default layout. Twins, duplicates and loops are drawn as in the default layout. Edges are drawn as polylines following waypoints for smooth curved routing.

```bash
# Circular layout from CLI
pypedigree tests/data/basic/three_gen.ped --circular -o three_gen_circular.png

# Combine with other options
pypedigree tests/data/complex/large_multi_union.ped --circular --dpi 300 -o large_multi_union_circular.png
```

**Python API**:

```python
from pypedigree import parse_pedigree, render_pedigree, PedigreeStyle

ped = parse_pedigree("tests/data/basic/three_gen.ped")
fig = render_pedigree(ped, circular=True, output="three_gen_circular.png")

# Tune ring spacing
style = PedigreeStyle(radial_gap=3.0, base_radius=2.5)
fig = render_pedigree(ped, style=style, circular=True, output="custom_circular.png")
```

**Builder API**:

```python
from pypedigree import PedigreeBuilder

builder = PedigreeBuilder("FAM1")
builder.add_individual("gf", sex="male")
builder.add_individual("gm", sex="female")
builder.add_individual("father", sex="male", father="gf", mother="gm")
builder.add_individual("mother", sex="female")
builder.add_individual("child", sex="male", father="father", mother="mother", phenotype="affected")
fig = builder.render(circular=True, output="builder_circular.png")
```

---

## Python API

All examples can also be rendered programmatically:

```python
from pypedigree import parse_pedigree, render_pedigree, PedigreeStyle

ped = parse_pedigree("tests/data/consanguinity/consanguinity.ped")
fig = render_pedigree(ped, style=PedigreeStyle(dpi=300), output="consanguinity.png")
```
