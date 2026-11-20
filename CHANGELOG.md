# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- pTDT deviation columns (#52): a column whose name starts with `pTDT` (`pTDT_autism`) holds each
  child's polygenic transmission deviation, (PGS_child − PGS_mid-parent) / SD(PGS_mid-parent). It
  is drawn with the PGS curve, always on the z scale, captioned "pTDT" instead of "PGS"; the
  legend has its own row ("pTDT deviation, (child − mid-parent) / SD, ±2 SD marked:" and its
  key) and the SVG tooltip reads "pTDT = +2.12", with no percentile. An individual's curve shows
  its pTDT deviations if it has any, otherwise its PGS, so the parents' PGS and the children's
  pTDT share a figure; an individual with both draws a `PedigreeDataWarning`. New
  `is_ptdt_column()`; new gallery example `tests/data/phenotype/ptdt.ped`. pypedigree does not
  compute the deviations (the issue's second step).
- A more compact legend for journal figures (#67). New style options: `legend_title` (`None`
  or `""` drops the "Legend" title), `legend_sex_symbols` (`false` drops the "= Male / =
  Female" row) and `legend_spacing` (a factor on the row height and the gaps, 0.6 to 3, e.g.
  0.7 for print; a row never gets shorter than its swatch and the captions under it). The
  figure is as much shorter as its legend. `legend: "minimal"` sets all three (no title, no
  sex symbols, 0.7, unless given) and draws the quadrant key with one symbol, the key's binary
  columns read from a caption under it ("Colored: affected · white: unaffected · grey:
  unknown") instead of a row of swatches, also with a single binary column. Every option
  defaults to the legend drawn before: no figure changes.
- `legend_max_width` style option: the widest the legend may be, in centimetres (`8.5`) or with
  a unit (`"85mm"`, `"8.5cm"`, `"3.3in"`), up to which the PGS row may grow before it wraps
  (#68).

- `gradient_bins` style option, for scores read against clinical cut-offs (#57): column →
  ascending cut-offs, e.g. `{"SRS": [60, 66, 76], "IQ": [70]}`. The band is filled with one of
  `len(cut-offs) + 1` steps of the column's color instead of a continuous shade (a value at a
  cut-off goes to the step above; a `gi_` column such as IQ, where low is atypical, bins the
  other way), and the legend shows the steps with their limits ("< 60", "60–65", "66–75",
  "≥ 76"). The value is still written in the band; the SVG tooltip adds its step ("SRS: 70
  (66–75)"). New gallery example: `tests/data/phenotype/gradient_bins.ped`.

- `gradient_bin_labels` style option (#81): names for the steps of `gradient_bins`, one per step
  from the lowest up, e.g. `{"SRS": [null, "mild", "moderate", "severe"]}` (`null` or `""` for an
  unnamed step; a `gi_` column lists them from its lowest step up too). The legend writes each
  name under its step's limits, the step as wide as its name; the SVG tooltip gives it with the
  step ("SRS: 81 (≥ 76, severe)"). The gallery's `gradient_bins` example names its steps.

- `legend_text` style option (#79): the legend's fixed wording by role (`key_note`,
  `shared_binary`, `affected`, `unaffected`, `unknown`, `pgs`, `ptdt`, `male`, `female`, `low`,
  `high`), for British spelling ("Coloured") or a figure in another language, in the PNG and
  SVG legends; a wrapped PGS row measures the new text, and a longer "male" text moves the
  female symbol along. An unknown role raises `PedigreeError` with a suggested spelling. The
  default wording is unchanged: it follows the documentation's spelling ("color", "grey").

### Changed

- With `legend: "minimal"`, the legend's frame hugs its content: no minimum width (#79). A
  minimal legend with one short row gets a frame as narrow as that row.
- A wrapped legend row keeps "±2 SD" on one line (#68 follow-up): its words are joined by a
  no-break space.
- With `gradient_range`, the legend's gradient bar writes the range's ends instead of "Low" /
  "High" (#57).

- The PGS legend row wraps instead of widening the legend (#68): when the caption and the
  score key are longer than the legend's widest other row (or than `legend_max_width`), the key
  goes on the lines below the caption, under its text, wrapped between scores (never inside a
  name), and a caption longer than that width wraps between words. The legend is that much
  taller. With four scores at 6 pt the legend of the issue's example is 8.3 cm wide instead of
  12.2 cm, so it fits a journal column. A row that fits, such as a single short score, is
  drawn on one line as before. PNG and SVG break at the same places (DejaVu Sans widths). Two
  goldens change (`multi_phenotype_pgs`, `circular_test_all`), and the gallery's
  `annotations.png`.

### Fixed

- The names of the legend's gradient, binned, categorical and binary rows start in one column
  again (#88). Since step names (#81) widen a binned row's steps, its name followed its own bar,
  so the names were staggered; the column is now past the widest binned bar, for every row, and
  the legend's width and the PGS row's wrapping (#68) follow it.
- The deceased line no longer reaches into the PGS box (#80): it overshoots the symbol by 0.2
  half-sizes at each end, while the box sat 0.15 beside it, so its end met the box's corner,
  where a low score's marker sits. Beside a struck-through symbol the box now keeps 0.12
  half-sizes clear of the line's end, on either side, and the layout counts the extra.
  Other symbols' boxes do not move.
- The right-hand partner's PGS curve and its values are drawn on their right, not on the line
  between the partners, which crossed the box, the curve, the dots and the values (#78). The
  layout already chose a side (`NodeLayout.pgs_side`); both renderers now follow it, the values
  going on the far side of the box. Someone with a partner on each side keeps the left. The
  layout keeps the box's room right of such a partner, so their right-hand neighbour's box does
  not meet it; an evaluated partner's box clears the asterisk; a proband's caption stays below
  a box on the right. Goldens `ptdt`, `legend_minimal`, `multi_phenotype_pgs` and
  `circular_test_all` change, and the gallery's `annotations.png`, `legend_minimal.png` and
  `ptdt.png`.
- The captions under the legend's last row (a gradient's ends, a binned row's limits, category
  names) no longer run past the legend's frame, nor the figure's edge at small sizes (#77): the
  panel gets a bottom margin for the part of them that hangs below the row. Legends whose last
  row has no captions keep their height. Four goldens change (`gradient_bins`,
  `multi_phenotype_inline`, `pheno_file`, `multi_phenotype_no_labels`: a gradient last), and
  the gallery's `gradient_bins.png`.
- Text enlarged by `min_font_size` keeps its gaps (#66). The notes and annotations below an ID
  label ("AFAB", "SB", `a_` columns) are spaced by their drawn size: each line starts at least
  the height of the line above plus 0.3 em below it, so they no longer read as one block with
  the label (or overlap it, when a pedigree is drawn into a small `ax`); the generation gap
  follows. The "PGS" caption keeps half a marker dot plus 0.3 em from the curve's box, so a
  dot on the box's edge (the family's lowest or highest score) no longer touches it. The PNG
  sets both gaps when the figure is drawn, at its scale, as the SVG does at its own. At the
  default sizes the fixed distances used before are larger, so figures are unchanged.

## [0.8.1] - 2026-11-19

### Added

- An individual affected in several outline (`o_`) columns gets one colored border per column
  (#46): the first on the symbol's outline, as before, the others around it, outwards, in
  column order, each one border's width and a white gap further out. Then every symbol is
  spaced and labelled for that many borders, so the labels of a generation stay level; the
  proband's arrow, the PGS curve, adoption brackets, the documented-evaluation asterisk and the
  deceased line go round the outer border. The legend adds "Several borders: one per column,
  inner first", and the SVG tooltip lists them ("Borders: WGS, WES"). Above three borders, a
  `PedigreeStyleWarning` says the symbols take much more room. Figures where nobody has more
  than one border are unchanged. New gallery example: `tests/data/phenotype/outline_rings.ped`.

### Fixed

- Only the first `o_` column an individual was affected in was drawn: the others did not show
  on the symbol, while the legend listed them all (#46).

## [0.8.0] - 2026-11-19

### Added

- `phenotype_labels` style option: the name a column is shown with, e.g.
  `{"SRS": "SRS T-score", "PGS_IQ": "IQ"}`, in the legend, the quadrant key, band text, the SVG
  tooltips and carrier rows (#56). The style's other options still use the column names.
- `legend` style option: `"auto"` (default), `"compact"` or `"full"` (#56).
- The PGS legend row says which color is which score: each name follows a dash of its color
  ("PGS distribution: -- autism -- ADHD"), instead of the names in parentheses (#56).
- `show_pgs_values` style option, `--show-pgs-values` on the command line and "Polygenic score
  values" in the web app: each PGS value is written beside its marker, in its color, so a
  printed figure says how high a score is, not only that it is high (#54). "+2.1" or "−0.4" on
  the z scale (a typographic minus), the rounded percentile, or the value to three significant
  digits on "minmax". The values go left of the curve's box, at their marker's height, spread
  apart when close; the layout widens the gaps by the widest value. Off by default.

### Changed

- Column names are shown with underscores as spaces ("intellectual_disability" ->
  "intellectual disability"), in the legend, band text and tooltips (#56); figures whose
  column names have none are unchanged.
- With the quadrant layout the legend is compact by default (`legend: "auto"`): the quadrant
  key gives each column's color, so two or more binary columns share one Affected /
  Unaffected / Unknown row, its first swatch striped with their colors, instead of one row
  each. With four phenotypes the legend loses three rows, which matters in a small journal
  figure. Gradient, categorical, outline, PGS and carrier rows are unchanged; `legend:
  "full"` restores one row per column. The golden images with quadrants are redrawn.

### Fixed

- The dots of two PGS markers with close values no longer hide each other: dots closer than
  their diameter are set side by side across the curve's box, the dashed lines staying at the
  exact heights (#54). The PNG decides at the finished figure's scale, as for the band text.

## [0.7.0] - 2026-11-19

### Added

- `min_font_size` style option, `--min-font-size PT` on the command line and "Smallest text" in
  the web app: no text is drawn smaller than this many points (#51). Notes, annotations, the
  text in bands and quadrants, "NA", the "PGS" captions and the legend are drawn at fractions of
  `font_size` (down to about half of it in the legend), so a figure meeting a journal's 6 pt
  rule needed much larger ID labels; with `min_font_size: 6` it meets it whatever `font_size`
  is. The layout measures the lines below a label at the size they are drawn, and the legend
  grows with its text. Without the option nothing changes.
- `PedigreeStyle.text_size(scale)`: the size of a text drawn at `scale` times `font_size`, never
  below `min_font_size`; both renderers take every text size from it.
- `band_text` style option: `"auto"` (default), `"name_value"`, `"value"` or `"none"` (#55).

### Changed

- Text in phenotype bands and quadrants no longer runs outside the symbol (#55). A gradient band
  wrote "name: value" whatever the room, so in the narrow top and bottom bands of a circle, with
  long column names, or in a figure drawn small for a journal, the text crossed the outline.
  Now each text is measured (DejaVu Sans Bold widths) against the width of the symbol at its
  height: "name: value" if it fits, else the value alone, else the value smaller, down to
  `min_font_size` or half `font_size`, else no text. The same goes for quadrant text and "NA".
  The PNG chooses when the figure is drawn, at its final scale (also into your own `ax`, or
  after a resize), the SVG at its own scale. `band_text: "name_value"` restores the previous
  text. Only the golden image with three-digit values in circle quadrants changed.

## [0.6.2] - 2026-11-19

### Added

- `--output-dir DIR`: one file per family in a single run, named after the family (`FAM1.png`,
  `FAM1.svg` with `--interactive`), in every format and with `--circular`, `--pheno`,
  `--annotations`, `--style` and `--family` (#49). Every file has the colors, column modes and
  legend of the whole file, so the figures of a cohort can be compared. Family IDs are made
  safe for file names (other characters become `_`, Windows device names get a `_`); two IDs
  giving the same name, ignoring case, get `_2`, `_3`... with a warning. A family that cannot
  be drawn is reported with its ID, the others are still written, and the exit code is 1.
- `render_families(pedigree, output_dir, ...)` does the same from Python and returns the files
  written; `PedigreeBatchError` lists the files written and the families that failed.
  `family_file_names(ids)` gives the file names.
- `Pedigree.subset(family_ids)`: a pedigree with only some families, keeping the phenotype and
  annotation columns of the whole file (#50). It keeps the whole pedigree as `source`, and
  `Pedigree.column_scope` and `Pedigree.individuals()` are new.

### Changed

- Several PGS columns get distinct marker colors when `pgs_colors` does not set them (#53).
  They were all red (`#e74c3c`), so the scores could not be told apart. A score named after a
  known phenotype (`PGS_autism`, `PGS_ADHD`, `PGS_ID`, `PGS_epilepsy`) takes that phenotype's
  color, the one its column has in the figure if there is one; the others take the first
  color unlike the other scores and the phenotype bands, starting from the red. A figure with
  one PGS column is unchanged.

### Fixed

- `--family` drew the family without its phenotype columns, PGS curves, annotations and
  legend, from column 6 only; so did the web app's family choice (#50). Both now keep the
  file's columns, so a family's colors and legend are those of the whole file.
- A family drawn alone (`--family`, the web app, `subset`) now has the column modes, category
  colors and carrier legend of the whole file. They were decided from the family's own
  values, while the parser normalises values over the whole file: an IQ column whose family
  had only the values 1 and 2 was drawn as a binary phenotype (#49).

## [0.6.1] - 2026-11-18

### Changed

- Columns without a color in `phenotype_colors` now get one from the palette even when other
  columns have one. They were drawn grey (`#999999`), so two such columns, two `o_` outlines
  say, could not be told apart.

### Fixed

- Two columns no longer get colors that look alike. With the four known phenotypes (autism,
  intellectual disability, ADHD, epilepsy) in a figure, the palette had only green and teal
  left for the other columns: a second `o_` outline was teal next to the first one's green,
  and a third repeated the green. When some colors are fixed (known phenotypes, or chosen by
  the user), a column whose palette color repeats one of them or looks like one takes the
  next color that does not, from the palette and three new colors (`EXTRA_PHENOTYPE_COLORS`:
  brown, light magenta, olive); past eight columns these come before the palette repeats.
  Figures with no known phenotype and no chosen color keep their colors, so the gallery and
  test images are unchanged.

## [0.6.0] - 2026-11-18

### Added

- Symbols are spaced by their labels' width, so long IDs do not run into each other (#41).
  Each label (ID, notes and annotations, these in their smaller type) is measured with the
  character widths of DejaVu Sans, at the scale of the SVG, where text takes the most room;
  where two neighbours' labels would overlap, those two are moved apart, by the room they
  need and no more (partners, siblings and neighbouring branches alike). Labels that fit
  change nothing: at the default 8 pt, IDs up to about 13 characters keep the usual gaps,
  and no gallery or test image moves. With `show_labels=False` (`--no-labels`) labels take
  no room. `compute_layout` takes `labels=True`; the circular layout is unchanged.

### Fixed

- The SVG's `viewBox` now includes labels wider than their symbol, so the outermost labels
  are no longer cut at the edge.

- `PedigreeBuilder`: a child given only a father or only a mother was drawn alone, joined to
  no one. The builder now gives it an unlabelled phantom for the other parent, as the parser
  has done since 0.4.4, so the child hangs from the couple (#40). `build()` and the parser
  now share one function for phantom parents.

### Deprecated

- `PedigreeStyle.min_node_distance` was never used by the layout (`sibling_gap` and
  `couple_gap` set the spacing). It is still accepted, and documented as having no effect;
  it will be removed in 1.0.

## [0.5.3] - 2026-11-17

### Fixed

- Input files (`.ped`, `.fam`, PLINK `.ped`, `--pheno`, `--annotations`) and style JSON files
  were read with the system's default encoding, cp1252 on most Windows systems, so non-ASCII
  IDs and phenotype names were garbled there. They are now read as UTF-8 on every platform; a
  byte order mark (Excel's "CSV UTF-8") is ignored, and a file that is not UTF-8 raises
  `PedigreeParseError` asking to save it as UTF-8 instead of a bare `UnicodeDecodeError`.
  `export_fam` writes UTF-8 with `\n` line endings on every platform.

### Changed

- CI runs the test suite on Windows too (Python 3.10–3.14). `.gitattributes` checks text files
  out with LF line endings everywhere.
- The notation images of `docs/pedigree_standard.md` (`docs/imgs/*.png`: sexes, affected,
  deceased, proband, relationship lines, twins, and the labelled diagram of the four kinds of
  line) are now drawn by pypedigree, from invented individuals, by
  `docs/imgs/make_notation_images.py`. They were cropped from a published nomenclature figure.
- The large synthetic fixture `complex/synesthesia.ped` is now `complex/large_multi_union.ped`,
  family `FAM3`, with the columns `b_ASD`, `b_Dyslexia` and `o_Sleep disorder` (were
  `FAM_SYN01`, `b_Synesthesia` and `o_Psychiatric comorbidities`). The data are unchanged and
  were already generated; the old labels suggested a study family. Its reference images, the
  gallery, the README and `pedigree_standard.md` follow.

## [0.5.2] - 2026-11-17

### Fixed

- PNG and PDF rendering grew with the square of the pedigree's size: matplotlib re-fitted the
  view to every artist after each one was added (a PGS curve's `fill_betweenx` triggered it),
  and the renderer counted the axes' artists before drawing every node, which only a rotated
  node needs. The view is now fitted once, at the end, and the artists counted only for
  rotated nodes. The figure is pixel-identical. On a synthetic pedigree with two phenotypes
  and a PGS, a PNG takes 5 s instead of 17 s for 1,000 people, 15 s instead of 52 s for 2,000,
  and 47 s instead of 264 s for 5,000.

## [0.5.1] - 2026-11-17

### Added

- Web app: a link can open the page on an example, by its file name without the extension
  (`?example=quadrants`). Clicking an example puts it in the address bar; opening a file of
  one's own takes it out.

### Changed

- Circular layout: an ID label no longer always goes radially outward, where the lines of the
  node's own descendants and the next ring often cross it. It takes the first direction, turning
  from outward to either side, where its text touches no line, shape or other label, and keeps
  a small gap in points from its shape. The labels are placed when the figure is drawn, at its
  final size, so they stay clear when a figure is resized after rendering. On the test
  fixtures at the default size, labels touching something drop from 15 of 53 to 0
  (`synesthesia.ped`) and from 49 of 79 to 2 (`multigenerational.ped`).

### Fixed

- Web app: the download buttons showed before anything was drawn, and stayed after an error,
  because a `display` rule overrode their `hidden` attribute.
- Web app: tab-separated files lined up badly in the text box when a column name was longer
  than the tab width (8). Tabs are now 12 characters wide. The "Four phenotypes in quadrants"
  example no longer forces the quadrant layout, which is the default since 0.5.0.

## [0.5.0] - 2026-11-17

### Changed

- **Phenotypes are drawn in quadrants by default.** `phenotype_layout` has a new value,
  `"auto"`, which is the default (and `--phenotype-layout auto` on the command line). It cuts
  each symbol in four, one quadrant per phenotype column or `cr_` carrier source, as
  `"quadrants"` does, up to four sources. Beyond four, it draws horizontal bands without a
  warning; `"quadrants"` still warns when it falls back. Figures with one to four phenotype
  columns therefore change. To keep the previous figures, set `"phenotype_layout": "bands"` in
  the style or pass `--phenotype-layout bands`. Pedigrees with no phenotype column, or with
  carrier sources only, are drawn as before. The web app offers the same three choices, with
  "auto" first.

### Fixed

- Warnings raised while drawing, such as `--phenotype-layout quadrants` falling back to bands
  with more than four phenotype columns, printed in the CLI as raw Python warnings with an
  internal path and a source line. Like those raised while parsing, they are now reported once
  as `Warning: <message>` on stderr, for PNG, PDF, circular and SVG output.

## [0.4.7] - 2026-11-17

### Added

- A web app, in `web/`: pypedigree runs in the browser through Pyodide, installed without its
  dependencies (parsing, layout and SVG need none). The user pastes, drops or opens a pedigree file,
  or picks one of eight synthetic examples. The family, bands or quadrants, the PGS scale, the
  carrier notation, labels, legend and any style option can be chosen. The page shows the
  interactive SVG and saves it, the same SVG without script, or a PNG made in the browser at 150,
  300 or 600 dpi. The file never leaves the computer. A help page has a tutorial and explains how
  to read the drawing. `web/build.py` assembles the site with the wheel of the checkout and the
  HTML documentation; `.github/workflows/pages.yml` publishes it on GitHub Pages for each release
  tag. `tests/test_web.py` tests the Python side, with matplotlib, numpy and click blocked.

### Fixed

- Drawn small, as in a journal figure (text large against the symbols), the proband's arrow
  shrank to its head and the "P" covered it. The arrow and the "P" are now sized from the label
  font (the arrow spans 1.4 em in x and in y; the "P" ends 0.2 em left of the arrow's start),
  in the PNG and the SVG, so the marker keeps its shape at any scale. At the default font a
  large pedigree's arrow is as long as before; a small pedigree drawn large gets a shorter
  arrow, in proportion to its text.
- With PGS curves, the "PGS" caption under the proband's (or consultand's) curve ran into the
  arrow and the "P". It is now written above that curve.

## [0.4.6] - 2026-11-17

### Fixed

- A loop between the descendants of two siblings (cousins marrying) was drawn with a
  duplicated individual when a third sibling married into a family placed later: that sibling
  was ordered between the two branches, which no longer met. Both branches now go just before
  the siblings who marry out, so the cousins are drawn side by side. On the test generators
  this removes one duplicate in every 100 (1,851 to 1,833 over 95 pedigrees); no fixture
  changes.

## [0.4.5] - 2026-11-17

### Fixed

- With `phenotype_layout: "quadrants"` and carrier sources but no phenotype column (a `carrier`
  column or `cr_` columns alone), the legend showed a quadrant key with the carrier in one
  quadrant, while the symbol was hatched whole. Nothing is cut into quadrants then: the symbol
  and the legend are drawn as with bands.
- With tight spacing or large text, the proband's "P" ran into the label of the individual on
  the left (for example `sibling_gap` 1.0 with `font_size` 18). The layout now keeps that
  neighbour, or the partner on the left, far enough for the "P" and their label
  (`symbols.marker_clearance`, using the longest label of the generation); a consultand's arrow
  gets the same room. With the default style no pedigree moves. Both renderers take the marker's
  geometry from `symbols.marker_points`.

### Added

- Zenodo DOI: [10.5281/zenodo.24289471](https://doi.org/10.5281/zenodo.24289471) for all versions
  (0.4.4: 10.5281/zenodo.24289472), in `CITATION.cff` and the README.

## [0.4.4] - 2026-11-16

### Added

- PGS on a fixed scale: `--pgs-scale zscore` / `percentile` (style
  `pgs_scale`) places each PGS value on the drawn standard normal, the same
  for every pedigree, instead of between the pedigree's lowest and highest
  score (`minmax`, still the default). The ±2 SD tails are marked, values
  beyond the box are clipped to its edge, the legend names the scale, and
  the SVG tooltip gives z and percentile ("z = +1.10, 86th percentile")
  (bourgeron-lab/pypedigree#7).
- Individuals with a row but never sampled: a `not_sampled` column (aliases
  `inferred`, `unsampled`, `dummy`; `1`, `yes` or `true`) or the builder's
  `inferred=True` draws "dummy" individuals, who only link real members, with
  a dashed outline like a phantom parent, while keeping their label, parents,
  phenotypes and their place in a `.fam` export; the SVG tooltip says "Not
  sampled (inferred)" (bourgeron-lab/pypedigree#6).
- Documented evaluation, records reviewed: an `evaluated` column (aliases
  `documented_evaluation`, `evaluation`) or the builder's `evaluated=` draws
  an asterisk beside the lower right of the symbol, as in Bennett et al.
  2022, Figure 2. Any text other than `1`/`yes`/`true` is the evaluation's
  result, written below the symbol without the "E" prefix the 2022 revision
  dropped (bourgeron-lab/pypedigree#11).
- Interactive SVG: each node carries `data-x`, `data-y` (the centre of its
  symbol) and `data-size`, in the coordinates of `#main-group`, so a page
  embedding the SVG can draw overlays (observed relatedness, links the
  pedigree does not show) in the same coordinate system through zoom and pan
  (bourgeron-lab/pypedigree#12).
- Hooks for host pages: a `css_class` column and `data_` columns
  (`data_qc` → `data-qc`), or the builder's `css_class=` and `data=`, add
  classes and `data-*` attributes to the individual's node in the interactive
  SVG, so a page can style or select individuals (dim deselected members,
  outline QC flags) without re-rendering. Duplicated individuals carry them
  too; the PNG ignores them. Names are validated (`PedigreeParseError` with
  the line, `PedigreeError` from the builder); `node`, `dimmed` and
  pypedigree's own `data-*` names are reserved (bourgeron-lab/pypedigree#9).
- A family holding only two parents, a man and a woman, with no child in the
  file (a trio whose child is missing) is drawn as their couple, father on the
  left, with a descent line to a "?", instead of two unrelated symbols, in
  the PNG, the SVG and the circular layout. The layout adds a placeholder
  child (`with_unknown_offspring`, `is_unknown_offspring`); the parsed
  pedigree and an exported `.fam` are unchanged.
- Twins from code: `PedigreeBuilder.add_individual(..., twin_of="S7",
  twin_type="MZ")` (a list for triplets) and `PedigreeBuilder.add_twins(["S7",
  "S8"], "MZ")`, recorded both ways as the `twins` column does. Types `MZ`,
  `DZ`, `TW` (default, zygosity unknown); an unknown type or a group of fewer
  than two raises `PedigreeError`, and a twin never added, twins in different
  families or with different parents give a `PedigreeDataWarning`
  (bourgeron-lab/pypedigree#10).
- Quadrant layout for phenotypes: `--phenotype-layout quadrants` (style
  `phenotype_layout: "quadrants"`) cuts each symbol in four, one quadrant per
  phenotype column or `cr_` carrier source, in the PNG, the interactive SVG and
  the circular layout. Autism, ADHD, intellectual disability and epilepsy take
  the top-left, top-right, bottom-left and bottom-right quadrants whatever the
  column order; other columns the free ones in order; `phenotype_quadrants`
  places any column. A quadrant shows the same codes as a band (color, white,
  grey, grey with "NA"); a gradient value is written alone. The legend starts
  with a key, a square and a circle cut in four with the column names around
  them. With more than four columns the figure falls back to bands with a
  `PedigreeStyleWarning`; an unknown layout or quadrant raises `PedigreeError`.
  Horizontal bands stay the default. Gallery example:
  `tests/data/phenotype/quadrants.ped` (synthetic).
- SVG for web pages and printed reports: `render_svg()` (and
  `PedigreeBuilder.render_svg()`, CLI `--interactive --no-script`) writes the
  same drawing as `render_interactive_svg()` without script or tooltip group;
  each individual keeps its `<title>` and `data-*` attributes. Both take
  `id_prefix=` (CLI `--id-prefix`), prepended to every id and every
  `url(#...)` / `href="#..."` reference, so several figures can share a page,
  and `responsive=True` (CLI `--responsive`), which writes `width="100%"`
  and no height, keeping the `viewBox`. An invalid prefix raises `ValueError`
  (bourgeron-lab/pypedigree#4).

### Changed

- `import pypedigree` no longer imports matplotlib (nor numpy):
  `render_pedigree` and `render_family` are imported on first use (PEP 562),
  and the CLI imports matplotlib only for PNG/PDF output. Parsing, layout and
  the interactive SVG never load it, so a web server avoids its import time
  and `pyplot`'s global state. Colors are read by a new `pypedigree.colors`
  module, with the same values as `matplotlib.colors.to_rgb` for hex, CSS4
  names, `tab:` colors, single letters, grey levels and RGB tuples; `xkcd:`
  names and `C0` shorthands are still accepted through matplotlib. The PNG
  and SVG outputs are unchanged (bourgeron-lab/pypedigree#8).
- Interactive SVG, inlined in a page: the root `<svg>` has the class
  `pypedigree`, and its CSS applies only under it (no more page-wide `.node`,
  `.edge` or `text` rules). The script acts only on its own root
  (`document.currentScript.closest('svg')`), not the page's first SVG, and
  listens on that root only (pointer capture for panning, no listener on the
  window). Opened on its own the wheel still zooms; inlined in a page only
  Ctrl/Cmd + wheel (or a trackpad pinch) zooms, so the page keeps scrolling
  (bourgeron-lab/pypedigree#4).

- Conventional colors for four phenotypes, recognised by column name (case,
  accents, prefix and aliases such as `ASD`, `TSA`, `ID`, `DI`, `TDAH`
  ignored): autism blue `#386eb2`, intellectual disability red `#d34436`, ADHD
  yellow `#e59635`, epilepsy purple `#5e287e`, in both layouts. They apply to
  any such column without a color in `phenotype_colors`, also when the user
  colored other columns (those used to fall back to grey). The palette skips
  its lookalike colors for the other columns, so figures with an autism, ADHD,
  ID or epilepsy column change colors; figures without one are unchanged.

### Fixed

- PGS curves, and the legend, were drawn only when the pedigree also had
  another phenotype column: a pedigree with PGS columns alone showed none
  (bourgeron-lab/pypedigree#7).
- A parent named in the father or mother column but with no row of their
  own left the child unconnected: the child was drawn alone, joined to
  neither parent. The parser now does what `PedigreeBuilder` does and makes
  that parent a phantom with the given ID (dashed outline, labelled, sex from
  the role), so the child hangs from both parents. The `PedigreeDataWarning`
  stays, once per missing parent, and says so: "father 'X' has no row in
  family 'F'; drawn as a phantom (dashed) parent". A `.fam` export still
  leaves the phantom out (bourgeron-lab/pypedigree#5).

## [0.4.3] - 2026-11-16

### Added

- `PedigreeStyle.show_na_labels` (default `True`) and the `--no-na-labels`
  CLI flag: with the option off, bands with no data (`-9`, `NA`, empty),
  `cr_` bands with no data and the centred "NA" of a symbol with no data at
  all are drawn grey without the "NA" text, in the PNG and the interactive
  SVG. They then look like an "unknown" (`0`) band.

## [0.4.2] - 2026-11-16

First release of the 0.4 series on PyPI. 0.4.0 and 0.4.1 were tagged on
GitHub but withdrawn before publication, so **0.4.2 ships everything listed
under [0.4.1] and [0.4.0]** plus the changes below. In short, compared
with 0.3.0:

- Standard pedigree notation (Bennett et al. 2008 and the 2022 revision):
  pregnancies and pregnancy losses, adoption, childlessness and infertility,
  gamete donors and surrogates, the consultand, and gender identity.
- Carrier status in the NSGC 2022 notation by default (`cr_` columns, one
  gene or condition each), with the legacy notation still available.
- A rewritten hierarchical layout: sibship lines never overlap, loops are
  drawn side by side where possible and with duplicated partners otherwise,
  every line is horizontal or vertical, and twins keep their spouses outside
  the twin connector.
- A circular layout built on the hierarchical one, with upright shapes.
- An interactive SVG that draws the same figure as the PNG, with tooltips,
  highlighting, zoom and pan.
- The PNG, its legend and the SVG take their phenotype decisions from one
  place, so they can no longer disagree; many rendering fixes follow.
- Python 3.10 to 3.14, typed (`py.typed`, `mypy --strict`), with warning
  categories callers can filter.

Figures change compared with 0.3.0 wherever a pedigree carries phenotypes;
see [0.4.0] for the breaking changes.

### Changed

- A missing value (`-9`, `NA`, `.`, empty) is drawn **grey**, the grey of
  `0` ("unknown"), with "NA" written in the band to tell the two apart; it
  was a white band. The same holds for a `cr_` carrier band with no data and
  for the single centred "NA" of a symbol whose bands all lack data, and a
  pregnancy-loss triangle whose first band has no data is grey. The "NA" text
  is darker (`#555555`, contrast 5:1 on the grey) so it stays legible.
  `NA_FILL` and `NA_TEXT_COLOR` in `pypedigree.phenotypes` hold the two
  colors. Figures with missing values change: the visual-regression goldens
  and reference images that contain them were regenerated.

### Documentation

- Every example now uses invented individuals and values, as the test
  fixtures have since 0.4.0: the Annotations section of the HTML/PDF
  documentation, the compact-header example in the README and
  `pedigree_standard.md`, and one row of the SVG tests. The fixture
  generator labels annotations with `CHD8`, `ARID1B`, `ADNP` and `SYNGAP1`,
  and `phenotype/annotations.ped` and `phenotype/test_all.ped` were
  regenerated with them.
- The project lives at <https://github.com/bourgeron-lab/pypedigree>;
  project URLs, the README badge, `CITATION.cff` and the documentation point
  there.
- `complex/synesthesia.ped` is generated with a random shape: founders 1
  and 2, thirteen children of whom one or two stay single, and 1 to 4
  children for each couple (55 individuals, 12 couples). It no longer keeps
  a fixed generation shape, founder IDs or proband.

## [0.4.1] - 2026-11-16 [NOT PUBLISHED]

Tagged on GitHub but withdrawn before it was published to PyPI; its
changes were first released in [0.4.2].

### Fixed

- A binary, gradient or unprefixed phenotype column with no value in the
  whole pedigree (every individual `-9`, `NA` or empty) drew grey "unknown"
  bands, the meaning of `0`, instead of white "NA" bands, the meaning of
  `-9`. It is now drawn like any other missing value. An individual is
  therefore drawn with bands even when no column holds a value for anyone
  (for example a `.pheno` file whose sample IDs match nobody in the
  pedigree): every band is "NA", where 0.4.0 fell back to the plain symbol
  filled from the affection column.
- A phenotype file (`--pheno`, `apply_phenotypes_to_pedigree`) whose sample
  IDs match no individual in the pedigree was applied silently. Unmatched IDs
  are now listed in a `PedigreeDataWarning` (`Warning: …` on the CLI), with
  their number and the first five, and the warning says when no ID matched
  at all, since every band is then drawn as "NA".
- The phenotype file of the `pheno_file` visual test used sample IDs 1-4
  against a pedigree of `father`, `mother` and `child`, so it matched nobody
  and the golden image tested an empty phenotype file. Its IDs now match, and
  the golden shows the binary and gradient bands it was meant to cover.

## [0.4.0] - 2026-11-13 [NOT PUBLISHED]

Tagged on GitHub but withdrawn before it was published to PyPI; its
changes were first released in [0.4.2].

Rendering-correctness release. Figures produced by this version differ from
0.3.0 wherever a pedigree carries phenotypes — the changes below are fixes,
but they change output, so regenerate any figure you keep under version
control.

### Breaking

- Carrier status follows the NSGC 2022 focused revision of the pedigree
  nomenclature (Bennett et al., J Genet Couns 2022;31(6):1238-48) by
  default. The symbol is divided into vertical sections, an affected one
  first, then one per carrier result, each with its own fill pattern defined
  in the legend; obligate carriers get the same pattern at half the density.
  Figures with a `carrier` column change. The previous half-filled shape and
  center dot remain available with `carrier_notation="legacy"` (style option,
  style file key, or `--carrier-notation legacy`).
- `render_pedigree()` now draws phenotype bands and a legend on a
  multi-phenotype pedigree where 0.3.0 drew plain black-and-white shapes.
  Nothing to change in calling code; the output is simply what the CLI
  already produced.
- An invalid color in a style raises `PedigreeError` instead of surfacing
  matplotlib's `ValueError`. Code with `except ValueError` around a render
  must catch `PedigreeError` (or `PedigreeError, ValueError`) instead.
- `PedigreeStyle.from_dict()` emits `PedigreeStyleWarning` for unknown keys.
  Under `-W error` this now raises; pass `warn_unknown=False` to opt out.
- The circular layout types a single-child union's edge as
  `EdgeType.DESCENT` rather than `EdgeType.INDIVIDUAL`. Rendering is
  unaffected; only callers inspecting `EdgeType` see the difference.
- Figures change even with unchanged input: legend rows for unprefixed
  continuous columns became gradient bars, categorical colors are assigned
  differently, and CLI figures with a column named exactly `PGS` are no
  longer spaced twice as wide as intended.
- `pypedigree.__version__` is read from installed package metadata. Running
  from an uninstalled source tree now reports `0.0.0.dev0`.
- A column named `gender` next to a `sex` column is read as the gender
  identity (see Added) instead of an extra phenotype column. Alone, `gender`
  is still the sex column.

### Added

- Standard notation beyond sex, affection status and death, following the
  NSGC nomenclature (Bennett et al. 2008 and its 2022 revision), in both the
  PNG and the interactive SVG, hierarchical and circular layouts. One header
  column each, case-insensitive codes:
  - `pregnancy`: `P` ongoing pregnancy ("P" inside the symbol), `SAB`
    miscarriage, `TOP` termination (struck through) and `ECT` ectopic
    pregnancy ("ECT" below), drawn as triangles, and `SB` stillbirth (struck
    through, "SB" below).
  - `adopted`: `in` or `out`. Brackets around the symbol; the line to
    adoptive parents is dashed.
  - `childless`: `choice` or `infertility`. A line below the symbol ending in
    one bar or two.
  - `art`: `D` gamete donor or `S` surrogate, written inside the symbol.
  - `gender_identity`: `man`, `woman` or `nonbinary`. The symbol follows the
    gender identity and the `sex` column becomes the sex assigned at birth,
    written below as AMAB or AFAB when the two differ.
  - `consultand`: `1` points an arrow without "P" at the symbol.

  The same values are `PedigreeBuilder.add_individual()` keywords
  (`pregnancy=`, `adopted=`, `childless=`, `art=`, `gender=`,
  `consultand=`), and `build_notation()` turns codes into the new
  `Notation` attached to every `Individual` and `NodeLayout`. New enums:
  `Pregnancy`, `Adoption`, `Childlessness`, `ReproductiveRole`, `Gender`;
  new edge type `EdgeType.CHILDLESS`. An unknown code is a
  `PedigreeParseError` naming the line. The decisions both renderers share
  live in the new `pypedigree.symbols` module. Pedigrees without these
  columns render exactly as before. Example fixture:
  `tests/data/nomenclature/nomenclature.ped`, also a visual-regression
  golden.
- The interactive SVG (`--interactive`, `render_interactive_svg()`) draws
  the same figure as the PNG: phenotype bands (binary, gradient, diverging,
  categorical), `o_` outline colors, PGS distribution curves, `a_`
  annotations and the legend, which `--no-legend` (or `show_legend=False`)
  hides. It used to draw shapes and clinical notation only, with phenotype
  values in the tooltips. Tooltips now also list categorical values, and
  "NA" rather than "nan" for missing data.
- `pypedigree.phenotypes`: the phenotype decisions both renderers share
  (effective display mode of each column, color and text of every band,
  rows and geometry of the legend). The matplotlib renderer's output is
  unchanged, pixel for pixel; the SVG can no longer diverge from it.
- `cr_` column prefix: carrier status of one gene or condition per column
  (`0` not a carrier, `1` carrier, `2` obligate carrier, `-9`/`NA` unknown).
  Each `cr_` column is its own patterned section and legend row, so a
  carrier of two conditions shows both, as NSGC 2022 recommends. With
  phenotype bands, each carrier source is an extra band below them.
  Example fixture: `tests/data/phenotype/carrier_genes.ped`.
- `PedigreeStyle.carrier_notation` (`"nsgc2022"` or `"legacy"`) and the
  matching `--carrier-notation` CLI option. Any other value raises
  `PedigreeError`.
- `resolve_style()`, `apply_pgs_spacing()` and `is_pgs_column()` are exported
  from the package for callers that need to know how a pedigree will be styled.
- `tests/test_phenotype_resolution.py`: regression tests for style resolution,
  per-column mode detection and categorical color assignment.
- `tests/test_cli.py::TestCLIWarnings`: coverage for CLI warning formatting.

- Warning categories so callers can filter pypedigree diagnostics
  independently of other libraries: `PedigreeWarning` (base),
  `PedigreeStyleWarning` and `PedigreeDataWarning`. For example
  `warnings.simplefilter("error", PedigreeDataWarning)` turns a
  sex-inconsistent parent into a hard failure.
- `PedigreeStyle.from_dict()` now reports unknown style options as a
  `PedigreeStyleWarning`, with a suggested spelling when one is close
  (`phenotype_color` -> "did you mean 'phenotype_colors'?"). Silently
  ignoring a typo was indistinguishable from the option having no effect.
  Pass `warn_unknown=False` to accept a dictionary carrying extra data.
- Invalid colors now raise `PedigreeError` naming the style option they came
  from, instead of surfacing a raw matplotlib `ValueError: Invalid RGBA
  argument`. Every color in the style is checked before any drawing starts,
  so all offending values are reported at once.
- `tests/test_layout_circular.py`: 169 tests for the circular layout, which
  previously had none — ring geometry, node spacing, rotation, edge
  composition, arc waypoints and rendering across 15 fixtures.
- `tests/test_legend.py`: legend panel geometry and row composition.

- `CITATION.cff` so the project is citable, with GitHub's "Cite this
  repository" button. Authors: Alexandre Mathieu, Freddy Cliquet.
- `.github/workflows/release.yml`: tag-driven PyPI publishing via trusted
  publishing (OIDC, no stored API token). The workflow refuses to build when
  the tag disagrees with the version in `pyproject.toml`, runs the test suite,
  checks metadata with `twine check --strict`, re-runs the tests from inside
  the extracted sdist, then publishes and creates the GitHub release.
- `docs/RELEASING.md`: release procedure and the one-time PyPI
  trusted-publishing setup (a manual step on pypi.org).
- `Alexandre Mathieu` added to the package authors, matching the commit history.

- `mypy --strict` is enforced in CI. The package ships `py.typed`, so the
  annotations are a public contract; they were never verified before. Fixing
  the 27 findings turned up one real annotation bug: `_detect_separator()` was
  annotated `-> str` but returns `None` to signal whitespace splitting.
- `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md` (Contributor Covenant 2.1), issue
  forms for bugs and feature requests, and a pull-request template.
- `tests/test_svg_renderer.py::TestSvgPhenotypes`: checks that every SVG band
  has the color the matplotlib renderer gives it, and that the legend, PGS
  curves, outline colors and annotations are drawn.

### Changed

- Fewer duplicated individuals in pedigrees with loops. Siblings are now
  ordered so the two branches of a loop meet, at every generation from the
  common ancestors down (previously only the partner themselves moved to
  the end of their sibship), and when lines still overlap, one loop per
  group of overlapping lines is broken before laying out again (previously
  every loop touching an overlapping line, at once). Over 24 synthetic
  clinical pedigrees with one cousin marriage, 1 needs a duplicate instead
  of 13; with three cousin marriages each, 35 duplicates instead of 57;
  12 % fewer over pedigrees with random mating. Layout takes about 6 ms for
  200 individuals with six loops. The reordering can put a sibling out of
  birth order. `TestLoopDuplicates` guards these counts.
- Circular layout rebuilt on the hierarchical one: each individual's
  horizontal position becomes an angle (counter-clockwise, centred at the
  bottom) and each generation a ring. It inherits the hierarchical
  guarantees whatever the number of founder couples: the sector model it
  replaces put duplicates among the founders, overlapping other nodes at the
  centre, and drew `consanguinity_first_cousins.ped` as a column instead of
  a circle. A small family now occupies only the arc it needs.
- Circular layout: shapes stay upright. A square turned to face outward read
  as the unknown-sex diamond, and bands changed direction around the circle.
  ID labels are set radially outside the shape, away from incoming lines
  (`NodeLayout.label_angle`).
- Hierarchical layout: horizontal placement rewritten. Each union's
  descendants now form a block; sibling blocks are packed against each
  other generation by generation, and every couple is centred over its
  children. The sweep, crossing-reduction and refinement passes, and the
  15-individual threshold that switched them on, are gone: small and large
  pedigrees follow the same rules. Figures change for 17 of the 40 test
  fixtures, and large tree pedigrees come out about a third narrower.
- Hierarchical layout: a loop (both partners descend from people in the
  pedigree, as in consanguinity) is drawn with the partners side by side,
  the second partner's sibship line stretching to reach them. Where that
  line would overlap another sibship, the partner is drawn twice instead —
  among their siblings and, as a duplicate joined by a dotted arc, next to
  their partner — as pedigree nomenclature allows. None of the test fixtures
  needs this; densely interbred pedigrees do (305 duplicates for 1,800
  individuals in the test suite's random-mating generator).
  `compute_layout()` returns those duplicates as extra nodes; the input
  family is not modified.

- Visual-regression golden images regenerated. Besides the descent-line fix,
  four multi-phenotype goldens still dated from March and showed the legend
  and gradient bugs fixed in this release; they passed only because the
  comparison tolerated a 2% mean pixel difference.
- Visual regression tests count changed pixels instead of averaging the
  difference over the image. A pixel counts as changed when a channel moves
  by more than 0.1, and a render fails above 0.01% changed pixels or on any
  size change (sizes up to 10% off used to be cropped and compared). Every
  past regression checked had passed the old 2% mean threshold, the smallest
  moving the mean by 0.0003; the new criterion catches them all, and a
  single missing proband arrow (175 pixels). Measured noise is 0 pixels,
  across platforms and on matplotlib 3.7 to 3.10. A failure saves the
  actual render and a diff image with the changed pixels in red.
- The PGS golden test no longer widens the gaps itself: `render_pedigree`
  does it, so the golden encoded doubled spacing.
- Ten new golden images: carrier notation (NSGC 2022, `cr_` columns,
  legacy), the circular layout (rings, twins, duplicates, bands), separate
  sibships, a second partnership and the uncle-niece loop.

- `numpy` is now a declared dependency (`numpy>=1.24`). The renderer imports
  it directly but relied on matplotlib to pull it in. The HTML documentation
  already listed it as a requirement.
- CI now also runs the suite on Python 3.10 with every dependency at its
  declared minimum version (matplotlib 3.7, click 8.0, numpy 1.24, pandas 1.5),
  so the lower bounds in `pyproject.toml` are tested rather than assumed.
- The `python -m pypedigree` test runs the current interpreter instead of
  `uv run`. Going through `uv run` re-synced the project's `.venv` to the
  lockfile in the middle of the test session, which removed pandas and pytest
  for every later test when the environment held other versions, and made the
  suite require uv.

- Reference images under `tests/data/` and `docs/imgs/examples/` regenerated
  for the legend and spacing fixes above. Several of them (`annotations.png`,
  `multi_phenotype_inline.png`, `test_all.png`, `large_family.png`,
  `real_world/multi_phenotype_pgs.png`, `synesthesia.png`) were already stale
  from an earlier legend revision and are now current.

- `__version__` is read from the installed package metadata instead of being
  duplicated in `__init__.py`, so `pyproject.toml` is the single source.
- The circular layout types the edge of a single-child union as
  `EdgeType.DESCENT` instead of `EdgeType.INDIVIDUAL`, matching what the
  hierarchical layout emits for the same structure. Rendering is unaffected
  (the renderer does not branch on edge type), but callers filtering on
  `EdgeType.DESCENT` now see parent-child links in both layouts.
- Sex-inconsistency warnings now carry `PedigreeDataWarning` rather than a
  bare `UserWarning`.
- Internal: the diverging-gradient color computation was duplicated between
  the bands and the legend bar; both now call `_diverging_color_at()`, so the
  two cannot drift apart.
- Internal: `_draw_legend()` went from 472 lines to about 90, with the panel
  geometry in a `_LegendMetrics` dataclass and one function per row type.
  Verified output-identical across 159 renders.

- Source distribution shrank from 4.5 MB to 328 kB. It carried the 2 MB
  documentation PDF, all reference and example images, `CLAUDE.md` and
  `uv.lock`. It now ships source, tests, the visual-regression golden images
  and the Markdown docs — enough to build and run the full suite, verified in
  CI on every release.
- `Development Status` classifier raised from `3 - Alpha` to `4 - Beta`.
- Project URLs, README badge and documentation now point at
  `github.com/bourgeron-lab/pypedigree`.

- **`requires-python` lowered from `>=3.12` to `>=3.10`.** The 3.12 floor was
  arbitrary: the full suite passes on 3.10 and 3.11, and the only 3.10-specific
  feature in the code is `zip(strict=True)`. This matters for bioinformatics
  clusters, many of which still ship 3.10 or 3.11. The documentation already
  claimed 3.10+, so `pyproject.toml` was the outlier.
- CI now tests Python 3.10, 3.11, 3.12, 3.13 and 3.14 on Linux and macOS (10
  jobs, up from 4). 3.14 is the version pinned in `.python-version` for
  development and had never been exercised in CI.
- `ruff` target version lowered to `py310` to match.

### Fixed

- The test suite passes from the source distribution with the newest
  matplotlib. matplotlib 3.11 places text about a pixel lower than 3.7-3.10,
  which the golden images were made with, so the visual comparison is now
  skipped outside that range instead of failing; and the check that every
  generated fixture has its image is skipped where the images are not
  shipped. A new CI job runs the suite with the newest release of every
  dependency, since `uv.lock` pins older ones.
- A duplicated individual (someone with several partners, or a partner drawn
  twice to break a loop) is drawn with the original's phenotype bands. It
  was a plain symbol, so in a figure with bands it lost every phenotype.
- The interactive SVG's default size keeps the drawing's proportions. Width
  and height were clamped separately (to 600-1600 and 400-1200 px), so a
  wide pedigree was letterboxed with white bands above and below. Small
  pedigrees are still at least 600 px wide.
- Circular layout: twin connectors and duplicate-individual arcs are drawn,
  as in the hierarchical layout. A union whose only children are one twin
  group has no sibship arc, only the twins' V.
- Band text (values, "NA") on a rotated shape is never drawn upside down.
- An affected carrier lost its affected fill: the shape was drawn white
  before the carrier half-fill. With the NSGC 2022 notation both are shown;
  the legacy notation is unchanged.
- The interactive SVG was drawn upside down, founders at the bottom: the
  layout's y axis points up and SVG's points down. Coordinates are now
  mirrored before rendering, duplicate arcs included.
- Interactive SVG of a file with several families: every family was drawn at
  the same place, on top of the others. They are now side by side, each
  titled, as in the PNG.
- Interactive SVG: element ids were built from sample IDs, so two families
  sharing an ID produced duplicate ids (one clip path served both), and IDs
  with spaces made invalid references. Ids are now generated.
- Interactive SVG: the proband arrow had no head (its marker was never
  defined); the tooltip background was a large black disc (corner radius 3
  layout units); the tooltip ignored zoom and pan and could leave the
  visible area.
- Interactive SVG: text was positioned with `dominant-baseline`, which
  librsvg and older Inkscape ignore; positions are now explicit offsets.
- Hierarchical layout: sibship lines of different couples overlapped, so the
  children of two couples hung from one continuous line and could not be
  told apart. It happened in one synthetic tree pedigree in ten between 8
  and 15 individuals, and in almost every pedigree of 50 or more (1,280
  overlapping pairs over 12 pedigrees of 100-200). Sibship lines no longer
  overlap; `TestSibshipSeparation` checks it on every fixture and on
  generated pedigrees up to 200 individuals, with and without loops.
- Hierarchical layout: the duplicate drawn for a second partnership stood
  on the founders' row whatever the generation of the individual it
  repeats. In `multiple_partners_complex.ped` the father appeared in
  generations I and II and his third child one generation too high.
- Hierarchical layout: a child pulled to a lower generation by their
  partner (`consanguinity_uncle_niece.ped`) was connected by a line that
  stopped one row above them. It now reaches their row.
- Legend panel sized from the font instead of the figure. It took a fixed
  share of the figure height, its swatches were sized on the figure width
  and its text in points: a wide pedigree got huge swatches and tiny text
  (about 3 pt once the 45-individual test pedigree is scaled to a printed
  page), a circular one a panel stretched over the whole square, a legend
  with four rows shrank its text to 4 pt, and raising `font_size` made the
  column captions overlap. One legend unit is now a fixed multiple of
  `font_size`: text, swatches and rows keep their proportions and the same
  physical size in every figure, the border hugs the entries instead of
  spanning the figure, and the figure grows downward to make room rather
  than squeezing the pedigree.
- Duplicate arcs bow by at most 0.4 generation gaps. Their curvature grew
  with their length, so an arc across a wide pedigree swept over several
  generations and, through matplotlib's autoscaling, shrank the whole
  drawing inside a mostly empty figure.
- Hierarchical layout: a single child that could not be centred under its
  parents was joined to them by a diagonal line. Pedigree nomenclature uses
  horizontal and vertical lines only; the child is now reached through a
  short horizontal segment at sibship level, as siblings are. This affected
  10 of the test fixtures, including `three_gen.ped` and both consanguinity
  examples. Node positions are unchanged — only the connecting lines are.
- Hierarchical layout: when a twin group was a union's only offspring and its
  convergence point was off the parents' midpoint, the descent line stopped
  beside it, leaving a gap. The same horizontal segment now joins them.
- Hierarchical layout: a married twin's partner could be drawn between the
  two twins (the side came from the order of the partnership, e.g. the wife
  of the left twin went to his right), so the twin connector spanned across
  the spouse. The partner of the leftmost twin of a group is now drawn on
  its left, and that of the rightmost twin on its right; a middle triplet,
  a twin's second partner, and a partner from a consanguineous loop keep the
  previous placement. `TestTwinSpouses` checks that no one is drawn between
  twins, on dedicated cases and on every fixture. No fixture or golden image
  changes.
- Header columns are matched by name, as documented, even when the header
  has no extra phenotype, annotation or notation column. Such headers were
  read by position, so a structural column outside the standard order was
  misread without a warning: a `Twins` column right after `Phenotype` was
  taken as the date of birth and the twins were lost, a `Deceased` column
  with no `Twins` column before it was parsed as twins, and reordered core
  columns (`IID FID Sex ...`) swapped sex, parents and IDs. Headers in the
  standard order parse exactly as before (every fixture is unchanged).
- Circular layout: ID labels were hidden by their own shape (13 of the 45
  individuals of the large test pedigree), and the proband "P" was drawn
  sideways with an arrow pointing at an unrotated corner, because shapes were
  rotated to face outward. Shapes now stay upright and labels sit radially
  outside them (see Changed). For nodes an API caller rotates, the label
  still moves below the rotated outline, with a white backing, and the
  proband marker stays upright, pointing at the rotated outline.

- Categorical phenotype colors (`c_` columns) are now resolved once per figure
  instead of being accumulated in a module-level cache. This fixes three defects:
  - Two categories could be assigned the same color when the column's own
    `phenotype_colors` entry also appeared in the categorical palette (the two
    palettes share six colors). In `prefix_modes.ped`, `ASD` and `ADHD` were
    both rendered green.
  - Categories leaked between figures: rendering family A then family B in one
    process listed A's categories in B's legend.
  - A second render with a different `phenotype_colors` silently reused the
    colors of the first, so changing a color in a notebook had no effect.
  - Category colors now also avoid the colors of other phenotype columns, so
    two bands on the same individual stay visually distinct.
- `render_pedigree()` now produces the same figure as the CLI. Palette
  auto-assignment lived only in `cli.py` and `PedigreeBuilder`, so
  `render_pedigree(parse_pedigree("multi_pheno.ped"))` silently rendered plain
  shapes with no bands and no legend. The logic moved to
  `resolve_style()` in `config.py`, which every entry point now shares.
- PGS spacing is no longer applied twice. Both the CLI and `render_pedigree()`
  widened `couple_gap`/`sibling_gap`, so figures from the CLI with a column
  named exactly `PGS` were spaced twice as far apart as intended.
  `apply_pgs_spacing()` is now applied by `render_pedigree()` only.
- `PedigreeBuilder.render()` no longer mutates the style passed to
  `set_style()`, so rendering the same builder twice no longer compounds
  PGS spacing or leaves auto-assigned colors behind.
- The legend now matches the bands for AUTO (unprefixed) phenotype columns.
  The legend had no AUTO branch and fell through to the binary
  Affected/Unaffected/Unknown swatches, while the bands were drawn as
  gradients. An unprefixed continuous column such as `IQ` showed gradient
  bands under a binary legend.
- The effective display mode of a column is decided once from all of its
  values instead of per individual, so a column can no longer mix binary and
  gradient bands across the same diagram (e.g. values `0,1,2,5,17`).
- The legend no longer reserves rows for `PGS`-prefixed columns as if they
  were bands, which left an oversized empty legend box.
- The Affected/Unaffected/Unknown legend header is only drawn when at least
  one band is actually binary.
- CLI validation warnings are reported as `Warning: <message>` on stderr
  instead of raw Python `UserWarning` output carrying internal file paths and
  an unrelated source line.
- Sex-consistency warnings are reported once per parent instead of once per
  child of that parent.
- `docs/examples_gallery.md`: `--style` takes a path to a JSON file; the
  inline-JSON example did not work.

### Removed

- `tests/visual/golden/multi_phenotype_inline.actual.png`: a visual-diff
  failure artifact committed by accident, sitting next to the golden images it
  could be mistaken for. `*.actual.png` is now gitignored.

- `NodeLayout.phenotype_display_modes`: a field neither layout ever populated
  and nothing ever read. The live field is `Pedigree.phenotype_display_modes`,
  which the renderer uses. Removing the dead twin stops it misleading readers.

### Documented

- A larger gallery example: `tests/data/complex/multigenerational.ped`, a
  synthetic 81-person, five-generation pedigree with a first-cousin
  marriage and a marriage between two branches, both drawn side by side
  without duplicates, plus bands, carriers, twins and standard notation. It
  is written by `generate_fixtures.py` and has a visual-regression golden.
  It replaces the 100-individual gallery image, which predated the current
  layout (crossing lines, unconnected spouses).
- The examples gallery no longer describes `synesthesia.ped` as data from a
  real study: like every fixture, it is synthetic.

- The interactive SVG does not support the circular layout. Stated in the
  README, the API reference and the HTML/PDF documentation.

### Security

- **Test fixtures are now fully synthetic.** Five fixtures previously held
  data derived from real participants: pedigree topologies plus diagnoses, IQ,
  SRS scores, polygenic scores and gene-level findings. De-identification alone
  was not enough — an accession code is a pseudonym rather than anonymisation,
  and a large pedigree carrying a rare phenotype can be identifying on its own
  — so the fixtures were regenerated from scratch:
  - `tests/data/generate_fixtures.py` builds them from a fixed seed, so they are
    reproducible and demonstrably fabricated.
  - Feature coverage is unchanged: the same column sets, display-mode prefixes,
    `-9` no-data values, comma decimal separators, generation shape (2 / 20 / 23
    individuals) and union count. No test assertion needed changing.
  - The git history was rewritten with `git filter-repo` so the two cohort
    accession codes appear in no commit, message or path.
  - `tests/data/real_world/` is renamed `tests/data/complex/`; the old name
    described provenance that no longer exists and would have misled the next
    reader into handling the files as real data.
  - `tests/data/README.md` records the provenance rule, and
    `tests/test_fixtures.py` enforces it: it fails on any accession-shaped
    identifier, any populated date-of-birth column, any family ID outside the
    synthetic set, and on a generated fixture that has drifted from the
    generator output.
- `CONTRIBUTING.md` and all three issue/PR templates now require test data to
  be synthetic or fully de-identified.

### Known limitations

- The interactive SVG has no circular layout: `--circular` with
  `--interactive` warns and draws the hierarchical layout.
- Gamete donors and surrogates are marked "D" and "S", but the diagonal
  lines that join them to a pregnancy in the NSGC figures are not drawn: a
  donor is a parent in the file and is joined to the other parent like any
  couple. Not covered at all: several individuals in one symbol ("n"),
  documented evaluation ("*"), the asymptomatic-carrier vertical line, and
  childlessness of a couple (the marker hangs from one individual).
- A pregnancy-loss triangle has no room for phenotype bands; in a figure
  with bands it takes the color of the individual's first band.

## [0.3.0] - 2026-04-03

### Added

- Annotation support:
  - Columns with `a_` prefix in .ped header are annotation columns (text displayed below individual's ID label)
  - Values rendered as italic grey text; missing values (0, -9, NA, empty) produce no annotation line
  - Generation gap auto-increases when annotations are present
  - External annotation file via `--annotations PATH` CLI option
  - `parse_annotation_file(path)` and `apply_annotations_to_pedigree()` API functions
  - `Individual.annotations: dict[str, str]` field on the data model
  - `NodeLayout.annotations: list[str]` field on the layout model
  - `Pedigree.annotation_names: list[str]` field
  - `PedigreeBuilder.add_individual(..., annotations=["line1", "line2"])` parameter
- Gradient inverted display mode:
  - New `gi_` prefix for `GRADIENT_INVERTED` mode
  - Color-to-white scale (high value = white, low value = full color)
  - Useful for metrics where low = bad (e.g., IQ: low IQ = more red)
- Multiple PGS columns:
  - Any column starting with `PGS` (case-insensitive) is now treated as a PGS column (e.g. `PGS_autism`, `PGS_IQ`, `PGS_ADHD`)
  - Single gaussian curve with multiple colored markers
  - `pgs_colors` style option for per-column marker colors
  - `couple_gap` and `sibling_gap` auto-increase when PGS columns are present
- New `PedigreeStyle` fields:
  - `pgs_colors: dict[str, str]` -- PGS column name to hex color for marker line
  - `gradient_range: dict[str, list[float]]` -- Custom min/max for gradient normalization
  - `gradient_diverging: dict[str, dict]` -- Tricolor diverging gradient configuration
  - `zero_pad_ids: int` -- Zero-pad sample IDs in labels to this width (0 = no padding)
- Comma decimal separator support in parser (European locale: `0,712` parsed as `0.712`)
- Gradient band values always displayed regardless of band color intensity
- AUTO mode now uses raw values (1/2) for binary detection instead of normalized values

## [0.2.2] - 2026-03-20

### Added

- Circular/radial layout mode:
  - `compute_layout_circular(family, style)` places founders at the center and descendants on concentric rings
  - `--circular` CLI flag to activate radial layout
  - `circular=True` parameter on `render_pedigree()` and `PedigreeBuilder.render()`
  - `radial_gap` and `base_radius` style parameters in `PedigreeStyle` for tuning ring spacing
  - `waypoints` field on `EdgeLayout` for curved/polyline edge routing
  - Waypoint-based polyline drawing in the renderer for non-orthogonal connections

### Changed

- Circular layout algorithm completely rewritten: sector-based radial model with per-lineage angular sectors, auto-scaling radius, and proper pedigree structure (DESCENT + SIBSHIP arcs + INDIVIDUAL lines instead of direct star-pattern lines).
- `base_radius` default changed from 2.0 to 3.0 in `PedigreeStyle`.
- Extra phenotype value semantics: `0` = unknown (grey band), `-9` = no data (white band with "NA" text). When all bands are NA for an individual, a single centered "NA" is shown instead of repeating per band.
- Legend layout is now fully adaptive: row height, swatch sizes, font sizes, and gaps all scale based on the number of phenotype rows (no more hardcoded sizes).
- Legend text simplified: removed "Top band:" / "Middle band:" prefixes, now shows just the phenotype name.
- Legend column headers ("Affected" / "Unaffected" / "Unknown") now shown once above all phenotype rows instead of repeated per row, preventing text overlap and reducing clutter.
- Legend text positioning fully dynamic: `text_x` (swatch rows) computed from actual swatch widths, `icon_text_x` (PGS, carrier, outline rows) computed from actual icon widths — both with proportional gap that scales with row count.
- Added automated legend overlap tests (`TestLegendNoOverlap`) that verify no text bounding boxes collide.
- `synesthesia.ped` updated: uses `b_ASD`, `b_Synesthesia`, `o_Psychiatric comorbidities` prefixed columns, IDs without leading spaces, and includes Deceased and Proband columns.

## [0.2.1] - 2026-03-19

### Added

- Smart header-based column mapping:
  - Parser now recognises column names (FID, ID, Father, Mother, Sex, ASD, etc.) and maps them to their roles
  - Extra columns not matching known structural names are treated as phenotype columns
  - Works with any column order (e.g., `FID ID Father Mother Sex ASD Synesthesia`)
  - No longer requires 10 fixed-position columns before phenotypes
  - Supports column name aliases (e.g., `Pat`/`Father`/`Dad`, `Pheno`/`Phenotype`/`Affected`)
  - First extra phenotype column used as primary phenotype when no explicit `Phenotype` column present
- Synesthesia pedigree test data (45 individuals, 2 phenotypes: ASD, Synesthesia)
- Standard pedigree encoding support for extra phenotypes (1=unaffected, 2=affected, 0/-9=unknown)
- Binary phenotype rendering: color=affected, white=unaffected, grey=unknown (no text by default)
- `show_phenotype_labels` style option to show YES/NO text inside binary bands
- `--show-phenotype-labels` CLI flag
- Prefix-based phenotype display modes in column headers:
  - `b_` → binary (color/white/grey bands)
  - `g_` → gradient (white→color continuous scale with value labels)
  - `o_` → outline (colored shape border when affected)
  - `c_` → categorical (distinct color per unique string value)
  - No prefix → auto-detect from data (backward compatible)
- `PhenotypeDisplayMode` enum exported from pypedigree

## [0.2.0] - 2026-03-19

### Added

- Carrier notation support:
  - Half-filled shapes for heterozygous carriers (`is_carrier`)
  - Center dot for obligate carriers (`is_obligate_carrier`)
  - `carrier` column in header-based .ped format (values: 0/1/2)
  - `carrier` and `obligate_carrier` parameters in `PedigreeBuilder.add_individual()`
  - Auto-generated legend entries for carrier symbols
- Interactive SVG export:
  - New `svg_renderer.py` module producing standalone SVG with embedded JavaScript
  - Tooltips on hover (sample_id, sex, phenotype, DOB, extra phenotypes)
  - Click to highlight connected relatives (parents, children, siblings)
  - Mouse wheel zoom and Ctrl+drag pan for large pedigrees
  - `--interactive` CLI flag
  - `render_interactive_svg()` public API function
  - `PedigreeBuilder.render_interactive_svg()` method
  - Carrier notation rendered in SVG (clip-path half-fill and center dot)

## [0.1.0] - 2026-03-16

### Added

- Initial release
- Standard 6-column pedigree file parsing (TSV, CSV, space-separated)
- Extended format support (DOB, twins, deceased, proband)
- Multi-phenotype rendering with horizontal color bands
  - Binary phenotypes (0/1): filled or empty bands with YES/NO labels
  - Continuous phenotypes: color gradient bands with value labels
- PGS (Polygenic Score) gaussian distribution curves per individual
- Inline phenotype parsing from columns 11+ (with header row)
- Separate `.pheno` TSV file support via `--pheno` CLI option
- Auto-generated phenotype legend with proper male/female symbols
- JSON-based style customization (`--style`)
- Multiple output formats: PNG, SVG, PDF
- Interactive display mode (`--show`)
- Family filtering (`--family`)
- Twins support: monozygotic, dizygotic, unknown zygosity
- Deceased markers, proband arrows, date of birth labels
- Consanguineous relationship rendering
- Missing/inferred parent handling (dashed shapes)
- Pandas DataFrame integration (`parse_pedigree_from_dataframe`)
- CLI options: `--outline-color`, `--no-labels`, `--no-legend`, `--dpi`
- Programmatic `PedigreeBuilder` API with fluent interface and Jupyter support
- Visual regression tests with golden image comparison (`--regen-golden`)
- Large pedigree layout (100-200+ individuals):
  - O(n) generation assignment with parent-child adjacency index
  - Barycenter heuristic for edge crossing reduction
  - Bidirectional overlap resolution (left-to-right + right-to-left pullback)
  - Bottom-up parent centering over children
  - Backward-compatible: new passes only activate for >15 individuals
- Auto-scaling figure dimensions based on layout bounds
- Comprehensive HTML and PDF documentation (`docs/`)
- PLINK bioinformatics integration:
  - `parse_fam_file()` for native PLINK `.fam` reading
  - `parse_plink_ped()` for PLINK `.ped` files (ignores genotype columns)
  - `export_fam()` for exporting pedigrees to `.fam` format (PLINK/KING interop)
  - CLI `--plink` flag and auto-detection of `.fam` extension
