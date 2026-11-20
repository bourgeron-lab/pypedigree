# CLAUDE.md - PyPedigree Development Guide

## Project Overview

pypedigree is a Python package that generates publication-quality pedigree (family tree) diagrams
from standard pedigree files used in genetics research.

Architecture follows a strict four-stage pipeline:

```
Parser -> Data Model -> Layout Engine -> Renderer
```

Each stage is independently testable with no circular dependencies.

## Quick Reference

```bash
# Install (development mode with all dev dependencies)
uv sync --extra dev

# Install with pandas support
uv sync --extra dev --extra pandas

# Run the CLI
uv run pypedigree tests/data/basic/simple_trio.ped -o output.png

# Run as module
uv run python -m pypedigree tests/data/basic/simple_trio.ped --show

# Run all tests
uv run pytest

# Run tests with coverage
uv run pytest --cov=pypedigree --cov-report=term-missing

# Lint
uv run ruff check src/ tests/ web/

# Format
uv run ruff format src/ tests/ web/
```

## Architecture

```
src/pypedigree/
  errors.py    -> Custom exceptions + warning categories (no dependencies)
  models.py    -> Data structures (depends on: errors)
  config.py    -> PedigreeStyle dataclass + style resolution (depends on: errors)
  colors.py    -> Color parsing as matplotlib.colors.to_rgb does, without matplotlib (no dependencies;
                  hands xkcd:/C0 colors to matplotlib only if installed)
  parser.py    -> File I/O + validation + pandas (depends on: models, errors; optional: pandas)
  layout.py    -> Coordinate math, incl. circular layout (depends on: models, config, symbols, errors)
  carriers.py  -> Carrier notation decisions shared by both renderers, no matplotlib (depends on: models, errors)
  symbols.py   -> Standard notation (pregnancies, adoption, childlessness, donors, gender, consultand):
                  decisions and geometry shared by layout and both renderers, no matplotlib (depends on: models)
  phenotypes.py -> Phenotype display decisions shared by both renderers: column modes, band colors/text,
                   legend rows and geometry, text fitted to bands; no matplotlib (depends on: models,
                   config, colors, carriers, symbols, errors)
  renderer.py      -> Matplotlib drawing (depends on: models, config, layout, carriers, phenotypes, symbols)
  svg_renderer.py  -> Interactive SVG rendering (depends on: models, config, layout, carriers, phenotypes, symbols)
  api.py           -> PedigreeBuilder fluent API (depends on: models, config, renderer, svg_renderer)
  batch.py         -> render_families / family_file_names: one file per family (#49) (depends on:
                      models, config, carriers, phenotypes; renderer or svg_renderer inside the function)
  cli.py       -> Click CLI (depends on: parser, config; renderer and matplotlib only for PNG/PDF)
  __init__.py  -> Public re-exports (includes compute_layout_circular); render_pedigree and
                  render_family are imported lazily through a module __getattr__ (PEP 562)
  __main__.py  -> python -m entry point
```

## Key Design Rules

1. **layout.py never imports matplotlib.** It produces pure numeric
   coordinates in FamilyLayout dataclasses. More generally, only `renderer.py`
   imports matplotlib (and numpy) at module level: `import pypedigree`, the
   parser, the layout, `phenotypes.py` and the SVG writer must not load it
   (#8). The CLI, `PedigreeBuilder.render()` and `_repr_png_()` import it
   inside the function. Colors go through `pypedigree.colors`, whose tables
   must stay equal to matplotlib's (`TestColorsMatchMatplotlib`);
   `tests/test_light_import.py` checks in a fresh interpreter that
   matplotlib, numpy and pandas stay out of `sys.modules`, and that SVG output
   works with matplotlib blocked.

2. **Hierarchical edges are orthogonal.** Every edge `compute_layout` emits is
   horizontal or vertical, as pedigree nomenclature requires; an offset single
   child hangs from a short sibship-level segment. `TestOrthogonalDescent`
   checks this on every fixture and on the large synthetic pedigrees. (The
   circular layout uses straight and curved lines by design.)

3. **Sibship lines never overlap.** `_HierarchicalPlacer` packs each union's
   descendants as a block against its neighbours' contours, so two couples'
   children never interleave. Loops are drawn side by side unless that
   stretches a sibship line over another one; then `_place_hierarchically`
   duplicates the partner and lays out again. To make that rare,
   `_children_group` orders each sibship so a loop's branches meet: a child
   whose line leads to a partner in a family placed later goes last, one whose
   line leads to a partner already drawn goes first, and the two ends of a loop
   closing inside the sibship (cousins marrying) go together just before the
   last ones, so a sibling marrying out never parts them. `loops_to_break` breaks
   one loop per tangle (group of overlapping lines) per pass, the one
   stretching the widest line, and every tangled loop only after half of
   `_MAX_LOOP_BREAKING_PASSES`. Don't add passes that move nodes individually
   after placement (the old sweep / barycenter / refinement passes broke
   this). `TestSibshipSeparation` checks it; `TestLoopDuplicates` guards the
   duplicate count. A partner normally takes the side given by the
   partnership's order, but the end twins of a twin group take theirs on the
   group's outer side (`_twin_outer_sides`), so no spouse sits between twins
   under the twin connector; loop partners keep their forced right side.
   `TestTwinSpouses` checks it.

   A proband's or consultand's marker sits at the lower left of the symbol,
   level with the labels. `_HierarchicalPlacer.marker_room` gives each of them
   the distance their left neighbour needs (`symbols.marker_clearance`, for the
   longest label of the generation, text converted with
   `symbols.UNITS_PER_POINT`, the SVG's scale); only those needing more than
   the sibling gap are listed. It is applied while packing, not after: as a
   left margin on the contour (`_Block.add(sid, x, gen, left, right)`, via
   `_HierarchicalPlacer._add`) and as a wider couple gap when they are the
   partner on the right (`_couple_gap`). At the default style it is never
   needed, so nothing moves.

   Labels get room the same way (#41). `_label_reach` gives each node half
   its widest line of text (`symbols.label_width`: ID, then notes and
   annotations at `SMALL_TEXT_SCALE`; character widths of DejaVu Sans in
   `symbols._CHAR_EM`, checked against the font by `tests/test_label_spacing.py`,
   converted with `UNITS_PER_POINT`) plus half of `LABEL_GAP_EM`, listing only
   those reaching beyond half the sibling gap. `_add` puts the excess on both
   sides of the contour (the left side takes the larger of it and the marker
   margin) and `_couple_gap` widens a couple to the sum of both reaches. Labels
   that fit change nothing: at 8 pt no fixture moves (IDs up to ~13
   characters). `compute_layout(..., labels=show_labels)`; the circular layout
   calls `_place_hierarchically(..., labels=False)` and is unchanged. The SVG's
   `_family_extent` also widens the `viewBox` by labels wider than their symbol.
   `min_node_distance` is deprecated: the layout never used it. The marker itself is sized in ems of the
   label font, never in layout units (`symbols.MARKER_ARROW_EM`,
   `MARKER_LETTER_GAP_EM`): the PNG draws the arrow and the right-aligned "P"
   as offsets in points from the tip (`annotate(..., textcoords="offset
   points")`), the SVG with `symbols.marker_points(tip, font_size * _PT)`. In
   layout units the arrow shrank to its head and the "P" covered it when a
   pedigree was drawn small (0.4.6). A curve's "PGS" caption goes above the curve
   for a proband or consultand, whose marker takes the place below
   (`phenotypes.pgs_caption_above`). `tests/test_proband_marker.py`.

   A family of exactly two founders, a man and a woman, with no union
   (parents whose child is missing) gets a placeholder child from
   `with_unknown_offspring()` before placement (`Individual.is_unknown_offspring`),
   so the couple, father left, and its descent line are laid out as usual; both
   renderers draw that node as a "?" (`symbols.unknown_offspring_anchor`) with no
   symbol, label or tooltip. The parsed `Family` is not modified.
   `tests/test_unknown_offspring.py` checks it.

4. **renderer.py never computes positions.** It only draws from
   pre-computed FamilyLayout objects. Keep drawing linear in the number of
   artists: `_render_layout` turns matplotlib's autoscaling off while it
   draws (otherwise each `fill_betweenx` re-fits the view to every artist)
   and `_finalize_axes` fits the view once; never call `len(ax.patches)` or
   the like per node (each call walks all artists). `TestDrawingCost`. It supports waypoint-based
   polyline edges for curved/non-orthogonal connections (circular layout).
   The one exception is where a circular-layout label goes around its node,
   decided at drawing time from the text's real size (`_RadialLabels`).

5. **models.py has no I/O.** It is pure data structures.

6. **All visual parameters live in PedigreeStyle.** No magic numbers
   anywhere in layout or renderer code.

7. **Type hints on all public functions.** Google-style docstrings.

8. **Style resolution is centralized in `config.py`.** `resolve_style()` fills
   in pedigree-dependent style defaults (palette colors) and is pure and
   idempotent — it never mutates its argument. Every entry point (CLI,
   `PedigreeBuilder`, `render_pedigree`, `render_interactive_svg`) calls it, so
   the API and the CLI always produce the same figure. `apply_pgs_spacing()` is
   *not* idempotent: each renderer's entry point (`render_pedigree()`,
   `render_interactive_svg()`) applies it exactly once, and nothing else does.

9. **Phenotype display decisions are resolved once per figure, in
   `phenotypes.py`, and shared by both renderers.** `PhenotypeRender` maps each
   column to its effective mode (`binary`/`gradient`/`categorical`/`outline`/
   `carrier`/`pgs`) and, for categorical columns, to its value -> color mapping.
   `band_cell()` gives each band's color and text, `legend_rows()` the legend's
   rows. With `phenotype_layout` "auto" (the default since 0.5.0) or "quadrants",
   `with_quadrants()` stores in
   `PhenotypeRender.quadrants` which quadrant each column or carrier source
   takes (`assign_quadrants()`: `phenotype_quadrants`, then the conventional
   quadrant of autism/ADHD/ID/epilepsy, then free ones in order; more than four
   columns fall back to bands, silently with "auto", with a warning with
   "quadrants"), `quadrant_box()`/`quadrant_center()`
   give the geometry (with an `up` sign for the SVG) and `quadrant_cell()` the
   fill and shorter text. With carrier sources but no band column,
   `with_quadrants()` assigns nothing: the symbol is drawn whole, as with bands,
   and the legend has no key. Conventional colors (`KNOWN_PHENOTYPE_COLORS`,
   `known_phenotype()`) are applied in `config.resolve_style()`. The matplotlib renderer and the SVG renderer only translate these into
   artists or elements, so the PNG, its legend and the interactive SVG can never
   disagree. Put any new display decision there, not in a renderer. Nothing is
   cached in module-level state between renders.

10. **No module-level mutable state in the renderer.** Anything derived from a
   pedigree or a style is computed per call and passed down. Caching it at
   module level leaks between renders in one process (Jupyter, batch loops,
   web services) and makes output depend on render order.

11. **Diagnostics use the `Pedigree*Warning` categories** from `errors.py`, never
   a bare `UserWarning`, so callers can filter or escalate them. The CLI
   runs parsing, style resolution and rendering inside `cli._reported_warnings()`,
   which captures them with `warnings.catch_warnings` and re-emits each message
   once as `Warning: <message>` on stderr, before "Saved". `TestRenderWarnings`
   checks it for every output, on a real stderr too.

12. **Standard notation is decided in `symbols.py`.** `models.Notation` holds an
    individual's pregnancy, adoption, childlessness, assisted-reproduction role,
    gender identity and consultand flag. `symbols.py` decides what they look like
    (shape sex from gender, AMAB/AFAB and "SB"/"ECT" notes, triangle and bracket
    geometry, the letter inside) with no matplotlib, and functions whose result
    is not symmetric top to bottom take an `up` sign so the SVG (y down) uses
    them too. The layout copies the decisions into `NodeLayout` (`sex`,
    `is_deceased`, `notes`, `marker_depth`) and emits the childlessness marker as
    `EdgeType.CHILDLESS` edges and an adopted-in child's line as dashed, so both
    renderers draw them without knowing about notation.

13. **The legend panel derives its own geometry, from the font.** `legend_metrics()`
    (in `phenotypes.py`) returns a `LegendMetrics` from the row count and
    `style.font_size` (`LegendMetrics.font`: more, with `min_font_size`, so its
    row text reaches it; captions through `caption_fs()`, never below it); one
    legend unit is a fixed multiple of that font size
    (`LEG_UNIT_X_EM`, `LEG_UNIT_Y_EM`), so text and swatches scale together and
    never with the figure. Both renderers use it: `_add_legend_axes()` grows the
    matplotlib figure downward by the panel height rather than taking a share of
    the pedigree's, and `_svg_legend()` places the panel below the pedigree.
    Each row type has its own drawing function; the quadrant key row is two
    rows tall (`LegendRow.slots`, `legend_slot_count()`, `legend_row_centers()`),
    is drawn before the Affected/Unaffected/Unknown captions, and takes its
    layout from `quadrant_key_geometry()`. With `legend` "compact" (the "auto"
    default with quadrants) two or more binary columns of the key become one
    `binary_shared` row (first swatch striped with their colors, #56). The PGS
    row's `key` (score name, color) is placed by `legend_key_items()`: a dash
    then the name, measured with `symbols.text_width`, after the label. Legend geometry constants
    (`LEG_*`) are not user-tunable style. `TestLegendPhysicalSize` checks it.

14. **Every text size comes from `PedigreeStyle.text_size(scale)`** (#51): `font_size`
    times the text's fraction, never below `min_font_size`. Fractions are named
    (`symbols.SMALL_TEXT_SCALE` notes/annotations, `BAND_TEXT_SCALE`, `NA_TEXT_SCALE`,
    `PGS_CAPTION_SCALE`, `phenotypes.QUADRANT_FONT_SCALE`, `LEG_TEXT_SCALE`); never write
    `style.font_size * 0.7` in a renderer. The layout measures labels with the same sizes
    (`symbols.label_width(label, below, size, below_size)`). Gaps around text follow its drawn size
    too (#66): the lines below a label start at `symbols.below_label_offsets` (each step the larger
    of `TEXT_LINE_STEP` and `TEXT_LINE_HEIGHT_EM` of the line above + `TEXT_LINE_GAP_EM` of its
    own), the "PGS" caption at `phenotypes.pgs_caption_offset` (the larger of 0.15 × half and half
    a dot + `PGS_CAPTION_GAP_EM`); the SVG converts with `_PT`, the PNG draws `renderer._OffsetText`,
    which moves itself at draw time at the figure's scale, and the layout's `_room_below_labels`
    uses the same offsets. At default sizes the fixed distances win, so nothing moves. `check_min_font_size()` runs
    with the other style checks. `tests/test_min_font_size.py` checks that no text of the
    PNG or SVG is below the minimum, on fixtures with every kind of text.

## File Format

See `docs/pedigree_standard.md` for the complete specification.

Standard 6-column format (TSV/CSV/space-separated):

| Col | Name | Values |
|-----|------|--------|
| 1 | family_id | Family identifier |
| 2 | sample_id | Individual identifier |
| 3 | paternal_id | Father's sample_id (0/NaN = unknown) |
| 4 | maternal_id | Mother's sample_id (0/NaN = unknown) |
| 5 | sex | 1=male, 2=female, 0/-9=unknown |
| 6 | phenotype | 1=unaffected, 2=affected, 0/-9=unknown |

Parent IDs with no row become phantoms keeping that ID (`_create_phantom_parents`,
which returns one warning message per missing parent); a child with one parent given
gets a `_phantom_*` parent with no label. `PedigreeBuilder.build()` calls the same function.

To keep some families, use `Pedigree.subset(ids)` (the CLI's `--family`, the web app,
`render_families`): it keeps `phenotype_names`, `phenotype_display_modes` and `annotation_names`,
and the whole pedigree as `source`. `Pedigree(families=...)` drops them, and the family is drawn
without bands, PGS curves or legend (#50). Both renderers resolve columns (`resolve_phenotype_render`:
AUTO modes, category colors; `resolve_carrier_scheme`) over `pedigree.column_scope` (`source`
or itself), since the parser normalised values over the whole file: a subset is drawn as in the
whole figure (#49). Decide anything else that depends on all values the same way.

`render_families` (`--output-dir`) checks the style once, then draws `subset([fid])` per family
into `family_file_names()` (safe characters, Windows device names, case-insensitive collisions
→ `_2` with a `PedigreeDataWarning`). A failing family is recorded and the others drawn; then
`PedigreeBatchError(written, failures)`; the CLI prints `Saved:` per file, `Error: family '<id>'`
per failure, exit 1. matplotlib figures are closed per family (`_save_figure`), only those it
opened. `tests/test_batch.py`.

Extended format adds:
- Column 7: date_of_birth (YYYYMMDD, YYYYMM, or YYYY)
- Column 8: twins (comma-separated TYPE:sample_id entries)
- Column 9: deceased (1=deceased, 0=alive)
- Column 10: proband (1=proband, 0=not proband)
- `css_class` column and `data_` prefixed columns: classes and `data-*` attributes on the
  individual's SVG node (`Individual.css_classes`, `.data_attributes`; validated by
  `parse_css_classes()` / `data_attribute_name()`; `node`, `dimmed` and pypedigree's own
  `data-*` names are reserved). Empty/./NA left out, 0 and -9 kept. Not drawn in the PNG.
- `not_sampled` column (aliases `inferred`, `unsampled`, `dummy`; 1/yes/true), by header name:
  `Individual.is_inferred`, a row never sampled (a "dummy"). Drawn dashed like a phantom
  (`NodeLayout.is_inferred`) but keeps label, phenotypes and `.fam` export; unlike
  `is_phantom`, it does not suppress bands. `tests/test_inferred.py`
- Columns with `a_` prefix: annotation columns (text displayed below individual's ID label, not phenotypes). Missing values (0, -9, NA, empty) produce no annotation line. Generation gap auto-increases when present.
- Standard notation columns, by header name (`_NOTATION_ROLES` in parser.py; codes
  case-insensitive, unknown code -> `PedigreeParseError` with the line):
  `pregnancy` (P/SAB/TOP/ECT/SB), `adopted` (in/out), `childless` (choice/infertility),
  `art` (D/S), `evaluated` (1/yes/true, or a result text written below; asterisk at the
  lower right, `symbols.evaluation_mark_position`), `gender_identity` (man/woman/nonbinary; `gender` means this only when
  a `sex` column is also present, otherwise it is the sex column), `consultand` (1).
  `build_notation()` accepts the same codes (used by `PedigreeBuilder`).
- Columns 11+: extra phenotypes (requires header row with column names)

Smart header-based column mapping:
- When a header is present, columns are matched by name (not position)
- Any column after the core 5 (FID, ID, Father, Mother, Sex) that doesn't match
  a known structural name is treated as a phenotype column
- Supports aliases: Pat/Father/Dad, Pheno/Phenotype/Affected, etc.
- `ColumnMapping` dataclass in parser.py handles the name→role mapping
- `_parse_individual_mapped()` uses ColumnMapping instead of fixed indices.
  `_parse_text` uses it whenever the header has extra phenotype, annotation or
  notation columns, or a structural column outside the standard order
  (`ColumnMapping.is_positional`); only a header in the standard order goes
  through the positional `_parse_individual()`

Multi-phenotype support:
- Inline phenotypes parsed from columns 11+ when header is present
- Separate .pheno TSV file via `--pheno` CLI option
- Binary phenotypes: 1=unaffected (white), 2=affected (color), 0=unknown (grey), -9=no data (grey with "NA" text: `NA_FILL`, the grey of 0, and `NA_TEXT_COLOR`; a `cr_` band with no data too; `show_na_labels=False` / `--no-na-labels` drops the text, leaving plain grey)
- When all bands are NA for an individual, a single centered "NA" is shown instead of repeating per band
- Continuous phenotypes → color gradient bands
- PGS scale: `pgs_scale` ("minmax" default, "zscore", "percentile") decides the position on the curve:
  `active_pgs(node, names, scale)`, `pgs_position(z) = 0.5 + z * PGS_CURVE_SIGMA` on fixed scales,
  `pgs_tail_positions()` (±2 SD dotted lines), `pgs_value_text()` for the SVG tooltip. Curves (and the
  legend) are drawn even with no other phenotype column. `tests/test_pgs_scale.py`
- PGS markers (#54): dots closer than `PGS_DOT_PT` are set side by side (`phenotypes.spread_dots`,
  `PGS_DOT_SPREAD` diameters apart); with `show_pgs_values` each value (`pgs_marker_text`:
  "+2.1"/U+2212 minus, percentile, 3 significant digits) is written left of the box, right-aligned,
  in the marker's color with a white halo (`PGS_VALUE_HALO_PT`; SVG `paint-order="stroke"`, PNG
  `patheffects.withStroke`), spread by `spread_labels` within the box height. The SVG decides at `_PT`
  (`_svg_pgs_values`); the PNG draws `renderer._PgsMarks`, an artist placing its dots and texts at
  draw time (its `get_children()` lets `findobj` see them). `pgs_values_width` (widest value +
  gap + `PGS_BOX_PAD`, at the SVG's scale) is passed to `apply_pgs_spacing(..., values_width)` and
  to the SVG viewBox. `tests/test_pgs_values.py`
- Several `o_` borders (#46): `phenotypes.outline_colors` gives one color per `o_` column the
  individual is affected in (column order); the first is the shape's outline (`OUTLINE_WIDTH_SCALE`
  × `shape_line_width`), the others lie `outline_ring_offsets` points further out (one border's
  width + `OUTLINE_RING_GAP_PT`), grown with `symbols.grown_half` (the diamond's vertices move
  √2/`DIAMOND_SCALE` times further). `with_outline_rings` stores the most borders one drawn
  individual has in `PhenotypeRender.outline_rings` (warns above `MAX_QUIET_OUTLINE_RINGS`);
  each renderer then applies `config.apply_outline_spacing(style, outline_reach(...))` once: both
  gaps widen by twice the reach and `label_offset` by the reach, for every symbol, so labels stay
  level. Per node, the PGS box, proband arrow, deceased line, brackets, asterisk and outside
  letter move by that node's own reach. The SVG draws `class="outline-ring"` shapes at `_PT`;
  the PNG draws `renderer._OutlineRings`, an artist sizing its patches in points at draw time
  (rotated with the node in the circular layout), and lifts the first border above the bands
  (zorder 5) so it is drawn whole, as wide as the others. Legend row `"outline_rings"`
  (`LEG_OUTLINE_RINGS_LABEL`), tooltip "Borders: …". Nothing changes when nobody has two.
  `tests/test_outline_rings.py`
- PGS columns: any column starting with `PGS` (case-insensitive) → gaussian distribution curve. Multiple PGS columns (e.g. `PGS_autism`, `PGS_IQ`) show colored markers on a single curve. Marker colors via `pgs_colors` style option; unset ones are filled by `config._assign_pgs_colors` when there are several PGS columns (#53): a known phenotype's color (its column's in the figure if any), else the first of `DEFAULT_PGS_COLOR` + palette + extras not `_alike` any other score or band. One PGS column keeps `DEFAULT_PGS_COLOR` (now in `config.py`). `couple_gap` and `sibling_gap` auto-increase when PGS present.
- Colors from `phenotype_colors` in style JSON, or auto-assigned from palette
  (`config.resolve_style`): by position among all columns, lookalikes of known
  phenotypes skipped; when any color is fixed (known phenotype or user's), a
  color equal or `_alike` (CIE76 < 20, or an `_ALIKE_PAIRS` pair such as green
  and teal) to one in use is replaced by the first distinct one of the palette +
  `EXTRA_PHENOTYPE_COLORS`. With none fixed the palette order is kept exactly (no
  figure changes). Columns are never left grey since 0.6.1.
- Comma decimal separator: parser accepts `0,712` as `0.712` (European locale)
- Gradient band values always written, fitted to their place (#55, `band_text`, default "auto"):
  `phenotypes.cell_texts` gives the candidates ("name: value", then `BandCell.short`, the value),
  `band_room` / `quadrant_room` the width inside the symbol at the text's height
  (`symbols.shape_half_width_at`), and `fit_text` picks the first that fits, else the last
  smaller down to `min_font_size` or `FIT_FLOOR_SCALE` x `font_size`, else None. Widths in
  DejaVu Sans Bold (`symbols._CHAR_EM_BOLD`, `text_width(bold=True)`). The SVG calls it with
  `_PT`; the PNG draws `renderer._FittedText`, a `Text` that fits itself in `draw` and
  `get_window_extent` from the axes' transform (works with a caller's `ax`, a resize, a tight
  bbox). Same for quadrant text and every "NA". `tests/test_band_text.py` checks no SVG band
  text leaves its symbol or quadrant.
- AUTO mode (no prefix) is resolved **once per column** from all of its values:
  any string value -> categorical, any numeric value outside 0/1/2 -> gradient,
  otherwise binary. Never per individual, so a column cannot mix band types.

Additional PedigreeStyle fields:
- `pgs_colors: dict[str, str]` — PGS column name → hex color for marker line (default: `#e74c3c` for one PGS column, distinct colors for several)
- `gradient_range: dict[str, list[float]]` — custom min/max for gradient normalization (e.g. `{"SRS": [40, 90]}`)
- `gradient_diverging: dict[str, dict]` — tricolor diverging gradient config (e.g. `{"IQ": {"min": 40, "mid": 100, "max": 160, "color_low": "#e74c3c", "color_high": "#2ecc71"}}`)
- `zero_pad_ids: int` — zero-pad IDs to this width (0 = no padding)
- `min_font_size: float | None` — smallest size of any text in points (None = no minimum, #51)
- `band_text: str` — "auto" (fit, default), "name_value", "value" or "none" (#55)
- `show_pgs_values: bool` — write PGS values beside their markers (#54)
- `phenotype_labels: dict[str, str]` — names shown for columns (#56); everything shown goes through
  `phenotypes.display_name()` (label, else `_` as spaces) or `pgs_trait_name()` (the PGS key), never
  the raw column name: legend, quadrant key, band text, SVG tooltips, `carriers.label(..., shown)`
- `legend: str` — "auto" (compact with quadrants), "compact" or "full" (#56); `check_legend()`
Warnings and errors:
- `PedigreeWarning` (base), `PedigreeStyleWarning`, `PedigreeDataWarning`
- Unknown keys in a style JSON warn with a suggested spelling (they are still ignored)
- Invalid colors raise `PedigreeError` naming the option, checked before any drawing

- Prefix-based display modes in column names:
  - `b_` → binary (color/white/grey/NA, no text by default)
  - `g_` → gradient (white→color continuous scale)
  - `gi_` → gradient inverted (color→white scale; high value = white, low = full color; for metrics where low = bad)
  - `o_` → outline (colored shape border, not a band; one border per `o_` column the individual is affected in, inner first)
  - `c_` → categorical (distinct color per unique string value)
  - `cr_` → carrier status of one gene or condition (0/1/2), drawn as a patterned section, not a band
  - No prefix → auto-detect from data

## CLI Options

```bash
pypedigree INPUT.ped [OPTIONS]

-o, --output PATH        Output file path
--format [png|svg|pdf]   Output format (auto-detected from extension)
--family ID              Render only this family ID
--output-dir DIR         One file per family in DIR (not with -o/--show)
--dpi INT                Output resolution (default: 150)
--show                   Display interactively instead of saving
--style PATH             JSON file with style overrides
--separator [tab|comma|space]  Force field separator
--pheno PATH             Separate .pheno TSV file with extra phenotype data
--annotations PATH       Separate TSV file with annotation text (displayed below IDs)
--outline-color COLOR    Shape/deceased outline color (e.g. 'grey', '#888888')
--carrier-notation [nsgc2022|legacy]  Carrier notation (default: nsgc2022)
--min-font-size PT       Smallest size of any text, in points (style `min_font_size`)
--no-labels              Hide sample ID labels
--no-legend              Hide the phenotype legend
--show-phenotype-labels  Show YES/NO text inside binary phenotype bands
--no-na-labels           Do not write "NA" in bands with no data (style `show_na_labels`)
--pgs-scale [minmax|zscore|percentile]  PGS position: min-max over the pedigree (default), z-score or percentile
--show-pgs-values        Write each PGS value beside its marker (style `show_pgs_values`)
--phenotype-layout [auto|quadrants|bands]  Quadrants up to four phenotypes, bands beyond (auto, default), or one of them
--interactive            Generate interactive SVG with tooltips and click highlighting
--no-script              With --interactive: same SVG, no script or tooltip (`render_svg`)
--id-prefix PREFIX       With --interactive: prefix every id in the SVG
--responsive             With --interactive: width="100%", no height
--circular               Use circular/radial layout (founders at center, descendants on concentric rings)
--plink                  Parse as PLINK format (ignore genotype columns after col 6)
```

Circular layout: `compute_layout_circular` reuses `_place_hierarchically` and
wraps its x order around the circle (`_circular_positions`), so it inherits the
hierarchical guarantees (no interleaved sibships, loops broken with duplicates).
Shapes stay upright; `NodeLayout.label_angle` is the direction of the ID label,
radially outward. The renderer places these labels in `_RadialLabels`, an
artist that runs when the figure is drawn: text is in points and the geometry in
layout units, so whether a label fits is only known at the final scale (a figure
may be resized after `render_pedigree`). Each label tries `_RADIAL_LABEL_TURNS`
(outward, then turning to either side, never inward) and takes the first whose
text touches no line, patch outline or text of the axes, nor an earlier label
(`_label_obstacles`, `_box_hits`); otherwise the least touching one. Notes and
annotations follow it in points. `TestRadialLabels`. Twin connectors hang from the sibship arc, and
duplicate arcs are emitted as in the hierarchical layout. The renderer still
supports `NodeLayout.rotation` for API callers: rotated labels sit below the
rotated outline, band text is never upside down, and the proband marker stays
upright.

Interactive SVG: `svg_renderer.py` draws the same figure as `renderer.py`
(bands, outline colors, PGS curves, annotations, carrier notation, legend) from
the same `phenotypes.py` / `carriers.py` / `symbols.py` decisions, plus per-individual tooltips,
click highlighting and zoom/pan. Families are laid side by side; element ids are
generated (`f<family>n<node>`), never derived from sample IDs. Each node group
carries `data-id`, `data-parents/children/siblings` and `data-x`, `data-y`,
`data-size` (symbol centre and size in `#main-group` coordinates, which zoom/pan
transform) for host overlays; `TestSvgNodePositions` checks them. Text is placed
with explicit baseline offsets (`_BASELINE_SHIFT`), not `dominant-baseline`,
which librsvg and older Inkscape ignore. Font sizes and line widths in points are
scaled by `_PT` layout units per point. No circular layout in the SVG.
Embedding (#4): `_wrap_svg` gives the root `class="pypedigree"` (`ROOT_CLASS`)
and scopes every CSS rule under it; `render_svg()` is the same drawing with
`_SvgWrap(interactive=False)` (no script, no tooltip group). `id_prefix` is
applied once, by `_prefix_ids` on the wrapped markup before the script is
appended: it rewrites ` id="`, `="url(#` and ` href="#` only, which text and
titles never contain because `html.escape` turns `"` into `&quot;`. So new ids
need no prefix at their call site, as long as they and their references use
those attribute forms. The script finds its own root (`document.currentScript
.closest('svg')`, else the last `svg.pypedigree`), looks elements up with
`[id="<prefix>..."]` inside it, never on `document`/`window`, and zooms on a
plain wheel only when the SVG is the document. `tests/test_svg_embedding.py`.

Validation warnings are printed to stderr as `Warning: <message>`.

## Testing

- All tests in `tests/` directory, mirroring `src/` structure
- Test data files in `tests/data/` with subdirectories: `basic/`, `consanguinity/`,
  `extended/`, `twins/`, `multi_partner/`, `phenotype/`, `plink/`, `format/`, `complex/`,
  `nomenclature/`
- **All fixtures are synthetic.** Never add pedigrees derived from real
  individuals, even pseudonymised — accession numbers and sample barcodes are
  pseudonyms, not anonymisation. See `tests/data/README.md`.
  - `tests/data/generate_fixtures.py` writes the generated ones from a fixed
    seed: `uv run python tests/data/generate_fixtures.py`. Among them,
    `complex/multigenerational.ped` (81 people, five generations, two loops) is
    the gallery's large example and a visual-regression golden
  - `tests/test_fixtures.py` asserts the files match the generator output and
    that no fixture carries an accession-shaped identifier or a populated DOB
- `conftest.py` sets matplotlib to Agg backend (no GUI windows)
- `tests/test_layout_large.py` — stress tests for large pedigrees (100-200 individuals)
- `conftest.py` generators: `generate_clinical_pedigree` (descent tree with outside
  partners, optional cousin marriages), `generate_realistic_pedigree` (every child
  marries into another branch: worst case for loops), deep and wide pedigrees
- `tests/visual/test_visual_regression.py` — comparison against golden images by
  changed-pixel count (a channel moving by > 0.1; fails above 0.01% of pixels or on
  any size change). Renders are pixel-reproducible across platforms and matplotlib
  3.7-3.10, so any failure is a real change; a mean over the image hid regressions.
  matplotlib 3.11 places text about a pixel lower, so outside `GOLDEN_MATPLOTLIB`
  (3.7 to <3.11) the comparison is skipped, not failed
  - Golden images in `tests/visual/golden/`; a failure writes `<name>.actual.png`
    and `<name>.diff.png` (changes in red) next to the golden
  - Regenerate: `uv run pytest tests/visual/ --regen-golden`
- `tests/test_api.py` — PedigreeBuilder programmatic API tests
- `tests/test_batch.py` — One file per family: file names, `render_families` (formats, failures,
  figures closed, same legend in every file), `--output-dir`
- `tests/test_carrier.py` — Carrier notation: NSGC 2022 sections and `cr_` columns, legacy half-fill / dot (model, parsing, rendering, API, CLI)
- `tests/test_light_import.py` — No matplotlib on import, parsing, layout, SVG or `--interactive` (fresh
  interpreter, matplotlib blocked); lazy `render_pedigree`; `pypedigree.colors` equal to `matplotlib.colors`
- `tests/test_svg_renderer.py` — Interactive SVG renderer (XML validity, nodes, tooltips, interactivity, CLI;
  `TestSvgPhenotypes` checks each band's color against `band_cell`, legend, PGS, outline, annotations;
  `TestSvgFixtures` renders every fixture and checks ids are unique; multi-family placement)
- `tests/test_svg_embedding.py` — SVG inlined in a page (#4): root class, scoped CSS, `id_prefix`
  (every id and reference), script scoped to its root, `render_svg()`, `responsive`, CLI options
- `tests/test_phenotype_resolution.py` — Style resolution (purity, idempotence),
  per-column AUTO mode detection, categorical color assignment, API/CLI parity,
  color validation
- `tests/test_layout_circular.py` — Circular layout: ring geometry, spacing,
  rotation, edge composition, arc waypoints, rendering across 15 fixtures
- `tests/test_legend.py` — Legend panel metrics, physical size and row composition
- `tests/test_legend_labels.py` — Names shown (`phenotype_labels`, `_` as spaces), compact legend
  with quadrants, the PGS key (#56)
- `tests/test_outline_rings.py` — One border per `o_` column (#46): order, spacing, labels level,
  PGS box, tooltip, legend row, PNG borders in points, warning above three
- `tests/test_fixtures.py` — Fixture provenance: reproducibility and identifier scan
- `tests/test_proband_marker.py` — The proband's "P" clear of the label on its left: room kept by
  the layout with tight spacing or large text, none at the default style, no overlap in the PNG
- `tests/test_web.py` — The web app: `web/pypedigree_web.py` (options, errors, warnings, without
  matplotlib/numpy/click), every example of `web/build.py`, the site build, ids used by `app.js`
- `tests/test_unknown_offspring.py` — A couple alone in its family: when the
  placeholder child is added, father left, "?" in PNG and SVG, no node or tooltip
- `tests/test_quadrants.py` — Quadrant layout: assignment (known names, explicit, free
  order, more than four), validation, geometry, PNG/SVG/circular drawing, carrier
  quadrant, legend key, CLI and style options
- `tests/test_nomenclature.py` — Standard notation: parsing and codes, `symbols.py`
  decisions and geometry, layout (shape from gender, dashed adopted-in line,
  childless edges, circular direction), PNG and SVG drawing, builder
- CI runs the suite on Linux, macOS and Windows, Python 3.10–3.14. Read and write text
  files with an explicit `encoding="utf-8"` (inputs: `parser._read_text`, UTF-8 with an
  optional BOM): the default is the locale's, cp1252 on Windows. `.gitattributes` checks
  text out with LF everywhere, so generated fixtures match byte for byte
- CI also runs the suite with the newest release of every dependency
  (`uv sync --upgrade`, job `test-latest-versions`), since `uv.lock` pins older ones
- CI also runs the suite with every dependency at its declared minimum version
  (`uv pip install --resolution lowest-direct`) on Python 3.10; raising a lower
  bound in `pyproject.toml` is a user-visible change
- Run: `uv run pytest -v`

## Visual Reference

Reference images for every visual element are in `docs/imgs/`, drawn by
pypedigree from invented individuals with
`uv run python docs/imgs/make_notation_images.py` (fixed scale, 0.8 inch per
layout unit, cropped to the content). Rerun it after a change to how a symbol
or line is drawn. The `diagram_informations.png` shows the labeled structure of a
complete pedigree diagram (relationship line, descent line,
sibship line, individual lines).

## Documentation

- `docs/pypedigree_documentation.html` — Full HTML documentation
- `docs/pypedigree_documentation.pdf` — PDF version of the documentation
- `docs/pedigree_standard.md` — File format specification
- `docs/api.md` — API reference for all public functions, classes, and types
- `docs/examples_gallery.md` — Gallery of example files with CLI commands
- `docs/RELEASING.md` — Release procedure and the one-time PyPI trusted-publishing setup
- `docs/imgs/` — Reference images for visual elements (`make_notation_images.py`)
- `docs/imgs/examples/` — Generated example images for documentation

## Web app

`web/` holds a static page that runs pypedigree in the browser with Pyodide
(version pinned in `web/build.py`, `PYODIDE_VERSION`), installed with
`micropip.install(wheel, deps=False)`: this relies on rule 1 (no matplotlib,
numpy or click for parsing, layout and SVG).
- `pypedigree_web.py` is the only Python the page runs: `render(text,
  options_json)` returns JSON (interactive SVG with `id_prefix="pv-"` and
  `responsive=True`, static SVG, families, warnings, or an error). Keep all
  pedigree logic in the package; this file only maps form options onto
  `PedigreeStyle` (`STYLE_CHOICES` must match the `<select>` options of
  `index.html`).
- `app.js` loads Pyodide and the wheel, shows the interactive SVG in a
  sandboxed `<iframe srcdoc>` (scripts set by `innerHTML` would not run), and
  makes the PNG from the static SVG on a canvas (its size is in CSS pixels, 96
  per inch). Every `$("id")` it uses must exist in `index.html`.
- `help.html`: tutorial, file format, how to read the drawing, privacy. Update it
  with user-visible changes, as for the HTML documentation.
- `build.py` writes `web/_site` (ignored): page files, the wheel of the
  checkout (`uv build --wheel`, or `--wheel PATH`), `config.json` (version,
  wheel, Pyodide, examples), the synthetic examples listed in `EXAMPLES` and
  the HTML documentation as `docs/`.
- `.github/workflows/pages.yml` publishes it on GitHub Pages on each `v*` tag
  (a release created by `release.yml` cannot trigger a workflow) and on
  `workflow_dispatch`.
- `tests/test_web.py` checks the glue (also with matplotlib, numpy and click
  blocked), every example, the build and the page's ids. `web/` is not in the
  sdist; the tests skip there. CI runs ruff on `web/` too.

## PLINK Integration

Native support for PLINK bioinformatics formats:

- `parse_fam_file(path)` — reads `.fam` files (6-column, no header)
- `parse_plink_ped(path)` — reads PLINK `.ped` files (ignores genotype columns after col 6)
- `export_fam(pedigree, path)` — writes `.fam` format (excludes phantoms/duplicates)
- CLI: `.fam` extension auto-detected; `--plink` flag for PLINK `.ped` files

## Pandas Integration

Optional pandas support via `parse_pedigree_from_dataframe()`:

```python
import pandas as pd
from pypedigree import parse_pedigree_from_dataframe

df = pd.read_csv("family.ped", sep="\t", header=None)
ped = parse_pedigree_from_dataframe(df)
```

Install with: `pip install pypedigree[pandas]`

## Adding a New Feature

1. If it affects data representation, update `models.py` first
2. If it affects file parsing, update `parser.py` with tests
3. If it affects positioning, update `layout.py` with coordinate tests
4. If it affects drawing, put the decision in `phenotypes.py` / `carriers.py` /
   `symbols.py` and draw it in both `renderer.py` and `svg_renderer.py`, with tests
5. Update `config.py` if new style parameters are needed
6. Update CLI if new flags are needed
7. Add example `.ped` file if the feature is user-visible
