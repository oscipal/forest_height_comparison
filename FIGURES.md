# Figure reference

What each output figure shows, and how the numbers behind it were computed.
This document is descriptive only — it says what is plotted, not what the
results mean.

Per-scene figures are written to `outputs/<site>/<scene>/` by `02_compare.py`
(via `bgt/viz.py`); pooled figures to `outputs/_combined/` by `03_combine.py`
(via `bgt/viz_combined.py`). Sections 2-9 cover the per-scene figures, section
10 the pooled ones.

---

## 1. Conventions shared by every figure

### Provenance caption

Every figure carries a two-line caption in muted text at the bottom left,
identifying exactly what was compared:

```
Site: Loundoungou   |   ALS: chm_lspikefree.tif (p90 per cell), acquired 19/03/2025 to 24/03/2025
BIOMASS scene: BIO_FP_FH__L2A_20260519T043740_20260531T043807_T_G01_M03_C___T033_F003   |   FP_FH__L2A, sensed 19 May 2026 to 31 May 2026
```

| Field | Source |
|---|---|
| Site name | `site` column of `03_processed_Loundoungou/summary_processing.csv` |
| ALS product and per-cell statistic | the `--chm` and `--primary-stat` settings in force for the run |
| ALS acquisition dates | `mindt` / `maxdt` in `summary_processing.csv` — the first and last day on which points were actually recorded |
| BIOMASS scene ID | the product name, parsed from the downloaded raster filename |
| BIOMASS sensing window | the two timestamps embedded in that name (`..._L2A_<start>_<end>_...`). An L2A forest height product is derived from a stack of passes, so both the first and last acquisition are given |

### The paired sample

Figures 1–6 and 8 (and every pooled figure) draw on the same set of *paired cells*, built once per
scene before any plotting:

1. The ALS CHM is read and masked (`mask_combined.tif`; pixels outside 0–100 m
   dropped).
2. It is warped once, with area-average resampling, onto a grid that is an
   exact integer refinement of the BIOMASS grid (1 m cells, matching the native
   ALS resolution). A second warp of a validity mask gives the
   valid-ALS **coverage fraction** of every fine cell.
3. Each BIOMASS cell is reduced from its own 93×93 block of fine cells. The
   mean is coverage-weighted; percentiles and the maximum are taken over the
   fine values directly.
4. A cell is kept only where BIOMASS reports a valid height, the BIOMASS
   quality value passes the filter (`--quality-max`, default 2.0), and at
   least 90 % of the cell area carries valid ALS.

The refinement factor follows the site: at ~93 m BIOMASS cells and a 1 m
target it is 93, so each BIOMASS cell is reduced from 93 x 93 = 8649 fine
cells.

The BIOMASS raster is never resampled — the ALS is always brought up to it.
Full detail in `README.md` and `bgt/coreg.py`.

### Sign convention

**Residual = BIOMASS − ALS** throughout. Positive means BIOMASS reports the
taller canopy.

### Colour

Colour is assigned by the job it does, never decoratively:

| Job | Encoding |
|---|---|
| Magnitude (heights, densities, counts) | one blue hue, light → dark |
| Polarity (differences, residuals) | diverging blue ↔ red with a neutral grey midpoint, symmetric about zero so grey means "no difference" |
| Identity (two distributions, several candidate statistics) | fixed categorical slots — blue, orange, aqua — in a fixed order, always with a legend |

---

## 2. `fig01_maps.png` — Canopy height on the BIOMASS grid

**Shows** three side-by-side maps of the ALS footprint, all on the BIOMASS cell
grid, in the BIOMASS CRS (longitude/latitude, degrees):

| Panel | Content |
|---|---|
| Left | the ALS reference — the chosen per-cell ALS statistic (default p90) |
| Middle | the BIOMASS L2A forest height for the same cells |
| Right | their difference, BIOMASS − ALS |

**How it is computed.** The paired-cell table is scattered back into a raster of
the BIOMASS window; cells that failed any filter stay blank. The difference
panel is a straight subtraction of the two panels to its left.

**Scaling.** The two height panels share one colour scale so they can be
compared directly; its limits are the 1st and 99th percentile of the two
panels' values pooled. The difference panel uses a diverging scale pinned to
zero at the grey midpoint, with symmetric limits at ±(98th percentile of the
absolute difference). All three are drawn with nearest-neighbour rendering, so
one screen block is one BIOMASS cell and the plot introduces no smoothing.

---

## 3. `fig02_scatter.png` — BIOMASS forest height vs ALS canopy height

**Shows** the cell-by-cell relationship between the two estimates, as a density
plot rather than dots, since the cells overplot heavily.

* **x** — ALS statistic per cell (m); **y** — BIOMASS forest height (m)
* **Hexagons** — count of paired cells falling in each hexagonal bin (42 bins
  across the range, empty bins left blank), shaded light → dark
* **Grey dashed line** — 1:1, i.e. perfect agreement
* **Orange line** — ordinary least squares fit of BIOMASS on ALS
  (`scipy.stats.linregress`)
* **Green line** — reduced major axis fit, slope = sign(r)·SD(BIOMASS)/SD(ALS),
  intercept = mean(BIOMASS) − slope·mean(ALS). RMA is symmetric in the two
  variables, unlike OLS, which treats x as error-free
* **Box** — n, bias, MAE, RMSE (absolute and as a % of the ALS mean), Pearson r,
  Lin's concordance correlation, and R² against the 1:1 line

Both axes span the same padded range and the aspect is locked to equal, so the
1:1 line sits at exactly 45°.

---

## 4. `fig03_residuals.png` — Residual structure across the height range

**Shows** how the residual behaves as canopy height increases.

* **Blue dots** — one per paired cell: residual (y) against the ALS statistic (x)
* **Grey dashed line** — zero residual
* **Orange line with error bars** — the mean residual within each ALS height
  bin, with whiskers at ±1 standard deviation, plotted at the bin centre

**Binning.** Bins are the fixed edges in `config.HEIGHT_BINS`
(0, 10, 15, 20, 25, 30, 35, 40, 45, 50, 100 m), left-closed. Bins holding fewer
than 3 cells are dropped rather than plotted with meaningless spread.

---

## 5. `fig04_distributions.png` — Height distributions of the two estimates

**Shows** the two marginal distributions over the paired cells, two ways.

* **Left — histograms.** ALS in blue, BIOMASS in orange, semi-transparent and
  overlaid. Both use one shared set of 40 bin edges computed over the two
  samples pooled, so the bars are directly comparable.
* **Right — empirical cumulative distributions.** Each sample sorted ascending
  and plotted against rank/n, giving the fraction of cells at or below each
  height. The annotation box reports the mean and standard deviation of each.

Because both curves are built from the *same* cells, any horizontal gap between
them is a difference in the estimates, not a difference in sampling.

---

## 6. `fig05_bland_altman.png` — Bland–Altman agreement

**Shows** the standard agreement plot for two methods measuring the same
quantity, where neither is assumed to be truth.

* **x** — the mean of the two estimates for the cell, (ALS + BIOMASS)/2
* **y** — their difference, BIOMASS − ALS
* **Orange line** — the mean difference (the bias), labelled with its value
* **Two grey dashed lines** — the 95 % limits of agreement, at
  mean ± 1.96 × SD of the differences (sample SD, n−1), each labelled

Plotting against the *mean* of the pair rather than against either one avoids
the regression artefact that appears when a difference is plotted against one of
its own terms.

---

## 7. `fig06_shift_search.png` — Co-registration shift search

Written only when the shift search runs (omitted under `--no-shift-search`).

**Shows** how well the two datasets agree at every candidate planimetric offset,
as two maps over the same offset grid:

| Panel | Content |
|---|---|
| Left | Pearson r between aggregated ALS and BIOMASS — darker is better |
| Right | RMSE between them — colour ramp reversed so darker is again better |

Both axes are the offset of the **ALS sampling window** in metres, east and
north; the origin is the nominal geolocation.

**How it is computed.** For each candidate offset the block origin is slid
across the fine grid and the ALS cell means are recomputed, then compared
against the BIOMASS heights. Offsets are tested every `--shift-step` metres
(default 4 m) across ±`--shift-search` metres (default 100 m). The block sums
come from a summed-area table, so every offset in the window is evaluated
exhaustively rather than by local optimisation. The ALS **mean** is used here
regardless of `--primary-stat`, because it is the statistic a summed-area table
can produce.

**The comparison set is held fixed across all offsets** — cells fully covered by
ALS at zero offset, then eroded by the search radius — so that no offset can
score well merely by pulling the window into better-sampled data.

**Markers.** `+` marks the best offset found, `o` marks no shift. The `+` is red
when the shift was applied and grey when it was rejected; the annotation box
states which, and why.

A search always returns some optimum, so four checks decide whether it is
applied. The first that fails blocks the shift:

| Check | Rejects |
|---|---|
| at least `--shift-min-cells` cells (default 150) in the search set | a search over thousands of offsets on a handful of cells |
| optimum within 60 % of the search radius | an objective that runs to the rim instead of peaking inside it |
| the r-optimum and the RMSE-optimum agree within 40 m | the two criteria pointing to different parts of the window |
| correlation gain at least `--shift-min-gain` (default 0.05) | a shift that buys nothing |

`--force-shift` applies the optimum regardless. Either way the offset that was
found is recorded in `summary_all_products.csv`, alongside the verdict.

---

## 8. `fig07_als_statistic.png` — Agreement by choice of ALS aggregate

**Shows** how the comparison changes with the choice of ALS statistic, as three
bar panels — RMSE, bias, and Pearson r — with one bar per candidate: mean, p50,
p90, p95, p99, max, in that fixed order. Each bar is labelled with its value.

**How it is computed.** Every candidate is aggregated from the *same* fine grid
over the *same* paired cells and compared against the *same* BIOMASS values; the
only thing that changes between bars is which statistic summarises the ALS
canopy surface within the cell. The panels are therefore directly comparable.
Each panel's subtitle states the direction of "better" (lower RMSE, zero bias,
higher r). The underlying numbers are in `metrics_by_als_stat.csv`.

---

## 9. `fig08_quality_classes.png` — Error and sample size by BIOMASS quality class

Written only when the scene ships a quality raster.

**Shows** how the error varies with the BIOMASS quality layer value:

| Panel | Content |
|---|---|
| Left | RMSE (blue) and bias (orange) as grouped bars, one pair per quality class |
| Right | the number of paired cells in each class |

The right panel exists so the left one can be read with its sample sizes in view
— a class holding a handful of cells produces an unstable RMSE.

**Classing.** If the quality layer holds 12 or fewer distinct values they are
used as classes directly. Otherwise — which is the case here, the layer being
continuous — the values are cut into 12 equal-width bins and the bin intervals
label the axis. Classes with fewer than 3 cells are dropped.

> **Note.** This figure is computed on the paired sample, which has *already*
> been quality-filtered by `--quality-max` (default 2.0). It therefore shows
> variation *within* the retained range only. To see the full quality range,
> including the cells the filter removes, run with `--quality-max none`.

---

## 10. Pooled figures (`outputs/_combined/`)

Written by `03_combine.py`, which concatenates every per-scene
`paired_cells.csv` into one sample. Colour here encodes **identity** -- which
site a cell came from -- using fixed categorical slots assigned in order and
never cycled.

**There is no pooled map.** The scenes sit on different grids in different
countries, so there is no shared space to draw one in. Every other per-scene
figure has a pooled counterpart.

Unless stated otherwise, fits, bias lines and metric boxes are computed on the
**pooled** sample, while colour splits the points by site. So the lines describe
all sites together and the colours show how each site sits relative to them.

| File | Relation to the per-scene figure |
|---|---|
| `fig02_scatter_combined.png` | as section 3, but points coloured by site rather than binned into hexagons, since the site is the more informative encoding once scenes are mixed. The 1:1, OLS and RMA lines are pooled, drawn in text ink so they are not mistaken for a site. The box reports the pooled metrics. |
| `fig03_residuals_combined.png` | as section 4, with one binned mean ± SD line **per site** rather than one overall, so a difference in residual slope between sites is visible. Points are coloured by site at low opacity behind them. |
| `fig04_distributions_combined.png` | left panel as section 5, pooled across everything. Right panel differs: one pair of cumulative curves per site, ALS solid and BIOMASS dashed in that site's colour, so the vertical gap between a matching pair is that site's shift. |
| `fig05_bland_altman_combined.png` | as section 6, with points coloured by site; the bias line and the two 95 % limits of agreement are pooled. |
| `fig07_als_statistic_combined.png` | as section 8, recomputed on the pooled sample. |
| `fig08_quality_classes_combined.png` | as section 9, recomputed on the pooled sample. Quality values are pooled across scenes before binning. |

Two figures exist only in the pooled set, since they compare scenes to each
other rather than pooling them:

### `fig09_per_scene.png` -- Agreement per scene

Three stacked bar panels -- bias, RMSE, Pearson r -- with one bar per BIOMASS
scene, coloured by site and labelled `<site>  <sensing start date>`. Scenes are
ordered by site, then by date. Each metric is computed on that scene's paired
cells alone, so the bars are independent of one another. Each bar carries its
value; the panel titles state the direction of "better".

### `fig10_per_site.png` -- Agreement per site

The same three metrics plus the paired-cell count, with one bar per site,
computed on all of that site's scenes pooled together. The fourth panel exists
so the other three can be read with their sample sizes in view.

---

## 11. Accompanying tables

Each figure has a machine-readable counterpart in the same folder.

Per scene, in `outputs/<site>/<scene>/`:

| File | Contents |
|---|---|
| `paired_cells.csv` | one row per paired cell — the raw material of every figure |
| `summary.csv` | the headline metrics for the scene |
| `metrics_by_als_stat.csv` | figure 8 |
| `metrics_by_height_bin.csv` | the height-stratified table behind figure 4 |
| `metrics_by_quality_class.csv` | figure 9 |
| `metrics_by_chm_product.csv` | the primary vs. secondary CHM cross-check |
| `../../summary_all_products.csv` | one row per scene, every site together, including the co-registration verdict |

Pooled, in `outputs/_combined/`:

| File | Contents |
|---|---|
| `paired_cells_all.csv` | every paired cell from every scene, with `site` and `scene` columns |
| `metrics_pooled.csv` | the single pooled metric row |
| `metrics_by_site.csv` | figure 10 |
| `metrics_by_scene.csv` | figure 9 |
| `metrics_by_als_stat.csv` | pooled figure 7 |
| `metrics_by_height_bin.csv` | pooled height-stratified table |
| `metrics_by_quality_class.csv` | pooled figure 8 |

Metric definitions are in `bgt/metrics.py`.
