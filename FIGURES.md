# Figure reference

What each output figure shows, and how the numbers behind it were computed.
This document is descriptive only — it says what is plotted, not what the
results mean.

Per-scene figures are written to `outputs/<site>/<scene>/` by `02_compare.py`
(via `bgt/viz.py`); pooled figures to `outputs/_combined/` by `03_combine.py`
(via `bgt/viz_combined.py`). Sections 2-12 cover the per-scene figures, section
13 the pooled ones, and section 16 the biomass figures `04_compare_agbd.py`
writes to `outputs_agbd/<site>/`.

---

## 1. Conventions shared by every figure

**Filenames.** Every file named in this document carries a run suffix on disk:
`_q2` for the headline run, `_q20` and `_allquality` for the two looser quality
filters — so `fig01_maps.png` is `fig01_maps_q2.png` in `outputs/`. The suffix is
left off below for readability; section 15 describes the three runs.

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

Both axes span **0–50 m** (`config.SCATTER_LIMITS_M`) and the aspect is locked to
equal, so the 1:1 line sits at exactly 45°. The range is fixed for every figure
in the repository — every scene, the pooled plots and the ETH grid — so any two
can be laid side by side and read against each other without checking their axes
first.

About 0.2 % of cells carry an ALS aggregate above 50 m (the tallest is 63 m on
the BIOMASS grid, 76 m on the ETH grid) and fall outside the drawn area. They are
still in every statistic: the metrics box, both fits and every table are computed
on the whole paired sample, not on the part of it that fits inside the axes.

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

The height axis uses the same fixed 0–50 m range as the scatters, so the panels
line up. The residual axis is left free, since its useful span differs by an
order of magnitude between sites.

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
found is recorded in `summary_all_products_q2.csv`, alongside the verdict.

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

## 10. The ETH comparison figures

Written only when a clipped ETH raster is available for the site (skipped under
`--no-eth`). The ETH global canopy height product (Lang et al. 2023) is a 10 m
canopy top height regressed from Sentinel-2 against GEDI, representative of
**2020**, delivered on a geographic grid as `uint8` whole metres with 255 as
no-data.

**The three-way cell set.** These figures use a stricter subset than the rest of
the run: the paired cells that *additionally* carry a finite ETH aggregate with
at least `--min-coverage` ETH coverage. So all three sources have data in every
cell plotted, and neither product is credited for covering ground the other one
misses. The cell count is therefore slightly below the paired count used
elsewhere.

**ETH is aggregated at the nominal geolocation**, never at the offset found for
the ALS — that offset describes where the ALS sits relative to BIOMASS and says
nothing about where a Sentinel-2 derived product sits. Otherwise it goes through
exactly the same machinery as the ALS: warped onto the same 1 m fine grid, then
block-reduced to the same BIOMASS cells, with the same per-cell statistics.

| File | Shows |
|---|---|
| `fig09_biomass_vs_eth.png` | the two spaceborne products against each other — x is the ETH aggregate, y BIOMASS FH. Same construction as section 3: hexbin density, 1:1, OLS and RMA, and a metrics box. Neither is a reference here, so read it as agreement between two estimates rather than as validation. |
| `fig10_eth_vs_als.png` | ETH against the ALS reference — x is the ALS aggregate, y the ETH aggregate. Directly comparable with `fig02_scatter.png`, which is the same plot for BIOMASS, so the two can be read side by side. |
| `fig11_eth_maps.png` | the ALS reference, ETH, and their difference, drawn on the **BIOMASS** grid so it lines up cell for cell with `fig01_maps.png`. Same scaling rules as section 2. |

Two artefacts to expect. ETH's whole-metre quantisation puts a floor of about
0.29 m on any RMSE against it and shows as horizontal banding in the scatters.
And the 2020 epoch means the ETH-to-ALS gap is larger than the BIOMASS-to-ALS
gap at every site — four years at Mbalmayo, five at Loundoungou and Luki.

---

## 11. `fig12_exclusions.png` — Where the BIOMASS cells go

**Shows** a funnel accounting for every cell of the analysis window: how many
were removed by each test, and how many survived into the paired sample.

| Bar | Removed because |
|---|---|
| outside ALS footprint | the cell has no valid ALS underneath it at all. The analysis window is a rectangle; a scan rarely is, so this is geometry rather than data loss |
| BIOMASS no-data | the FH raster carries its no-data value (−9999) there |
| BIOMASS height out of range | the FH value is outside 0–100 m |
| BIOMASS quality above threshold | the quality layer exceeds `--quality-max` (default 2.0) |
| ALS coverage below minimum | some ALS is present, but less than `--min-coverage` (default 90 %) of the cell's area |
| ALS statistic unavailable | no usable ALS aggregate despite sufficient coverage — a safety net that should normally read zero |
| kept (paired) | survived every test; this is the sample every other figure uses |

**How it is computed.** The tests are applied **in the order listed**, and each
cell is attributed to the *first* test it fails. That makes the bars additive:
they sum exactly to the size of the analysis window. The run asserts this
reconciles with the pairing mask before writing anything, so the chart cannot
drift from the sample it describes.

The ALS-footprint test comes first deliberately. Most of what it removes was
never a candidate, and putting it first keeps it from inflating the other
categories.

**Labels.** Each bar is annotated with its count and, in brackets, the share of
the cells *still remaining* when that test was applied — which is what says
whether a test did a little or a lot of work. The "kept" bar instead reports its
share of the whole window. Bars that read zero are kept rather than hidden:
"this was not the cause" is a real answer.

---

## 12. `fig13_mask_map.png` — Why each BIOMASS cell was kept or dropped

**Shows** the same accounting in space: the analysis window drawn cell by cell,
each coloured by the test that removed it (or by "kept"). This is where a bare
count becomes diagnostic — an exclusion concentrated at the footprint edge means
something different from one scattered through the interior.

**How it is computed.** Exactly the codes behind figure 11, drawn on the BIOMASS
grid with nearest-neighbour rendering, so one screen block is one cell. The
legend carries each class's cell count and its share of the window.

**Colour.** The four hues carrying real exclusion reasons — blue, orange, aqua,
violet — were validated as an all-pairs set on this surface (worst CVD ΔE 9.2,
worst normal-vision ΔE 16.3). "Kept" and "outside ALS footprint" take recessive
neutrals rather than a fifth and sixth hue: they are the two states a reader is
not asked to discriminate between, and holding the palette to four hues is what
keeps the four that matter apart. Aqua sits below 3:1 against the surface, so
the legend carries counts and the same figures appear in figure 11 and in
`exclusion_counts.csv` — identity is never colour alone.

---

## 13. Pooled figures (`outputs/_combined/`)

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
| `fig11_exclusions_by_scene.png` | as section 11, one stacked bar per scene. Stacked to 100 % rather than to counts, because the analysis windows differ in size between sites and the question is what proportion of each scene is lost to what. Each bar is annotated with its absolute window size; the counts themselves are in `exclusion_counts_all.csv`. |

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

## 14. Accompanying tables

Each figure has a machine-readable counterpart in the same folder.

Per scene, in `outputs/<site>/<scene>/`:

| File | Contents |
|---|---|
| `paired_cells.csv` | one row per paired cell — the raw material of every figure |
| `metrics_three_way.csv` | the BIOMASS, ETH and ALS comparisons over the common three-way cell set |
| `metrics_by_als_stat.csv` | figure 8 |
| `metrics_by_height_bin.csv` | the height-stratified table behind figure 4 |
| `metrics_by_quality_class.csv` | figure 9 |
| `metrics_by_chm_product.csv` | the primary vs. secondary CHM cross-check |
| `exclusion_counts.csv` | figures 11 and 12 -- one row of cell counts per exclusion reason |
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
| `exclusion_counts_all.csv` | figure 11 of the pooled set -- one row per scene, absolute counts |

Metric definitions are in `bgt/metrics.py`.

---

## 15. The three runs (`_q2`, `_q20`, `_allquality`)

Every figure and table above exists three times, in the same folders,
distinguished by a suffix inserted before the extension that names the BIOMASS
quality filter the run used. Filenames carry no unsuffixed form — a figure
always says which filter produced it:

| Suffix | Run | BIOMASS quality filter |
|---|---|---|
| `_q2` | `02_compare.py` | `q ≤ 2` — the headline results |
| `_q20` | `02_compare.py --quality-max 20` | `q ≤ 20` |
| `_allquality` | `02_compare.py --quality-max none` | none |

So `fig02_scatter_q2.png`, `fig02_scatter_q20.png` and
`fig02_scatter_allquality.png` sit side by side, as do their metric tables,
`paired_cells_*.csv`, `exclusion_counts_*.csv` and the pooled
`outputs/_combined/*` files. The suffix is derived from `--quality-max`, so a
run cannot overwrite another run's outputs; `03_combine.py --suffix <s>` pools
the matching per-scene tables, so a pooled figure never mixes two runs.

Everything else about the three runs is identical — same co-registration, same
coverage rule, same ALS statistic — and no shift was applied in any of them, so
a difference between two of these figures is the quality filter and nothing
else. Read them together: the default set is the headline result, and the two
looser sets are what says whether the filter is doing real work. The axis range
is the same fixed 0–50 m in all three, so a variant's scatter is directly
comparable with its default counterpart.

---

## 16. The biomass figures (`outputs_agbd/<site>/`)

Written by `04_compare_agbd.py` via `viz.plot_relation`. These compare the ALS
**above-ground biomass density** map against a P-band observable, so the two
axes carry different units. That changes what the figure can say and how it is
built: there is no 1:1 line, no bias and no RMSE against the reference, because
nothing here is an estimate of anything else. What the figures show is the
*association* — and, in the binned median, its shape.

Both use the same construction:

* **x** — ALS above-ground biomass density per 50 m cell (t/ha), band `AGBD` of
  `data/AGBD_<site>.tif`
* **Hexagons** — count of paired cells per bin (44 across the range), light → dark
* **Orange line** — ordinary least squares fit of the observable on AGBD
* **Green line and band** — median of the observable in each 25 t/ha bin, with
  the interquartile range shaded. This is the element that shows saturation: a
  relation can hold a respectable correlation while its median flattens out
* **Box** — n, Pearson r and its R², Spearman ρ, the scatter about the fitted
  line, and how many cells fall outside the drawn y range

The y range is the 0.5th to 99.5th percentile of the observable, rounded
outward. A handful of very dark cells — open water, a clearing — otherwise span
more of the axis than the whole population does. Every statistic is computed on
the complete paired sample, and the box states how many cells the axes cut off.

**The paired sample** is every cell that carries an AGBD estimate, has at most
`--max-masked` (default 25 %) of its ALS area rejected by `mask_combined`, and
has a value in the layer being compared.

### `scatter_tomo_hv_<h>m.png` — Tomographic HV power vs biomass

* **y** — normalised tomographic HV power at *h* metres above the ground, in the
  layer's own 8-bit DN, from `data/Tomo/Tomo_hv_norm_<h>.tif`

The tomographic layers carry no geotransform: they are written on exactly the
grid of the AGBD raster, and are paired with it **by pixel index**, which the
loader enforces by rejecting any layer whose shape differs. They write a literal
0 outside the reconstructed area — the same wedge along one edge at every
height — and that 0 is read as no-data, not as zero returned power.

### `scatter_l1b_<pol>.png` — BIOMASS L1b backscatter vs biomass

* **y** — γ⁰ (or σ⁰) of one polarimetric channel of a single L1b DGM scene, in dB

The caption names the scene, its acquisition date, orbit direction and WRS grid
code. The L1b product is a detected **ground-range** image with no geotransform,
so it is geolocated from its own annotation and LUT: image row → azimuth time,
image column → ground range → slant-range time, then latitude, longitude and the
beta-to-gamma factor are interpolated from the LUT at that pair of times. Pixels
are beta-nought *amplitude*; squaring gives beta-nought power and the LUT factor
converts it.

The pixels landing in a cell are averaged **in linear power**, and only the cell
mean is turned into dB — averaging dB would bias the mean low. At about 4-5 L1b
pixels per 50 m cell the speckle in a single scene is barely reduced, which is
most of the scatter in the figure.

That grid is also saved as a raster, on the AGBD grid and CRS, for mapping or
for re-checking the geolocation with a shift test. It goes to
`data/l1b/<scene>/<site>_<pol>_db.tif`, beside the annotation and LUT it came
from: `outputs_agbd/` holds only figures and tables.

### Tables

`paired_cells_agbd_<h>m.csv` is the paired sample itself, one row per cell —
its row and column on the AGBD grid, the AGBD value, the masked fraction, the
tomographic DN, and, when the L1b half ran, the backscatter in dB and how many
L1b pixels landed in the cell. `summary_agbd_<h>m.csv` is one row per
observable: n, Pearson r, Spearman ρ, both fit coefficients and R². Both are
named after the tomographic height of the run that wrote them, so a run at a
different height cannot overwrite another's tables. `l1b_scenes.csv` lists every
L1b scene that fully covers the site, which is what `--list-scenes` prints.
