"""Central configuration for the ALS-vs-BIOMASS forest height comparison.

Everything that a user is likely to want to change lives here. Both entry-point
scripts (``01_download_biomass.py`` and ``02_compare.py``) read their defaults
from this module and expose them as command line overrides.
"""

from pathlib import Path

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #

ROOT = Path(__file__).resolve().parent

#: Glob matching every processed ALS site folder in the repository.
ALS_DIR_GLOB = "03_processed_*"


def site_dirs(root: Path = ROOT) -> list[Path]:
    """Every processed ALS site folder, in a stable order.

    Sites are discovered rather than listed, so dropping a new
    ``03_processed_<site>`` folder into the repository is enough to include it
    in the next run.
    """
    return sorted(p for p in root.glob(ALS_DIR_GLOB) if p.is_dir())


def site_name(als_dir: Path) -> str:
    """Short site label derived from the folder name."""
    return Path(als_dir).name.replace("03_processed_", "")


#: Default site, used only as a fallback for functions called without one.
ALS_DIR = ROOT / "03_processed_Loundoungou"

#: Where downloaded BIOMASS products are stored.
BIOMASS_DIR = ROOT / "data" / "biomass"

#: Where figures, tables and logs are written.
OUTPUT_DIR = ROOT / "outputs"

#: Sub-folder of OUTPUT_DIR holding the pooled, all-scene results.
COMBINED_DIR_NAME = "_combined"


# --------------------------------------------------------------------------- #
# ALS reference product
# --------------------------------------------------------------------------- #
#
# Product choice (see README.md for the full rationale):
#
#   chm_lspikefree.tif  -- locally adaptive spike-free CHM. The GCA pipeline
#                          documentation (_INFO_products.csv) explicitly
#                          recommends this layer "for comparative analyses
#                          across time and space" because it is robust against
#                          differences in pulse density and instrumentation.
#                          That is exactly the situation here: we compare an
#                          ALS scan against a spaceborne product.
#
#   chm_highest.tif     -- rejected. Documented as "not recommended for
#                          comparative analyses"; its height estimate is
#                          density dependent.
#   chm_tin.tif         -- rejected as the primary layer. Pit/spike-ridden and
#                          interpretable as mean interception height rather
#                          than top-of-canopy height. Kept as an optional
#                          cross-check via ALS_CHM_SECONDARY.
#
#: Primary ALS canopy height model used for the comparison.
ALS_CHM = "chm_lspikefree.tif"

#: Optional second CHM, aggregated alongside the primary one for a sensitivity
#: check. Set to ``None`` to skip.
ALS_CHM_SECONDARY = "chm_tin.tif"

#: Quality mask applied to the ALS CHM. ``mask_combined.tif`` bundles the cloud,
#: no-ground and pulse-density<2 masks, and the pipeline documentation states it
#: "should always be applied".
ALS_MASK = "mask_combined.tif"

#: Additional, stricter masks. They are documented as optionally useful but may
#: discard valid terrain, so they are off by default.
ALS_MASK_EXTRA = []  # e.g. ["mask_pd04.tif", "mask_steep.tif"]

#: Outline of the scan in WGS84, used to build the BIOMASS search box.
ALS_OUTLINE_WGS84 = "outline_WGS84.shp"

#: Physically implausible ALS canopy heights are dropped before aggregation.
ALS_HEIGHT_MIN = 0.0
ALS_HEIGHT_MAX = 100.0


# --------------------------------------------------------------------------- #
# BIOMASS product
# --------------------------------------------------------------------------- #

#: ESA MAAP STAC catalogue.
MAAP_STAC_URL = "https://catalog.maap.eo.esa.int/catalogue"

#: Keycloak endpoint that trades an offline token for a short-lived access token.
MAAP_TOKEN_URL = (
    "https://iam.maap.eo.esa.int/realms/esa-maap/protocol/openid-connect/token"
)
MAAP_CLIENT_ID = "offline-token"
MAAP_CLIENT_SECRET = "p1eL7uonXs6MDxtGbgKdPVRAmnGxHpVE"  # public, from ESA docs

#: STAC collection and product type of the BIOMASS Level-2a forest height product.
BIOMASS_COLLECTION = "BiomassLevel2a"
BIOMASS_PRODUCT_TYPE = "FP_FH__L2A"

#: Assets pulled for every matching item. ``enclosure_i_fh_tiff`` is the height
#: raster, ``enclosure_i_quality_tiff`` its per-pixel quality layer, and
#: ``enclosure_xml`` the annotation file holding the product metadata.
BIOMASS_ASSETS = ["enclosure_i_fh_tiff", "enclosure_i_quality_tiff", "enclosure_xml"]

#: Buffer (degrees) added around the ALS outline when searching the catalogue.
SEARCH_BUFFER_DEG = 0.05

#: Plausible range for BIOMASS forest height values; anything outside is treated
#: as no-data. The L2A FH product carries heights in metres.
FH_MIN = 0.0
FH_MAX = 100.0

#: Quality-layer filter: pixels are kept where ``quality <= FH_QUALITY_MAX``.
#:
#: The L2A quality layer is continuous, not a class code, and the annotation
#: files label it only as "Forest height quality" without units.
#:
#: The threshold of 2.0 was derived at Loundoungou, where the distribution is
#: sharply bimodal: a main population at or below ~1.14 and, after a clean empty
#: gap, a tail from ~8 to 100 whose cells carry a -15 to -40 m bias against ALS.
#: There, 2.0 sits inside the gap and its exact value does not matter.
#:
#: **The value distribution is not bimodal everywhere.** At Luki and Mbalmayo it
#: is continuous (median 0.3-4.1, 90th percentile 31-50), so there is no empty
#: gap for the threshold to sit in, and it removes 15-60 % of the in-footprint
#: cells against 0-3 % at Loundoungou.
#:
#: The threshold still holds up there, but on outcome rather than on the shape
#: of the histogram: at Mbalmayo the rejected cells sit under the same ALS
#: canopy height as the kept ones (30.3 vs 30.6 m) while BIOMASS reports ~12 m
#: against ~25 m, i.e. they are retrieval failures rather than short forest.
#: See the exclusion section of README.md, and fig12/fig13 per scene.
#:
#: Set to ``None`` to disable; the comparison always reports metrics stratified
#: by quality so the choice stays visible and checkable.
FH_QUALITY_MAX = 2.0


# --------------------------------------------------------------------------- #
# ETH global canopy height (Lang et al. 2023)
# --------------------------------------------------------------------------- #
#
# A second spaceborne product, for context: 10 m canopy top height regressed
# from Sentinel-2 with GEDI as the training target, representative of 2020.
# Published under CC BY 4.0; cite Lang, Jetz, Schindler & Wegner, Nature Ecology
# & Evolution 7, 1778-1789 (2023), doi:10.1038/s41559-023-02206-6.
#
# The tiles are cloud-optimised GeoTIFFs, so the download step reads only the
# window covering each ALS footprint rather than the whole ~400 MB tile.

#: Base of the ETH libdrive share holding the 3-degree COG tiles.
ETH_BASE_URL = (
    "https://libdrive.ethz.ch/index.php/s/cO8or7iOe5dT2Rt/download"
    "?path=%2F3deg_cogs&files="
)

#: Tile file name; ``{tile}`` is the south-west corner, e.g. ``N00E015``.
ETH_TILE_TEMPLATE = "ETH_GlobalCanopyHeight_10m_2020_{tile}_Map.tif"

#: The companion predictive-standard-deviation tile.
ETH_SD_TEMPLATE = "ETH_GlobalCanopyHeight_10m_2020_{tile}_Map_SD.tif"

#: Edge length of one tile, in degrees.
ETH_TILE_DEG = 3

#: Where the clipped ETH windows are stored.
ETH_DIR = ROOT / "data" / "eth"

#: Degrees added around the ALS outline when clipping.
ETH_BUFFER_DEG = 0.02

#: Where the ETH-vs-ALS results go. Deliberately a separate tree: 03_combine.py
#: pools every ``paired_cells.csv`` it finds under OUTPUT_DIR, and these pairs
#: live on a different grid, so they must not be swept into that sample.
ETH_OUTPUT_DIR = ROOT / "outputs_eth"

#: The product is delivered as uint8 with 255 as no-data, so heights are whole
#: metres in 0-254. Values above this are treated as no-data.
ETH_HEIGHT_MAX = 254.0

# --------------------------------------------------------------------------- #
# Co-registration and aggregation
# --------------------------------------------------------------------------- #

#: The BIOMASS grid is the target grid: the coarse product is never resampled.
#: The ALS CHM is first warped onto a grid that is an exact integer refinement
#: of the BIOMASS grid, then block-reduced. This is the target edge length of a
#: fine cell, in metres.
#:
#: 1 m keeps the ALS at its native resolution through the intermediate step, so
#: the per-cell percentiles are taken over the canopy surface as measured rather
#: than over pre-averaged blocks. At a ~93 m BIOMASS cell the refinement factor
#: is 93 and each cell is reduced from 93 x 93 = 8649 fine cells.
#:
#: The cost is memory: the fine grid scales with the inverse square of this
#: value, so 1 m needs about four times what 2 m did. Raise it if a run does not
#: fit; 2 m biases the ALS reference low by roughly 0.15 m (measured at
#: Loundoungou), which is small against an RMSE of several metres but is a real
#: smoothing of the extremes.
FINE_CELL_M = 1.0

#: Maximum planimetric shift tested during co-registration, in metres.
SHIFT_SEARCH_M = 100.0

#: Step of the shift search, in metres. Must be >= FINE_CELL_M; it is rounded to
#: a whole number of fine cells.
SHIFT_STEP_M = 4.0

#: Criterion optimised by the shift search: "pearson" (maximise) or "rmse"
#: (minimise). Pearson correlation is the standard choice because it is
#: insensitive to an overall height bias between the two sensors.
SHIFT_CRITERION = "pearson"

#: ALS statistic used while searching for the best shift. The mean is cheap to
#: evaluate for every candidate offset.
SHIFT_ALS_STAT = "mean"

#: Minimum gain in Pearson correlation over the unshifted comparison before a
#: shift is applied. A trustworthy co-registration produces a clear interior
#: peak; a marginal or edge-of-window "optimum" means the objective surface is
#: tracking a large-scale trend rather than a geolocation offset, and applying
#: it would fit noise. Such solutions are reported but not applied (override
#: with --force-shift).
SHIFT_MIN_GAIN = 0.05

#: A shift is applied only if its optimum lies within this fraction of the
#: search radius. A real geolocation error smaller than the search window peaks
#: near the middle of it; an objective surface that keeps improving out to the
#: rim is tracking a large-scale trend, and its "optimum" is wherever the window
#: happens to stop.
SHIFT_MAX_RIM_FRAC = 0.6

#: Minimum number of cells in the co-registration set before a shift may be
#: applied. The search tests thousands of candidate offsets, so on a small
#: sample the best of them is easily a chance fit; below this the result is
#: reported but not applied.
SHIFT_MIN_CELLS = 150

#: The correlation-based and RMSE-based optima must agree to within this
#: distance (m). When maximising r and minimising RMSE point to different parts
#: of the search window, there is no consistent alignment to be had.
SHIFT_MAX_DISAGREE_M = 40.0

#: A BIOMASS cell enters the comparison only if at least this fraction of its
#: area is covered by valid, unmasked ALS canopy height.
MIN_ALS_COVERAGE = 0.90

#: ALS statistics computed per BIOMASS cell. The BIOMASS FH product estimates
#: *top canopy height*, so the appropriate ALS analogue is not obvious a priori;
#: the comparison evaluates all of these and reports which one matches best.
ALS_STATS = ["mean", "p50", "p90", "p95", "p99", "max"]

#: The statistic used for the headline figures and tables. p90 is chosen from
#: the evidence rather than by convention: over this site it minimises both RMSE
#: and bias against BIOMASS FH among all candidates (see
#: outputs/*/metrics_by_als_stat.csv). p95 -- the usual default for "top height"
#: -- comes second.
PRIMARY_ALS_STAT = "p90"


# --------------------------------------------------------------------------- #
# Reporting
# --------------------------------------------------------------------------- #

#: Edges of the ALS height bins used for the stratified error table.
HEIGHT_BINS = [0, 10, 15, 20, 25, 30, 35, 40, 45, 50, 100]

#: Figure resolution.
FIG_DPI = 200
