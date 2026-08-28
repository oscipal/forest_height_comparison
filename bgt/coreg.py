"""Co-registration of the 1 m ALS canopy height model onto the BIOMASS grid.

Design decision: **the BIOMASS grid is the target grid and is never resampled.**
Interpolating a coarse product onto a fine grid would invent detail it does not
have and would correlate neighbouring residuals. Instead the ALS CHM is
aggregated *up* to the BIOMASS cells:

1. Build a fine grid that is an exact integer refinement (factor ``n``) of the
   BIOMASS grid over the ALS footprint, padded so that shifted sampling windows
   stay inside it.
2. Warp the ALS CHM once onto that fine grid with area-average resampling,
   carrying a companion coverage grid that records the valid ALS area fraction
   of every fine cell.
3. Block-reduce ``n x n`` fine cells per BIOMASS cell. Because the fine grid is
   an exact refinement, every block corresponds to exactly one BIOMASS pixel --
   no interpolation, no partial cells.
4. Search a planimetric shift by sliding the block origin over the fine grid and
   maximising the agreement between the aggregated ALS and BIOMASS heights. The
   search is exhaustive but cheap: block sums come from an integral image.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import rasterio
from affine import Affine
from rasterio.crs import CRS
from rasterio.warp import Resampling, reproject, transform_bounds
from scipy import ndimage

import config
from bgt.als import AlsChm

_EARTH_M_PER_DEG = 111_320.0


# --------------------------------------------------------------------------- #
# Why a cell is excluded
# --------------------------------------------------------------------------- #
#
# Every cell of the analysis window carries exactly one code. The codes are
# assigned in the order listed, so a cell that fails several tests is attributed
# to the first one it fails -- which makes the counts additive and lets them be
# read as a funnel from "every cell in the window" down to "paired".
#
# The ALS footprint test comes first deliberately: the window is a rectangle
# around a scan that is rarely rectangular, so most of what it excludes was
# never a candidate rather than a loss of usable data.

KEPT = 0
ALS_ABSENT = 1
FH_NODATA = 2
FH_OUT_OF_RANGE = 3
FH_QUALITY = 4
ALS_PARTIAL = 5
ALS_STAT_MISSING = 6

#: Ordered so that iteration yields the funnel from top to bottom.
EXCLUSION_ORDER = [
    ALS_ABSENT,
    FH_NODATA,
    FH_OUT_OF_RANGE,
    FH_QUALITY,
    ALS_PARTIAL,
    ALS_STAT_MISSING,
    KEPT,
]

EXCLUSION_LABELS = {
    KEPT: "kept (paired)",
    ALS_ABSENT: "outside ALS footprint",
    FH_NODATA: "BIOMASS no-data",
    FH_OUT_OF_RANGE: "BIOMASS height out of range",
    FH_QUALITY: "BIOMASS quality above threshold",
    ALS_PARTIAL: "ALS coverage below minimum",
    ALS_STAT_MISSING: "ALS statistic unavailable",
}

#: Labels keyed by short name, in funnel order -- iteration order is what the
#: exclusion chart uses to lay out its bars.
EXCLUSION_LABELS_BY_KEY: dict[str, str] = {}

#: Short keys for the CSV columns.
EXCLUSION_KEYS = {
    KEPT: "kept",
    ALS_ABSENT: "als_absent",
    FH_NODATA: "fh_nodata",
    FH_OUT_OF_RANGE: "fh_out_of_range",
    FH_QUALITY: "fh_quality_rejected",
    ALS_PARTIAL: "als_coverage_low",
    ALS_STAT_MISSING: "als_stat_missing",
}

EXCLUSION_LABELS_BY_KEY.update(
    {EXCLUSION_KEYS[code]: EXCLUSION_LABELS[code] for code in EXCLUSION_ORDER}
)


# --------------------------------------------------------------------------- #
# BIOMASS side
# --------------------------------------------------------------------------- #


@dataclass
class BiomassFh:
    """A BIOMASS L2A forest height window plus its quality layer."""

    height: np.ndarray
    quality: np.ndarray | None
    reject: np.ndarray  #: per-cell exclusion code, BIOMASS-side reasons only
    transform: Affine
    crs: CRS
    item_id: str
    source: Path
    meta: dict[str, Any] = field(default_factory=dict)

    @property
    def shape(self) -> tuple[int, int]:
        return self.height.shape  # type: ignore[return-value]

    def cell_size_m(self) -> tuple[float, float]:
        """Approximate cell size in metres (handles geographic CRSs)."""
        return _pixel_size_m(self.transform, self.crs, self.shape)


def _pixel_size_m(
    transform: Affine, crs: CRS, shape: tuple[int, int]
) -> tuple[float, float]:
    if crs.is_geographic:
        centre_lat = transform.f + transform.e * shape[0] / 2.0
        x_m = abs(transform.a) * _EARTH_M_PER_DEG * np.cos(np.deg2rad(centre_lat))
        y_m = abs(transform.e) * _EARTH_M_PER_DEG
        return float(x_m), float(y_m)
    return float(abs(transform.a)), float(abs(transform.e))


def load_fh(
    fh_path: Path,
    quality_path: Path | None,
    footprint_bounds: tuple[float, float, float, float],
    footprint_crs: CRS,
    pad_m: float,
    fh_min: float = config.FH_MIN,
    fh_max: float = config.FH_MAX,
    quality_max: float | None = config.FH_QUALITY_MAX,
) -> BiomassFh | None:
    """Read the BIOMASS FH raster over the ALS footprint.

    Returns ``None`` when the product does not actually cover the footprint with
    valid data -- BIOMASS scene bounding boxes are generous, so a scene can be
    listed as overlapping while its usable swath is elsewhere.
    """
    with rasterio.open(fh_path) as src:
        left, bottom, right, top = transform_bounds(
            footprint_crs, src.crs, *footprint_bounds, densify_pts=51
        )
        px_x, px_y = _pixel_size_m(src.transform, src.crs, src.shape)
        pad_x = pad_m / px_x * abs(src.transform.a)
        pad_y = pad_m / px_y * abs(src.transform.e)

        window = rasterio.windows.from_bounds(
            left - pad_x,
            bottom - pad_y,
            right + pad_x,
            top + pad_y,
            transform=src.transform,
        )
        window = window.round_offsets(op="floor").round_lengths(op="ceil")
        full = rasterio.windows.Window(0, 0, src.width, src.height)
        # A scene can be catalogued as overlapping the footprint and still not
        # reach it: the frames are large and their bounding boxes generous. That
        # leaves the two windows disjoint, which rasterio raises on rather than
        # returning an empty window, so it is tested for first.
        if not rasterio.windows.intersect([window, full]):
            return None
        window = window.intersection(full)
        if window.width < 2 or window.height < 2:
            return None

        height = src.read(1, window=window, masked=True).filled(np.nan).astype(np.float64)
        transform = src.window_transform(window)
        crs = src.crs
        meta = {
            "fh_path": str(fh_path),
            "fh_crs": str(crs),
            "fh_dtype": src.dtypes[0],
            "fh_nodata": src.nodata,
            "fh_full_shape": src.shape,
            "fh_pixel_m_x": px_x,
            "fh_pixel_m_y": px_y,
        }

    # Record why each BIOMASS cell drops out, before the value is discarded.
    reject = np.where(np.isfinite(height), KEPT, FH_NODATA).astype(np.uint8)
    out_of_range = np.isfinite(height) & ((height < fh_min) | (height > fh_max))
    reject[out_of_range] = FH_OUT_OF_RANGE
    height[out_of_range] = np.nan

    quality = None
    if quality_path is not None:
        # Opened rather than stat'ed first: BIOMASS product names are long enough
        # that a nested path can exceed the Windows MAX_PATH limit, under which
        # Path.is_file() reports False for a file GDAL opens without trouble.
        try:
            # The FH and quality rasters ship on the same product grid, so the
            # same window applies to both.
            with rasterio.open(quality_path) as src:
                if src.shape != meta["fh_full_shape"]:
                    raise ValueError(
                        "The FH and quality rasters are not on the same grid; "
                        "cannot align them by window."
                    )
                quality = (
                    src.read(1, window=window, masked=True)
                    .filled(np.nan)
                    .astype(np.float64)
                )
            meta["quality_path"] = str(quality_path)
        except rasterio.errors.RasterioIOError as exc:
            print(f"  quality layer unreadable, continuing without it: {exc}")

    if quality is not None and quality_max is not None:
        # `~(q <= max)` also catches a nan quality value, which is the intended
        # conservative reading: unknown quality is not passing quality.
        failed = (reject == KEPT) & ~(quality <= quality_max)
        reject[failed] = FH_QUALITY
        height[failed] = np.nan

    if not np.isfinite(height).any():
        return None

    return BiomassFh(
        height=height,
        quality=quality,
        reject=reject,
        transform=transform,
        crs=crs,
        item_id=Path(fh_path).parent.parent.name,
        source=Path(fh_path),
        meta=meta,
    )


def exclusion_grid(
    fh: BiomassFh,
    coverage: np.ndarray,
    als_stat: np.ndarray,
    min_coverage: float = config.MIN_ALS_COVERAGE,
) -> np.ndarray:
    """One exclusion code per cell, combining the BIOMASS and ALS side tests.

    Codes are assigned in ``EXCLUSION_ORDER``: each test only sees the cells
    that survived the tests before it, so the resulting counts sum to the size
    of the analysis window and read as a funnel.
    """
    reason = np.full(fh.shape, KEPT, dtype=np.uint8)

    # 1. No ALS at all: the window is a rectangle, the scan usually is not.
    reason[coverage <= 0.0] = ALS_ABSENT

    # 2. BIOMASS-side reasons, for cells that do have ALS underneath them.
    from_fh = (reason == KEPT) & (fh.reject != KEPT)
    reason[from_fh] = fh.reject[from_fh]

    # 3. Partial ALS coverage: some, but not enough to summarise the cell.
    reason[(reason == KEPT) & (coverage < min_coverage)] = ALS_PARTIAL

    # 4. Anything left without a usable ALS statistic (rare; belt and braces).
    reason[(reason == KEPT) & ~np.isfinite(als_stat)] = ALS_STAT_MISSING
    return reason


def exclusion_counts(reason: np.ndarray) -> dict[str, int]:
    """Cell counts per exclusion code, keyed by short name, funnel-ordered."""
    return {
        EXCLUSION_KEYS[code]: int((reason == code).sum())
        for code in EXCLUSION_ORDER
    }


# --------------------------------------------------------------------------- #
# Fine grid and warping
# --------------------------------------------------------------------------- #


@dataclass
class FineGrid:
    """The ALS CHM resampled onto an integer refinement of the BIOMASS grid."""

    values: np.ndarray  #: area-averaged ALS height per fine cell (nan where empty)
    coverage: np.ndarray  #: valid ALS area fraction per fine cell, 0..1
    transform: Affine
    crs: CRS
    refine: int  #: fine cells per BIOMASS cell edge
    pad: int  #: BIOMASS cells of padding on each side
    fine_m: tuple[float, float]  #: fine cell size in metres (x, y)


def build_fine_grid(
    chm: AlsChm,
    fh: BiomassFh,
    fine_cell_m: float = config.FINE_CELL_M,
    shift_search_m: float = config.SHIFT_SEARCH_M,
) -> FineGrid:
    """Warp ``chm`` onto a padded integer refinement of the BIOMASS grid."""
    px_x, px_y = fh.cell_size_m()
    refine = max(1, int(np.ceil(min(px_x, px_y) / fine_cell_m)))
    fine_x_m, fine_y_m = px_x / refine, px_y / refine

    pad = int(np.ceil(shift_search_m / min(px_x, px_y))) + 1
    h, w = fh.shape
    fine_shape = ((h + 2 * pad) * refine, (w + 2 * pad) * refine)
    fine_transform = fh.transform * Affine.translation(-pad, -pad) * Affine.scale(
        1.0 / refine
    )

    print(
        f"  fine grid: refine={refine} ({fine_x_m:.2f} x {fine_y_m:.2f} m cells), "
        f"pad={pad} BIOMASS cells, shape={fine_shape}"
    )

    values = np.full(fine_shape, np.nan, dtype=np.float32)
    reproject(
        source=chm.height,
        destination=values,
        src_transform=chm.transform,
        src_crs=chm.crs,
        src_nodata=np.nan,
        dst_transform=fine_transform,
        dst_crs=fh.crs,
        dst_nodata=np.nan,
        resampling=Resampling.average,
        num_threads=4,
    )

    valid_src = np.isfinite(chm.height).astype(np.float32)
    coverage = np.zeros(fine_shape, dtype=np.float32)
    reproject(
        source=valid_src,
        destination=coverage,
        src_transform=chm.transform,
        src_crs=chm.crs,
        src_nodata=None,
        dst_transform=fine_transform,
        dst_crs=fh.crs,
        dst_nodata=0.0,
        init_dest_nodata=True,
        resampling=Resampling.average,
        num_threads=4,
    )

    # A fine cell with no valid ALS contribution must not pretend to have one.
    coverage[~np.isfinite(values)] = 0.0
    values[coverage <= 0.0] = np.nan

    return FineGrid(
        values=values,
        coverage=coverage,
        transform=fine_transform,
        crs=fh.crs,
        refine=refine,
        pad=pad,
        fine_m=(fine_x_m, fine_y_m),
    )


# --------------------------------------------------------------------------- #
# Block aggregation
# --------------------------------------------------------------------------- #


def _integral(a: np.ndarray) -> np.ndarray:
    """Summed-area table with a zero row and column prepended."""
    out = np.zeros((a.shape[0] + 1, a.shape[1] + 1), dtype=np.float64)
    np.cumsum(np.cumsum(a, axis=0), axis=1, out=out[1:, 1:])
    return out


def _block_sums(integral: np.ndarray, starts_y: np.ndarray, starts_x: np.ndarray,
                n: int) -> np.ndarray:
    y0, y1 = starts_y, starts_y + n
    x0, x1 = starts_x, starts_x + n
    return (
        integral[np.ix_(y1, x1)]
        - integral[np.ix_(y0, x1)]
        - integral[np.ix_(y1, x0)]
        + integral[np.ix_(y0, x0)]
    )


class BlockAggregator:
    """Aggregates the fine ALS grid into BIOMASS cells for any sampling offset.

    An offset of ``(oy, ox)`` fine cells means the ALS window for a BIOMASS cell
    is taken ``ox`` fine cells to the east and ``oy`` fine cells to the south of
    its nominal position.
    """

    def __init__(self, grid: FineGrid, fh_shape: tuple[int, int]):
        self.grid = grid
        self.fh_shape = fh_shape
        n, pad = grid.refine, grid.pad

        # The summed-area tables are what make the shift search affordable, but
        # they are float64 and the size of the fine grid -- at a 1 m fine cell
        # that is hundreds of megabytes each. An aggregator that is only ever
        # asked for one offset (the secondary CHM, ETH) never needs them, so
        # they are built on first use rather than in the constructor.
        self._int_w: np.ndarray | None = None
        self._int_vw: np.ndarray | None = None

        h, w = fh_shape
        self._base_y = pad * n + n * np.arange(h)
        self._base_x = pad * n + n * np.arange(w)
        self._n = n
        self.max_offset = pad * n  # keeps every shifted window inside the grid

    def _integrals(self) -> tuple[np.ndarray, np.ndarray]:
        if self._int_w is None:
            weights = np.nan_to_num(self.grid.coverage, nan=0.0)
            weighted = np.nan_to_num(self.grid.values, nan=0.0) * weights
            self._int_w = _integral(weights)
            self._int_vw = _integral(weighted)
        return self._int_w, self._int_vw  # type: ignore[return-value]

    def release(self) -> None:
        """Drop the cached summed-area tables once the shift search is done."""
        self._int_w = None
        self._int_vw = None

    def mean_and_coverage(self, oy: int, ox: int) -> tuple[np.ndarray, np.ndarray]:
        """Coverage-weighted mean height and area coverage per BIOMASS cell.

        Uses the summed-area tables, which pay for themselves across the
        thousands of offsets the shift search evaluates.
        """
        if abs(oy) > self.max_offset or abs(ox) > self.max_offset:
            raise ValueError(f"offset ({oy}, {ox}) exceeds the padded fine grid")
        int_w, int_vw = self._integrals()
        sy, sx = self._base_y + oy, self._base_x + ox
        w_sum = _block_sums(int_w, sy, sx, self._n)
        vw_sum = _block_sums(int_vw, sy, sx, self._n)
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = np.where(w_sum > 0, vw_sum / w_sum, np.nan)
        coverage = w_sum / (self._n * self._n)
        return mean, coverage

    def _mean_and_coverage_direct(
        self, oy: int, ox: int
    ) -> tuple[np.ndarray, np.ndarray]:
        """The same quantities for a single offset, without a summed-area table.

        Reduces the blocks directly. Slower per call, but it avoids allocating
        two float64 copies of the fine grid for an aggregator that will only
        ever be asked about one offset.

        It is also the more accurate of the two: summing within a block in
        float64 is exact, whereas the summed-area tables accumulate in float32
        across the whole grid and then difference large numbers, which costs
        about 1e-4 m on these grids. That is irrelevant to the shift search,
        which only ranks correlations, but the statistics that get reported come
        through here, so they take the exact route.
        """
        n, (h, w) = self._n, self.fh_shape
        sy, sx = self._base_y[0] + oy, self._base_x[0] + ox
        weights = np.nan_to_num(
            self.grid.coverage[sy : sy + h * n, sx : sx + w * n], nan=0.0
        ).astype(np.float64)
        values = np.nan_to_num(
            self.grid.values[sy : sy + h * n, sx : sx + w * n], nan=0.0
        ).astype(np.float64)
        w_sum = weights.reshape(h, n, w, n).sum(axis=(1, 3))
        vw_sum = (values * weights).reshape(h, n, w, n).sum(axis=(1, 3))
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = np.where(w_sum > 0, vw_sum / w_sum, np.nan)
        return mean, w_sum / (n * n)

    def full_stats(self, oy: int, ox: int, stats: list[str]) -> dict[str, np.ndarray]:
        """All requested per-cell statistics at one offset.

        Percentiles are taken over the fine cells of the block. Fine cells are
        equal-area by construction, so this is an unweighted percentile of an
        area-regular sample of the canopy surface.
        """
        n, h, w = self._n, *self.fh_shape
        sy, sx = self._base_y[0] + oy, self._base_x[0] + ox
        block = self.grid.values[sy : sy + h * n, sx : sx + w * n]
        block = block.reshape(h, n, w, n).transpose(0, 2, 1, 3).reshape(h, w, n * n)

        # A single offset: reduce directly rather than paying for the tables.
        mean, coverage = (
            self.mean_and_coverage(oy, ox)
            if self._int_w is not None
            else self._mean_and_coverage_direct(oy, ox)
        )
        out: dict[str, np.ndarray] = {"coverage": coverage}
        finite = np.isfinite(block)
        out["n_fine"] = finite.sum(axis=2).astype(np.int32)

        percentiles = sorted(
            {int(s[1:]) for s in stats if s.startswith("p") and s[1:].isdigit()}
        )
        # Cells outside the ALS footprint are legitimately all-NaN; numpy's
        # warning about them carries no information here.
        with np.errstate(invalid="ignore"), warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="All-NaN", category=RuntimeWarning)
            warnings.filterwarnings("ignore", message="Mean of empty slice",
                                    category=RuntimeWarning)
            warnings.filterwarnings("ignore", message="Degrees of freedom",
                                    category=RuntimeWarning)
            if percentiles:
                pvals = np.nanpercentile(block, percentiles, axis=2)
                for q, arr in zip(percentiles, pvals):
                    out[f"p{q}"] = arr
            if "mean" in stats:
                out["mean"] = mean
            if "max" in stats:
                out["max"] = np.nanmax(np.where(finite, block, -np.inf), axis=2)
                out["max"][out["n_fine"] == 0] = np.nan
            if "std" in stats:
                out["std"] = np.nanstd(block, axis=2)
        return out


# --------------------------------------------------------------------------- #
# Shift search
# --------------------------------------------------------------------------- #


@dataclass
class ShiftResult:
    offset_yx: tuple[int, int]
    shift_east_m: float
    shift_north_m: float
    offsets_y: np.ndarray
    offsets_x: np.ndarray
    fine_m: tuple[float, float]  #: fine cell size (x, y) in metres
    pearson: np.ndarray
    rmse: np.ndarray
    criterion: str
    n_cells: int
    pearson_at_zero: float
    rmse_at_zero: float
    pearson_at_best: float
    rmse_at_best: float
    accepted: bool = True  #: whether the optimum was applied
    reason: str = ""  #: why it was or was not applied

    @property
    def applied_offset(self) -> tuple[int, int]:
        """The offset actually used: the optimum, or no shift if it was rejected."""
        return self.offset_yx if self.accepted else (0, 0)


def find_shift(
    aggregator: BlockAggregator,
    fh_height: np.ndarray,
    search_m: float = config.SHIFT_SEARCH_M,
    step_m: float = config.SHIFT_STEP_M,
    criterion: str = config.SHIFT_CRITERION,
    min_coverage: float = config.MIN_ALS_COVERAGE,
    min_gain: float = config.SHIFT_MIN_GAIN,
    max_rim_frac: float = config.SHIFT_MAX_RIM_FRAC,
    min_cells: int = config.SHIFT_MIN_CELLS,
    max_disagree_m: float = config.SHIFT_MAX_DISAGREE_M,
    force: bool = False,
) -> ShiftResult:
    """Exhaustively search the planimetric shift that best aligns ALS and BIOMASS.

    The comparison set is fixed across all candidate offsets: cells that are
    fully covered by ALS at zero offset and far enough from the footprint edge
    that no candidate window can run off it. Without this, offsets that pull the
    window towards the dense interior of the scan would win for the wrong reason.
    """
    grid = aggregator.grid
    fine_x_m, fine_y_m = grid.fine_m

    max_off_x = min(int(round(search_m / fine_x_m)), aggregator.max_offset)
    max_off_y = min(int(round(search_m / fine_y_m)), aggregator.max_offset)
    step_x = max(1, int(round(step_m / fine_x_m)))
    step_y = max(1, int(round(step_m / fine_y_m)))

    offsets_x = np.arange(-max_off_x, max_off_x + 1, step_x)
    offsets_y = np.arange(-max_off_y, max_off_y + 1, step_y)
    if 0 not in offsets_x:
        offsets_x = np.sort(np.append(offsets_x, 0))
    if 0 not in offsets_y:
        offsets_y = np.sort(np.append(offsets_y, 0))

    _, cov0 = aggregator.mean_and_coverage(0, 0)
    core = np.isfinite(fh_height) & (cov0 >= min_coverage)
    fh_cell_m = min(fine_x_m, fine_y_m) * grid.refine
    erode_cells = int(np.ceil(search_m / fh_cell_m))
    if erode_cells > 0:
        core = ndimage.binary_erosion(
            core, structure=np.ones((3, 3), bool), iterations=erode_cells,
            border_value=0,
        )

    n_cells = int(core.sum())
    if n_cells < 20:
        raise ValueError(
            f"Only {n_cells} BIOMASS cells survive the co-registration mask "
            "(need >= 20). The ALS footprint is probably too small relative to "
            "the BIOMASS cell size or to SHIFT_SEARCH_M."
        )
    print(
        f"  shift search: {len(offsets_y)} x {len(offsets_x)} offsets over "
        f"+/-{search_m:.0f} m, on {n_cells} fully covered cells"
    )

    target = fh_height[core]
    pearson = np.full((len(offsets_y), len(offsets_x)), np.nan)
    rmse = np.full_like(pearson, np.nan)

    for i, oy in enumerate(offsets_y):
        for j, ox in enumerate(offsets_x):
            mean, _ = aggregator.mean_and_coverage(int(oy), int(ox))
            sample = mean[core]
            ok = np.isfinite(sample)
            if ok.sum() < 20:
                continue
            a, b = sample[ok], target[ok]
            if a.std() == 0 or b.std() == 0:
                continue
            pearson[i, j] = float(np.corrcoef(a, b)[0, 1])
            rmse[i, j] = float(np.sqrt(np.mean((a - b) ** 2)))

    if criterion == "pearson":
        flat = np.nanargmax(pearson)
    elif criterion == "rmse":
        flat = np.nanargmin(rmse)
    else:
        raise ValueError(f"unknown shift criterion: {criterion!r}")
    bi, bj = np.unravel_index(flat, pearson.shape)
    oy, ox = int(offsets_y[bi]), int(offsets_x[bj])

    zi = int(np.where(offsets_y == 0)[0][0])
    zj = int(np.where(offsets_x == 0)[0][0])

    # A trustworthy co-registration shows a genuine, well-supported *peak*.
    # Each check below rules out a different way of being fooled; the first one
    # that fails blocks the shift, which is then reported but not applied.
    ri, rj = np.unravel_index(np.nanargmin(rmse), rmse.shape)
    disagree_m = float(
        np.hypot(
            (offsets_y[bi] - offsets_y[ri]) * fine_y_m,
            (offsets_x[bj] - offsets_x[rj]) * fine_x_m,
        )
    )
    rim_frac = max(
        abs(offsets_y[bi]) / max(max_off_y, 1),
        abs(offsets_x[bj]) / max(max_off_x, 1),
    )
    gain = float(pearson[bi, bj] - pearson[zi, zj])

    accepted, reason = True, "well-supported interior peak, accepted"
    if n_cells < min_cells:
        accepted = False
        reason = (
            f"only {n_cells} cells in the search set (minimum {min_cells}): the "
            "best of thousands of candidate offsets is not distinguishable from "
            "a chance fit, so no shift is applied"
        )
    elif rim_frac > max_rim_frac:
        accepted = False
        reason = (
            f"optimum sits at {rim_frac:.0%} of the search radius (limit "
            f"{max_rim_frac:.0%}): the objective runs to the rim rather than "
            "peaking inside it, so no shift is applied"
        )
    elif disagree_m > max_disagree_m:
        accepted = False
        reason = (
            f"the correlation and RMSE optima disagree by {disagree_m:.0f} m "
            f"(limit {max_disagree_m:.0f} m): no consistent alignment, so no "
            "shift is applied"
        )
    elif gain < min_gain:
        accepted = False
        reason = (
            f"correlation gain {gain:+.3f} is below the {min_gain:.3f} "
            "threshold: no shift is applied"
        )
    if force:
        accepted, reason = True, reason + " (overridden by --force-shift)"

    result = ShiftResult(
        offset_yx=(oy, ox),
        shift_east_m=ox * fine_x_m,
        shift_north_m=-oy * fine_y_m,
        offsets_y=offsets_y,
        offsets_x=offsets_x,
        fine_m=(fine_x_m, fine_y_m),
        pearson=pearson,
        rmse=rmse,
        criterion=criterion,
        n_cells=n_cells,
        pearson_at_zero=float(pearson[zi, zj]),
        rmse_at_zero=float(rmse[zi, zj]),
        pearson_at_best=float(pearson[bi, bj]),
        rmse_at_best=float(rmse[bi, bj]),
        accepted=accepted,
        reason=reason,
    )
    print(
        f"  best offset: east {result.shift_east_m:+.1f} m, "
        f"north {result.shift_north_m:+.1f} m  "
        f"(r {result.pearson_at_zero:.3f} -> {result.pearson_at_best:.3f}, "
        f"RMSE {result.rmse_at_zero:.2f} -> {result.rmse_at_best:.2f} m)"
    )
    print(f"  {'ACCEPTED' if accepted else 'REJECTED'}: {reason}")
    return result
