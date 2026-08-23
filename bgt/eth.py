"""The ETH global canopy height product (Lang et al. 2023) as a comparison target.

10 m canopy top height regressed from Sentinel-2 with GEDI as the training
reference, representative of the year 2020. ``01b_download_eth.py`` clips it to
each ALS footprint; this module reads those clips.

Two properties of the product shape everything downstream:

* it is delivered on a **geographic** grid (EPSG:4326), so its cells are not
  square in metres -- roughly 9.26 m east-west and 9.28 m north-south at these
  latitudes. ``coreg._pixel_size_m`` already handles that, so the same fine-grid
  and block-aggregation machinery works with it as a target grid;
* it is **uint8 with 255 as no-data**, so heights are whole metres in 0-254.
  The quantisation puts a floor of about 0.29 m on any RMSE computed against it,
  which is far below the errors of interest but visible as banding in scatters.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import rasterio
from rasterio.crs import CRS
from rasterio.transform import Affine
from rasterio.warp import transform_bounds

import config
from bgt.coreg import _pixel_size_m


@dataclass
class EthChm:
    """An ETH canopy height window, plus its predictive standard deviation.

    Carries the same attributes the co-registration code expects of a target
    grid (``height``, ``transform``, ``crs``, ``shape``, ``cell_size_m``), so it
    can be used either as the grid to aggregate onto or as a source to aggregate
    up from.
    """

    height: np.ndarray
    sd: np.ndarray | None
    transform: Affine
    crs: CRS
    source: Path
    meta: dict[str, Any] = field(default_factory=dict)
    name: str = "ETH canopy height 10 m (2020)"

    @property
    def shape(self) -> tuple[int, int]:
        return self.height.shape  # type: ignore[return-value]

    def cell_size_m(self) -> tuple[float, float]:
        return _pixel_size_m(self.transform, self.crs, self.shape)


def clip_path(site: str, eth_dir: Path = config.ETH_DIR, sd: bool = False) -> Path:
    """Where ``01b_download_eth.py`` puts the clip for ``site``."""
    stem = "eth_canopy_height_sd" if sd else "eth_canopy_height"
    return eth_dir / site / f"{stem}_{site}.tif"


def load_eth(
    site: str,
    footprint_bounds: tuple[float, float, float, float],
    footprint_crs: CRS,
    pad_m: float = 0.0,
    eth_dir: Path = config.ETH_DIR,
    height_min: float = 0.0,
    height_max: float = config.ETH_HEIGHT_MAX,
    with_sd: bool = True,
) -> EthChm | None:
    """Read the ETH clip for ``site`` over the ALS footprint.

    Returns ``None`` when no clip exists or it holds no valid pixel over the
    footprint, so a missing site degrades to "skip the ETH comparison" rather
    than failing the run.
    """
    path = clip_path(site, eth_dir)
    if not path.is_file():
        print(f"  no ETH clip for {site} at {path} -- run 01b_download_eth.py")
        return None

    with rasterio.open(path) as src:
        left, bottom, right, top = transform_bounds(
            footprint_crs, src.crs, *footprint_bounds, densify_pts=51
        )
        px_x, px_y = _pixel_size_m(src.transform, src.crs, src.shape)
        pad_x = pad_m / px_x * abs(src.transform.a)
        pad_y = pad_m / px_y * abs(src.transform.e)

        window = rasterio.windows.from_bounds(
            left - pad_x, bottom - pad_y, right + pad_x, top + pad_y,
            transform=src.transform,
        )
        window = window.round_offsets(op="floor").round_lengths(op="ceil")
        window = window.intersection(
            rasterio.windows.Window(0, 0, src.width, src.height)
        )
        if window.width < 2 or window.height < 2:
            return None

        # uint8 cannot hold nan, so widen before filling the mask.
        height = (
            src.read(1, window=window, masked=True).astype(np.float64).filled(np.nan)
        )
        transform = src.window_transform(window)
        crs = src.crs
        full_shape = src.shape
        meta = {
            "eth_path": str(path),
            "eth_crs": str(crs),
            "eth_dtype": src.dtypes[0],
            "eth_nodata": src.nodata,
            "eth_pixel_m_x": px_x,
            "eth_pixel_m_y": px_y,
            "eth_tiles": src.tags().get("tiles", ""),
        }

    outside = np.isfinite(height) & ((height < height_min) | (height > height_max))
    height[outside] = np.nan

    sd = None
    sd_path = clip_path(site, eth_dir, sd=True)
    if with_sd and sd_path.is_file():
        with rasterio.open(sd_path) as src:
            # The two clips are written from the same window, so the same window
            # applies to both -- but only if they really are on one grid.
            if src.shape == full_shape:
                sd = (
                    src.read(1, window=window, masked=True)
                    .astype(np.float64)
                    .filled(np.nan)
                )
                meta["eth_sd_path"] = str(sd_path)
            else:
                print("  ETH SD clip is not on the height grid, ignored")

    if not np.isfinite(height).any():
        print(f"  ETH clip for {site} has no valid pixel over the footprint")
        return None

    n = int(np.isfinite(height).sum())
    print(f"  ETH window {height.shape} cells, ~{px_x:.1f} x {px_y:.1f} m, "
          f"{n:,} valid cells")

    return EthChm(height=height, sd=sd, transform=transform, crs=crs,
                  source=path, meta=meta)
