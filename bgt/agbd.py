"""The 50 m AGBD reference map and the tomographic HV layers that share its grid.

``data/AGBD_<site>.tif`` is a ten-band raster produced from the ALS: band 2 is
the above-ground biomass density in t/ha, band 10 the fraction of the cell that
the ALS quality masks reject. ``data/Tomo/Tomo_hv_norm_<h>.tif`` holds the
normalised tomographic HV power at height ``h`` above the ground.

The tomographic layers carry no geotransform -- they are written on exactly the
same row/column grid as the AGBD raster, so they are paired by pixel index and
the AGBD file supplies the geometry for both.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import rasterio
from affine import Affine
from rasterio.crs import CRS

import config


@dataclass(frozen=True)
class AgbdGrid:
    """The AGBD map together with the grid every other layer is paired on."""

    agbd: np.ndarray          #: t/ha, nan where the map has no estimate
    masked_fraction: np.ndarray
    transform: Affine
    crs: CRS
    path: Path

    @property
    def shape(self) -> tuple[int, int]:
        return self.agbd.shape

    def valid(self, max_masked: float) -> np.ndarray:
        """Cells with an AGBD estimate and little enough ALS masking."""
        return np.isfinite(self.agbd) & (self.masked_fraction <= max_masked)


def agbd_path(site: str) -> Path:
    return config.AGBD_DIR / config.AGBD_TEMPLATE.format(site=site)


def load_agbd(site: str) -> AgbdGrid:
    """Read the AGBD band and the mask band of the site's reference map."""
    path = agbd_path(site)
    if not path.is_file():
        raise FileNotFoundError(f"No AGBD map for {site}: {path}")

    with rasterio.open(path) as src:
        names = list(src.descriptions)
        agbd = src.read(_band(names, config.AGBD_BAND_NAME), masked=True)
        mask = src.read(_band(names, config.AGBD_MASK_BAND_NAME), masked=True)
        return AgbdGrid(
            agbd=agbd.filled(np.nan).astype(float),
            # A cell the mask band says nothing about is treated as unmasked;
            # its AGBD value is what decides whether it is used.
            masked_fraction=mask.filled(0.0).astype(float),
            transform=src.transform,
            crs=src.crs,
            path=path,
        )


def _band(names: list[str | None], wanted: str) -> int:
    """1-based index of the band called ``wanted``."""
    for index, name in enumerate(names, start=1):
        if name == wanted:
            return index
    raise KeyError(f"No band named {wanted!r}; the file has {names}")


def tomo_path(height_m: int) -> Path:
    return config.TOMO_DIR / config.TOMO_TEMPLATE.format(height=height_m)


def tomo_heights() -> list[int]:
    """Heights, in metres, for which a tomographic layer exists."""
    prefix, suffix = config.TOMO_TEMPLATE.split("{height}")
    heights = []
    for path in config.TOMO_DIR.glob(config.TOMO_TEMPLATE.format(height="*")):
        heights.append(int(path.name[len(prefix):len(path.name) - len(suffix)]))
    return sorted(heights)


def load_tomo(height_m: int, shape: tuple[int, int]) -> np.ndarray:
    """Read the tomographic HV layer at ``height_m`` as float, nan where no data.

    ``shape`` is the AGBD grid the layer is paired on; a layer that does not
    match it cannot be paired by pixel index and is rejected rather than
    silently resampled.
    """
    path = tomo_path(height_m)
    if not path.is_file():
        available = tomo_heights()
        raise FileNotFoundError(
            f"No tomographic layer at {height_m} m: {path}. Available: {available}"
        )

    with rasterio.open(path) as src:
        if src.shape != shape:
            raise ValueError(
                f"{path.name} is {src.shape}, but the AGBD grid is {shape}; the "
                "two are paired by pixel index and must agree"
            )
        layer = src.read(1).astype(float)
        if src.nodata is not None:
            layer[layer == src.nodata] = np.nan
    if config.TOMO_NODATA_DN is not None:
        layer[layer == config.TOMO_NODATA_DN] = np.nan
    return layer
