"""The Meta global canopy height map v2 (DINOv3) as a comparison target.

~1.2 m canopy height predicted from Maxar very-high-resolution imagery,
published by Meta on AWS (``s3://dataforgood-fb-data/forests/v2/global/
dinov3_global_chm_v2_ml3/``). ``01c_download_meta.py`` clips it to each ALS
footprint; this module locates the tiles and reads those clips.

Three properties of the product shape everything downstream:

* it is tiled by **zoom-10 quadkey** in Web Mercator (EPSG:3857), so a tile
  name is computed from a coordinate rather than looked up;
* it is **uint8 in whole metres with no no-data value**. Zero is a real
  prediction -- open water, a clearing -- not a gap, and is kept;
* it is **not a single epoch**. Each tile is a mosaic of imagery from different
  dates, recorded per polygon in ``metadata/<quadkey>.geojson``. The clip step
  writes the share of the footprint under each date, since that date, not the
  product's release, is what the ALS time gap has to be measured against.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import rasterio
from affine import Affine
from rasterio.crs import CRS

import config
from bgt.als import AlsChm


def quadkey(lon: float, lat: float, zoom: int = config.META_TILE_ZOOM) -> str:
    """Quadkey of the Web Mercator tile containing ``(lon, lat)``."""
    n = 2 ** zoom
    x = int((lon + 180.0) / 360.0 * n)
    s = math.sin(math.radians(lat))
    y = int((0.5 - math.log((1 + s) / (1 - s)) / (4 * math.pi)) * n)
    x, y = min(max(x, 0), n - 1), min(max(y, 0), n - 1)
    digits = []
    for i in range(zoom, 0, -1):
        mask = 1 << (i - 1)
        digits.append(str((1 if x & mask else 0) + (2 if y & mask else 0)))
    return "".join(digits)


def tiles_for_bbox(bbox: tuple[float, float, float, float],
                   zoom: int = config.META_TILE_ZOOM) -> list[str]:
    """Every tile touched by ``(minx, miny, maxx, maxy)`` in lon/lat.

    Sampled on a lattice finer than one tile, so a box spanning several tiles
    picks up all of them, not only its corners'.
    """
    minx, miny, maxx, maxy = bbox
    step = 360.0 / 2 ** zoom / 4
    xs = np.append(np.arange(minx, maxx, step), maxx)
    ys = np.append(np.arange(miny, maxy, step), maxy)
    return sorted({quadkey(x, y, zoom) for x in xs for y in ys})


def tile_url(tile: str) -> str:
    return f"{config.META_BASE_URL}/chm/{tile}.tif"


def metadata_url(tile: str) -> str:
    return f"{config.META_BASE_URL}/metadata/{tile}.geojson"


def clip_path(site: str, meta_dir: Path = config.META_DIR) -> Path:
    """Where ``01c_download_meta.py`` puts the height clip for ``site``."""
    return meta_dir / site / f"meta_chm_{site}.tif"


def dates_path(site: str, meta_dir: Path = config.META_DIR) -> Path:
    """Where ``01c_download_meta.py`` puts the imagery dates for ``site``."""
    return meta_dir / site / f"meta_imagery_dates_{site}.csv"


def load_clip(site: str, meta_dir: Path = config.META_DIR) -> AlsChm | None:
    """Read the Meta clip for ``site`` as a canopy height model.

    It comes back in the same container as the ALS CHM, because at ~1.2 m it is
    one: the fine-grid warp and block reduction in :mod:`bgt.coreg` then treat
    both identically, and every per-cell statistic is computed the same way on
    either side. Returns ``None`` when the site has no clip.
    """
    path = clip_path(site, meta_dir)
    if not path.is_file():
        print(f"  no Meta clip for {site} at {path} -- run 01c_download_meta.py")
        return None
    with rasterio.open(path) as src:
        # uint8 cannot hold nan; nothing is masked, since the product declares
        # no no-data and a 0 is a real prediction.
        height = src.read(1).astype(np.float32)
        return AlsChm(height=height, transform=src.transform, crs=src.crs,
                      name=path.name)


@dataclass(frozen=True)
class TargetGrid:
    """A grid to aggregate both CHMs onto.

    Carries the attributes :func:`bgt.coreg.build_fine_grid` reads from its
    target (``shape``, ``transform``, ``crs``, ``cell_size_m``).
    """

    transform: Affine
    crs: CRS
    shape: tuple[int, int]

    def cell_size_m(self) -> tuple[float, float]:
        return abs(self.transform.a), abs(self.transform.e)


def footprint_grid(chm: AlsChm, cell_m: float) -> TargetGrid:
    """A north-up ``cell_m`` grid over the ALS footprint, in the ALS CRS.

    Its origin is the ALS raster's own, so with a whole-metre ``cell_m`` the
    1 m fine grid built on it coincides with the ALS pixels and the ALS is
    copied onto it rather than resampled.
    """
    left, bottom, right, top = rasterio.transform.array_bounds(*chm.shape,
                                                               chm.transform)
    shape = (int(np.ceil((top - bottom) / cell_m)),
             int(np.ceil((right - left) / cell_m)))
    return TargetGrid(Affine(cell_m, 0.0, left, 0.0, -cell_m, top), chm.crs, shape)
