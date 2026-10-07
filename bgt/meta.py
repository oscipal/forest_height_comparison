"""The Meta global canopy height map v2 (DINOv3) as a comparison target.

~1.2 m canopy height predicted from Maxar very-high-resolution imagery,
published by Meta on AWS (``s3://dataforgood-fb-data/forests/v2/global/
dinov3_global_chm_v2_ml3/``). ``01c_download_meta.py`` clips it to each ALS
footprint; this module locates the tiles and names the clips.

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
from pathlib import Path

import numpy as np

import config


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

