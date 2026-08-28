"""Loading and masking of the airborne laser scanning (ALS) reference products."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
import pandas as pd
import numpy as np
import rasterio
from rasterio.crs import CRS

import config


@dataclass
class AlsChm:
    """A masked ALS canopy height model held in memory.

    ``height`` is float32 with ``nan`` wherever the pixel is missing, masked out
    by the quality masks, or outside the plausible height range.
    """

    height: np.ndarray
    transform: rasterio.Affine
    crs: CRS
    name: str

    @property
    def shape(self) -> tuple[int, int]:
        return self.height.shape  # type: ignore[return-value]

    @property
    def res(self) -> float:
        return abs(self.transform.a)

    @property
    def valid_fraction(self) -> float:
        return float(np.isfinite(self.height).mean())


def load_chm(
    als_dir: Path = config.ALS_DIR,
    chm_name: str = config.ALS_CHM,
    mask_name: str | None = config.ALS_MASK,
    extra_masks: list[str] | None = None,
    height_min: float = config.ALS_HEIGHT_MIN,
    height_max: float = config.ALS_HEIGHT_MAX,
) -> AlsChm:
    """Read a CHM and apply the GCA quality masks.

    In the GCA pipeline a mask layer marks *invalid* pixels by setting them to
    no-data, so a pixel is kept only where the mask raster is finite.
    """
    chm_path = als_dir / chm_name
    if not chm_path.is_file():
        raise FileNotFoundError(f"ALS CHM not found: {chm_path}")

    with rasterio.open(chm_path) as src:
        height = src.read(1, masked=True).filled(np.nan).astype(np.float32)
        transform, crs = src.transform, src.crs

    n_raw = int(np.isfinite(height).sum())

    for name in [m for m in [mask_name] if m] + list(extra_masks or []):
        mask_path = als_dir / name
        if not mask_path.is_file():
            raise FileNotFoundError(f"ALS mask not found: {mask_path}")
        with rasterio.open(mask_path) as src:
            if (src.shape != height.shape) or (src.transform != transform):
                raise ValueError(
                    f"{name} is not on the same grid as {chm_name}; "
                    "the GCA products are expected to share one grid."
                )
            mask = src.read(1, masked=True).filled(np.nan)
        height[~np.isfinite(mask)] = np.nan

    outside = np.isfinite(height) & ((height < height_min) | (height > height_max))
    height[outside] = np.nan

    n_kept = int(np.isfinite(height).sum())
    print(
        f"  {chm_name}: {n_kept:,} valid pixels of {height.size:,} "
        f"({100.0 * n_kept / height.size:.1f}%); masking removed "
        f"{n_raw - n_kept:,} pixels ({int(outside.sum()):,} out of range)"
    )

    return AlsChm(height=height, transform=transform, crs=crs, name=chm_name)


def load_outline(
    als_dir: Path = config.ALS_DIR,
    name: str = config.ALS_OUTLINE_WGS84,
) -> gpd.GeoDataFrame:
    """Read the scan outline shapefile (WGS84).

    Some sites ship an attribute table whose encoding disagrees with its .cpg
    file -- Mbalmayo declares UTF-8 but stores Latin-1 accents in the contact
    field -- which makes GDAL fail on the whole layer. Only the geometry is
    needed here, so fall back to reading that alone rather than failing on
    bytes that have no bearing on the outline.
    """
    path = als_dir / name
    try:
        return gpd.read_file(path)
    except UnicodeDecodeError:
        print(f"  {path.name}: undecodable attribute table, reading geometry only")
        return gpd.read_file(path, columns=[])


def search_bbox(
    als_dir: Path = config.ALS_DIR,
    buffer_deg: float = config.SEARCH_BUFFER_DEG,
) -> tuple[float, float, float, float]:
    """Bounding box of the ALS scan in WGS84, padded by ``buffer_deg``."""
    outline = load_outline(als_dir).to_crs("EPSG:4326")
    minx, miny, maxx, maxy = outline.total_bounds
    return (
        float(minx - buffer_deg),
        float(miny - buffer_deg),
        float(maxx + buffer_deg),
        float(maxy + buffer_deg),
    )


def load_summary(als_dir: Path = config.ALS_DIR) -> dict[str, str]:
    """Read the one-row ``summary_processing.csv`` shipped with the ALS products.

    Supplies the provenance stamped on every figure: which site was flown, when,
    and with what.
    """
    path = als_dir / "summary_processing.csv"
    if not path.is_file():
        return {}
    try:
        table = pd.read_csv(path, dtype=str, encoding="utf-8")
    except UnicodeDecodeError:
        table = pd.read_csv(path, dtype=str, encoding="latin-1")
    table = table.fillna("")
    return {} if table.empty else table.iloc[0].to_dict()


def acquisition_label(als_dir: Path = config.ALS_DIR) -> tuple[str, str]:
    """``(site name, acquisition date range)`` from the processing summary.

    ``mindt``/``maxdt`` are the first and last day on which points were actually
    recorded, which is the honest acquisition window (``acq_mindt``/``acq_maxdt``
    are the planned campaign dates).
    """
    summary = load_summary(als_dir)
    # Not every summary fills in the site name (Luki2025 leaves it blank), so
    # fall back to the product folder, which always carries it. A site flown
    # twice writes the same name into both summaries (Ipassa), which would pool
    # two overlapping scans into one sample, so there the folder label -- which
    # carries the acquisition -- wins over the summary.
    if not Path(als_dir).name.startswith("03_processed_"):
        site = config.site_name(als_dir)
    else:
        site = (
            summary.get("site")
            or summary.get("dataset")
            or Path(als_dir).name.replace("03_processed_", "")
        )
    start, end = summary.get("mindt", ""), summary.get("maxdt", "")
    if start and end:
        dates = start if start == end else f"{start} to {end}"
    else:
        dates = "date unknown"
    return site, dates
