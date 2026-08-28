"""Normalise a non-standard ALS product folder into the GCA layout.

Most sites arrive as Global Canopy Atlas output and are read directly. Some do
not: Amacayacu ships site-suffixed filenames, its masks stacked as bands of one
raster with the opposite polarity (1 marks a *flagged* pixel, where GCA marks an
invalid pixel with no-data), and no scan outline at all.

This step writes the missing GCA-named files beside the originals, so that every
later step sees one layout and needs no per-site special case. Nothing is
overwritten and no original is modified; re-running is a no-op unless
``--overwrite`` is given.

Usage:
    python 00_normalize_als.py --site Amacayacu
    python 00_normalize_als.py --site Amacayacu --overwrite
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio
from rasterio import features
from shapely.geometry import shape
from shapely.ops import unary_union

import config

#: Mask bands worth materialising, and the GCA filename each one takes.
MASK_BANDS = {
    "mask_combined": "mask_combined.tif",
    "mask_pd02": "mask_pd02.tif",
    "mask_pd04": "mask_pd04.tif",
    "mask_noground": "mask_noground.tif",
    "mask_unstabledtm": "mask_unstabledtm.tif",
}


def _single_band(src, band: int, data: np.ndarray, out_path: Path) -> None:
    profile = src.profile | {"count": 1, "dtype": "float32", "nodata": float("nan"),
                             "compress": "deflate"}
    with rasterio.open(out_path, "w", **profile) as dst:
        dst.write(data.astype(np.float32), 1)
        dst.set_band_description(1, out_path.stem)


def normalise_masks(als_dir: Path, stacked: Path, overwrite: bool) -> list[Path]:
    """Split a multiband mask raster into GCA-polarity single-band masks.

    The source marks a flagged pixel with 1 and a clean one with 0. The GCA
    convention -- the one `bgt.als.load_chm` implements -- is the opposite: a
    mask carries data where the pixel is usable and no-data where it is not. So
    the band is inverted, not merely copied: 1 becomes ``nan`` and 0 becomes 1.
    """
    written = []
    with rasterio.open(stacked) as src:
        names = list(src.descriptions)
        for name, filename in MASK_BANDS.items():
            if name not in names:
                print(f"  {filename}: no '{name}' band in {stacked.name}, skipped")
                continue
            out_path = als_dir / filename
            if out_path.is_file() and not overwrite:
                print(f"  {filename}: exists, kept")
                continue
            band = names.index(name) + 1
            flagged = src.read(band)
            keep = np.where(flagged == 0, np.float32(1.0), np.float32(np.nan))
            # Pixels the source left as no-data were never scanned; they carry
            # no ALS either way, so they stay no-data here too.
            keep[~np.isfinite(flagged)] = np.nan
            _single_band(src, band, keep, out_path)
            share = 100.0 * float(np.isfinite(keep).mean())
            print(f"  {filename}: written, {share:.1f}% of the grid usable")
            written.append(out_path)
    return written


def derive_outline(als_dir: Path, chm_path: Path, overwrite: bool) -> Path | None:
    """Trace the scan outline from the CHM's valid pixels.

    The BIOMASS catalogue search and the ETH clip both need a footprint polygon.
    Where the products ship none, the footprint of the canopy height model is
    the honest substitute: it is exactly the ground the comparison can use.
    The traced boundary is simplified to the raster's own resolution, which
    drops the pixel staircase without moving the edge more than one cell.
    """
    out_wgs84 = als_dir / config.ALS_OUTLINE_WGS84
    if out_wgs84.is_file() and not overwrite:
        print(f"  {out_wgs84.name}: exists, kept")
        return out_wgs84

    with rasterio.open(chm_path) as src:
        valid = np.isfinite(src.read(1)).astype(np.uint8)
        transform, crs, res = src.transform, src.crs, abs(src.transform.a)

    shapes = [shape(geom) for geom, value
              in features.shapes(valid, mask=valid.astype(bool), transform=transform)
              if value == 1]
    if not shapes:
        print("  no valid CHM pixels; cannot derive an outline")
        return None

    polygon = unary_union(shapes).simplify(res)
    local = gpd.GeoDataFrame({"site": [config.site_name(als_dir)]},
                             geometry=[polygon], crs=crs)
    local.to_file(als_dir / "outline_localCRS.shp")
    local.to_crs("EPSG:4326").to_file(out_wgs84)
    area_km2 = float(local.geometry.area.iloc[0]) / 1e6
    print(f"  {out_wgs84.name}: traced from {chm_path.name}, {area_km2:.1f} km2")
    return out_wgs84


def normalise_site(als_dir: Path, overwrite: bool) -> bool:
    site = als_dir.name.replace("03_processed_", "")
    print(f"\n=== {site} ===")

    chm_path = als_dir / config.ALS_CHM
    if not chm_path.is_file():
        suffixed = als_dir / f"{Path(config.ALS_CHM).stem}_{site}.tif"
        if not suffixed.is_file():
            print(f"  no {config.ALS_CHM} and no {suffixed.name}; nothing to do")
            return False
        if overwrite or not chm_path.is_file():
            chm_path.write_bytes(suffixed.read_bytes())
            print(f"  {config.ALS_CHM}: copied from {suffixed.name}")

    stacked = next((p for p in als_dir.glob("masks*.tif")), None)
    if (als_dir / config.ALS_MASK).is_file() and not overwrite:
        print(f"  {config.ALS_MASK}: exists, kept")
    elif stacked is not None:
        normalise_masks(als_dir, stacked, overwrite)
    else:
        print(f"  no {config.ALS_MASK} and no stacked masks*.tif found")

    derive_outline(als_dir, chm_path, overwrite)
    return True


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--site", action="append", default=None,
                   help="site folder to normalise (repeatable); default: every "
                        "03_processed_* folder that is not already GCA-shaped")
    p.add_argument("--overwrite", action="store_true",
                   help="rewrite files this step has already produced")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    candidates = sorted(p for p in config.ROOT.glob(config.ALS_DIR_GLOB) if p.is_dir())
    if args.site:
        wanted = {s.lower() for s in args.site}
        candidates = [p for p in candidates
                      if p.name.replace("03_processed_", "").lower() in wanted]
    if not candidates:
        print("No matching site folders.")
        return 1

    done = sum(normalise_site(p, args.overwrite) for p in candidates)
    print(f"\nNormalised {done} site folder(s)  ({datetime.now():%Y-%m-%d %H:%M})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
