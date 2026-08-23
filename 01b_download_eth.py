"""Step 1b -- clip the ETH global canopy height product to each ALS site.

A second spaceborne height product to compare against: Lang et al. (2023), 10 m
canopy top height from Sentinel-2 trained on GEDI, representative of 2020.

The published tiles are 3 x 3 degree cloud-optimised GeoTIFFs of ~400 MB each.
Because they are COGs, there is no need to fetch a whole tile: this script reads
only the window covering the ALS outline straight over HTTP and writes that clip
to disk, which is a few hundred kilobytes per site. No token and no account are
needed -- the data is CC BY 4.0.

Clips land in ``data/eth/<site>/``, alongside the predictive standard deviation
that the authors publish with the height.

Examples
--------
    python 01b_download_eth.py --list-only
    python 01b_download_eth.py
    python 01b_download_eth.py --site Luki2025 --no-sd
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import rasterio
from rasterio.merge import merge

import config
from bgt import als


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-dir", type=Path, default=config.ETH_DIR,
                   help="destination for the clips (default: %(default)s)")
    p.add_argument("--site", action="append", default=None,
                   help="restrict to this site (repeatable); default: every "
                        "03_processed_* folder found")
    p.add_argument("--buffer", type=float, default=config.ETH_BUFFER_DEG,
                   help="degrees added around the ALS outline (default: %(default)s)")
    p.add_argument("--no-sd", action="store_true",
                   help="skip the standard-deviation layer")
    p.add_argument("--list-only", action="store_true",
                   help="report the tiles each site needs, download nothing")
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args(argv)


def tile_name(lat: float, lon: float, step: int = config.ETH_TILE_DEG) -> str:
    """Name of the tile containing ``(lat, lon)``.

    Tiles are named after their south-west corner, so the corner is the
    coordinate floored to the tile grid: 17.05 E, 2.37 N lies in ``N00E015``.
    """
    lat0 = int(math.floor(lat / step) * step)
    lon0 = int(math.floor(lon / step) * step)
    ns = "N" if lat0 >= 0 else "S"
    ew = "E" if lon0 >= 0 else "W"
    return f"{ns}{abs(lat0):02d}{ew}{abs(lon0):03d}"


def tiles_for_bbox(bbox: tuple[float, float, float, float],
                   step: int = config.ETH_TILE_DEG) -> list[str]:
    """Every tile touched by ``(minx, miny, maxx, maxy)``, in a stable order."""
    minx, miny, maxx, maxy = bbox
    lon0 = math.floor(minx / step) * step
    lat0 = math.floor(miny / step) * step
    names = []
    lat = lat0
    while lat <= maxy:
        lon = lon0
        while lon <= maxx:
            names.append(tile_name(lat, lon, step))
            lon += step
        lat += step
    return sorted(dict.fromkeys(names))


def tile_url(tile: str, sd: bool = False) -> str:
    template = config.ETH_SD_TEMPLATE if sd else config.ETH_TILE_TEMPLATE
    return config.ETH_BASE_URL + template.format(tile=tile)


def clip_tiles(tiles: list[str], bbox: tuple[float, float, float, float],
               out_path: Path, sd: bool = False) -> Path | None:
    """Read ``bbox`` out of the remote tiles and write it to ``out_path``.

    Opening a COG over ``/vsicurl/`` fetches the header and only the blocks the
    window actually touches, so the transfer is proportional to the footprint
    rather than to the tile.
    """
    sources = []
    try:
        for tile in tiles:
            url = f"/vsicurl/{tile_url(tile, sd=sd)}"
            try:
                sources.append(rasterio.open(url))
            except rasterio.RasterioIOError as exc:
                print(f"    tile {tile} unreadable, skipped: {exc}")
        if not sources:
            return None

        # merge() handles the single-tile case too, and clips to bounds for us.
        data, transform = merge(sources, bounds=bbox)
        profile = sources[0].profile
    finally:
        for src in sources:
            src.close()

    valid = int((data[0] != profile.get("nodata", 255)).sum())
    if valid == 0:
        print("    the clip is entirely no-data -- check the site coordinates")

    profile.update(driver="GTiff", height=data.shape[1], width=data.shape[2],
                   count=1, transform=transform, compress="deflate",
                   tiled=True, blockxsize=256, blockysize=256)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(out_path, "w", **profile) as dst:
        dst.write(data[0], 1)
        dst.update_tags(
            source="ETH global canopy height 10 m, 2020 (Lang et al. 2023)",
            doi="10.1038/s41559-023-02206-6",
            licence="CC BY 4.0",
            tiles=",".join(tiles),
        )
    size_kb = out_path.stat().st_size / 1024
    print(f"    wrote {out_path.name}  {data.shape[2]} x {data.shape[1]} px, "
          f"{size_kb:,.0f} kB, {valid:,} valid pixels")
    return out_path


def select_sites(requested: list[str] | None) -> list[Path]:
    dirs = config.site_dirs()
    if not requested:
        return dirs
    by_name = {config.site_name(d): d for d in dirs}
    missing = [s for s in requested if s not in by_name]
    if missing:
        raise SystemExit(f"unknown site(s): {', '.join(missing)}; "
                         f"available: {', '.join(by_name)}")
    return [by_name[s] for s in requested]


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    site_dirs = select_sites(args.site)
    if not site_dirs:
        print("No 03_processed_* site folders found.")
        return 1

    for als_dir in site_dirs:
        site = config.site_name(als_dir)
        print(f"\n=== {site} ===")
        bbox = als.search_bbox(als_dir, buffer_deg=args.buffer)
        tiles = tiles_for_bbox(bbox)
        print(f"  outline bbox {tuple(round(v, 4) for v in bbox)}")
        print(f"  tile(s): {', '.join(tiles)}")
        if args.list_only:
            for tile in tiles:
                print(f"    {tile_url(tile)}")
            continue

        out_dir = args.out_dir / site
        wanted = [(out_dir / f"eth_canopy_height_{site}.tif", False)]
        if not args.no_sd:
            wanted.append((out_dir / f"eth_canopy_height_sd_{site}.tif", True))

        for path, is_sd in wanted:
            if path.is_file() and not args.overwrite:
                print(f"    {path.name} exists, skipped (--overwrite to refetch)")
                continue
            clip_tiles(tiles, bbox, path, sd=is_sd)

    if not args.list_only:
        print(f"\nClips are in {args.out_dir}")
        print("Next: python 02_compare.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
